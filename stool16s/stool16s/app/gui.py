"""Windows desktop window for the stool 16S classifier (PySide6)."""
import html
import threading
import traceback
from pathlib import Path

from PySide6.QtCore import QObject, Qt, Signal
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QFileDialog, QHBoxLayout, QHeaderView, QLabel,
    QLineEdit, QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget,
)

import pipeline

APP_NAME = "Stool 16S Classifier"
FILE_FILTER = "Sequence files (*.fasta *.fa *.fna *.fas *.fastq *.fq *.txt);;All files (*)"


class _Signals(QObject):
    status = Signal(str)
    done = Signal(object)
    failed = Signal(str, str)


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(APP_NAME)
        self.resize(920, 680)
        self.setAcceptDrops(True)
        self.result = None
        self.busy = False

        self.signals = _Signals()
        self.signals.status.connect(self._on_status)
        self.signals.done.connect(self.show_result)
        self.signals.failed.connect(self._on_failed)

        root = QVBoxLayout(self)
        root.setContentsMargins(16, 16, 16, 16)
        root.setSpacing(10)

        title = QLabel(f"<h2 style='margin:0'>{APP_NAME}</h2>")
        hint = QLabel("Choose your stool sample file (FASTA or FASTQ), or drag it onto this window, "
                      "then press Run.")
        hint.setWordWrap(True)
        root.addWidget(title)
        root.addWidget(hint)

        row = QHBoxLayout()
        self.path_edit = QLineEdit()
        self.path_edit.setPlaceholderText("No file chosen")
        self.path_edit.setReadOnly(True)
        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.clicked.connect(self.browse)
        row.addWidget(self.path_edit, 1)
        row.addWidget(self.browse_btn)
        root.addLayout(row)

        row2 = QHBoxLayout()
        self.run_btn = QPushButton("Run")
        self.run_btn.setMinimumHeight(34)
        self.run_btn.setStyleSheet("font-weight: bold;")
        self.run_btn.clicked.connect(self.start_run)
        self.save_btn = QPushButton("Save results as CSV...")
        self.save_btn.setMinimumHeight(34)
        self.save_btn.setEnabled(False)
        self.save_btn.clicked.connect(self.save_csv)
        row2.addWidget(self.run_btn, 1)
        row2.addWidget(self.save_btn, 1)
        root.addLayout(row2)

        self.progress = QProgressBar()
        self.progress.setRange(0, 0)  # busy animation
        self.progress.setVisible(False)
        self.status_label = QLabel("Ready.")
        root.addWidget(self.progress)
        root.addWidget(self.status_label)

        self.summary = QLabel("")
        self.summary.setTextFormat(Qt.RichText)
        self.summary.setWordWrap(True)
        root.addWidget(self.summary)

        self.table = QTableWidget(0, 0)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        root.addWidget(self.table, 1)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setPlaceholderText("Details (cleaning and classification log) appear here.")
        self.log.setMaximumHeight(130)
        root.addWidget(self.log)

    # --- choosing a file -----------------------------------------------------------
    def browse(self):
        start = str(Path(self.path_edit.text()).parent) if self.path_edit.text() else ""
        path, _ = QFileDialog.getOpenFileName(self, "Choose sample file", start, FILE_FILTER)
        if path:
            self.path_edit.setText(path)

    def dragEnterEvent(self, event):
        if not self.busy and event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        urls = event.mimeData().urls()
        if urls and urls[0].isLocalFile():
            self.path_edit.setText(urls[0].toLocalFile())

    # --- running -------------------------------------------------------------------
    def _set_busy(self, busy):
        self.busy = busy
        self.run_btn.setEnabled(not busy)
        self.browse_btn.setEnabled(not busy)
        self.save_btn.setEnabled(not busy and self.result is not None)
        self.progress.setVisible(busy)

    def start_run(self):
        path = self.path_edit.text().strip()
        if not path or not Path(path).is_file():
            QMessageBox.warning(self, APP_NAME, "Please choose a sample file first.")
            return
        self.result = None
        self.table.setRowCount(0)
        self.table.setColumnCount(0)
        self.summary.setText("")
        self.log.clear()
        self._set_busy(True)
        self.status_label.setText("Starting...")

        def work():
            try:
                res = pipeline.run_pipeline(path, status=self.signals.status.emit)
            except pipeline.PipelineError as e:
                self.signals.failed.emit(str(e), "")
            except Exception:
                self.signals.failed.emit("Something unexpected went wrong.", traceback.format_exc())
            else:
                self.signals.done.emit(res)

        threading.Thread(target=work, daemon=True).start()

    def _on_status(self, text):
        self.status_label.setText(text)

    def _on_failed(self, message, details):
        self._set_busy(False)
        self.status_label.setText("Failed.")
        self.log.setPlainText((message + "\n\n" + details).strip())
        QMessageBox.critical(self, APP_NAME, message)

    # --- showing results -----------------------------------------------------------
    def show_result(self, res):
        self.result = res
        self._set_busy(False)
        self.status_label.setText(f"Done. {res.total} sequences classified from {res.sample_name}.")

        lines = [f"<b>{res.total}</b> sequences classified:"]
        for label, n in res.counts:
            pct = 100 * n / res.total if res.total else 0
            lines.append(f"&nbsp;&nbsp;{html.escape(label)}: <b>{n}</b> ({pct:.1f}%)")
        self.summary.setText("<br>".join(lines))

        self.table.setColumnCount(len(res.header))
        self.table.setHorizontalHeaderLabels(res.header)
        self.table.setRowCount(len(res.rows))
        grey = QBrush(QColor("#888888"))
        label_col = res.header.index("label")
        for r, row in enumerate(res.rows):
            for c, value in enumerate(row):
                item = QTableWidgetItem(value)
                if row[label_col] == "Unknown":
                    item.setForeground(grey)
                self.table.setItem(r, c, item)
        self.table.resizeColumnsToContents()
        self.table.horizontalHeader().setStretchLastSection(True)

        self.log.setPlainText(("Cleaning:\n" + res.clean_log + "\nClassification:\n" + res.classify_log).strip())
        self.save_btn.setEnabled(True)

    def save_csv(self):
        if self.result is None:
            return
        default = str(Path(self.path_edit.text()).with_name(Path(self.path_edit.text()).stem + "_results.csv"))
        path, _ = QFileDialog.getSaveFileName(self, "Save results", default, "CSV files (*.csv)")
        if not path:
            return
        if not path.lower().endswith(".csv"):
            path += ".csv"
        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                f.write(self.result.csv_text)
        except OSError as e:
            QMessageBox.critical(
                self, APP_NAME,
                f"Could not save the file:\n{e}\n\nIf it is open in Excel, close it and try again.")
            return
        self.status_label.setText(f"Saved: {path}")


def run():
    app = QApplication.instance() or QApplication([])
    app.setApplicationName(APP_NAME)
    win = MainWindow()
    win.show()
    return app.exec()
