"""Entry point.

  Stool16S.exe                    opens the window
  Stool16S.exe --selftest out.txt runs a built-in check and writes the report to out.txt
                                  (used by the GitHub Actions build; exit code 0 = pass)
"""
import sys
import traceback


def selftest(report_path):
    lines = []
    ok = False
    try:
        import os
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        import pipeline

        sample = pipeline.resource_dir() / "data" / "selftest_sample.fasta"
        res = pipeline.run_pipeline(sample)
        by_id = {r[0]: r for r in res.rows}
        lines.append(f"Classified {res.total} sequences: {res.counts}")

        randoms = [r for r in res.rows if r[0].startswith("random_")]
        assert randoms, "no random test sequences were classified"
        assert all(r[1] == "Unknown" for r in randoms), "a random sequence got a bacterial label"

        ref = by_id.get("Faecalibacterium__prausnitzii__OP782675")
        assert ref is not None, "the reference test sequence is missing from the results"
        assert ref[1] == "Faecalibacterium", f"wrong label for reference test: {ref}"
        assert ref[2] == "Faecalibacterium prausnitzii", f"wrong species call: {ref}"
        assert float(ref[4]) >= 99.9, f"identity too low: {ref}"
        lines.append("Pipeline check passed (random = Unknown, reference = F. prausnitzii).")

        # Make sure the window itself can start and show results (catches missing Qt files).
        import gui
        from PySide6.QtWidgets import QApplication
        app = QApplication.instance() or QApplication([])
        win = gui.MainWindow()
        win.show_result(res)
        assert win.table.rowCount() == len(res.rows), "results table has the wrong number of rows"
        lines.append("Window check passed.")
        ok = True
    except Exception:
        lines.append("SELFTEST FAILED:\n" + traceback.format_exc())
    lines.append("RESULT: " + ("PASS" if ok else "FAIL"))
    text = "\n".join(lines) + "\n"
    if report_path:
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(text)
    else:
        print(text)
    return 0 if ok else 1


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == "--selftest":
        return selftest(sys.argv[2] if len(sys.argv) > 2 else None)
    import gui
    return gui.run()


if __name__ == "__main__":
    sys.exit(main())
