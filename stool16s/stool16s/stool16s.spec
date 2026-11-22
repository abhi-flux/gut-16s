# -*- mode: python ; coding: utf-8 -*-
# PyInstaller build recipe. Build with:  pyinstaller --noconfirm --clean stool16s.spec
import os

datas = [("data", "data")]
if os.path.isdir("blast_bin"):
    datas.append(("blast_bin", "blast_bin"))

a = Analysis(
    ["app/main.py"],
    pathex=["app"],
    datas=datas,
    hiddenimports=["Bio.Entrez", "Bio.SeqIO", "classify", "clean_sample", "pipeline", "gui"],
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="Stool16S",
    console=False,  # no black terminal window
)
coll = COLLECT(exe, a.binaries, a.datas, name="Stool16S")
