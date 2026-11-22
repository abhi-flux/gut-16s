# Stool 16S Classifier (Windows app)

Windows desktop app around the original `classify.py` and `clean_sample.py` pipeline
(both included **unchanged** in `app/`). The user picks a stool sample file (FASTA or
FASTQ), presses Run, and gets a table of results that can be saved as CSV.

## For the user
1. Unzip the download and open the `Stool16S` folder.
2. Double-click `Stool16S.exe` (keep the whole folder together).
3. Choose or drag in your sample file, press **Run**, then **Save results as CSV...**.

Windows may show a blue "Windows protected your PC" box because the app is not code-signed.
Click **More info**, then **Run anyway**.

No Python, BLAST or internet connection is needed. The 197 reference sequences and the 21
targets are built in.

## For the developer
- Get the app: push this folder to a GitHub repository (branch `main`). Open the **Actions**
  tab, wait for "Build Windows app" to turn green, open the run, and download the
  **Stool16S-windows** artifact.
- The build downloads BLAST+ from NCBI, builds the app with PyInstaller, then starts the
  built `.exe` in `--selftest` mode. The run fails if random sequences are not "Unknown"
  or if the real *F. prausnitzii* reference is not labelled correctly.
- Run from source (needs Python, `pip install -r requirements.txt`, and BLAST+ on PATH):
  `python app/main.py`

| Path | What it is |
|---|---|
| `app/classify.py`, `app/clean_sample.py` | Original scripts, unchanged |
| `app/pipeline.py` | Runs the two scripts in a temp folder and returns the results |
| `app/gui.py` | The window |
| `app/main.py` | Entry point and `--selftest` |
| `data/` | `references.fasta`, `targets.txt`, self-test sample |
| `blast_bin/` | Filled by the build with Windows BLAST+ |
| `stool16s.spec` | PyInstaller recipe |

## Known limits (carried over, to fix in the next round)
- No references for *Eubacterium*, *Prevotella*, *Campylobacter*, so those can never be labelled.
- *Escherichia* and *Shigella* cannot be told apart by 16S.
- Reference set contains only the targets (no outgroup), thresholds are not validated.
- Counts are sequences per label, not bacterial abundance.
