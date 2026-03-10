# flb_utils

Portable Python utilities for reading REKLAB `.flb` binary files, exporting
trials to CSV, and generating signal plots. Drop the `utils/` folder into any
project that needs to process `.flb` recordings.

---

## Contents

| File | Description |
|------|-------------|
| `flb_reader.py` | Low-level binary parser — `read_flb()` and `flb_to_csv()` |
| `main.py` | CLI entry point: reads a `.flb` file, exports CSVs, saves plots |
| `requirements.txt` | Python package dependencies |

---

## Requirements

Python 3.8+ is required.

```bash
pip install -r requirements.txt
```

| Package | Purpose |
|---------|---------|
| `numpy` | Array operations and binary parsing |
| `pandas` | Trial data as DataFrames, CSV export |
| `matplotlib` | Signal plots |

---

## Quick start

```bash
# Minimal usage — exports all CSVs + plots trial 1 + overlay
python main.py YES_281124.flb

# Plot a specific trial (e.g. trial 11)
python main.py YES_281124.flb --trial 11

# Plot trial 11, skip CSV export
python main.py YES_281124.flb --trial 11 --no_csv

# Export CSVs only, skip all plots
python main.py YES_281124.flb --no_plots

# Custom output directory
python main.py YES_281124.flb --output_dir ./results

# IES01 subject (different raw channel order)
python main.py YES_281124.flb --subject_id IES01 --trial 5

# Reduce overlay plot density (default downsample = 10)
python main.py YES_281124.flb --downsample 20 --no_csv
```

### All options

| Option | Default | Description |
|--------|---------|-------------|
| `flb_file` | *(required)* | Path to the `.flb` file |
| `--trial N` | `1` | Trial number to plot (1-based) |
| `--output_dir DIR` | `<name>_output/` | Root output folder |
| `--subject_id ID` | `default` | Use `IES01` for that subject's channel order |
| `--downsample N` | `10` | Keep every N-th sample in the overlay plot |
| `--no_csv` | off | Skip CSV export |
| `--no_plots` | off | Skip plot generation |

---

## Output structure

```
<output_dir>/
    csv/
        <name>_trial001.csv
        <name>_trial002.csv
        ...
    plots/
        trial01_all_channels.png   ← all 6 channels for Trial 1
        all_trials_overlay.png     ← all trials overlaid per channel
```

---

## Using as a library

```python
from flb_reader import read_flb, flb_to_csv

# Read trials into a list of DataFrames
trials = read_flb('recording.flb')
df = trials[0]          # columns: time, pos, tq, gm, gl, sol, ta
print(df.head())

# Export all trials to CSV
csv_paths = flb_to_csv('recording.flb', output_dir='./csv')
```

---

## Channel layout

| Column | Signal | Units |
|--------|--------|-------|
| `time` | Time axis | s |
| `pos`  | Ankle angle | rad |
| `tq`   | Torque | Nm |
| `gm`   | Gastrocnemius Medial EMG | mV |
| `gl`   | Gastrocnemius Lateral EMG | mV |
| `sol`  | Soleus EMG | mV |
| `ta`   | Tibialis Anterior EMG | mV |

> **Note:** Subject `IES01` has a different raw channel order. Pass
> `--subject_id IES01` (CLI) or `subject_id='IES01'` (API) to remap correctly.

---

## FLB format reference

The `.flb` format was reverse-engineered from `flbio.m` / `flb2mat.m`
(R. Kearney, McGill University). Supported versions: 2, 3, and 4.
See `flb_reader.py` for full format documentation.
