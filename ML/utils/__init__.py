# utils package – FLB file reader and signal processing helpers
from .flb_reader import read_flb, flb_to_csv
from .emg_envelope import (extract_envelope, compute_mvc_max,
                            process_trials, normalize_with_max, detect_emg_columns)
from .dataset_builder import (build_dataset, classify_trials,
                               get_passive_torque_map, downsample_trial, build_windows)
