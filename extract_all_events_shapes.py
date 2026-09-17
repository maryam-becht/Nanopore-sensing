# ============================================================
# EXTRACT ALL NANPORE EVENT SHAPES
# ============================================================

print("==========================================")
print("STARTING EVENT SHAPE EXTRACTION")
print("==========================================")

from pathlib import Path
import numpy as np
import pandas as pd
import pyabf


# ============================================================
# CHANGE ONLY THIS PATH
# ============================================================

EXCEL_PATH = Path(
    r"C:\Users\mbech\Documents\Peptides\asynpY125\asynpY125_P3_nanopore_R_data.xlsx"
)

# ============================================================
# SETTINGS
# ============================================================

# Amount of baseline shown before/after each event
PADDING_BEFORE_MS = 1.0
PADDING_AFTER_MS = 1.0


# ============================================================
# CHECK EXCEL
# ============================================================

print("\n1. Checking Excel file...")

if not EXCEL_PATH.exists():
    raise FileNotFoundError(
        f"\nExcel file not found:\n{EXCEL_PATH}\n"
        "\nChange EXCEL_PATH at the beginning of the script."
    )

print("Excel found:")
print(EXCEL_PATH)


# ============================================================
# READ EVENTS
# ============================================================

print("\n2. Reading events...")

events = pd.read_excel(
    EXCEL_PATH,
    sheet_name="events"
)

print(f"{len(events)} filtered events found.")

required_columns = [
    "event_id",
    "file_name",
    "channel",
    "sweep",
    "start_idx_in_sweep",
    "end_idx_in_sweep",
    "baseline_local_pA",
    "dwell_total_ms",
    "I_over_I0_pct"
]

missing = [
    col for col in required_columns
    if col not in events.columns
]

if missing:
    raise RuntimeError(
        f"Missing columns in events sheet: {missing}"
    )


# ============================================================
# READ ABF PATH FROM EXCEL
# ============================================================

print("\n3. Looking for ABF path...")

discovery = pd.read_excel(
    EXCEL_PATH,
    sheet_name="discovery_audit"
)

if "abf_file" not in discovery.columns:
    raise RuntimeError(
        "'abf_file' missing from discovery_audit."
    )

if "abf_path" not in discovery.columns:
    raise RuntimeError(
        "'abf_path' missing from discovery_audit."
    )

abf_mapping = {}

for _, row in discovery.iterrows():

    if pd.isna(row["abf_file"]) or pd.isna(row["abf_path"]):
        continue

    filename = str(row["abf_file"])
    path = Path(str(row["abf_path"]))

    abf_mapping[filename] = path


# ============================================================
# CHECK ABF FILES
# ============================================================

for filename in events["file_name"].unique():

    filename = str(filename)

    if filename not in abf_mapping:
        raise FileNotFoundError(
            f"No ABF path found for {filename}"
        )

    path = abf_mapping[filename]

    print(f"\nABF expected:")
    print(path)

    if not path.exists():

        print("\nWARNING")
        print("The ABF path stored in Excel is not accessible.")
        print("You need to change the ABF location.")

        # Ask for folder manually
        folder_text = input(
            "\nPaste the folder containing the ABF file "
            "and press Enter:\n"
        ).strip().strip('"')

        new_path = Path(folder_text) / filename

        if not new_path.exists():
            raise FileNotFoundError(
                f"\nABF still not found:\n{new_path}"
            )

        abf_mapping[filename] = new_path

    print("ABF found!")


# ============================================================
# FUNCTION: CONVERT CURRENT TO pA
# ============================================================

def convert_to_pA(signal, unit):

    signal = np.asarray(
        signal,
        dtype=np.float64
    )

    unit = str(unit).strip().lower()

    if unit == "pa":
        return signal

    elif unit == "na":
        return signal * 1000

    elif unit in ["ua", "µa"]:
        return signal * 1_000_000

    elif unit == "ma":
        return signal * 1_000_000_000

    elif unit == "a":
        return signal * 1_000_000_000_000

    else:
        print(
            f"WARNING: unknown current unit '{unit}'."
        )
        return signal


# ============================================================
# EXTRACT EVENT SHAPES
# ============================================================

print("\n4. Extracting event shapes...")

all_event_shapes = []

groups = events.groupby(
    [
        "file_name",
        "channel",
        "sweep"
    ],
    sort=False
)

for (
    file_name,
    channel,
    sweep
), group_events in groups:

    file_name = str(file_name)
    channel = int(channel)
    sweep = int(sweep)

    abf_path = abf_mapping[file_name]

    print("\n------------------------------------------")
    print(f"ABF: {file_name}")
    print(f"Channel: {channel}")
    print(f"Sweep: {sweep}")
    print(f"Events: {len(group_events)}")
    print("------------------------------------------")

    # Load ABF
    abf = pyabf.ABF(
        str(abf_path)
    )

    # Select correct sweep + channel
    abf.setSweep(
        sweepNumber=sweep,
        channel=channel
    )

    fs = float(
        abf.dataRate
    )

    print(
        f"Sampling rate: {fs:.0f} Hz"
    )

    # Actual current trace
    current_pA = convert_to_pA(
        abf.sweepY,
        abf.sweepUnitsY
    )

    # Convert padding from ms to samples
    pad_before = int(
        round(
            PADDING_BEFORE_MS *
            fs /
            1000
        )
    )

    pad_after = int(
        round(
            PADDING_AFTER_MS *
            fs /
            1000
        )
    )

    # ----------------------------------------
    # EACH EVENT
    # ----------------------------------------

    for count, (_, event) in enumerate(
        group_events.iterrows(),
        start=1
    ):

        event_id = int(
            event["event_id"]
        )

        start_idx = int(
            event["start_idx_in_sweep"]
        )

        end_idx = int(
            event["end_idx_in_sweep"]
        )

        # Workbook indices are inclusive
        end_exclusive = end_idx + 1

        # Add baseline before and after
        left = max(
            0,
            start_idx - pad_before
        )

        right = min(
            len(current_pA),
            end_exclusive + pad_after
        )

        sample_indices = np.arange(
            left,
            right
        )

        # THIS IS THE REAL ABF CURRENT
        event_current = current_pA[
            left:right
        ]

        # Time relative to event onset
        time_relative_ms = (
            sample_indices - start_idx
        ) / fs * 1000

        # Which samples belong to event itself
        inside_event = (
            (sample_indices >= start_idx) &
            (sample_indices < end_exclusive)
        )

        baseline = float(
            event["baseline_local_pA"]
        )

        # Instantaneous I/I0
        if (
            np.isfinite(baseline)
            and baseline != 0
        ):

            trace_ratio = (
                event_current /
                baseline *
                100
            )

        else:

            trace_ratio = np.full(
                len(event_current),
                np.nan
            )

        event_df = pd.DataFrame({

            "event_id":
                event_id,

            "file_name":
                file_name,

            "channel":
                channel,

            "sweep":
                sweep,

            "sample_index":
                sample_indices,

            "time_relative_ms":
                time_relative_ms,

            "current_pA":
                event_current,

            "baseline_pA":
                baseline,

            "I_over_I0_trace_pct":
                trace_ratio,

            "inside_event":
                inside_event.astype(int),

            "dwell_total_ms":
                float(
                    event["dwell_total_ms"]
                ),

            "event_I_over_I0_pct":
                float(
                    event["I_over_I0_pct"]
                )
        })

        all_event_shapes.append(
            event_df
        )

        # Progress
        if (
            count % 50 == 0
            or count == len(group_events)
        ):
            print(
                f"Extracted "
                f"{count}/{len(group_events)} events"
            )


# ============================================================
# COMBINE
# ============================================================

print("\n5. Combining events...")

event_shapes = pd.concat(
    all_event_shapes,
    ignore_index=True
)

event_shapes = event_shapes.sort_values(
    [
        "event_id",
        "time_relative_ms"
    ]
)


# ============================================================
# SAVE CSV
# ============================================================

OUTPUT_PATH = EXCEL_PATH.with_name(
    EXCEL_PATH.stem +
    "_event_shapes.csv"
)

event_shapes.to_csv(
    OUTPUT_PATH,
    index=False
)


# ============================================================
# FINAL CHECK
# ============================================================

n_events_exported = (
    event_shapes["event_id"]
    .nunique()
)

print("\n==========================================")
print("DONE!")
print("==========================================")

print(
    f"\nEvents exported: "
    f"{n_events_exported}"
)

print(
    f"\nTotal trace samples exported: "
    f"{len(event_shapes)}"
)

print(
    "\nOutput file:"
)

print(
    OUTPUT_PATH
)

print(
    "\nYou can now use this CSV "
    "in the R Shiny Event Explorer."
)