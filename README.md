# EDF Alignment Tools

A collection of Python utilities for working with EDF/EDF+ (European Data Format) files, focused on viewing, time alignment, and lag estimation between recordings.

## Overview

This toolkit provides three main utilities:

1. **edf_viewer.py** - Visualize EDF files and compare multiple recordings
2. **estimate_lag.py** - Estimate time lag between two EDF files using cross-correlation
3. **apply_lag_correction.py** - Apply time lag corrections to EDF files

These tools are particularly useful when working with simultaneous recordings from multiple devices that may have clock drift or synchronization offsets.

## Installation

### Requirements

- Python 3.11+
- Dependencies listed in `requirements.txt`

### Setup

```bash
# Clone or download this repository
cd edf-tools

# Create a virtual environment (recommended)
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Tools

### 1. EDF Viewer

View and compare EDF recordings with support for multiple files, time alignment, and annotations.

**Basic usage:**

```bash
# View a single EDF file
python edf_viewer.py recording.edf

# View specific time window (30 seconds starting at 60s)
python edf_viewer.py recording.edf --offset 60 --duration 30

# View specific channels only
python edf_viewer.py recording.edf --channels F7,F8

# Compare two recordings
python edf_viewer.py file1.edf file2.edf --duration 60

# Compare with time shift (shift file2 backward by 15.5 seconds)
python edf_viewer.py file1.edf file2.edf --shift -15.5 --duration 30
```

**Features:**
- Load and display all channels from EDF/EDF+ files
- Multi-file comparison with automatic time alignment
- Annotation display with color-coded markers
- Customizable time windows and channel selection
- Support for different sampling rates (automatic resampling for display)

**Options:**
- `--info` - Show file metadata without plotting
- `--duration N` - Plot N seconds of data
- `--offset N` - Start at N seconds into the recording
- `--channels CH1,CH2,...` - Plot only specified channels
- `--shift SECONDS` - Apply time shift for alignment (for multi-file mode)
- `--annotations-all-channels` - Show annotations on all subplots
- `--title "TEXT"` - Custom plot title

### 2. Estimate Lag

Estimate the time offset between two EDF recordings using cross-correlation analysis.

**Basic usage:**

```bash
# Estimate lag using channel F7
python estimate_lag.py file1.edf file2.edf --channel F7

# Use different channels from each file
python estimate_lag.py file1.edf file2.edf --channel1 F7 --channel2 F8

# Show diagnostic plots
python estimate_lag.py file1.edf file2.edf --channel F7 --plot

# Analyze correlation stability over time (detect drift)
python estimate_lag.py file1.edf file2.edf --channel F7 --stability-plot
```

**Output interpretation:**

The tool reports the lag in seconds:
- **Negative lag**: File 2 is ahead of File 1 (events occur earlier)
- **Positive lag**: File 2 is behind File 1 (events occur later)

Example output:
```
Estimated lag: -15.234 s
  • File 2 is 15.234 seconds AHEAD of File 1
  • To align File 2 with File 1, shift by -15.234 seconds
```

**Options:**
- `--channel NAME` - Use same channel from both files
- `--channel1 NAME` / `--channel2 NAME` - Use different channels
- `--duration N` - Use only first N seconds
- `--max-lag N` - Maximum lag to search in seconds (default: 120)
- `--plot` - Show cross-correlation plots
- `--stability-plot` - Analyze correlation quality and drift over time
- `--no-bootstrap` - Disable confidence interval calculation (faster)
- `--quiet` - Output only the lag value (for scripting)
- `--output FILE.json` - Save results to JSON

### 3. Apply Lag Correction

Apply a time offset correction to an EDF file by adjusting the start time and trimming samples.

**Basic usage:**

```bash
# Apply correction (typically from estimate_lag.py output)
python apply_lag_correction.py input.edf -15.234 -o corrected.edf

# Verbose output showing details
python apply_lag_correction.py input.edf -15.234 -o corrected.edf --verbose
```

**How it works:**

1. Adjusts the EDF start time metadata by the specified lag
2. Trims signal samples to align with whole-second boundaries (EDF limitation)
3. Adjusts annotations to maintain correct absolute timing

**Note:** EDF format only supports 1-second resolution for start times. Fractional seconds are handled by trimming the appropriate number of samples from the beginning of each signal.

**Options:**
- `-o, --output FILE` - Output file path (required)
- `-v, --verbose` - Show detailed processing information

## Typical Workflow

When working with two synchronized recordings that may have time offsets:

```bash
# Step 1: Estimate the lag between recordings
python estimate_lag.py device1.edf device2.edf --channel F7 --plot

# Output: Estimated lag: -15.234 s

# Step 2: Verify the alignment visually
python edf_viewer.py device1.edf device2.edf --shift -15.234 --duration 60

# Step 3: Apply correction to create aligned file
python apply_lag_correction.py device2.edf -15.234 -o device2_aligned.edf

# Step 4: Verify the corrected alignment
python edf_viewer.py device1.edf device2_aligned.edf --duration 60
```

## Understanding Cross-Correlation Lag Estimation

The `estimate_lag.py` tool uses cross-correlation to find the time offset between two signals:

1. **Load signals**: Extracts specified channels from both EDF files
2. **Find overlap**: Determines the overlapping time period based on file timestamps
3. **Preprocess**: Applies z-score normalization (data should be pre-filtered before analysis)
4. **Correlate**: Computes cross-correlation and finds the lag with maximum correlation
5. **Refine**: Uses parabolic interpolation for sub-sample accuracy
6. **Validate**: Optionally computes bootstrap confidence intervals and stability analysis

**Key parameters:**

- **Max lag search**: Limits the search range (default ±120s) for efficiency
- **Bootstrap**: Estimates confidence intervals by resampling signal windows
- **Stability analysis**: Detects time-varying drift by analyzing correlation in sliding windows

**Quality metrics:**

- **Peak correlation**: Higher is better (>0.6 good, >0.8 excellent)
- **Confidence interval**: Narrower intervals indicate more reliable estimates
- **Lag variability**: Low standard deviation across windows indicates stable alignment
- **Drift rate**: Should be near zero for constant time offset

## File Format Details

**EDF/EDF+ format:**
- European Data Format for multi-channel time series
- Supports multiple channels with different sampling rates
- Start time resolution: 1 second (fractional seconds handled by sample alignment)
- Can include annotations (EDF+) for marking events

**Supported features:**
- Multi-channel signals with independent sampling rates
- Annotations (stimulation markers, events, etc.)
- Physical units and calibration metadata
- Patient and recording identification fields

## Troubleshooting

**Low correlation (<0.3)**
- Signals may not be from the same session
- Try different channels
- Ensure signals are properly pre-filtered (e.g., bandpass 0.5-35 Hz for EEG)
- Check that files actually overlap in time

**Bootstrap confidence interval fails**
- Recording too short (need >90 seconds typically)
- Use `--no-bootstrap` to skip CI calculation

**Stability analysis shows high drift**
- Time offset is not constant (clock drift)
- Consider using shorter analysis segments
- May need more sophisticated drift correction

**File has no overlapping time**
- Check file start times with `--info` mode
- Files may be from different sessions
- Use `edf_viewer.py file.edf --info` to inspect timestamps

## License

This project is provided as-is for research and educational purposes.
