# Quick Start Guide

## Installation

```bash
cd edf-tools
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Quick Examples

### 1. View an EDF file

```bash
python edf_viewer.py your_recording.edf
```

### 2. Compare two recordings

```bash
# View both files overlaid
python edf_viewer.py recording1.edf recording2.edf --duration 60

# Estimate the time lag between them
python estimate_lag.py recording1.edf recording2.edf --channel F7

# Example output: -15.234 s (recording2 is ahead)

# Visualize with the correction applied
python edf_viewer.py recording1.edf recording2.edf --shift -15.234 --duration 60
```

### 3. Create time-corrected file

```bash
# Apply the estimated lag to create aligned file
python set_start_time.py recording2.edf --adjust-by -15.234 -o recording2_aligned.edf --verbose

# Verify the correction
python edf_viewer.py recording1.edf recording2_aligned.edf --duration 60
```

## Common Options

**edf_viewer.py**
- `--duration 60` - Show 60 seconds
- `--offset 30` - Start at 30 seconds
- `--channels F7,F8` - Show only these channels
- `--shift -15.2` - Shift file(s) by this many seconds

**estimate_lag.py**
- `--channel F7` - Use F7 channel from both files
- `--plot` - Show correlation plots
- `--stability-plot` - Analyze drift over time
- `--quiet` - Output only the lag value

**set_start_time.py**
- `--adjust-by SECONDS` - Adjust by lag offset
- `--set-time DATETIME` - Set absolute start time
- `-o output.edf` - Output file (required)
- `--verbose` - Show details

## Troubleshooting

**"Channel not found"**
- Run with `--info` to see available channels
- Use `--channel1` and `--channel2` for different channel names

**"No time overlap"**
- Files are from different sessions
- Check timestamps with `--info` mode

**Low correlation warning**
- Ensure data is properly pre-filtered (e.g., bandpass 0.5-35 Hz for EEG)
- Try different channels

## Need Help?

See the full README.md for detailed documentation and advanced usage.
