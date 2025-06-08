import pandas as pd
import re
import sys
from pathlib import Path

def read_txt_to_df(input_file):
    """
    Reads a space-aligned TXT file and returns a pandas DataFrame.
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()
    # Skip comment and blank lines
    data_lines = [line.strip() for line in lines if line.strip() and not line.startswith("#")]

    # Define the column headers
    columns = [
        "target_name", "accession_1", "query_name", "accession_2",
        "E-value_full", "score_full", "bias_full",
        "E-value_best", "score_best", "bias_best",
        "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description"
    ]

    # Parse lines by splitting on 2+ spaces
    rows = []
    for line in data_lines:
        parts = re.split(r'\s{2,}', line)
        if len(parts) < len(columns):
            parts += [''] * (len(columns) - len(parts))
        rows.append(parts[:len(columns)])

    return pd.DataFrame(rows, columns=columns)

def convert_txts_to_excel(output_file, input_files):
    """
    Writes multiple DataFrames (from TXT files) to an Excel workbook, each on its own sheet.
    """
    with pd.ExcelWriter(output_file) as writer:
        for txt_file in input_files:
            df = read_txt_to_df(txt_file)
            # Use file stem as sheet name (max 31 chars)
            sheet_name = Path(txt_file).stem[:31]
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    print(f"✅ Excel workbook saved to: {output_file}")

def main():
    """
    Usage:
      python multi_txt_to_excel_converter.py output.xlsx
      python multi_txt_to_excel_converter.py output.xlsx path/to/directory
      python multi_txt_to_excel_converter.py output.xlsx file1.txt file2.txt ...
    """
    if len(sys.argv) == 2:
        # No inputs specified: use all .txt in current directory
        output_excel = sys.argv[1]
        input_txts = sorted(str(p) for p in Path.cwd().glob("*.txt"))
    elif len(sys.argv) == 3 and Path(sys.argv[2]).is_dir():
        # Directory specified
        output_excel = sys.argv[1]
        input_txts = sorted(str(p) for p in Path(sys.argv[2]).glob("*.txt"))
    elif len(sys.argv) >= 3:
        # Explicit list of files
        output_excel = sys.argv[1]
        input_txts = sys.argv[2:]
    else:
        print(main.__doc__)
        sys.exit(1)

    if not input_txts:
        print("❌ No .txt files found to process.")
        sys.exit(1)

    convert_txts_to_excel(output_excel, input_txts)

if __name__ == "__main__":
    main()
