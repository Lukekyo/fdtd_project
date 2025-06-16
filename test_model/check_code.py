import os
import re
import pandas as pd

# 修改為你自己電腦上 repo 的路徑
repo_path = "C:\Pierre\Coding\Python\\fdtd\\fdtd"
print("Scanning path:", repo_path)
patterns = {
    "float_dtype": re.compile(r"dtype\s*=\s*(np\.)?(float32|float64)"),
    "real_access": re.compile(r"\.real"),
    "imag_access": re.compile(r"\.imag"),
    "np_zeros": re.compile(r"np\.zeros\s*\(.*\)"),
    "np_empty": re.compile(r"np\.empty\s*\(.*\)"),
    "bd_zeros": re.compile(r"bd\.zeros\s*\(.*\)"),
    "bd_empty": re.compile(r"bd\.empty\s*\(.*\)")
}

findings = {key: [] for key in patterns}

for root, dirs, files in os.walk(repo_path):
    for file in files:
        if file.endswith(".py"):
            filepath = os.path.join(root, file)
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    for key, pattern in patterns.items():
                        if pattern.search(line):
                            findings[key].append((filepath, i + 1, line.strip()))

# 輸出成 CSV 檔案
for key, records in findings.items():
    if records:
        df = pd.DataFrame(records, columns=["File", "Line", "Code"])
        df.to_csv(f"scan_{key}.csv", index=False)
        print(f"Exported scan_{key}.csv with {len(df)} findings.")