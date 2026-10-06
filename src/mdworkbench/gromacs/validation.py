from pathlib import Path
import math
import re


def validate_gro(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) < 4:
        raise ValueError("Truncated GRO")
    count = int(lines[1])
    if count <= 0 or len(lines) != count + 3:
        raise ValueError("GRO atom count mismatch")
    for line in lines[2:2 + count]:
        values = [float(line[a:b]) for a, b in ((20, 28), (28, 36), (36, 44))]
        if not all(math.isfinite(x) for x in values):
            raise ValueError("Nonfinite GRO position")
    box = [float(x) for x in lines[-1].split()]
    if len(box) not in (3, 9) or not all(math.isfinite(x) for x in box) or any(x <= 0 for x in box[:3]):
        raise ValueError("GRO requires a finite positive periodic box")
    return count


def mdp_values(path):
    values = {}
    for line in Path(path).read_text().splitlines():
        line = line.split(";", 1)[0].strip()
        if "=" in line:
            key, value = line.split("=", 1)
            key = key.strip().lower().replace("_", "-")
            if key in values:
                raise ValueError(f"Duplicate MDP key: {key}")
            values[key] = value.strip()
    return values


INCLUDE = re.compile(r'^\s*#include\s+["<]([^">]+)[">]', re.MULTILINE)
