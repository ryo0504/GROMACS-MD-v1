"""Export a portable version inventory and the non-Conda pip supplement only."""
from importlib import metadata
from pathlib import Path
import json
import sys

output = Path(sys.argv[1])
platform = sys.argv[2]
prefix = Path(sys.prefix)
owned = set()
for file in (prefix / "conda-meta").glob("*.json"):
    for name in json.loads(file.read_text()).get("files", []):
        if ".dist-info/" in name:
            owned.add(name.split(".dist-info/", 1)[0] + ".dist-info")
inventory = {}
supplement = []
for dist in metadata.distributions():
    name = dist.metadata.get("Name")
    if not name:
        continue
    inventory[name] = dist.version
    # locate METADATA through the public distribution file list.
    meta = next((p for p in dist.files or [] if str(p).endswith(".dist-info/METADATA")), None)
    if meta is None or name == "mdworkbench":
        continue
    directory = Path(dist.locate_file(meta)).parent
    try:
        relative = directory.relative_to(prefix).as_posix()
    except ValueError:
        continue
    if relative not in owned:
        supplement.append(f"{name}=={dist.version}")
(output / f"{platform}-packages.json").write_text(json.dumps(dict(sorted(inventory.items())), indent=2) + "\n")
(output / f"{platform}-pip.txt").write_text(
    "# Host-specific pip supplement. Recreate the matching Conda lock FIRST.\n"
    "# Install with --no-deps; never apply across incompatible PyTorch backends.\n" +
    "\n".join(sorted(supplement, key=str.lower)) + "\n")
