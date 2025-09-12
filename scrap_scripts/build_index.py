import pickle
from pathlib import Path
from collections import defaultdict


def build_file_index(root_dir: Path) -> dict[str, list[Path]]:
    index = defaultdict(list)
    for path in root_dir.rglob("*"):
        if path.is_file():
            index[path.name].append(path.resolve())
    return index


# build index
root = Path("/run/media/james/T7/globdb_r226")
file_index = build_file_index(root)

# write index as pkl file
with open("./globd_r226_index.pkl", "wb") as f:
    print("Writing pickle file...")
    pickle.dump(file_index, f)
    print("Done writing pickle file!")
