# using: script_name.py path_to_json [path_to_table]

import json
import sys
import os

import pandas as pd

# path_to_json = "../data/overall_mutations_2.json"


def read_json(path_to_json: str) -> list:
    assert os.path.exists(path_to_json), f"file {path_to_json} doesn't exist"

    with open(path_to_json) as fin:
        substitutions = json.load(fin)
    return substitutions


def converter(substitutions=None, path_to_json=None) -> pd.DataFrame:
    """convert json file with substitutions to csv-table

    params:
        - substitutions: list of substitutions
        - path_to_json: path to json file with substitutions

    substitutions structure:
    [
        [parent, child, pair_substs[
            pos(0-based), parent_nucl, child_cnucl
        ]]
    ]
    """
    assert substitutions is not None or path_to_json is not None, "pass something"
    substitutions = substitutions or read_json(path_to_json)

    scollection = []
    for parent, child, pair_substs in substitutions:
        for pos, pnuc, cnuc, pnuc_con, cnuc_con in pair_substs:
            scollection.append((pos, pnuc, cnuc, pnuc_con, cnuc_con, parent, child))

    cols = ["pos", "parent_nucl", "child_nucl", "parent_nucl_context", 
            "child_nucl_context", "parent_node", "child_node"]
    df = pd.DataFrame(scollection, columns=cols)
    return df


def main():
    if len(sys.argv) == 1:
        print("pass arguments", file=sys.stderr)
        print(
            f"using: {sys.argv[0]} path_to_json [path_to_table]",
            file=sys.stderr)
        sys.exit()

    elif len(sys.argv) == 2:
        path_to_json = sys.argv[1]
        assert path_to_json.endswith('.json'), "filename must end with .json"
        path_to_table = path_to_json.replace("json", "csv")

    elif len(sys.argv) == 3:
        path_to_json = sys.argv[1]
        path_to_table = sys.argv[2]

    print(f"path_to_json is {path_to_json}", file=sys.stderr)
    print(f"path_to_table is {path_to_table}", file=sys.stderr)

    df = converter(path_to_json=path_to_json)
    df.to_csv(path_to_table, index=None)
    print("done", file=sys.stderr)


if __name__ == "__main__":
    main()
