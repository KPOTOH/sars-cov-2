import click
import pandas as pd
from ete3 import PhyloTree
import tqdm

from src.utils import node_parent


def read_tree(path: str, frmt=1) -> PhyloTree:
    tree = PhyloTree(path, format=frmt)
    return tree


def add_edge_level_to_table(substitutions: pd.DataFrame, tree_path: str):
    """ запускать, только если нужно пересчитать уровни для замещений """
    tree = read_tree(tree_path)
    edge_levels = dict()

    for leaf in tqdm.tqdm(tree.iter_leaves(), total=55000):
        cur_level = 1
        parent = leaf
        while parent != tree:
            parent = node_parent(leaf)
            pair = (parent.name, leaf.name)
            edge_levels[pair] = cur_level

            leaf = parent
            cur_level += 1

    assert len(edge_levels) == len(tree.get_descendants())  # E = V - 1

    substitutions = substitutions.copy()
    substitutions['edge_level'] = substitutions.apply(
        lambda r: edge_levels[(r.parent_node, r.child_node)], axis=1)

    return substitutions


def add_id_col(df: pd.DataFrame, col_name="id"):
    """ add index column at first position """
    df = df.copy()
    columns = list(df.columns)
    df[col_name] = df.index
    new_columns = [col_name] + columns
    df = df[new_columns]
    return df


@click.command("featirizer", help="add some features to mutations table")
@click.option("--input-table", required=True, help="path to input csv")
@click.option("--out-table", required=True, help="path to output mutations csv")
@click.option("--tree", required=True, help="path to tree in newick format")
def main(input_table: str, out_table: str, tree: str):
    df = pd.read_csv(input_table)
    df = add_edge_level_to_table(df, tree)
    df = add_id_col(df)
    df.to_csv(out_table, index=None)


if __name__ == "__main__":
    main()
