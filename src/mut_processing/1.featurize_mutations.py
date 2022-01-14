import click
import pandas as pd
from ete3 import PhyloTree
import tqdm
from queue import Queue

from utils import node_parent


def read_tree(path: str, frmt=1) -> PhyloTree:
    tree = PhyloTree(path, format=frmt)
    return tree


def _add_dist2terminals(substitutions: pd.DataFrame, tree: PhyloTree):
    """ говно тупого говна """
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


def add_dist2root(substitutions: pd.DataFrame, tree: PhyloTree) -> pd.DataFrame:
    """ 
    calculate distance to root from each node and assign to mutations in table;
    mutation has distance of child node of edge
    """
    phylo_dist = dict()
    topology_dist = dict()
    phylo_dist[tree.name] = 0
    topology_dist[tree.name] = 0

    for node in tqdm.tqdm(tree.iter_descendants(), total=len(tree.get_descendants())):
        pd = tree.get_distance(node)
        td = int(tree.get_distance(node, topology_only=True))
       
        phylo_dist[node.name] = pd
        topology_dist[node.name] = td

    assert len(phylo_dist) == len(tree.get_descendants()) + 1
    assert len(topology_dist) == len(tree.get_descendants()) + 1

    substitutions = substitutions.copy()
    # substitutions['phylo_dist'] = substitutions.child_node.map(phylo_dist)
    substitutions['topology_dist'] = substitutions.child_node.map(topology_dist)
    return substitutions


def add_id_col(df: pd.DataFrame, col_name="mut_id"):
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
    dnd = read_tree(tree)
    df = add_dist2root(df, dnd)
    df = add_id_col(df)
    df.to_csv(out_table, index=None)


if __name__ == "__main__":
    main()
