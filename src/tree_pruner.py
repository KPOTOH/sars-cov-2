"""
Чтобы удалить ноду в бифуркации (или метелке) и сохранить действительное расстояние
до оставшегося единственного (или нет) узла без промежуточной ноды, надо обратиться
по ссылке к этой ноде и вызвать метод `.delete(preserve_branch_length=True)`

"""

from random import randint
import sys

from ete3 import PhyloTree

DEFAULT_TREE_FORMAT = 1


def node_parent(node):
    return next(node.iter_ancestors())


def filter_parents(nodes: set) -> list:
    """ return only nodes that are really pre_leaves

    won't show all "dublicates" and because didn't use"""
    pre_leaves = []
    for node in nodes:
        if len(node) == len(node.get_children()):
            pre_leaves.append(node)
    return pre_leaves


def is_binary_tree(tree: PhyloTree):
    nchildrens_set = set()
    for node in tree.iter_descendants():
        nchildrens_set.add(len(node.get_children()))
    if nchildrens_set == {0, 2}:
        return True
    return False


def umbreallas2singles(tree: PhyloTree):
    """delete all umbreallas; automatically resolve polytomies
    (maybe not in other cases)
    """
    leaves = tree.get_leaves()
    leaves_parents = set(map(node_parent, leaves))

    for parent_node in leaves_parents:
        cur_leaves = []
        for potential_leaf in parent_node.get_children():
            if len(potential_leaf) == 1:
                cur_leaves.append(potential_leaf)
        if len(cur_leaves) == 1:
            continue

        accounted_dists = set()
        for leaf in cur_leaves:
            d = leaf.dist
            if d in accounted_dists:
                leaf.delete(preserve_branch_length=True)
            else:
                accounted_dists.add(d)
    return tree


def is_approx_equal(n1, n2):
    if n1 == n2:
        return True
    max_, min_ = max(n1, n2), min(n1, n2)
    cutoff = .25  # max_ >4 times more than min_ => is not equal
    if min_ / max_ > cutoff:
        return True
    return False


def random_deletion_in_bifurcations(tree: PhyloTree):
    """ delete random nodes in biburcations """
    if not is_binary_tree(tree):
        raise ValueError("tree must be binary")

    leaves = tree.get_leaves()
    leaves_parents = set(map(node_parent, leaves))
    pre_leaves = filter_parents(leaves_parents)

    for node in pre_leaves:
        childrens = node.get_children()
        assert len(childrens) == 2, "tree is not binary"

        if is_approx_equal(childrens[0].dist, childrens[1].dist):
            idx = randint(0, 1)
            childrens[idx].delete(preserve_branch_length=True)
    return tree


def pruning(newick, max_tree_len=55_000):
    tree = PhyloTree(newick, format=DEFAULT_TREE_FORMAT)
    umbreallas2singles(tree)
    i = 0
    while len(tree) > max_tree_len:
        i += 1
        random_deletion_in_bifurcations(tree)  # but it's single object
        if i > 10:
            break
    return tree, i


def main():
    if len(sys.argv) == 1:
        sys.stderr.write(
            "where arguments?!\nusage: ./tree_pruner.py path_to_newick [path_to_pruned_tree]\n")
        return

    path_to_newick = sys.argv[1]
    path_to_pruned_tree = sys.argv[2] if len(
        sys.argv) > 2 else path_to_newick + '.pruned'

    sys.stderr.write(
        f"Run command: ./tree_pruner.py {path_to_newick} {path_to_pruned_tree}\n"
    )
    # tree = PhyloTree(path_to_newick, format=DEFAULT_TREE_FORMAT)
    # pruned_tree = umbreallas2singles(tree)
    pruned_tree, n_iterations = pruning(path_to_newick)
    sys.stderr.write(
        f"{len(pruned_tree)} nodes in pruned tree after umbreallas deletion "
        f"and {n_iterations} iterations of random bifurcations deletion\n"
    )
    pruned_tree.write(outfile=path_to_pruned_tree, format=DEFAULT_TREE_FORMAT)
    sys.stderr.write("Job done.\n")


if __name__ == "__main__":
    main()
