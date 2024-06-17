import numpy as np
import pandas as pd
from ast import literal_eval
from ete3 import Tree

from Solvers.FastME.pharser_newik.newwik_handler import compute_newick
from Solvers.solver import Solver

df = pd.read_csv('../results/resultsrdpii_K2P.csv.csv')

def get_ete_tree(code, n_taxa):
    solver = Solver()
    solver.n_taxa = n_taxa
    solver.m = solver.n_taxa * 2 - 2

    adj_mat = solver.initial_adj_mat().astype(int)
    for step in range(3, solver.n_taxa):
        idx = np.array(np.nonzero(np.triu(adj_mat))).T[code[step - 3]]
        adj_mat = solver.add_node(adj_mat, idx, step, solver.n_taxa)

    T = solver.get_tau(adj_mat)
    return compute_newick(T)

def compute_robinson_distance(code1, code2, n_taxa):
    tree1 = Tree(get_ete_tree(code1, n_taxa))
    tree2 = Tree(get_ete_tree(code2, n_taxa))
    res = tree1.robinson_foulds(tree2, unrooted_trees=True)
    return res[0]

tree_fast_me = literal_eval(df[df.run == 1].fast_traj.iloc[0])
tree_rand_me = literal_eval(df[df.run == 1].rand_traj.iloc[0])
tree_phyloes_me = literal_eval(df[df.run == 1].phyloes_tj.iloc[0])


print(compute_robinson_distance(tree_fast_me, tree_phyloes_me, 100))
print(compute_robinson_distance(tree_fast_me, tree_rand_me, 100))
print(compute_robinson_distance(tree_rand_me, tree_phyloes_me, 100))






