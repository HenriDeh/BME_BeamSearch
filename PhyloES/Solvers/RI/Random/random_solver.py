import random

import numpy as np
from Solvers.solver import Solver


class RandomSolver(Solver):

    def __init__(self, d, sorted_d=False):
        super(RandomSolver, self).__init__(d, sorted_d)

    def solve(self, start=3, adj_mat=None):
        adj_mat = self.initial_adj_mat() if adj_mat is None else adj_mat
        for i in range(start, self.n_taxa):
            idxs_list = np.array(np.nonzero(np.triu(adj_mat))).T
            idxs = random.choice(idxs_list)
            adj_mat = self.add_node(adj_mat, idxs, i, self.n_taxa)

        self.solution = adj_mat
        self.obj_val = self.compute_obj_val_from_adj_mat(adj_mat, self.d, self.n_taxa)
        self.T = self.get_tau(self.solution)


