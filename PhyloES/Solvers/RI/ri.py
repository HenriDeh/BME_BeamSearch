import os
from typing import List

import numpy as np

from Solvers.PhyloES.PhyloEsUtils.utils import random_trees_generator
from Solvers.RI.Random.random_solver import RandomSolver
from Solvers.BNNI_BSPR.fast_cpp import FastCpp
from Solvers.solver import Solver


class RI(Solver):
    def __init__(self, d, parallel=False, labels: List[str] = None):
        super().__init__(d, labels=labels)
        self.d_np = self.d.astype(np.double)
        self.random_solver = RandomSolver(d)
        self.parallel = parallel
        self.counter = 0
        self.best_iteration = None
        self.numProcs = os.cpu_count()
        self.better_solutions = []

        self.fast_me_solver = FastCpp()
        self.nni_counter, self.spr_counter = 0, 0

    def solve(self, iterations, ob_init_val=10, sol_init=None, count_solutions=False):
        best_val, best_sol = 10 ** 5, sol_init
        self.better_solutions.append(sol_init)
        tot_iter = 0
        left = iterations
        while left > 0:
            batch = min([self.numProcs, left])
            init_mats = self.initial_adj_mat(self.device, batch)
            left -= self.numProcs
            adj_mats = random_trees_generator(3, init_mats, self.n_taxa)
            mats = adj_mats.to('cpu').numpy().astype(dtype=np.int32)
            adjs, objs, nni_counts, spr_counts = \
                self.fast_me_solver.run_parallel(self.d_np, mats, self.n_taxa, self.m, batch)
            self.fast_me_solver.free_result_memory()
            tot_iter += objs.shape[0]
            self.nni_counter += nni_counts
            self.spr_counter += spr_counts

            best_val_idx = np.argmin(objs)
            if objs[best_val_idx] < best_val:
                best_val, best_sol = objs[best_val_idx], adjs[best_val_idx]
                # self.best_iteration = i
        self.solution = best_sol
        self.obj_val = best_val
        self.T = self.get_tau(self.solution)
