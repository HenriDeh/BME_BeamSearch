import random
from typing import Union, List, Tuple

import numpy as np
import torch
from Solvers.PhyloES.PhyloEsUtils.utils import random_trees_generator_and_objs, adjust_matrices
from Solvers.BNNI_BSPR.fast_cpp import FastCpp
from Solvers.solver import Solver


class PhyloES(Solver):
    def __init__(self, d, population_size: Union[int, List[Tuple]] = 16, max_iterations=1000, replace=True,
                 max_non_improve_iter=None, min_tol=1e-16, labels: List[str] = None):

        """

        d: np.array(n,n), the distance matrix

        population_size: int or  list of tuples, A single value specifies the fixed

        tree population size. A list of tuples of the kind [(iter1, size1), (iter2, size2),...] indicates the size of the population up to
        the given iteration

        max_iterations: int, maximum number of iterations

        replace: bool, whether to apply individual replacement

        max_non_improve_iter: int, if not None it stops the algorithm after a maximum number of iterations without improvement is reached

        min_tol: float, minimum numerical tolerance to compare the fitness of the best and the worst individuals of the population. If the
        difference between the two individuals is lower than the min_tol the algorithm stops.

        """

        super().__init__(d, labels=labels)

        self.d_np = self.d.astype(np.double)
        self.d = torch.tensor(self.d, device=self.device)
        self.replace = replace

        self.counter = 0
        self.best_iteration = None
        self.better_solutions = []
        self.max_iterations = max_iterations
        self.min_tol = min_tol
        self.n_trees = None
        self.fast_time = 0
        self.fast_me_solver = FastCpp()
        self.nni_counter, self.spr_counter = [], []
        self.adaptive = True if type(population_size) == list else False
        if self.adaptive:
            self.population_size = sorted(population_size, key=lambda t: t[0])
        else:
            self.population_size = population_size
        self.max_non_improve_iter = max_non_improve_iter if max_non_improve_iter is not None else max_iterations * population_size
        self.stop_criterion = None

        self.obj_vals = None

        self.best_vals = []
        self.worse_vals = []
        self.iterations = None

    def solve(self):

        """
        run the algorithm
        """

        iteration = 0
        population_size = self.set_batch(iteration) if self.adaptive else self.population_size

        init_mats = self.initial_adj_mat(self.device, population_size)
        obj_vals, adj_mats = random_trees_generator_and_objs(3, self.d, init_mats, self.n_taxa, self.powers,
                                                             self.device)
        obj_vals, adj_mats = self.run_bnni_spr(adj_mats, population_size)
        best = torch.argmin(obj_vals)
        trajectories = self.tree_encoding(adj_mats)

        self.obj_val, self.solution = obj_vals[best], adj_mats[best]

        tj = torch.zeros((2 * population_size, self.n_taxa - 3), device=self.device, dtype=torch.long)
        objs = torch.ones(2 * population_size, device=self.device, dtype=torch.float64) * 1000

        tj[population_size: 2 * population_size] = trajectories
        objs[population_size: 2 * population_size] = obj_vals
        objs, sorted_idxs = torch.sort(objs)
        tj = tj[sorted_idxs]

        not_improved_counter = 0
        iteration = 1
        combs = 2
        while combs > 1 and iteration < self.max_iterations and objs[population_size - 1] - objs[0] > self.min_tol:

            population_size = self.set_batch(iteration) if self.adaptive else self.population_size
            objs, tj = objs[:2 * population_size], tj[:2 * population_size]

            init_mats = self.initial_adj_mat(self.device, population_size)
            combs, adj_mats, objs, tj = self.distribution_policy(init_mats, tj, objs, population_size)
            self.best_vals.append(objs[0].item())
            self.worse_vals.append(objs[population_size - 1].item())

            obj_vals, adj_mats = self.run_bnni_spr(adj_mats, population_size)
            best = torch.argmin(obj_vals)
            trajectories = self.tree_encoding(adj_mats)

            if obj_vals[best] < self.obj_val:
                self.obj_val, self.solution = obj_vals[best], adj_mats[best]
                not_improved_counter = 0
            else:
                not_improved_counter += 1

            tj[population_size: 2 * population_size] = trajectories
            objs[population_size: 2 * population_size] = obj_vals

            objs, sorted_idxs = torch.sort(objs)
            tj = tj[sorted_idxs]

            iteration += 1

        self.n_trees = iteration * population_size

        self.stop_criterion = 'convergence' if combs == 1 else \
            ('min_tol' if objs[population_size - 1] - objs[0] < self.min_tol else
             ('max_iterations' if iteration == self.max_iterations else 'local_plateau'))
        self.iterations = iteration
        self.T = self.get_tau(self.solution.to('cpu'))
        self.d = self.d.to('cpu')
        self.obj_val = self.compute_obj()
        self.solution = self.solution.to('cpu').numpy()

    def set_batch(self, iteration):
        for i in range(len(self.population_size) - 1):
            if self.population_size[i][0] <= iteration < self.population_size[i + 1][0]:
                return self.population_size[i][1]
        return self.population_size[-1][1]

    def tree_encoding(self, trees):
        sols = trees
        trajectories = self.decoding(sols)
        return trajectories

    def decoding(self, adj_mats):
        last_inserted_taxa = self.n_taxa - 1
        n_internals = self.m - self.n_taxa

        adj_mats = adjust_matrices(adj_mats, last_inserted_taxa, n_internals, self.n_taxa)
        adj_mats = adj_mats.unsqueeze(1).repeat(1, 2, 1, 1)
        reversed_idxs = torch.tensor([[i, i - 1] for i in range(1, adj_mats.shape[0] * 2, 2)],
                                     device=adj_mats.device).flatten()
        trajectories = torch.zeros((adj_mats.shape[0], self.n_taxa - 3), device=self.device, dtype=torch.int)

        last_inserted_taxa = self.n_taxa - 1
        for step in range(self.m - 1, self.n_taxa, -1):
            adj_mats[:, 1, :, last_inserted_taxa] = adj_mats[:, 1, last_inserted_taxa, :] = 0
            idxs = torch.nonzero(adj_mats[:, 1, step])
            idxs = torch.column_stack([idxs, idxs[:, 1][reversed_idxs]])

            adj_mats[idxs[:, 0], 1, idxs[:, 1], idxs[:, 2]] = adj_mats[idxs[:, 0], 1, idxs[:, 2], idxs[:, 1]] = 1
            adj_mats[:, 1, :, step] = adj_mats[:, 1, step, :] = 0

            k = (last_inserted_taxa - 1) * 2 - 1
            all_non_zeros = torch.nonzero(torch.triu(adj_mats[:, 1, :, :])).view(adj_mats.shape[0], k, 3)
            chosen = idxs[range(0, adj_mats.shape[0] * 2, 2)].repeat_interleave(k, dim=0).view(adj_mats.shape[0], k, 3)
            tj = torch.argwhere((all_non_zeros == chosen).prod(dim=-1))
            trajectories[:, last_inserted_taxa - 3] = tj[:, 1]

            last_inserted_taxa -= 1
        return trajectories

    @staticmethod
    def permute(adj_mats, step, idx):
        adj_mats[:, step, :] += adj_mats[:, idx, :]
        adj_mats[:, idx, :] = adj_mats[:, step, :] - adj_mats[:, idx, :]
        adj_mats[:, step, :] -= adj_mats[:, idx, :]

        adj_mats[:, :, step] += adj_mats[:, :, idx]
        adj_mats[:, :, idx] = adj_mats[:, :, step] - adj_mats[:, :, idx]
        adj_mats[:, :, step] -= adj_mats[:, :, idx]
        return adj_mats

    def distribution_policy(self, adj_mats, trajectories, obj_vals, population_size):
        if self.replace:
            equals = [obj_vals[i] == obj_vals[i + 1] for i in range(population_size - 1)]
            idx = population_size - 2
            while idx > 0 and equals[idx]:
                idx -= 1
            if 0 < idx < population_size - 2:
                # print("replacement")
                trajectories[idx + 1] = trajectories[0]
        tj = trajectories[:population_size]

        combs = 1
        for step in range(3, self.n_taxa):
            c = torch.unique(tj[:, step - 3], return_counts=True)
            combs *= c[1].shape[0]
            idxs = random.choices(tj[:, step - 3], k=population_size)
            idxs_list = torch.nonzero(torch.triu(adj_mats)).reshape(population_size, -1, 3)
            idxs_list = idxs_list[[range(population_size)], idxs, :]
            idxs_list = (idxs_list[:, :, 0], idxs_list[:, :, 1], idxs_list[:, :, 2])
            adj_mats = self.add_nodes(adj_mats, idxs_list, new_node_idx=step, n=self.n_taxa)

        return combs, adj_mats, obj_vals, trajectories

    def run_bnni_spr(self, adj_mats, population_size):
        mats = adj_mats.to('cpu').numpy().astype(dtype=np.int32)
        adjs, objs, nni_counts, spr_counts = \
            self.fast_me_solver.run_parallel(self.d_np, mats, self.n_taxa, self.m, population_size)
        self.nni_counter += [nni_counts]
        self.spr_counter += [spr_counts]

        objs_tensor = torch.tensor(objs, device=self.device)
        adjs_tensor = torch.tensor(adjs, device=self.device)
        self.fast_me_solver.free_result_memory()

        return objs_tensor, adjs_tensor
