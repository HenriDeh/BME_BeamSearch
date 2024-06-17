import os
import random

import numpy as np
import pandas as pd
import torch

from Solvers.FastME.fast_me import FastMeSolver
from Solvers.PhyloES.phyloes import PhyloES

from Solvers.RI.ri import RI

torch.set_printoptions(precision=16)

torch.manual_seed(0)
np.random.seed(0)
random.seed(0)
seed = 0

max_iter = 10_000
pop_size = [(0, 64), (5, 32), (25, 16)]

files = []
for (dir_path, dir_names, file_names) in os.walk('Data_/benchmarks/matrices/experiment_mats'):
    files.extend(file_names)

files = sorted(files)
print(files, '\n')


for file in files:
    print("##########################")
    print(f'Working on datsaset: {file}')
    print(f'Algorithm params: \n> population size (n_iter, pop_size): {pop_size}\n> max iterations: {max_iter}\n')
    data_set_list = []
    run_list, pop_size_list, max_iter_list, stop_list, n_better_solution = [], [], [], [], []

    result_list, best_list, worse_list, iterations = [], [], [], []

    run_per_problem = 10

    d = np.abs(np.around(np.loadtxt('Data_/benchmarks/matrices/experiment_mats/' + file), 10))

    for run in range(run_per_problem):
        print('\n##> Run', run)
        result_run = [d.shape[0], file]
        run_list.append(run)
        pop_size_list.append(pop_size)
        max_iter_list.append(max_iter)

        random.seed(seed)
        phyloes = PhyloES(d, population_size=pop_size, max_iterations=max_iter, replace=True)
        phyloes.solve_timed()
        print("phyloes\ttime:", phyloes.time, '\tobj:', phyloes.obj_val, '   n_trees:',
              phyloes.n_trees, '  stop criterion:', phyloes.stop_criterion)
        stop_list.append(phyloes.stop_criterion)
        best_list.append(phyloes.best_vals)
        worse_list.append(phyloes.worse_vals)
        iterations.append(phyloes.iterations)

        rand_fast = RI(d, parallel=False)
        rand_fast.solve_timed(phyloes.n_trees)
        print("rand_fa\ttime:", rand_fast.time, '\tobj:', rand_fast.obj_val, '   n_trees:',
              phyloes.n_trees)

        fast = FastMeSolver(d, bme=True, nni=True, digits=17, bootrstap=True, post_processing=True, triangular_inequality=False,
                            logs=False)
        fast.solve_all_flags()
        print("fast\ttime:",  fast.time, '\tobj:', fast.obj_val)
        result_run += [fast.obj_val, rand_fast.obj_val,  phyloes.obj_val, fast.time, rand_fast.time,
                       phyloes.time, fast.method, rand_fast.nni_counter, rand_fast.spr_counter, phyloes.nni_counter, phyloes.spr_counter]
        result_list.append(result_run)

        seed += 1

    df = pd.DataFrame(result_list, columns=['Taxa', 'Problem', 'fast_obj', 'random_obj', 'phyloes_obj', 'fast_time',
                                                'rand_time', 'phyloes_time', 'fastMe_init', 'rf_nni', 'rf_spr', 'p_nni', 'p_spr'])
    df['run'] = run_list
    df['pop_size'] = pop_size_list
    df['max_iter'] = max_iter_list
    df['stop'] = stop_list
    df['iterations'] = iterations
    df['best_list'] = best_list
    df['worst_list'] = worse_list


# df.to_csv('results/results' + file + '.csv', index_label=False, index=False)




