import networkx as nx
import numpy as np
import torch
import random

from matplotlib import pyplot as plt

from Solvers.solver import Solver


def random_trees_generator(start, adj_mats, n_taxa):
    batch_size = adj_mats.shape[0]
    for step in range(start, n_taxa):
        choices = 3 + (step - 3) * 2
        idxs_list = torch.nonzero(torch.triu(adj_mats)).reshape(batch_size, -1, 3)
        rand_idxs = random.choices(range(choices), k=batch_size)
        idxs_list = idxs_list[[range(batch_size)], rand_idxs, :]
        idxs_list = (idxs_list[:, :, 0], idxs_list[:, :, 1], idxs_list[:, :, 2])
        adj_mats = Solver.add_nodes(adj_mats, idxs_list, new_node_idx=step, n=n_taxa)
    return adj_mats


def random_trees_generator_and_objs(start, d, adj_mats, n_taxa, powers, device):
    adj_mats = random_trees_generator(start, adj_mats, n_taxa)
    obj_vals = Solver.compute_obj_val_batch(adj_mats, d, powers, n_taxa, device)
    return obj_vals, adj_mats


def adjust_matrices(adj_mat, last_inserted_taxa, n_internals, n_taxa):
    adj_mat = adj_mat.unsqueeze(1).repeat(1, 2, 1, 1)

    # reorder matrix according to encoding
    for step in range(last_inserted_taxa + n_internals, n_taxa, -1):
        not_in_position = torch.argwhere(adj_mat[:, 1, step, last_inserted_taxa] == 0).squeeze(1)

        if len(not_in_position) > 0:
            idxs = torch.argwhere(adj_mat[not_in_position, 1, last_inserted_taxa, :] == 1)

            adj_mat = permute(adj_mat, step, not_in_position, idxs)

        adj_mat[:, 1, :, last_inserted_taxa] = adj_mat[:, 1, last_inserted_taxa, :] = 0

        i = torch.nonzero(adj_mat[:, 1, step, :])
        idx = i[:, 1]

        idxs = idx.reshape(-1, 2).T
        batch_idxs = range(adj_mat.shape[0])
        adj_mat[batch_idxs, 1, idxs[0], idxs[1]] = adj_mat[batch_idxs, 1, idxs[1], idxs[0]] = 1
        adj_mat[:, 1, :, step] = adj_mat[:, 1, step, :] = 0

        last_inserted_taxa -= 1
    adj_mat = adj_mat[:, 0, :, :]
    return adj_mat


def permute(adj_mats, step, not_in_position, idx):
    adj_mats[not_in_position, :, step, :] += adj_mats[not_in_position, :, idx[:, 1], :]
    adj_mats[not_in_position, :, idx[:, 1], :] = adj_mats[not_in_position, :, step, :] - adj_mats[not_in_position, :, idx[:, 1], :]
    adj_mats[not_in_position, :, step, :] -= adj_mats[not_in_position, :, idx[:, 1], :]

    adj_mats[not_in_position, :, :, step] += adj_mats[not_in_position, :, :, idx[:, 1]]
    adj_mats[not_in_position, :, :, idx[:, 1]] = adj_mats[not_in_position, :, :, step] - adj_mats[not_in_position, :, :, idx[:, 1]]
    adj_mats[not_in_position, :, :, step] -= adj_mats[not_in_position, :, :, idx[:, 1]]
    return adj_mats


