import string

import numpy as np

from Solvers.PhyloES.phyloes import PhyloES


'''simple run from numpy array distance matrix'''

# distance matrix as numpy array
d = np.loadtxt('Data_/benchmarks/matrices/experiment_mats/100_rdpii_F81.csv')

# reducing d size to 30x30 for the example purposes
d = d[:30, :30]

# init phyloes
pes = PhyloES(d)

# run the algorithm
pes.solve()

# get solution tree length
print(pes.obj_val)

# get solution tree adj mat
print(pes.solution)

# get networkx graph
nx_graph = pes.get_nx_graph()


'''example for plotting solutions'''

d = d[:20, :20]

# optional taxa labels: the order of the labels must be consistent with the rows of the input matrix d
labels = [char for char in string.ascii_uppercase[:20]]

# init phyloes
pes = PhyloES(d, labels=labels)

# run the algorithm
pes.solve()

# plot phylogeny
pes.plot_phylogeny(filename=None, node_size=300)


'''all params example'''

pes = PhyloES(d, population_size=16, max_iterations=1000, replace=True,
              max_non_improve_iter=10, min_tol=1e-16, labels=labels)

# full params explanation
print(help(pes))

