import os

import numpy as np

files = []
for (dir_path, dir_names, file_names) in os.walk('Data_/benchmarks/matrices/full_mats/'):
    files.extend(file_names)

files = sorted(files)
files = files[:4] + files[8:12]

for file in files:

    d = np.abs(np.around(np.loadtxt('Data_/benchmarks/matrices/full_mats/' + file), 10))

    if file[0] == 'r':
        d_ = d[:100, :100]
        np.savetxt('Data_/benchmarks/matrices/experiment_mats/100_' + file + '.csv', d_)
        d_ = d[:200, :200]
        np.savetxt('Data_/benchmarks/matrices/experiment_mats/200_' + file + '.csv', d_)
    else:
        d_ = d[:300, :300]
        np.savetxt('Data_/benchmarks/matrices/experiment_mats/200_' + file + '.csv', d_)