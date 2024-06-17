import csv
import os

import numpy as np
from os import walk


class DistanceMaker:

    def __init__(self):
        self.fastme_path = 'Data_/benchmarks/data_sets/'
        self.rdpii_path = 'Data_/benchmarks/data_sets/RDPII_218.phy'
        self.zilla_path = 'Data_/benchmarks/data_sets/500_ZILLA.phy'
        self.fastme_dist_folder = 'Data_/benchmarks/fastme_distances/'
        self.matrices_folder = 'Data_/benchmarks/matrices/full_mats/'

        self.distances = {
            'p_distance': 'p',
            'RY_symmetric': 'Y',
            'RY': 'R',
            'JC69': 'J',
            'K2P': 'K',
            'F81': '1',
            'F84': '4',
            'TN93': 'T'

        }

    def compute_matrices(self):
        for d in self.distances.keys():
            os.system(self.fastme_path + "fastme_for_distance -i " + self.rdpii_path + ' --remove_gap'
                      + ' --dna=' + self.distances[d] + ' -O ' + self.fastme_dist_folder + 'rdpii_' + d)

            os.system(self.fastme_path + "fastme_for_distance -i " + self.zilla_path + ' --remove_gap'
                      + ' --dna=' + self.distances[d] + ' -O ' + self.fastme_dist_folder + 'zilla_' + d)

    def make_np_csv(self):

        # list to store files name
        files = []
        for (dir_path, dir_names, file_names) in walk(self.fastme_dist_folder):
            files.extend(file_names)

        for file in files:
            with open(self.fastme_dist_folder + file, newline='') as csvfile:
                spamreader = csv.reader(csvfile, quotechar='|')
                i = 0
                mat = []
                for row in spamreader:
                    if i > 0 and len(row) > 0:
                        string_row = row[0][15:].split('    ')
                        mat.append([float(num) for num in string_row])

                    i += 1

            mat = np.array(mat)
            np.savetxt(self.matrices_folder + file + '.csv', mat)


dist_maker = DistanceMaker()

dist_maker.compute_matrices()

dist_maker.make_np_csv()



