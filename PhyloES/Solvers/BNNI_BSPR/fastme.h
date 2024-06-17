//
//  Copyright 2002-2007 Rick Desper, Olivier Gascuel
//  Copyright 2007-2014 Olivier Gascuel, Stephane Guindon, Vincent Lefort
//
//  This file is part of FastME.
//
//  FastME is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  FastME is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastME.  If not, see <http://www.gnu.org/licenses/>
//


#ifndef FASTME_H_
#define FASTME_H_

#include "graph.h"
#include "bNNI.h"
#include "SPR.h"
// #include "distance.h"
#include "initialiser.h"

void run(double **d, int **init_adj, int* solution_mat, int n_taxa, int m, double & obj_val, int & nni_count, int & spr_count);
	
tree *ImproveTree (tree *T0, double **D, double **A,
	int *nniCount, int *sprCount);

// void explainedVariance (double **D, tree *T, int n, int precision,
// 	int input_type);

#endif /*FASTME_H_*/

