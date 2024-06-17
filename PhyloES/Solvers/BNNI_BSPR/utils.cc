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


#include "utils.h"


/*********************************************************/

int *initZeroArray (int l)
{
	int *x = new int[l];
	int i;

	for (i=0; i<l; i++)
		x[i] = 0;

	return x;
}

/*********************************************************/

int *initOneArray (int l)
{
	int i;

	int *x = new int[l];

	for (i=0; i<l; i++)
		x[i] = 1;

	return x;
}

/*********************************************************/

double **initDoubleMatrix (int d)
{
	int i,j;
	double **A = new double*[d];


	for (i=0; i<d; i++)
	{
		A[i] = new double[d];
		for (j=0; j<d; j++)
			A[i][j] = 0.0;
	}

	return A;
}

/*********************************************************/

void fillZeroMatrix (double ***A, int d)
{
	int i,j;
	
	for (i=0; i<d; i++)
	{
		for (j=0; j<d; j++)
			(*A)[i][j] = 0.0;
	}
	
	return;
}


void deleteMatrix (double **D, int size)
{
	int i=0;

	if (nullptr != D)
	{
		for (i=0; i<size; i++)
			if (nullptr != D[i])
				delete[] D[i];
		delete[] D;
	}

	return;
}