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


#include "SPR.h"
#include <iostream>

/*********************************************************/

void zero3DMatrix (double ***X, int h, int l, int w)
{
	int i, j, k;

	for (i=0; i<h; i++)
		for (j=0; j<l; j++)
			for (k=0; k<w; k++)
				X[i][j][k] = 0.0;

	return;
}

/*********************************************************/

void findTableMin (int *imin, int *jmin, int *kmin, int n, double ***X, double *min)
{
	int i, j, k;

	for (i=0; i<2; i++)
		for (j=0; j<n; j++)
			for (k=0; k<n; k++)
			{
				if (X[i][j][k] < *min)
				{
					*min = X[i][j][k];
					*imin = i;
					*jmin = j;
					*kmin = k;
				}
			}

	return;
}

/*********************************************************/

void SPR(tree* T, double** D, double **A, int *count)
{
	int i, j, k, countMinWeight;
	node *v;
	edge *e, *f;
	double ***swapWeights;
	double treeWeightBefore, treeWeightMin;
	double swapValue = 0.0;

	double eps = 0.0000000000000000001; //FLT_EPSILON
	
	// std::cout<<"Assign BME weights and re-weight tree"<<std::endl;
	
	swapWeights = new double** [2];
	makeBMEAveragesTable (T, D, A);
	assignBMEWeights (T, A);
	weighTree (T);
	treeWeightBefore = treeWeightMin = T->weight;
	countMinWeight = 0;

	// std::cout<<"Before SPR: tree length is "<< treeWeightBefore<<std::endl;


	for (i=0; i<2; i++)
		swapWeights[i] = initDoubleMatrix (T->size);

	do
	{
		swapValue = 0.0;
		zero3DMatrix (swapWeights, 2, T->size, T->size);
		i = j = k = 0;

		for (e = depthFirstTraverse (T, nullptr); nullptr != e; e = depthFirstTraverse (T, e))
			assignSPRWeights (e->head, A, swapWeights);

		findTableMin (&i, &j, &k, T->size, swapWeights, &swapValue);
		swapValue = swapWeights[i][j][k];


		if (swapValue < -eps)
		{
			// std::cout<<"New tree length should be "<< T->weight + 0.25 * swapValue<<std::endl;

			v = indexedNode (T, j);
			f = indexedEdge (T, k);


			
				// if ((nullptr == f->head->index) || (strlen (f->head->label) == 0))
				// {
				// 	if ((nullptr == f->tail->label) || (strlen (f->tail->label) == 0))
				// 		Debug ( (char*)"Swapping tree below '%s' to split edge '%s' with internal head and tail.",
				// 			v->parentEdge->label, f->label);

				// 	else
				// 		Debug ( (char*)"Swapping tree below '%s' to split edge '%s' with internal head and tail '%s'.",
				// 			v->parentEdge->label, f->label, f->tail->label);
				// }
				// else
				// {
				// 	if ((nullptr == f->tail->label) || (strlen (f->tail->label) == 0))
				// 		Debug ( (char*)"Swapping tree below '%s' to split edge '%s' with head '%s' and internal tail.",
				// 			v->parentEdge->label, f->label, f->head->label, f->tail->label);

				// 	else
				// 		Debug ( (char*)"Swapping tree below '%s' to split edge '%s' with head '%s' and tail '%s'.",
				// 			v->parentEdge->label, f->label, f->head->label, f->tail->label);
				// }
			

			SPRTopShift (v, f, 2-i);
			makeBMEAveragesTable (T, D, A);
			assignBMEWeights (T, A);
			weighTree (T);
			(*count)++;
			
			
			// std::cout<<"SPR "<<*count<<": new tree length is "<<T->weight<<std::endl;
			// std::cout<<"min roba "<<eps<<std::endl;

			// If SPR lead to the same minimum tree weight twice, then stop
			if (T->weight < treeWeightMin - eps)
			{
				treeWeightMin = T->weight;
				countMinWeight = 0;
			}
			else
			{
				// New tree weight equals minimum tree weight
				if (T->weight < treeWeightMin + eps)
					countMinWeight++;
			}

		}
	} while (swapValue < -eps && countMinWeight < 1);


	for (i=0; i<2; i++)
		deleteMatrix (swapWeights[i], T->size);

	delete[] swapWeights;

	return;
}

/*********************************************************/

/* Assigns values to array swapWeights
 * swapWeights[0][j][k] will be the value of removing the tree below the
 * edge whose head node has index j and reattaching it to split the edge
 * whose head has the index k
 * swapWeights[1][j][k] will be the value of removing the tree above the
 * edge whose head node has index j and reattaching it to split the edge
 * whose head has the index k */
void assignSPRWeights (node *vtest, double **A, double ***swapWeights)
{
	edge *etest, *left, *right, *sib, *par;

	etest = vtest->parentEdge;
	left = vtest->leftEdge;
	right = vtest->rightEdge;
	par = etest->tail->parentEdge;
	sib = siblingEdge (etest);

	if (nullptr != par)
		assignDownWeightsUp (par, vtest, sib->head, nullptr, nullptr, 0.0, 1.0, A, swapWeights);

	if (nullptr != sib)
		assignDownWeightsSkew (sib, vtest, sib->tail, nullptr, nullptr, 0.0, 1.0, A, swapWeights);

	/* Assigns values for moving subtree rooted at vtest, starting with
	 * edge parental to tail of edge parental to vtest */

	if (nullptr != left)
	{
		assignUpWeights (left, vtest, right->head, nullptr, nullptr, 0.0, 1.0, A, swapWeights);
		assignUpWeights (right, vtest, left->head, nullptr, nullptr, 0.0, 1.0, A, swapWeights);
	}

	return;
}

/*********************************************************/

/* Recall NNI formula: change in tree length from AB|CD split to AC|BD
 * split is proportional to D_AC + D_BD - D_AB - D_CD
 * In our case B is the tree being moved (below vtest), A is the tree
 * backwards below back, but with the vtest subtree removed, C is the
 * sibling tree of back and D is the tree above etest
 * Use va to denote the root of the sibling tree to B in the original
 * tree.
 * Please excuse the multiple uses of the same letters: A,D, etc. */
void assignDownWeightsUp (edge *etest, node *vtest, node *va, edge *back,
	node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights)
{
	edge *par, *sib, *skew;
	double D_AC, D_BD, D_AB, D_CD;
	par = etest->tail->parentEdge;
	skew = siblingEdge(etest);
	if (nullptr == back)	// first recursive call
	{
		if (nullptr == par)
			return;

		else	// start the process of assigning weights recursively
		{
			assignDownWeightsUp (par, vtest, va, etest, va,
						A[va->index][vtest->index], 0.5, A, swapWeights);
			assignDownWeightsSkew (skew, vtest, va, etest, va,
						A[va->index][vtest->index], 0.5, A, swapWeights);
		}
	}
	else	// second or later recursive call
	{
		sib = siblingEdge (back);
		D_BD = A[vtest->index][etest->head->index];		// straightforward
		D_CD = A[sib->head->index][etest->head->index];	// this one too
		D_AC = A[sib->head->index][back->head->index]
			+ coeff*(A[sib->head->index][va->index]
			- A[sib->head->index][vtest->index]);
		D_AB = 0.5 * (oldD_AB + A[vtest->index][cprev->index]);
		swapWeights[0][vtest->index][etest->head->index] = swapWeights[0][vtest->index][back->head->index]
			+ (D_AC + D_BD - D_AB - D_CD);

		if (nullptr != par)
		{
			assignDownWeightsUp (par, vtest, va, etest, sib->head,
								D_AB, 0.5 * coeff, A, swapWeights);
			assignDownWeightsSkew (skew, vtest, va, etest, sib->head,
								D_AB, 0.5 * coeff, A, swapWeights);
		}
	}

	return;
}

/*********************************************************/

void assignDownWeightsSkew (edge *etest, node *vtest, node *va, edge *back,
	node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights)
{
	/* same general idea as assignDownWeights,
	 * except needing to keep track of things a bit differently */

	edge *par, *left, *right;
	/* par here = sib before
	 * left, right here = par, skew before */

	double D_AB, D_CD, D_AC, D_BD;
	/* B is subtree being moved - below vtest
	 * A is subtree remaining fixed - below va, unioned with all trees already passed by B
	 * C is subtree being passed by B, in this case above par
	 * D is subtree below etest, fixed on other side */

	par = etest->tail->parentEdge;
	left = etest->head->leftEdge;
	right = etest->head->rightEdge;
	if (nullptr == back)
	{
		if (nullptr == left)
			return;

		else
		{
			assignDownWeightsDown (left, vtest, va, etest, etest->tail,
					A[vtest->index][etest->tail->index], 0.5, A, swapWeights);
			assignDownWeightsDown (right, vtest, va, etest, etest->tail,
					A[vtest->index][etest->tail->index], 0.5, A, swapWeights);
		}
	}
	else
	{
		D_BD = A[vtest->index][etest->head->index];
		D_CD = A[par->head->index][etest->head->index];
		D_AC = A[back->head->index][par->head->index] +
			coeff * (A[va->index][par->head->index] -
			A[vtest->index][par->head->index]);
		D_AB = 0.5 * (oldD_AB + A[vtest->index][cprev->index]);
		swapWeights[0][vtest->index][etest->head->index] =
			swapWeights[0][vtest->index][back->head->index] +
			(D_AC + D_BD - D_AB - D_CD);

		if (nullptr != left)
		{
			assignDownWeightsDown (left, vtest, va, etest, etest->tail,
				D_AB, 0.5 * coeff, A, swapWeights);
			assignDownWeightsDown (right, vtest, va, etest, etest->tail,
				D_AB, 0.5 * coeff, A, swapWeights);
		}
	}

	return;
}

/*********************************************************/

void assignDownWeightsDown (edge *etest, node *vtest, node *va, edge *back,
	node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights)
{
	/* again the same general idea */

	edge *sib, *left, *right;
	/* sib here = par in assignDownWeightsSkew
	 * rest is the same as assignDownWeightsSkew */

	double D_AB, D_CD, D_AC, D_BD;
	/* B is below vtest, A is below va unioned with all trees already passed by B
	 * C is subtree being passed - below sib
	 * D is tree below etest */

	sib = siblingEdge (etest);
	left = etest->head->leftEdge;
	right = etest->head->rightEdge;
	D_BD = A[vtest->index][etest->head->index];
	D_CD = A[sib->head->index][etest->head->index];
	D_AC = A[sib->head->index][back->head->index] +
		coeff * (A[sib->head->index][va->index] -
		A[sib->head->index][vtest->index]);
	D_AB = 0.5 * (oldD_AB + A[vtest->index][cprev->index]);
	swapWeights[0][vtest->index][etest->head->index] =
		swapWeights[0][vtest->index][back->head->index] +
		( D_AC + D_BD - D_AB - D_CD);

	if (nullptr != left)
	{
		assignDownWeightsDown (left, vtest, va, etest, sib->head, D_AB,
			0.5 * coeff, A, swapWeights);
		assignDownWeightsDown (right, vtest, va, etest, sib->head, D_AB,
			0.5 * coeff, A, swapWeights);
	}

	return;
}

/*********************************************************/

void assignUpWeights (edge *etest, node *vtest, node *va, edge *back,
	node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights)
{
	/* SPR performed on tree above vtest...
	 * same idea as above, with appropriate selections of edges and nodes */

	edge *sib, *left, *right;
	double D_AB, D_CD, D_AC, D_BD;
	/* B is above vtest, A is other tree below vtest unioned with trees in path to vtest
	 * sib is tree C being passed by B
	 * D is tree below etest */

	sib = siblingEdge (etest);
	left = etest->head->leftEdge;
	right = etest->head->rightEdge;
	if (nullptr == back)	/* first recursive call */
	{
		if (nullptr == left)
			return;

		else	/* start the process of assigning weights recursively */
		{
			assignUpWeights (left, vtest, va, etest, va,
				A[va->index][vtest->index], 0.5, A, swapWeights);
			assignUpWeights (right, vtest, va, etest, va,
				A[va->index][vtest->index], 0.5, A, swapWeights);
		}
	}
	else	/* second or later recursive call */
	{
		D_BD = A[vtest->index][etest->head->index];
		D_CD = A[sib->head->index][etest->head->index];
		D_AC = A[back->head->index][sib->head->index] +
			coeff * (A[va->index][sib->head->index] -
			A[vtest->index][sib->head->index]);
		D_AB = 0.5 * (oldD_AB + A[vtest->index][cprev->index]);
		swapWeights[1][vtest->index][etest->head->index] =
			swapWeights[1][vtest->index][back->head->index] +
			(D_AC + D_BD - D_AB - D_CD);

		if (nullptr != left)
		{
			assignUpWeights (left, vtest, va, etest, sib->head, D_AB,
				0.5 * coeff, A, swapWeights);
			assignUpWeights (right, vtest, va, etest, sib->head, D_AB,
				0.5 * coeff, A, swapWeights);
		}
	}

	return;
}

/*********************************************************/

/* Starting with edge u above edges p, d
 * removes p, d from tree, u connects to d->head to compensate */
void pruneSubtree (edge *p, edge *u, edge *d)
{
	p->tail->parentEdge = nullptr;	/* remove p subtree */
	u->head = d->head;
	d->head->parentEdge = u;	/* u connected to d->head */
	d->head = nullptr;				/* d removed from tree */

	return;
}

/*********************************************************/

/* Splits edge e to make it parental to p,d.
 * d is parental to what previously was below e */
void SPRsplitEdge (edge *e, edge *p, edge *d)
{
	d->head = e->head;
	e->head = p->tail;
	p->tail->parentEdge = e;
	d->head->parentEdge = d;

	return;
}

/*********************************************************/

/* Topological shift function
 * removes subtree rooted at v and re-inserts to split e */
void SPRDownShift (node *v, edge *e)
{
	edge *vup, *vdown, *vpar;

	vpar = v->parentEdge;
	vdown = siblingEdge(vpar);
	vup = vpar->tail->parentEdge;

	/* topological shift */
	pruneSubtree (vpar, vup, vdown);

	/* removes v subtree and vdown, extends vup */
	SPRsplitEdge (e, vpar, vdown);		/* splits e to make e sibling
										 * edge to vpar, both below vup */

	return;
}

/*********************************************************/

/* An inelegant iterative version */
void SPRUpShift (node *vmove, edge *esplit)
{
	edge *f;
	edge **EPath;
	edge **sib;
	node **v;
	int i;
	int pathLength;

	for (f=esplit->tail->parentEdge,pathLength=1; f->tail != vmove; f=f->tail->parentEdge)
		pathLength++;
	/* count number of edges to vmove
	 * note that pathLength > 0 will hold */

	EPath = new edge*[pathLength];
	v = new node* [pathLength];
	sib = new edge* [pathLength + 1];
	/* there are pathLength + 1 side trees, one at the head and tail of each edge in the path */

	sib[pathLength] = siblingEdge (esplit);
	i = pathLength;
	f = esplit->tail->parentEdge;
	while (i > 0)
	{
		i--;
		EPath[i] = f;
		sib[i] = siblingEdge (f);
		v[i] = f->head;
		f = f->tail->parentEdge;
	}
	/* indexed so head of Epath is v[i], tail is v[i-1] and sibling edge is sib[i]
	 * need to assign head, tail of each edge in path
	 * as well as have proper values for the left and right fields */

	if (esplit == esplit->tail->leftEdge)
	{
		vmove->leftEdge = esplit;
		vmove->rightEdge = EPath[pathLength-1];
	}
	else
	{
		vmove->rightEdge = esplit;
		vmove->leftEdge = EPath[pathLength-1];
	}
	esplit->tail = vmove;
	/* espilt->head remains unchanged
	 * vmove has the proper fields for left, right, and parentEdge */

	for(i=0; i<(pathLength-1); i++)
		EPath[i]->tail = v[i+1];

	/* this bit flips the orientation along the path
	 * the tail of Epath[i] is now v[i+1] instead of v[i-1] */

	EPath[pathLength-1]->tail = vmove;

	for (i=1; i<pathLength; i++)
	{
		if (sib[i+1] == v[i]->leftEdge)
			v[i]->rightEdge = EPath[i-1];

		else
			v[i]->leftEdge = EPath[i-1];
	}

	if (sib[1] == v[0]->leftEdge)
		v[0]->rightEdge = sib[0];

	else
		v[0]->leftEdge = sib[0];

	sib[0]->tail = v[0];

	delete[] EPath;
	delete[] v;
	delete[] sib;

	return;
}

/*********************************************************/

void SPRTopShift (node *vmove, edge *esplit, int UpOrDown)
{
	if (2 == UpOrDown) // DOWN = 2
		SPRDownShift (vmove, esplit);

	else
		SPRUpShift (vmove, esplit);

	return;
}
