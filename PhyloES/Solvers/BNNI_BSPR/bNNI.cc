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


#include "bNNI.h"
#include <iostream>


/*********************************************************/

void bNNIRetestEdge (int *p, int *q, edge *e,tree *T, double **avgDistArray,
	double *weights, int *location, int *possibleSwaps)
{
	int tloc;

	tloc = location[e->head->index+1];
	location[e->head->index+1] = bNNIEdgeTest (e, T, avgDistArray, weights + e->head->index+1);

	if (0 == location[e->head->index+1])  // NONE = 0
	{
		if (0 != tloc)
			popHeap (p, q, weights, (*possibleSwaps)--, q[e->head->index+1]);
	}
	else
	{
		if (0 == tloc)
			pushHeap (p, q, weights, (*possibleSwaps)++, q[e->head->index+1]);

		else
			reHeapElement (p, q, weights, *possibleSwaps, q[e->head->index+1]);
	}

	return;
}

/*********************************************************/

void bNNItopSwitch (edge *e, int direction, double **A)
{
	edge *down, *swap, *fixed;
	node *u, *v;

	down = siblingEdge (e);
	u = e->tail;
	v = e->head;

	if (3 == direction) // LEFT = 3
	{
		swap = e->head->leftEdge;
		fixed = e->head->rightEdge;
		v->leftEdge = down;
	}
	else
	{
		swap = e->head->rightEdge;
		fixed = e->head->leftEdge;
		v->rightEdge = down;
	}

	swap->tail = u;
	down->tail = v;

	if (e->tail->leftEdge == e)
		u->rightEdge = swap;

	else
		u->leftEdge = swap;

	bNNIupdateAverages (A, v, e->tail->parentEdge, down, swap, fixed);

	return;
}

/*********************************************************/

void bNNI (tree *T, double **avgDistArray, int *count)
{
	edge *e;
	edge **edgeArray;
	int *p, *location, *q;
	int i;
	int possibleSwaps;
	double *weights;
	
	
	p = initPerm (T->size+1);
	q = initPerm (T->size+1);
	edgeArray = new edge*[T->size + 1];
	weights = new double [T->size + 1];
	location = new int [T->size + 1];
	for (i=0; i<T->size+1; i++)
	{
		weights[i] = 0.0;
		location[i] = 0; // 0 = 0
	}

	assignBMEWeights (T, avgDistArray);
	weighTree (T);

	// std::cout<<"Before NNI: tree length is "<< T->weight<<std::endl;


	e = findBottomLeft (T->root->leftEdge);
	while (nullptr != e)
	{
		edgeArray[e->head->index+1] = e;
		location[e->head->index+1] = bNNIEdgeTest (e, T, avgDistArray,
			weights + e->head->index + 1);
		e = depthFirstTraverse (T,e);
	}

	possibleSwaps = makeThreshHeap (p, q, weights, T->size+1,0.0);
	permInverse (p, q, T->size+1);

	/* We put the negative values of weights into a heap, indexed by p
	 * with the minimum value pointed to by p[1]
	 * p[i] is index (in edgeArray) of edge with i-th position in the
	 * heap, q[j] is the position of edge j in the heap */
	
	while (weights[p[1]] < -DBL_EPSILON)
	{
		(*count)++;
		T->weight = T->weight + weights[p[1]];
		
		// std::cout<<"NNI: "<< *count<<" new tree length is "<<T->weight<<std::endl;
		

		bNNItopSwitch (edgeArray[p[1]], location[p[1]], avgDistArray);
		location[p[1]] = 0;
		weights[p[1]] = 0.0;	//after the bNNI, this edge is in optimal configuration
		popHeap (p, q, weights, possibleSwaps--, 1);

		/* but we must retest the other edges of T */
		/* CHANGE 2/28/2003 expanding retesting to _all_ edges of T */

		e = depthFirstTraverse (T, nullptr);
		while (nullptr != e)
		{
			bNNIRetestEdge (p, q, e, T, avgDistArray, weights, location, &possibleSwaps);
			e = depthFirstTraverse (T, e);
		}
	}
	// std::cout<<std::endl;

	delete[] p;
	delete[] q;
	delete[] location;
	delete[] edgeArray;
	delete[] weights;
	assignBMEWeights (T, avgDistArray);

	return;
}

/* This function is the meat of the average distance matrix recalculation.
 * Idea is: we are looking at the subtree rooted at rootEdge. The subtree
 * rooted at closer is closer to rootEdge after the NNI, while the subtree
 * rooted at further is further to rootEdge after the NNI. direction tells
 * the direction of the NNI with respect to rootEdge */

/*********************************************************/

void updateSubTreeAfterNNI (double **A, node *v, edge *rootEdge,
	node *closer, node *further, double dcoeff, int direction)
{
	edge *sib;

	switch (direction)
	{
		case 1 :	/* UP = 1 rootEdge is below the center edge of the NNI
					 * recursive calls to subtrees, if necessary */
			if (nullptr != rootEdge->head->leftEdge)
				updateSubTreeAfterNNI (A, v, rootEdge->head->leftEdge,
					closer, further, 0.5 * dcoeff, 1);

			if (nullptr != rootEdge->head->rightEdge)
				updateSubTreeAfterNNI (A, v, rootEdge->head->rightEdge,
					closer, further, 0.5 * dcoeff, 1);

			updatePair (A, rootEdge, rootEdge, closer, further, dcoeff, 1);
			sib = siblingEdge (v->parentEdge);
			A[rootEdge->head->index][v->index] =
			A[v->index][rootEdge->head->index] =
				0.5 * A[rootEdge->head->index][sib->head->index] +
				0.5 * A[rootEdge->head->index][v->parentEdge->tail->index];
			break;

		case 2 :	/* DOWN =2 rootEdge is above the center edge of the NNI */
			sib = siblingEdge (rootEdge);
			if (nullptr != sib)
				updateSubTreeAfterNNI (A, v, sib, closer, further,
					0.5 * dcoeff, 5);

			if (nullptr != rootEdge->tail->parentEdge)
				updateSubTreeAfterNNI (A, v, rootEdge->tail->parentEdge,
					closer, further, 0.5 * dcoeff, 2);

			updatePair (A, rootEdge, rootEdge, closer, further, dcoeff, 2);
			A[rootEdge->head->index][v->index] =
			A[v->index][rootEdge->head->index] =
				0.5 * A[rootEdge->head->index][v->leftEdge->head->index] +
				0.5 * A[rootEdge->head->index][v->rightEdge->head->index];
			break;

		case 5 :	/* rootEdge is in subtree skew to v */
			if (nullptr != rootEdge->head->leftEdge)
				updateSubTreeAfterNNI (A, v, rootEdge->head->leftEdge,
					closer, further, 0.5 * dcoeff, 5);

			if (nullptr != rootEdge->head->rightEdge)
				updateSubTreeAfterNNI (A, v, rootEdge->head->rightEdge,
					closer, further, 0.5 * dcoeff, 5);

			updatePair (A, rootEdge, rootEdge, closer, further, dcoeff, 1);
			A[rootEdge->head->index][v->index] =
			A[v->index][rootEdge->head->index] =
				0.5 * A[rootEdge->head->index][v->leftEdge->head->index] +
				0.5 * A[rootEdge->head->index][v->rightEdge->head->index];
			break;
	}

	return;
}

/* swapping across edge whose head is v */

/*********************************************************/

void bNNIupdateAverages (double **A, node *v, edge *par, edge *skew,
	edge *swap, edge *fixed)
{
	A[v->index][v->index] = 0.25 *
		(A[fixed->head->index][par->head->index] +
		A[fixed->head->index][swap->head->index] +
		A[skew->head->index][par->head->index] +
		A[skew->head->index][swap->head->index]);

	updateSubTreeAfterNNI (A, v, fixed, skew->head, swap->head, 0.25, 1);
	updateSubTreeAfterNNI (A, v, par, swap->head, skew->head, 0.25, 2);
	updateSubTreeAfterNNI (A, v, skew, fixed->head, par->head, 0.25, 1);
	updateSubTreeAfterNNI (A, v, swap, par->head, fixed->head, 0.25, 5);

	return;
}

/*********************************************************/

double wf5 (double D_AD, double D_BC, double D_AC, double D_BD,
	double D_AB, double D_CD)
{
	double weight;

	weight = 0.25 * (D_AC + D_BD + D_AD + D_BC) + 0.5 * (D_AB + D_CD);

	return (weight);
}

/*********************************************************/

int bNNIEdgeTest (edge *e, tree *T, double **A, double *weight)
{
	edge *f;
	double D_LR, D_LU, D_LD, D_RD, D_RU, D_DU;
	double w1, w2, w0;

	if ((leaf(e->tail)) || (leaf(e->head)))
		return (0);

	f = siblingEdge (e);

	D_LR = A[e->head->leftEdge->head->index][e->head->rightEdge->head->index];
	D_LU = A[e->head->leftEdge->head->index][e->tail->index];
	D_LD = A[e->head->leftEdge->head->index][f->head->index];
	D_RU = A[e->head->rightEdge->head->index][e->tail->index];
	D_RD = A[e->head->rightEdge->head->index][f->head->index];
	D_DU = A[e->tail->index][f->head->index];

	w0 = wf5 (D_RU, D_LD, D_LU, D_RD, D_DU, D_LR);	// weight of current config
	w1 = wf5 (D_RU, D_LD, D_DU, D_LR, D_LU, D_RD);	// weight with L<->D switch
	w2 = wf5 (D_DU, D_LR, D_LU, D_RD, D_RU, D_LD);	// weight with R<->D switch

	if (w0 <= w1)
	{
		if (w0 <= w2)	// w0 <= w1,w2
		{
			*weight = 0.0;
			return (0);
		}
		else			// w2 < w0 <= w1
		{
			*weight = w2 - w0;
			// std::cout<< "Possible swap across "<<e->head-> index<< " Weight dropping by"<< e->tail->index<< " "<< w0 - w2<< std::endl;
			// std::cout<<"New tree length should be "<< T->weight + w2 - w0<<std::endl;
			
			return 4; // right
		}
	}
	else if (w2 <= w1)	// w2 <= w1 < w0
	{
		*weight = w2 - w0;
		
		// std::cout<< "Possible swap across "<<e->head-> index<< " Weight dropping by"<< e->tail->index<< " "<< w0 - w2<< std::endl;
		// std::cout<<"New tree length should be "<< T->weight + w2 - w0<<std::endl;
		return 4;
	}
	else				// w1 < w2, w0
	{
		*weight = w1 - w0;
		// std::cout<< "Possible swap across "<<e->head-> index<< " Weight dropping by"<< e->tail->index<< " "<< w0 - w1<< std::endl;
		// std::cout<<"New tree length should be "<< T->weight + w1 - w0<<std::endl;
		return 3; // LEFT = 3
	}
}

/* limitedFillTableUp fills all the entries in D associated with
 * e->head, f->head and those edges g->head above e->head, working
 * recursively and stopping when trigger is reached */

/*********************************************************/

void limitedFillTableUp (edge *e, edge *f, double **A, edge *trigger)
{
	edge *g, *h;
	g = f->tail->parentEdge;
	if (f != trigger)
		limitedFillTableUp (e, g, A, trigger);

	h = siblingEdge (f);
	A[e->head->index][f->head->index] =
	A[f->head->index][e->head->index] =
	0.5 * (A[e->head->index][g->head->index] + A[e->head->index][h->head->index]);

	return;
}

