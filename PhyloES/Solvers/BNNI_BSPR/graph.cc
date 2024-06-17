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


#include "graph.h"
#include <iostream>

/*********************************************************/

bool leaf (node *v)
{
	int count = 0;

	if (nullptr != v->parentEdge)
		count++;

	if (nullptr != v->leftEdge)
		count++;

	if (nullptr != v->rightEdge)
		count++;

	if (count > 1)
		return (false);

	return (true);
}

/*********************************************************/

set *addToSet (node *v, set *X)
{
	if (nullptr == X)
	{
		X = new set;
		X->firstNode = v;
		X->secondNode = nullptr;
	}
	else if (nullptr == X->firstNode)
		X->firstNode = v;

	else
		X->secondNode = addToSet (v, X->secondNode);

	return (X);
}

/*********************************************************/

set *copySet (set *X)
{
	set *ret = new set;

	if (nullptr != X)
	{
		ret->firstNode = copyNode (X->firstNode);
		ret->secondNode = copySet (X->secondNode);
	}

	return ret;
}

/*********************************************************/

node *makeNode (int index)
{
	node *newNode= new node;		/* points to new node added to the graph */
	newNode->index = index;
	newNode->index2 = -1;
	newNode->parentEdge = nullptr;
	newNode->leftEdge = nullptr;
	newNode->rightEdge = nullptr;
	/* all fields have been initialized */

	return newNode;
}

/*********************************************************/

edge *makeEdge (node *tail, node *head, double weight)
{
	edge *newEdge = new edge;
	newEdge->tail = tail;
	newEdge->head = head;
	newEdge->distance = weight;
	newEdge->totalweight = 0.0;

	return newEdge;
}

/*********************************************************/

tree *newTree ()
{
	tree *T = new tree;

	T->root = nullptr;
	T->size = 0;
	T->weight = -1;

	return T;
}

/*********************************************************/

/* frees subtree below e->head, recursively, then frees e->head and e */
void deleteSubTree (edge *e)
{
	node *v;
	edge *e1, *e2;

	v = e->head;
	e1 = v->leftEdge;
	if (nullptr != e1)
		deleteSubTree(e1);

	e2 = v->rightEdge;
	if (nullptr != e2)
		deleteSubTree(e2);

	if (nullptr != v)
		deleteNode(v);

	if (nullptr != e)
		delete e;

	return;
}

/*********************************************************/

void deleteTree (tree *T)
{
	if (nullptr != T->root->leftEdge)
		deleteSubTree (T->root->leftEdge);

	if (nullptr != T->root)
		deleteNode(T->root);
	

	if (nullptr != T) delete T;
		T = nullptr;

	return;
}

/*********************************************************/

void deleteSet (set *S)
{
	if (nullptr != S)
	{
		deleteNode (S->firstNode);
		deleteSet (S->secondNode);
		delete S;
	}

	return;
}

/*********************************************************/

void deleteNode (node *n)
{
	if (nullptr != n)
		delete n;

	return;
}

/*********************************************************/

/* copyNode returns a copy of v which has all of the fields identical
 * to those of v, except the node pointer fields */
node *copyNode (node *v)
{
	node *w;

	w = makeNode (v->index);
	w->index2 = v->index2;

	return (w);
}

/*********************************************************/

/* copyEdge calls makeEdge to make a copy of a given edge
 * does not copy all fields */
edge *copyEdge (edge *e)
{
	edge *newEdge;

	newEdge = makeEdge (e->tail, e->head, e->distance);
	newEdge->topsize = e->topsize;
	newEdge->bottomsize = e->bottomsize;
	e->totalweight = 1.0;
	newEdge->totalweight = -1.0;

	return (newEdge);
}

/*********************************************************/

edge *siblingEdge (edge *e)
{
	if (e == e->tail->leftEdge)
		return (e->tail->rightEdge);

	else
		return (e->tail->leftEdge);
}

/*********************************************************/

void updateSizes (edge *e, int direction)
{
	edge *f;

	switch (direction)
	{
		case 1 :
			f = e->head->leftEdge;
			if (nullptr != f)
				updateSizes (f, 1);

			f = e->head->rightEdge;
			if (nullptr != f)
				updateSizes (f, 1);

			e->topsize++;
			break;
		case 2 :
			f = siblingEdge (e);
			if (nullptr != f)
				updateSizes (f, 1);

			f = e->tail->parentEdge;
			if (nullptr != f)
				updateSizes (f, 2);

			e->bottomsize++;
			break;
	}

	return;
}

/*********************************************************/

node *copySubtree (node *v)
{
	node *newNode;
	newNode = copyNode (v);

	if (nullptr != v->leftEdge)
	{
		newNode->leftEdge = copyEdge (v->leftEdge);
		newNode->leftEdge->tail = newNode;
		newNode->leftEdge->head = copySubtree (v->leftEdge->head);
		newNode->leftEdge->head->parentEdge = newNode->leftEdge;
	}
	if (nullptr != v->rightEdge)
	{
		newNode->rightEdge = copyEdge (v->rightEdge);
		newNode->rightEdge->tail = newNode;
		newNode->rightEdge->head = copySubtree (v->rightEdge->head);
		newNode->rightEdge->head->parentEdge = newNode->rightEdge;
	}

	return (newNode);
}

/*********************************************************/

tree *copyTree (tree *T)
{
	tree *Tnew;
	node *n1, *n2, *n3;
	edge *e1, *e2;

	n1 = copyNode (T->root);
	Tnew = newTree ();
	Tnew->root = n1;
	if (nullptr != T->root->leftEdge)
	{
		e1 = copyEdge (T->root->leftEdge);
		n1->leftEdge = e1;
		n2 = copySubtree (e1->head);
		e1->head = n2;
		e1->tail = n1;
		n2->parentEdge = e1;
	}
	if (nullptr != T->root->rightEdge)
	{
		e2 = copyEdge (T->root->rightEdge);
		n1->rightEdge = e2;
		n3 = copySubtree (e2->head);
		e2->tail = n1;
		e2->head = n3;
		n3->parentEdge = e2;
	}

	Tnew->size = T->size;
	Tnew->weight = T->weight;

	return (Tnew);
}

/*********************************************************/

void weighTree (tree *T)
{
	edge *e;

	T->weight = 0;
	for (e=depthFirstTraverse(T,nullptr); nullptr!=e; e=depthFirstTraverse(T,e))
		T->weight += e->distance;

	return;
}

/*********************************************************/

edge *findEdge (tree *T, edge *e)
{
	edge *f;

	// for (f=depthFirstTraverse(T,nullptr); nullptr!=f; f=depthFirstTraverse(T,f))
	// 	if (0 == strcmp (e->label, f->label))
	// 		return (f);

	// Exit ( (char*)"Cannot find edge %s with tail %s and head %s", e->label, e->tail->label, e->head->label);

	return nullptr;
}

/*********************************************************/

node *indexedNode (tree *T, int i)
{
	edge *e;

	for (e=depthFirstTraverse(T,nullptr); nullptr!=e; e=depthFirstTraverse(T,e))
		if (i == e->head->index)
			return (e->head);

	if (i == T->root->index)
		return (T->root);

	return (nullptr);
}

/*********************************************************/

edge *indexedEdge (tree *T, int i)
{
	edge *e;

	for (e=depthFirstTraverse(T,nullptr); nullptr!=e; e=depthFirstTraverse(T,e))
		if (i == e->head->index)
			return (e);

	return (nullptr);
}




