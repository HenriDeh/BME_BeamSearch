#include <iostream>
#include "fastme.h"



void run(double **d, int **init_adj, int* solution_mat, int n_taxa, int m, double & obj_val, int & nni_count, int & spr_count)
{

	time_t t_beg, t_end;

	tree *T = nullptr;

	int **sparse_A;
	sparse_A = get_sparse_A(init_adj, m);


	T = adj_to_tree(sparse_A, n_taxa, m);
	for (int i = 0; i < m; i++)
		delete[] sparse_A[i];
	delete[] sparse_A;

	int i = 0;

	int numSpecies = n_taxa;
	int setCounter = 0;

	int nniCount, sprCount, repCounter, repToPrint, printedRep;

	int seqLength = -1;

	// Random numbers for bootstraps
	int **rnd = nullptr;

	double **A;
	A = nullptr;

	nniCount = sprCount = repCounter = 0;
	setCounter++;

	A = initDoubleMatrix(m);

	T = ImproveTree(T, d, A, &nniCount, &sprCount);

	tree_to_adj_mat(T, solution_mat);

	// explainedVariance (D, T, numSpecies, options->precision, options->input_type, options->fpO_stat_file);

	deleteMatrix(A, m);

	obj_val = T -> weight;
	nni_count = nniCount;
	spr_count = sprCount;

	deleteTree(T);

}

tree *ImproveTree(tree *T0, double **D, double **A,
				  int *nniCount, int *sprCount)
{
	// T0 = ME
	// T1 = ME + NNI
	// T2 = ME + SPR
	// T3 = ME + NNI + SPR
	tree *T1, *T2, *T3;

	T1 = T2 = T3 = nullptr;

	// Print tree length
	weighTree(T0);

	// std::cout << "Tree length is " << T0->weight << std::endl;

	// std::cout << "NNI" << std::endl;

	T1 = copyTree(T0);

	makeBMEAveragesTable(T1, D, A);

	bNNI(T1, A, nniCount);
	assignBMEWeights(T1, A);

	// std::cout << "Performed NNI(s) " << *nniCount << std::endl;

	// std::cout << "SPR" << std::endl;

	T2 = copyTree(T0);

	makeBMEAveragesTable(T2, D, A);
	SPR(T2, D, A, sprCount);
	assignBMEWeights(T2, A);

	// std::cout << "Performed SPR(s) " << *sprCount << std::endl;

	if (nullptr != T1)
		weighTree(T1);

	if (nullptr != T2)
		weighTree(T2);

	if (nullptr != T1)
	{
		if (T0->weight > T1->weight)
		{
			deleteTree(T0);
			T0 = T1; // using T0 as the place to store the minimum evolution tree in all cases
			T1 = nullptr;
		}
	}
	else if (nullptr != T1)
		deleteTree(T1);

	if (nullptr != T2)
	{
		if (T0->weight > T2->weight)
		{
			deleteTree(T0);
			T0 = T2;
			T2 = nullptr;
		}
	}
	else if (nullptr != T2)
		deleteTree(T2);

	return T0;
}