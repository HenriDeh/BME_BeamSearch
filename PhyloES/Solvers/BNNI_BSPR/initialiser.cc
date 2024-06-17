#include "initialiser.h"
#include <iostream>

void make_set(int n_taxa, set *S) {
	node *v = nullptr;
	for(int i=0; i< n_taxa; i++) {
		char c= i + '0';
		v = makeNode (-1);
		v->index2 = i;
		S = addToSet (v, S);
	}

}


void recursion_print (edge *e)
{	
	if (nullptr!= e) {
	std::cout<<"node "<<e->head->index2<<" "<< e->head->index<<std::endl;
	recursion_print(e->head->leftEdge);
	recursion_print(e->head->rightEdge);
	}	
}


void printTree (tree *T)
{	
	// if(T-> root -> leftEdge != nullptr) printf("here we go");
	printf("size %i \n", T->size);
	printf("node %i %i\n" ,T->root->index2, T->root->index);
	recursion_print(T->root->leftEdge);
	recursion_print(T->root->rightEdge);

}

int recursion_fill_adj (int size, int from, int internal_idx, int* A, edge *e){
	if (nullptr!= e) {

		if(e->head->index2 < 0){
			A[from* size + internal_idx] = 1;
			A[internal_idx*size + from] = 1;	
			from = internal_idx;
			internal_idx ++;	
		}
		else{
			A[from * size + e->head->index2] = 1;
			A[e->head->index2* size + from] = 1;
			from = e->head->index2;
		}

	
	internal_idx = recursion_fill_adj(size, from, internal_idx, A, e->head->leftEdge);
	internal_idx = recursion_fill_adj(size, from, internal_idx, A, e->head->rightEdge);
	
	}	
	return internal_idx;

}

void tree_to_adj_mat(tree *T, int* solution_mat) {
	int size = T->size;
    
    for(int i = 0;i<size;i++) {
		for(int j=0; j<size; j++) solution_mat[i*size + j] = 0;
    }
	int internal_idx = (T-> size + 2) / 2;
	internal_idx = recursion_fill_adj(size, T->root->index, internal_idx, solution_mat, T->root->leftEdge);
	recursion_fill_adj(size, T->root->index, internal_idx, solution_mat, T->root->rightEdge);

	// std::cout<<"solution"<<std::endl;
	for(int i = 0;i<size;i++) {
		for(int j = 0; j < size; j++) {
			// std::cout<< A[i][j] <<" ";
		}
		// std::cout<<std::endl;
	}
	
}

int** get_sparse_A(int** A, int n_nodes) {
    int** sparse_A;
    sparse_A = new int*[n_nodes]; 
    int k = 0;

    for(int i = 0; i<n_nodes; i++) {
        sparse_A[i] = new int[3];
        sparse_A[i][0] = -1; sparse_A[i][1] = -1; sparse_A[i][2] = -1; 

        for(int j = 0; j< n_nodes; j++){
			
            if(A[i][j]==1) {
                sparse_A[i][k] = j;
                k++;
            }
        }
        k = 0; 
    }


	return sparse_A;
}



void adj_to_tree_recursion(int mat_index, int parent_mat_index, int *node_index, int *leaf_index, int n_taxa, int tree_size, int** sparse_A, node* parent, int is_left) {
	
	
	if(*node_index + 1 < tree_size){
		// printf("%i %i %i\n", *node_index, is_left, mat_index);
		*node_index = *node_index + 1;
		// printf("%i p\n", *node_index);
		node* child;
		child = makeNode(*node_index);
		edge* new_edge;
		new_edge = makeEdge (parent, child, 1); // weight 1

		child -> parentEdge = new_edge;

		if(is_left == 0) parent -> leftEdge = new_edge;
		else parent -> rightEdge = new_edge;
		

		if (sparse_A[mat_index][1] == -1) { // is a leaf
			*leaf_index = mat_index;
			child -> index2 = mat_index;
			}
		else{
			is_left = 0;
			for(int i=0; i< 3; i++){
				if(sparse_A[mat_index][i] != parent_mat_index)
					{
					adj_to_tree_recursion(sparse_A[mat_index][i], mat_index, node_index, leaf_index, n_taxa, tree_size, sparse_A, child, is_left);
					is_left ++;
					}
			}
		}

	}
	
}

tree* adj_to_tree(int** sparse_A, int n_taxa, int m) {
	tree *T;
	node *root;

	T = newTree();
	T -> size = m;
    int node_index = 0;
	root = makeNode(0);
	root->index2 = 0;
	T -> root = root;
	int leaf_index = 1;
	adj_to_tree_recursion(sparse_A[0][0], 0, &node_index, &leaf_index, n_taxa, T->size, sparse_A, root, 0);
	// printf("%i root %i\n", root-> leftEdge -> head->index, root-> leftEdge -> tail->index);
	return T;
}