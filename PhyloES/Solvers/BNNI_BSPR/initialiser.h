
#include "graph.h"


void make_set(int n_taxa, set *S);
void recursion_print (edge *e);
void printTree (tree *T);
int recursion_fill_adj (int size, int from, int internal_idx, int* A, edge *e);
void tree_to_adj_mat(tree *T, int* solution_mat);
int** get_sparse_A(int** A, int n_nodes);
void adj_to_tree_recursion(int mat_index, int parent_mat_index, 
                    int *node_index, int *leaf_index, int n_taxa, int tree_size, int** sparse_A, node* parent, int is_left);
tree* adj_to_tree(int** sparse_A, int n_taxa, int m);