#ifndef TREES_HYBRID_H_
#define TREES_HYBRID_H_

#include <math.h>
#include <stdlib.h>

typedef struct tree_hybrid {
	unsigned int size_tree;
	unsigned int size_vector;
	double* tree;
} tree_hybrid;

// Constructors
tree_hybrid *tree_hybrid_create(unsigned int size_tree, unsigned int size_vector);

// Destructor
void tree_hybrid_destroy(tree_hybrid *tree);

// Getters and setters
int tree_hybrid_get_safe(tree_hybrid *tree, double *result, unsigned int time, unsigned int row, unsigned int column);

double tree_hybrid_get(tree_hybrid *tree, unsigned int time, unsigned int row, unsigned int column);

int tree_hybrid_set_safe(tree_hybrid *tree, double value, unsigned int time, unsigned int row, unsigned int column);

double* tree_hybrid_get_vector(tree_hybrid *tree, unsigned int time, unsigned int column);

// Pretty-printers

int tree_hybrid_print_slice(tree_hybrid *tree, unsigned int time);

#endif
