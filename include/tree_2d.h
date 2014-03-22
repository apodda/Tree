#ifndef TREES_2D_H_
#define TREES_2D_H_

#include <math.h>
#include <stdlib.h>

typedef struct tree_2d {
	unsigned int size;
	double* tree;
} tree_2d;

// Constructors
tree_2d *tree_2d_create(unsigned int size);

tree_2d *tree_2d_create_standard(unsigned int size, double h, double start);

tree_2d *tree_2d_create_fun(unsigned int size, double h, double start, double (*fun)(unsigned int, unsigned int, double));

// Destructor
void tree_2d_destroy(tree_2d *tree);

// Getters and setters
int tree_2d_get_safe(tree_2d *tree, double *result, unsigned int row, unsigned int column);

double tree_2d_get(tree_2d *tree, unsigned int row, unsigned int column);

int tree_2d_set_safe(tree_2d *tree, double value, unsigned int row, unsigned int column);

// Print
int tree_2d_print(tree_2d *tree);
#endif
