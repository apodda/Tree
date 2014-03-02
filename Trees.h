#include <math.h>
#include <stdlib.h>

typedef struct tree_2d {
	unsigned int rows;
	unsigned int columns;
	double* tree;
} tree_2d;

tree_2d create_tree_2d(unsigned int rows, unsigned int columns)
{
	struct tree_2d *tree = malloc(sizeof(*tree));
	
	tree->rows = rows;
	tree->columns = columns;
	tree->tree = malloc(sizeof(double) * (rows) * (columns));
	
	return *tree;
}

int tree_2d_get(double *result, unsigned int row, unsigned int column, tree_2d *tree)
{
	// Returns 0 on correct execution, 1 if an index is out of bonds
	if(row < tree->rows && column < tree->columns){
		*result = tree->tree[row + column * tree->rows];
		return 0;
	} else {
		return 1;
	}
}

int tree_2d_set(double value, unsigned int row, unsigned int column, tree_2d *tree)
{
	// Returns 0 on correct execution, 1 if an index is out of bonds
	if(row < tree->rows && column < tree->columns){
		tree->tree[row + column * tree->rows] = value;
		return 0;
	} else {
		return 1;
	}
}
