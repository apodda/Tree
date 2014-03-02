#include <math.h>
#include <stdlib.h>

typedef struct tree_2d {
	unsigned int size;
	double* tree;
} tree_2d;

tree_2d tree_2d_create(unsigned int size)
{
	struct tree_2d *tree = malloc(sizeof(*tree));
	
	tree->size = size;
	tree->tree = malloc(sizeof(double) * (size) * (size + 1) / 2);
	
	return *tree;
}

tree_2d tree_2d_create_standard(unsigned int size, double h, double start)
{
	int row = 0;
	int column = 0;
	struct tree_2d tree = tree_2d_create(size);
	
	// We *don't* use tree_2d_set for speed reasons (?)
	for (column = 0; column < tree.size; column += 1)
	{
		for (row = 0; row <= column; row += 1)
		{
			tree.tree[column * (column + 1) / 2 + row] = start + (2 * column - row) * sqrt(h);
		}
	}
	
	return tree;
}

int tree_2d_get(tree_2d *tree, double *result, unsigned int row, unsigned int column)
{
	// Returns 0 on correct execution, 1 if an index is out of bonds
	if(row < tree->size && row <= column) {
		// We use column-wise ordering:
		// 0 1 3 6 ...
		//   2 4 7 ...
		//     5 8 ...
		//       9 ...
		// The (linear) index of the first element of each column is
		// obtained with the summation formula 
		// 0 + 1 + 2 + ... + i = i * (i+1) * 0.5
		*result = tree->tree[column * (column + 1) / 2 + row];
		return 0;
	} else {
		return 1;
	}
}

int tree_2d_set(tree_2d *tree, double value, unsigned int row, unsigned int column)
{
	// Returns 0 on correct execution, 1 if an index is out of bonds
	if(row < tree->size && row <= column) {
		tree->tree[column * (column + 1) / 2 + row] = value;
		return 0;
	} else {
		return 1;
	}
}
