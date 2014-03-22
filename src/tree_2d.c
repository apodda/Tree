#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <assert.h>
#include "../include/tree_2d.h"
#include "../include/macros.h"

// Constructors and destructors

tree_2d *tree_2d_create(unsigned int size)
{
	tree_2d *tree = malloc(sizeof(*tree));
	// TODO Error checking?
	
	tree->size = size;
	tree->tree = malloc(sizeof(double) * (size) * (size + 1) / 2);
	
	return tree;
}

void tree_2d_destroy(tree_2d *tree)
{
	free(tree->tree);
	free(tree);
}

tree_2d *tree_2d_create_standard(unsigned int size, double h, double start)
{
	int row = 0;
	int time = 0;
	double tmp;
	
	struct tree_2d *tree = tree_2d_create(size);
	
	for (time = 0; time < tree->size; time += 1)
	{
		for (row = 0; row <= time; row += 1)
		{
			tmp = start + (2 * row - time) * sqrt(h);
			tree_2d_set_safe(tree, tmp, time, row);
		}
	}
	
	return tree;
}

tree_2d *tree_2d_create_map(unsigned int size, double h, double start, double (*fun)(unsigned int, unsigned int, double))
{
	int row = 0;
	int time = 0;
	double tmp;
	
	struct tree_2d *tree = tree_2d_create(size);
	
	for (time = 0; time < tree->size; time += 1)
	{
		for (row = 0; row <= time; row += 1)
		{
			tmp = start + fun(time, row, h);
			tree_2d_set_safe(tree, tmp, time, row);
		}
	}
	
	return tree;
}

// Getters and setters
// FIXME should these be macros?
static int _check_boundaries(tree_2d *tree, unsigned int time, unsigned int row)
{
	// Returns 0 on correct execution, 1 if an index is out of bonds
	if(time < tree->size && row <= time){
		return 1;
	}
	else {
		PRINTF("Index %d %d out of bounds %d %d\n", time, row, tree->size, time);
		return 0;
	}
}

static unsigned int _linear_index(unsigned int time, unsigned int row)
{
	// We use column-major ordering. The columns reprent time:
	//       9 ...
	//     5 8 ...
	//   2 4 7 ...
	// 0 1 4 6 ...
	// The (linear) index of the first element of each time is
	// obtained with the summation formula 
	// 0 + 1 + 2 + ... + i = i * (i+1) * 0.5
	return time * (time + 1) / 2 + row;
}


int tree_2d_get_safe(tree_2d *tree, double *result, unsigned int time, unsigned int row)
{
	if(_check_boundaries(tree, time, row)) {
		*result = tree->tree[_linear_index(time, row)];
		return 0;
	} else {
		return 1;
	}
}

double tree_2d_get(tree_2d *tree, unsigned int time, unsigned int row)
{
	return tree->tree[_linear_index(time, row)];
}

int tree_2d_set_safe(tree_2d *tree, double value, unsigned int time, unsigned int row)
{
	// Returns 0 on correct execution, 1 if an index is out of bonds
	if(_check_boundaries(tree, time, row)) {
		tree->tree[_linear_index(time, row)] = value;
		return 0;
	} else {
		return 1;
	}
}

int tree_2d_print(tree_2d *tree)
{
	double tmp;
	int time, row;
	
	for(time = 0; time < tree->size; time++) {
		for(row = 0; row <= time; row++) {
			if (!tree_2d_get_safe(tree, &tmp, time, row)) {
				printf("%f ", tmp);
			} else {
				assert(0);
			}
		}
		printf("\n");
	}
	return 0;
}
