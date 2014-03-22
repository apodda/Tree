#include <math.h>
#include <stdlib.h>
// #include <stdarg.h> // varargs
#include <stdio.h> // printf
#include <assert.h>
#include "../include/tree_hybrid.h"
#include "../include/macros.h"

tree_hybrid *tree_hybrid_create(unsigned int size_tree, unsigned int size_vector)
{
	tree_hybrid *tree = malloc(sizeof(*tree));
	// TODO Error checking?
	
	tree->size_tree = size_tree;
	tree->size_vector = size_vector;
	tree->tree = malloc(sizeof(double) * (size_tree) * (size_tree + 1) * size_vector / 2);
	
	return tree;
}

void tree_hybrid_destroy(tree_hybrid *tree)
{
	free(tree->tree);
	free(tree);
}

// FIXME should these be macros?
static unsigned int _linear_index(tree_hybrid *tree, unsigned int time, unsigned int row, unsigned int space)
{
	return (time * (time + 1) / 2 + row) * tree->size_vector + space;
}

static int _check_boundaries(tree_hybrid *tree, unsigned int time, unsigned int row, unsigned int space)
{
	if (time < tree->size_tree && row <= time && space < tree->size_vector) {
		return 1;
	} else {
		PRINTF("Index %d %d %d out of bounds %d %d %d\n", time, row, space, tree->size_tree, time, tree->size_vector);
		return 0;
	}
}

int tree_hybrid_get_safe(tree_hybrid *tree, double *result, unsigned int time, unsigned int row, unsigned int space)
{
	// Returns 0 on correct execution, 1 if an index is out of bonds
	if(_check_boundaries(tree, time, row, space)) {
		*result = tree->tree[_linear_index(tree, time, row, space)];
		return 0;
	} else {
		return 1;
	}
}

double tree_hybrid_get(tree_hybrid *tree, unsigned int time, unsigned int row, unsigned int space)
{
	return tree->tree[_linear_index(tree, time, row, space)];
}

double* tree_hybrid_get_vector(tree_hybrid *tree, unsigned int time, unsigned int row)
{
	return tree->tree + _linear_index(tree, time, row, 0);
}

int tree_hybrid_set_safe(tree_hybrid *tree, double value, unsigned int time, unsigned int row, unsigned int space)
{
	// Returns 0 on correct execution, 1 if an index is out of bonds
	if(_check_boundaries(tree, time, row, space)) {
		tree->tree[_linear_index(tree, time, row, space)] = value;
		return 0;
	} else {
		return 1;
	}
}

int tree_hybrid_print_slice(tree_hybrid *tree, unsigned int time)
{
	double tmp;
	int row, space;
	
	if(time >= tree->size_tree) {
		return 1;
	}
	
	for(row = 0; row <= time; row++) {
		for(space = 0; space < tree->size_vector; space++) {
			if (!tree_hybrid_get_safe(tree, &tmp, time, row, space)) {
				printf("%f ", tmp);
			} else {
				assert(0);
			}
		}
		printf("\n");
	}
	return 0;
}
