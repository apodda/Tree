#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "../include/tree_2d.h"
#include "../include/tree_hybrid.h"
#include "../include/macros.h"

void print_vector(double *vec, int size)
{
	int i = 0;
	for(i = 0; i < size; i++)
	{
		printf("%f ", vec[i]);
	}
	printf("\n");
}

typedef struct heston{
	double interest_rate;
	double dividend; // delta
	double kappa; // mean_reversion_rate
	double theta; // long run variance
	double volvol; // sigma_V
	double ro;
} heston;

// Returns the biggest l such that S_ij + mu * h >= S_i+1,l
/*int s_jump_index(int i, int j, double h, double sigma_Y, double r)
{
	double tmp = log(1 + h * r) / sigma_Y;
	
	return (int) floor(j + 1 - tmp);
}

// Returns the biggest l such that S_ik + mu * h >= S_i+1,l
int r_jump_index(int i, int k, double h, double sigma_V, double theta, double R_0)
{
	double tmp;
	
	tmp = R_0 + (2 * k - i) * sqrt(h);
	tmp = (1 / h - k) * tmp * tmp + theta * k / (sigma_V * sigma_V);
	tmp = sqrt(tmp) - R_0 / sqrt(h) + i + 1;
	
	return (int) floor(tmp / 2);
}*/

// Naive implementation of
// max(0 <= j_star <= j | S_ij + h * mu >= S_{i+1,j_star})
unsigned int jump_down_index(unsigned int i, unsigned int j, double mu, double h, tree_2d * tree)
{
	int j_star;
	double step = tree_2d_get(tree, i, j) + mu * h;
	double tmp;
	
	//return j;
	
	for(j_star = j; j_star >= 0; --j_star){
		tmp = tree_2d_get(tree, i+1, j_star);
		
		if(step >= tmp){
			return j_star;
		}
	}
	return 0;
}

// Naive implementation of
// min(j+1 <= j_star <= i+1 | S_ij + h * mu <== S_{i+1,j_star})
unsigned int jump_up_index(unsigned int i, unsigned int j, double mu, double h, tree_2d * tree)
{
	int j_star;
	double step = tree_2d_get(tree, i, j) + mu * h;
	double tmp;
	
	//return j+1;
	
	for(j_star = j+1; j_star <= i+1; ++j_star){
		tmp = tree_2d_get(tree, i+1, j_star);
		
		if(step <= tmp){
			return j_star;
		}
	}
	return i+1;
}

double jump_probability(unsigned int i,unsigned  int k, double h, double mu, tree_2d *tree)
{
	// TODO Error handling (ie, i >= m, k >= N...)
	
	double probability;
	int kd = jump_down_index(i, k, mu, h, tree);
	int ku = jump_up_index(i, k, mu, h, tree);
	
	probability = mu * h + tree_2d_get(tree, k, i) - tree_2d_get(tree, kd, i+1);
	probability = probability / (tree_2d_get(tree, ku, i+1) - tree_2d_get(tree, kd, i+1));
	probability = fmax(probability, 0.0);
	probability = fmin(probability, 1.0);
	
	return probability;
}

// LU factorization **without** pivoting. Will happily crash should pivoting be needed
void lu_factorization(double *matrix, unsigned int size, unsigned int stride)
{
	int row, column;
	while(size > 2){
		// For column-major indexing:
		// matrix(row, column) = matrix[row + stride * column]	
		for(row = 1; row <= size; row++){
			matrix[row] *= matrix[row] / matrix[1];
		}
		for(row = 1; row <= size; row++){
			for (column = 1; column <= size; column++){
				// Update the lower-right submatrix, subtracting
				// the product of the first column (minus the
				// first element), and the first row (minus the
				// first element)
				matrix[row + stride*column] -= matrix[row] * matrix[stride * column];
			}
		}
		// Repeat the same algorithm on the lower-right submatrix
		size--;
		matrix = matrix + 1 + stride;
	}
}

void lu_solve(double *matrix, double *solution, double *known_term, unsigned int size)
{
	int row, column;
	double tmp;
	// Solve Lx' = b
	for(row = 0; row < size; row++){
		tmp = known_term[row];
		for(column = 0; column < row; column++){
			// For column-major indexing:
			// matrix(row, column) = matrix[row + size * column]
			tmp -= matrix[row + size * column] * solution[column];
		}
		// We know that l_ii = 1
		known_term[row] = tmp;
	}
	// Solve Ux = x'
	for(row = size - 1; row >= 0; row--){
		tmp = known_term[row];
		for(column = row+1; column < size; column++)
		{
			tmp -= matrix[row + size * column] * solution[column];
		}
		solution[row] = tmp / matrix[row + size * row];
	}
}

// Algorithm for solving tridiagonal systems. Will happily crash on a zero pivot
/*void thomas_algorithm(double *sub_diag, double *diag, double *super_diag, double *result, double *known_term, unsigned int size)
{
	int i;
	double *gamma = malloc(sizeof(double) * size);
	double beta = diag[0];
	
	result[0] = known_term[0] / diag[0];
	for(i = 1; i < size; i++) {
		gamma[i] = super_diag[i-1] / beta;
		beta = diag[i] - sub_diag[i-1] * gamma[i];
		result[i] = (known_term[i] - sub_diag[i-1] * result[i-1]) / beta;
	}
	
	for(i = size - 2; i >= 0; i--) {
		result[i] -= gamma[i+1] * result[i+1];
	}
	free(gamma);
}*/

void thomas_algorithm(double *sub_diag, double *diag, double *super_diag, double *result, double *known_term, unsigned int size)
{
	int i;
	double *gamma = malloc(sizeof(double) * size);
	double *delta = malloc(sizeof(double) * size);

	gamma[0] = super_diag[0] / diag[0];
	delta[0] = known_term[0] / diag[0];
	
	for (i = 1; i < size-1; i++)
	{
		gamma[i] = super_diag[i] / (diag[i] - gamma[i - 1] * sub_diag[i - 1]);
		delta[i] = (known_term[i] - delta[i-1] * sub_diag[i - 1]) / (diag[i] - gamma[i-1] * sub_diag[i - 1]);
	}
	
	delta[size - 1] = (known_term[size - 1] - delta[size - 2] * sub_diag[size - 2]) / (diag[size - 1] - gamma[size - 2] * sub_diag[size - 2]);
	
	result[size - 1] = delta[size - 1];
	for (i = size - 2; i > -1; i--)
	{
		result[i] = delta[i] - gamma[i] * result[i + 1];
	}
	free(delta);
	free(gamma);
}

/*void thomas_algorithm_neumann(double convection, double diffusion, double *result, double *known_term, unsigned int size)
{
	int i;
	double a, b, c;
	double beta, gamma;

	
	b = 1 + 2 * diffusion;
	c = - 2 * diffusion;
	
	// i = 0
	beta = b;
	result[0] = known_term[0] / b;
	
	// i = 1
	gamma = c / beta;
	beta = b - a * result[0];
	result[1] = (known_term[1] - a * result[0]) / beta;
	
	c = - convection - diffusion;
	for (i = 2; i < size - 1; i++)
	{
		gamma = c / beta;
		beta = b - a * result[i-1];
		result[i] = (known_term[i] - a * result[i-1]) / beta;
	}
	// i = size -1
	a = - 2 * diffusion;
	gamma = c / beta;
	beta = b - a * result[0];
	result[1] = (known_term[1] - a * result[0]) / beta;
	
	//
}*/

void explicit_step(double alpha, double beta, double *result, double *known_term, unsigned int size)
{
	int i;
	double alpha_abs = fabs(alpha);
	double alpha_down = (alpha < 0 ? alpha_abs : 0);
	double alpha_up = (alpha > 0 ? alpha_abs : 0);
	
	double diag = 1 - 2 * beta - 2 * alpha_abs;
	double sub_diag = beta + 2 * alpha_down;
	double super_diag = beta + 2 * alpha_up;
	
	double tmp = 2 * beta + 2 * alpha_abs;
	
	// The first and the last row are special, due to boundary conditions
	result[0] = diag * known_term[0] + tmp * known_term[1];

	for (i = 1; i < size - 1; i++)
	{
		result[i] = known_term[i - 1] * sub_diag +
			known_term[i] * diag + known_term[i+1] * super_diag;
	}
	
	result[size-1] = known_term[size-2] * tmp + known_term[size-1] * diag;
}

void initialize_tree_V(tree_2d *tree, double start, double sigma, double h)
{
	int n, k;
	double tmp = 0;
	
	for (n = 0; n < tree->size; n++)
	{
		for (k = 0; k <= n; k++)
		{
			tmp = sqrt(start) + (sigma * 0.5) * (double) (2 * k - n) * sqrt(h);
			tmp = fmax(tmp, 0.0);
			tmp = tmp * tmp;
			
			if(tree_2d_set_safe(tree, tmp, n, k)) {
				abort();
			}
		}
	}
}

// The result is the matrix associated with the spatial discretization of a 1D 
// advection-diffusion problem with centered differences, and Neumann conditions
void initialize_matrix(double *sub_diag, double *diag, double *super_diag, int size, double alpha, double beta)
{
	int i;
	
	diag[0] = 1 + 2 * beta;
	super_diag[0] = - 2 * beta;
	
	for (i = 1; i < size - 1; i++)
	{
		sub_diag[i-1] = alpha - beta;
		diag[i] = 1 + 2 * beta;
		super_diag[i] = - alpha - beta;
	}
	
	diag[size-1] = 1 + 2 * beta;
	sub_diag[size-2] = - 2 * beta;
	
	//print_vector(sub_diag, size);
	//print_vector(diag, size);
	//print_vector(super_diag, size);
}

double payoff(double S, double strike)
{
	return fmax(strike - S, 0.0);
}

// TODO remove/polish
void tridiagsolver(double a[], double b[], double c[], double bind[], int size) {
    int n = size, i = 0;
    c[0] /= b[0];
    bind[0] /= b[0];
    for (i=1; i<n; i++) {
        c[i] /= b[i]-a[i]*c[i-1];
        bind[i] = (bind[i]-a[i]*bind[i-1])/(b[i]-a[i]*c[i-1]);
    }
    bind[n] = (bind[n]-a[n]*bind[n-1])/(b[n]-a[n]*c[n-1]);
    for (i=n-1; i>=0; i--) {
        bind[i] -= c[i]*bind[i+1];
    }
}

int main (int argc, char *argv[])
{
	unsigned int N = 100;
	int i, j, k;
	unsigned int ku, kd;
	double eps = 0.0001;
	
	double tmp_payoff;
	double tmp_y;
	double tmp_v;
	double p_up, p_down;

	// Model parameters
	// Fixme move those into a function;
	double ro = -0.5;
	double r = log(1.1);
	double kappa = 2;
	double theta = 0.1;
	double strike = 100;
	double delta = 0; // dividend
	double maturity = 1.0;
	
	double sigma_V = 0.04;
	
	double start_S = 100;
	double start_V = 0.1;
	
	// Derived parameters
	double h = maturity / (double) N;
	// FIXME choose a better grid for y
	int M = N;
	//double boundary = 10;
	double delta_y = h;
	
	double start_Y = log(start_S) - ro / sigma_V * start_V;
	
	double mu_Y, mu_V;
	
	double *row_up, *row_down;
	
	tree_2d *tree_V = tree_2d_create(N+1);
	initialize_tree_V(tree_V, start_V, sigma_V, h);
	
	tree_hybrid *tree_P = tree_hybrid_create(N + 1, 2 * M + 1);
	double alpha, beta;
	double *tmp_vector = malloc(sizeof(double) * tree_P->size_vector);
	double *result;
	double v;
	
	double *super_diag = malloc(sizeof(double) * (tree_P->size_vector - 1));
	double *diag = malloc(sizeof(double) * tree_P->size_vector);
	double *sub_diag = malloc(sizeof(double) * (tree_P->size_vector - 1));
	
	// FIXME Ensure tree_P->size_tree == tree_V->size
	for(j = 0; j < tree_P->size_tree; j++)
	{
		for(k = 0; k < tree_P->size_vector; k++)
		{
			tmp_y = start_Y + delta_y * (k - M);
			//printf("%f ", tmp_y);
			if(tree_2d_get_safe(tree_V, &tmp_v, tree_P->size_tree-1, j)) {
				PRINTF("Error: ");
				// TODO Handle index-out-of-bonds errors
			}
			
			tmp_payoff = payoff(exp(tmp_y + ro / sigma_V * tmp_v), strike);
			if(tree_hybrid_set_safe(tree_P, tmp_payoff, tree_P->size_tree-1, j, k)) {
				PRINTF("Error: ");
				// TODO Handle index-out-of-bonds errors
			}
		}
		//printf("\n");
	}
	
	for(i = tree_P->size_tree-2; i > -1; i--) {
		// Scorri v (ie le colonne):
		//   trova ku, kd
		//     fai un passo di differenze finite con dato iniziale i+1,ku,:
		//     idem con i+1,kd,:
		//   assegna P(i, :, v) = p_kd * diff. finite kd + p_ku * diff. finite ku
		for(j = 0; j <= i; j++) {
			if(tree_2d_get_safe(tree_V, &v, i, j)) {
				abort();
			}
			//PRINTF("---\n");
			//PRINTF("The interest rate at index %d, %d is %lf \n", i, j, v);
			
			// FIXME move all into a struct...
			mu_Y = r - delta - 0.5 * v - ro / sigma_V * kappa * (theta - v);
			mu_V = kappa * (theta - v);
			
			alpha = h * 0.5 / delta_y * mu_Y;
			beta = h * 0.5 / (delta_y * delta_y) * (1 - ro * ro) * v;
			//printf("mu_Y = %f, mu-V = %f\n", mu_Y, mu_V);
			
			ku = jump_up_index(i, j, mu_V, h, tree_V);
			row_up = tree_hybrid_get_vector(tree_P, i+1, ku);
			p_up = jump_probability(i, j, h, mu_V, tree_V);
			
			/*if(ku != j + 1) {
				printf("Jumping high at index %d, %d: ku = %d\n", i, j, ku);
			}*/
			
			//PRINTF("The up jump index at index %d, %d is %d \n", i, j, ku);
			//PRINTF("The up jump probability at index %d, %d is %lf \n", i, j, p_up);
			
			kd = jump_down_index(i, j, mu_V, h, tree_V);
			row_down = tree_hybrid_get_vector(tree_P, i+1, kd);
			p_down = 1 - p_up;
			
			/*if(kd != j) {
				printf("Jumping low at index %d, %d: kd = %d\n", i, j, kd);
			}*/
			
			//PRINTF("The down jump index at index %d, %d is %d \n", i, j, kd);
			//PRINTF("The down jump probability at index %d, %d is %lf \n", i, j, p_down);
			
			// FIXME Check indexing
			result = tree_hybrid_get_vector(tree_P, i, j);
			
			//tmp_vector = result;
			
			for(k = 0; k < tree_P->size_vector; k++) {
				tmp_vector[k] = exp(-r * h) * (p_up * row_up[k] + p_down * row_down[k]);
			}
			
			initialize_matrix(sub_diag, diag, super_diag, tree_P->size_vector, alpha, beta);
			/*diag[0] = 1 + 2 * beta;
			super_diag[0] = -2 * beta;
			for (k=0; k<=tree_P->size_vector - 2; k++) {
				sub_diag[k] = alpha - beta;
				diag[k] = 1 + 2 * beta;
				super_diag[k] = -alpha - beta;
			}
			sub_diag[tree_P->size_vector-1] = -2 * beta;
			diag[tree_P->size_vector-1] = 1 + 2*beta;*/
			
			if (v > eps) {
				// A u_n = u_n+1
				// FIXME Check size and matrix initialization
				thomas_algorithm(sub_diag, diag, super_diag, result, tmp_vector, tree_P->size_vector);
				//tridiagsolver(sub_diag, diag, super_diag, tmp_vector, tree_P->size_vector-1);
				
			} else {
				// u_n = A u_n+1
				explicit_step(alpha, beta, result, tmp_vector, tree_P->size_vector);
				//abort(); // FIXME Remove !!!
			}
			
			/*for(k = 0; k < tree_P->size_vector; k++) {
				printf("%f ", result[k]);
			}
			printf("\n");cd
			abort();*/
		}
	}
	
	/*for(i = tree_P->size_tree-1; i >= 0; i--) {
		printf("--- %d ---\n", i);
		tree_hybrid_print_slice(tree_P, i);
	}*/
	
	double tmp;
	if(!tree_hybrid_get_safe(tree_P, &tmp, 0, 0, M)) {
		printf("The price is %f\n", tmp);
	} else {
		printf("Error: Index out of bonds\n");
	}
	
	//tree_hybrid_print_slice(tree_P, 0);
	
	free(sub_diag);
	free(diag);
	free(super_diag);
	
	return 0;
}
