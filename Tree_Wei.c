#include <math.h>
#include <stdlib.h>

// Returns the biggest l such that Y_ik + mu * h >= Y_i+1,l
int jump_index(int i, int k, double h, double mu)
{
	double tmp = (mu * sqrt(h) + 1) / 2;

	if (tmp > 0) {
		return k + (int) floor(tmp);
	} else {
		return k - (int) -floor(-tmp);
	}
}

double jump_probability(int i, int k, double h, double mu, double **tree, int m, int N)
{
	// TODO Error handling (ie, i >= m, k >= N...)
	
	double probability;
	int kd = jump_index(i, k, h, mu);
	int ku = jump_index(i, k, h, mu);
	
	// C indexes double arrays row-first, so k goes before i
	probability = mu * h + tree[k][i] - tree[kd][i+1];
	probability = probability / (tree[ku][i+1] - tree[kd][i+1]);
	probability = fmax(probability, 0.0);
	probability = fmin(probability, 1.0);
	
	return probability;
}

double interest_rate(int i, int k, double sigma_r, double **tree)
{
	// TODO Error handling
	double interest_rate = tree[k][i] * tree[k][i] * sigma_r * sigma_r * 0.25;
	
	return fmax(interest_rate, 0.0);
}

double asset_price(int i, int j, int k, double sigma_s, double ro, double** tree_y, double** tree_R)
{
	return exp(sigma_s * (sqrt(1 - ro * ro) * tree_y[j][i] + ro * tree_R[k][i]));
}

int main (int argc, char *argv[])
{
	int m = 10;
	int N = m;
	int i, j, k;
	
	double maturity = 1;
	double h = maturity / m;
	
	double sigma_s = 100;
	double sigma_r = 0.5;
	double ro = 0.1;
	
	double start_s = 100;
	double start_r = 0.5;
	
	double y_tree[m][m];
	double r_tree[m][m];
	
	double start_y = (log(start_s) / sigma_s 
			- 2 * ro * sqrt(start_r) / sigma_r) / sqrt(1 - ro * ro);
	double start_R = 2 * sqrt(start_r) / sigma_r;
	
	for (i = 0; i < N; i += 1) {
		for (j = 0; j < i; j += 1) {
			// We're building two different trees of the same size
			// C uses row-first ordering...
			r_tree[j][i] = start_R + (2 * j - i) * sqrt(h);
			y_tree[j][i] = start_y + (2 * j - i) * sqrt(h);
		}
	}
	
	return 0;
}
