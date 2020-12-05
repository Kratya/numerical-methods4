#include "slae.h"

SLAE::SLAE(int size)
{
	n = size;
	A = new vectorD[n];
	for (int i = 0; i < n; ++i)
	{
		A[i] = new double[n];
	}
	b = new double[n];
	x = new double[n];
}

SLAE::~SLAE()
{
	for (int i = 0; i < n; ++i)
	{
		delete[] A[i];
	}
	delete[] A;
	delete[] x;
	delete[] b;
}

double Scal(vectorD a, vectorD b, int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += a[i] * b[i];
	}
	return sum;
}

double Norm(vectorD a, int n)
{
	return sqrt(Scal(a, a, n));
}

bool GausSolver(SLAE* S)
{
	int n = S->n;
	vectorD* matrix = S->A;
	vectorD b = S->b;
	vectorD x = S->x;
	double max;
	const double eps_gauss = 1e-14;
	for (int k = 0; k < n; ++k)
	{
		max = fabs(matrix[k][k]);
		int ind = k;
		for (int i = k + 1; i < n; i++)
		{
			if (fabs(matrix[i][k]) > max)
			{
				max = fabs(matrix[i][k]);
				ind = i;
			}
		}
		{
			double* temp = matrix[k];
			matrix[k] = matrix[ind];
			matrix[ind] = temp;
		}
		{
			double temp = b[k];
			b[k] = b[ind];
			b[ind] = temp;
		}
		for (int i = k; i < n; i++) {
			double temp = matrix[i][k];
			if (fabs(temp) > eps_gauss) {
				for (int j = 0; j < n; j++) {
					matrix[i][j] = matrix[i][j] / temp;
				}
				b[i] = b[i] / temp;
				if (i != k) {
					for (int j = 0; j < n; j++) {
						matrix[i][j] = matrix[i][j] - matrix[k][j];
					}
					b[i] = b[i] - b[k];
				}
			}
		}
	}
	if (fabs(b[n - 1]) > eps_gauss && matrix[n - 1][n - 1] <= eps_gauss)
	{
		return false;
	}
	for (int k = n - 1; k >= 0; --k)
	{
		x[k] = b[k];
		for (int i = 0; i < k; i++) {
			b[i] = b[i] - matrix[i][k] * x[k];
		}
	}
	return true;
}
