#include "slae.h"

using namespace std;

SLAE::SLAE(int size)
{
	n = size;
	A.resize(n);
	for (int i = 0; i < n; ++i)
	{
		A[i].resize(n);
	}
	b.resize(n);
	x.resize(n);
}

SLAE::~SLAE()
{
}

double Scal(vector<double> &a, vector<double> &b, int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++) 
	{
		sum += a[i] * b[i];
	}
	return sum;
}

double Norm(vector<double> &a, int n)
{
	return sqrt(Scal(a, a, n));
}

bool GausSolver(SLAE* S)
{
	int n = S->n;
	vector<vector<double>> matrix = S->A;
	vector<double> b = S->b;
	vector<double> x = S->x;
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
			vector<double> temp = matrix[k];
			matrix[k] = matrix[ind];
			matrix[ind] = temp;
		}
		{
			double temp = b[k];
			b[k] = b[ind];
			b[ind] = temp;
		}
		for (int i = k; i < n; i++) 
		{
			double temp = matrix[i][k];
			if (fabs(temp) > eps_gauss) 
			{
				for (int j = 0; j < n; j++)
				{
					matrix[i][j] = matrix[i][j] / temp;
				}
				b[i] = b[i] / temp;
				if (i != k) {
					for (int j = 0; j < n; j++)
					{
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
		for (int i = 0; i < k; i++)
		{
			b[i] = b[i] - matrix[i][k] * x[k];
		}
	}
	S->A = matrix;
	S->b = b;
	S->x = x;
	return true;
}
