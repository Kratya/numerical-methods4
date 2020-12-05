#include "SNU.h"

SNU::SNU(int lm, int ln, func* system, function <void(vectorD*, vectorD)> jacob, vectorD lx, double lerr, double lmaxIt, double lbetaMin)
{
	mode = 0;
	sys = system;
	Jacobian = jacob;
	m = lm;
	n = ln;
	err = lerr;
	maxIt = lmaxIt;
	betaMin = lbetaMin;
	S = new SLAE(n);
	rx = lx;
	Jac = new vectorD[m];
	for (int i = 0; i < m; ++i)
	{
		Jac[i] = new double[n];
	}
	index = new int[m];
	y = new double[m];
	F = new double[m];
}

SNU::SNU(int lm, int ln, func* system, vectorD lx, double lerr, double lmaxIt, double lbetaMin, double ldx)
{
	mode = 1;
	sys = system;
	Jacobian = {};
	m = lm;
	n = ln;
	err = lerr;
	maxIt = lmaxIt;
	betaMin = lbetaMin;
	dx = ldx;
	S = new SLAE(n);
	rx = lx;
	Jac = new vectorD[m];
	for (int i = 0; i < m; ++i)
	{
		Jac[i] = new double[n];
	}
	index = new int[m];
	y = new double[m];
	F = new double[m];
	chx = new double[m];
	chf = new double[m];
}

void SNU::updateF(vectorD x)
{
	for (int i = 0; i < m; ++i)
	{
		F[i] = -sys[i](x);
	}
}

void SNU::updateJAC(vectorD x)
{
	if (mode == 0)
	{
		Jacobian(Jac, x);
	}
	else
	{
		for (int i = 0; i < n; ++i) {
			for (int k = 0; k < n; ++k)
			{
				chx[k] = x[k];
			}
			chx[i] += dx;
			for (int k = 0; k < m; ++k)
			{
				chf[k] = -sys[k](chx);
			}
			for (int j = 0; j < m; ++j) {
				Jac[j][i] = (F[j] - chf[j]) / dx;
			}
		}
	}
	for (int i = 0; i < m; ++i)
	{
		index[i] = i;
	}
}

void SNU::makeSlau()
{
	double temp;
	vectorD* matrix = S->A;
	vectorD b = S->b;
	int minInd;
	int count = m - n + 1;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			matrix[i][j] = 0;
		}
	}
	sort(index, index + m, [&](const int a, const int b) {return fabs(F[a]) < fabs(F[b]); });
	minInd = index[0];
	for (int i = 0; i < count; ++i) {
		if (index[i] < minInd)
			minInd = index[i];
	}
	for (int i = 0; i < minInd; ++i) {
		b[i] = F[i];
		for (int im = 0; im < n; ++im)
			matrix[i][im] = Jac[i][im];
	}
	b[minInd] = 0;
	for (int i = 0; i < count; ++i) {
		b[minInd] += F[index[i]] * F[index[i]];
		for (int jm = 0; jm < n; ++jm)
		{
			matrix[minInd][jm] += 2 * F[index[i]] * Jac[index[i]][jm];
		}
	}
	int j = minInd + 1;
	for (int i = j; i < m; ++i) {
		bool flag = true;
		for (int k = 0; k < count; ++k)
		{
			if (i == index[k])
			{
				flag = false;
				break;
			}
		}
		if (flag)
		{
			b[j] = F[i];
			for (int jm = 0; jm < n; ++jm)
			{
				matrix[j][jm] = Jac[i][jm];
			}
			++j;
		}
	}
}

void SNU::calcSolv()
{
	double beta = 1;
	double norm, norm2;
	vectorD x = S->x;
	updateF(rx);
	updateJAC(rx);
	norm = Norm(F, n);
	int iter = 0;
	while (norm > err && iter < maxIt)
	{
		makeSlau();
		if (!GausSolver(S))
		{
			printf_s("SOLVE NOT FOUND\n");
			return;
		}
		beta = 2;
		do {
			beta /= 2;
			for (int i = 0; i < n; ++i)
			{
				y[i] = rx[i] + beta * x[i];
			}
			updateF(y);
			norm2 = Norm(F, n);
		} while (norm2 >= norm && beta > betaMin);
		for (int i = 0; i < n; ++i)
			rx[i] += beta * x[i];
		updateJAC(rx);
		makeSlau();
		norm = Norm(F, n);
		iter++;
		printf_s("%d:  norm = %.2le  b = %.2le \t", iter, norm, beta);
		for (int i = 0; i < n; ++i)
		{
			printf_s("%lf\t", rx[i]);
		}
		printf_s("\n");
		if (beta <= betaMin)
		{
			return;
		}
	}
}

SNU::~SNU()
{
	delete S;
	delete[] rx;
	delete[] F;
	delete[] index;
	for (int i = 0; i < m; ++i)
	{
		delete[] Jac[i];
	}
	delete[] chx;
	delete[] chf;
	delete[] Jac;
}
