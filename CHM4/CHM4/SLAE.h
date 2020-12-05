#pragma once
#include <cmath>
typedef double* vectorD;

class SLAE
{
public:
	int n;
	vectorD* A;
	vectorD x;
	vectorD b;
	SLAE(int size);
	~SLAE();
};
double Scal(vectorD a, vectorD b, int n);
double Norm(vectorD a, int n);
bool GausSolver(SLAE* S);
