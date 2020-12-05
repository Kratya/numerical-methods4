#pragma once
#include "SLAE.h"
#include <string>
#include <cstdlib>
#include <vector>
#include <functional>
#include <algorithm>
using namespace std;

typedef function <double(vectorD)> func;
class SNU
{
private:
	SLAE* S;
	int n, m;
	int mode;
	vectorD rx, F, y, chx, chf;
	func* sys;
	function <void(vectorD*, vectorD)> Jacobian;
	double err, maxIt, betaMin, dx;
	vectorD* Jac;
	int* index;
	void updateF(vectorD x);
	void updateJAC(vectorD x);
	void makeSlau();
public:
	SNU(int m, int n, func* system, function <void(vectorD*, vectorD)> jacob, vectorD x, double err, double maxIt, double betaMin);
	SNU(int m, int n, func* system, vectorD x, double err, double maxIt, double betaMin, double dx);
	void calcSolv();
	~SNU();
};
