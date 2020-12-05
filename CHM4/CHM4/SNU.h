#pragma once
#include "SLAE.h"
#include <string>
#include <cstdlib>
#include <vector>
#include <functional>
#include <algorithm>
using namespace std;

typedef function <double(vector<double> &xt)> func;

class SNU
{
private:
	SLAE* S;
	int n, m;
	int mode;
	vector<double> rx, F, y, chx, chf;
	func* sys;
	function <void(vector<vector<double>>&, vector<double>&)> Jacobian;
	double err, maxIt, betaMin, dx;
	vector<vector<double>> Jac;
	int* index;
	void updateF(vector<double> &x);
	void updateJAC(vector<double> &x);
	void makeSlau();
public:
	SNU(int m, int n, func* system, function <void(vector<vector<double>>&, vector<double>&)> jacob,
		vector<double> &x, double err, double maxIt, double betaMin);
	SNU(int m, int n, func* system, vector<double> &x, double err, double maxIt, double betaMin, double dx);
	void calcSolv();
	~SNU();
};
