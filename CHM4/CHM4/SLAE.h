#pragma once
#include <cmath>
#include <vector>

using namespace std;

class SLAE
{
public:
	int n;
	vector<vector<double>> A;
	vector<double> x;
	vector<double> b;
	SLAE(int size);
	~SLAE();
};
double Scal(vector<double> &a, vector<double> &b, int n);
double Norm(vector<double> &a, int n);
bool GausSolver(SLAE* S);
