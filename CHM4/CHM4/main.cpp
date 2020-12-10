#include <iostream>
#include <fstream>
#include "SNU.h"
#include "tests.h"
using namespace std;

int main()
{
	setlocale(LC_ALL, "Russian");
	funcV funcs = initAllFunc();
	ifstream params;
	params.open("config.txt", ios_base::in);
	if (!params.is_open())
	{
		cout << "Cannot open file" << endl;
		return 1;
	}
	double eps, blim;
	int maxiter;
	params >> eps >> maxiter >> blim;
	vector<double> x;
	x.resize(g_n);
	for (int i = 0; i < g_n; ++i)
	{
		params >> x[i];
	}
	SNU* S = new SNU(g_m, g_n, funcs, jacob, x, eps, maxiter, blim);
	//SNU* S = new SNU(g_m, g_n, funcs, x, eps, maxiter, blim, 0.1);
	S->calcSolv();
	delete S;
}
