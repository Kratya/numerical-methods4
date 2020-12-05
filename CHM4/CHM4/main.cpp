#include "SNU.h"
#include "tests.h"

int main()
{
	func* funcs = initAllFunc();

	FILE* params;
	fopen_s(&params, "config.txt", "r");
	if (params == 0)
	{
		return 1;
	}
	double eps, blim;
	int maxiter;
	fscanf_s(params, "%lf %d %lf", &eps, &maxiter, &blim);
	vector<double> x;
	x.resize(g_n);
	for (int i = 0; i < g_n; ++i)
	{
		fscanf_s(params, "%lf", &x[i]);
	}
	SNU* S = new SNU(g_m, g_n, funcs, jacob, x, eps, maxiter, blim);
	//SNU* S = new SNU(g_m, g_n, funcs, x, eps, maxiter, blim, 0.1);
	S->calcSolv();
	delete S;
}
