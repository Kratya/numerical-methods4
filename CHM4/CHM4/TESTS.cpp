#include "tests.h"

#ifdef test1
int g_m = 2; int g_n = 2;
double f1(vector<double> &x)
{
	return x[0] * x[0] + x[1] * x[1] - 1;
}
double f2(vector<double> &x)
{
	return (x[0] - 4) * (x[0] - 4) + x[1] * x[1] - 4;
}
void jacob(vector<vector<double>> &jac, vector<double>&x)
{
	jac[0][0] = 2 * x[0]; jac[0][1] = 2 * x[1];
	jac[1][0] = 2 * x[0] - 8; jac[1][1] = 2 * x[1];
}
funcV initAllFunc()
{
	funcV funcs(g_m);
	funcs[0] = f1;
	funcs[1] = f2;
	return funcs;
}
#endif

#ifdef test2
int g_m = 2; int g_n = 2;
double f1(vector<double> &x)
{
	return x[0] * x[0] + x[1] * x[1] - 1;
}
double f2(vector<double> &x)
{
	return (x[0] - 2) * (x[0] - 2) + x[1] * x[1] - 1;
}
void jacob(vector<vector<double>>& jac, vector<double>& x)
{
	jac[0][0] = 2 * x[0]; jac[0][1] = 2 * x[1];
	jac[1][0] = 2 * x[0] - 4; jac[1][1] = 2 * x[1];
}
funcV initAllFunc()
{
	funcV funcs(g_m);
	funcs[0] = f1;
	funcs[1] = f2;
	return funcs;
}
#endif

#ifdef test3
int g_m = 2; int g_n = 2;
double f1(vector<double> &x)
{
	return x[0] * x[0] + x[1] * x[1] - 1;
}
double f2(vector<double> &x)
{
	return (x[0] - 1) * (x[0] - 1) + x[1] * x[1] - 2;
}
void jacob(vector<vector<double>> &jac, vector<double> &x)
{
	jac[0][0] = 2 * x[0]; jac[0][1] = 2 * x[1];
	jac[1][0] = 2 * x[0] - 2; jac[1][1] = 2 * x[1];
}
funcV initAllFunc()
{
	funcV funcs(g_m);
	funcs[0] = f1;
	funcs[1] = f2;
	return funcs;
}
#endif

#ifdef test4
int g_m = 3; int g_n = 2;
double f1(vector<double> &x)
{
	return x[0] * x[0] + x[1] * x[1] - 1;
}
double f2(vector<double> &x)
{
	return (x[0] - 1) * (x[0] - 1) + x[1] * x[1] - 2;
}
double f3(vector<double> &x)
{
	return x[0] / 2 + 2 - x[1];
}
void jacob(vector<vector<double>> &jac, vector<double> &x)
{
	jac[0][0] = 2 * x[0]; jac[0][1] = 2 * x[1];
	jac[1][0] = 2 * x[0] - 2; jac[1][1] = 2 * x[1];
	jac[2][0] = 1. / 2; jac[2][1] = -1;
}
funcV initAllFunc()
{
	funcV funcs(g_m);
	funcs[0] = f1;
	funcs[1] = f2;
	funcs[2] = f3;
	return funcs;
}
#endif

#ifdef test5
int g_m = 3; int g_n = 2;
double f1(vector<double> &x)
{
	return -x[0] - x[1];
}
double f2(vector<double> &x)
{
	return x[0] + 1 - x[1];
}
double f3(vector<double> &x)
{
	return x[0] / 2 + 2 - x[1];
}
void jacob(vector<vector<double>> &jac, vector<double> &x)
{
	jac[0][0] = -1; jac[0][1] = -1;
	jac[1][0] = 1; jac[1][1] = -1;
	jac[2][0] = 1. / 2; jac[2][1] = -1;
}
funcV initAllFunc()
{
	funcV funcs(g_m);
	funcs[0] = f1;
	funcs[1] = f2;
	funcs[2] = f3;
	return funcs;
}
#endif

#ifdef test6
int g_m = 3; int g_n = 2;
double f1(vector<double> &x)
{
	return -x[0] - x[1];
}
double f2(vector<double> &x)
{
	return x[0] + 1 - x[1];
}
double f3(vector<double> &x)
{
	return 2 * x[0] / 2 + 4 - 2 * x[1];
}
void jacob(vector<vector<double>> &jac, vector<double> x)
{
	jac[0][0] = -1; jac[0][1] = -1;
	jac[1][0] = 1; jac[1][1] = -1;
	jac[2][0] = 2. / 2; jac[2][1] = -2;
}
funcV initAllFunc()
{
	funcV funcs(g_m);
	funcs[0] = f1;
	funcs[1] = f2;
	funcs[2] = f3;
	return funcs;
}
#endif

#ifdef test7
int g_m = 3; int g_n = 2;
double f1(vector<double> &x)
{
	return -x[0] - x[1];
}
double f2(vector<double> &x)
{
	return x[0] + 1 - x[1];
}
double f3(vector<double> &x)
{
	return 20 * x[0] / 2 + 40 - 20 * x[1];
}
void jacob(vector<vector<double>> &jac, vector<double> &x)
{
	jac[0][0] = -1; jac[0][1] = -1;
	jac[1][0] = 1; jac[1][1] = -1;
	jac[2][0] = 20. / 2; jac[2][1] = -20;
}
funcV initAllFunc()
{
	funcV funcs(g_m);
	funcs[0] = f1;
	funcs[1] = f2;
	funcs[2] = f3;
	return funcs;
}
#endif

#ifdef test8
int g_m = 2; int g_n = 2;
double f1(vector<double> &x)
{
	return sin(x[0]) - x[1];
}
double f2(vector<double> &x)
{
	return x[0] / 4 - 0.5 - x[1];
}
void jacob(vector<vector<double>> &jac, vector<double> &x)
{
	jac[0][0] = cos(x[0]); jac[0][1] = -1;
	jac[1][0] = 1. / 4; jac[1][1] = -1;
}
funcV initAllFunc()
{
	funcV funcs(g_m);
	funcs[0] = f1;
	funcs[1] = f2;
	return funcs;
}
#endif


#ifdef test9
int g_m = 3; int g_n = 2;
double f1(vector<double> &x)
{
	return x[0] * x[0] - x[1];
}
double f2(vector<double> &x)
{
	return x[0] - x[1];
}
double f3(vector<double> &x)
{
	return -x[0] + 2 - x[1];
}
void jacob(vector<vector<double>> &jac, vector<double> &x)
{
	jac[0][0] = 2 * x[0]; jac[0][1] = -1;
	jac[1][0] = 1; jac[1][1] = -1;
	jac[2][0] = -1; jac[2][1] = -1;
}
funcV initAllFunc()
{
	funcV funcs(g_m);
	funcs[0] = f1;
	funcs[1] = f2;
	funcs[2] = f3;
	return funcs;
}
#endif

