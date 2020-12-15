#include "Head.h"

double specf::fañtorial(double x)
{
	if (x == 0.)
		return 1.;
	return x * fañtorial(x - 1.);
}

double specf::SecondDerivative(std::function<double(double)> f, double point)
{
	double h = 1.;
	return ((f(point + h) - 2 * f(point) + f(point - h)) / (h*h));
}

bool specf::next_combination(std::vector<int>& vec, int n)
{
	int m = vec.size();
	for (int i = m - 1; i >= 0; --i)
		if (vec[i] < n - m + i)
		{
			++vec[i];
			for (int j = i + 1; j < m; ++j)
				vec[j] = vec[j - 1] + 1;
			return true;
		}
	return false;
}

void specf::Difference(double* x, double* y, double* z,int n)
{
	for (int i = 0; i < n; i++)
		z[i] = x[i] - y[i];
}

double specf::NormInf(double* x, int n)
{
	int indexMax = 0;
	for (int index = 0; index < n; index++)
	{
		if (abs(x[indexMax]) < abs(x[index]))
		{
			indexMax = index;
		}
	}
	return abs(x[indexMax]);
}

//Òåñòîâûå ôóíêöèè
double testf::fun0(double x)
{
	return  x;
}

double testf::fun1(double x)
{
	return x * x;
}

double testf::fun2(double x)
{
	return 1 / (1 + x * x);
}

double testf::fun3(double x)
{
	return 1 / atan(1 + 10 * x * x);
}

double testf::fun4(double x)
{
	double R1 = pow(4 * pow(x, 3) + 2 * x*x - 4 * x + 2, sqrt(2));
	double R2 = asin(1 / (5 + x - x * x));
	return R1 + R2 - 5;
}

double testf::fun5(double x)
{
	double R1 = sin((pow(x, 4) + pow(x, 3) - 3 * x + 3 - pow(30, 1 / 3)) / 2);
	double R2 = tanh((4 * sqrt(3)*pow(x, 3) - 2 * x - 6 * sqrt(2) + 1) / (-2 * sqrt(3)*pow(x, 3) + x + 3 * sqrt(2)));
	return R1 + R2 + 1.2;
}

double testf::fun6(double x)
{
	return exp(-x);
}

double testf::fun7(double x)
{
	return 1;
}

double testf::fun8(double x)
{
	return 1 / (1 + 25 * x * x);
}

double testf::fun9(double x)
{
	return sin(M_PI*x);
}