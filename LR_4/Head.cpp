#include "Head.h"

double specf::fa�torial(double x)
{
	if (x == 0.)
		return 1.;
	return x * fa�torial(x - 1.);
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

//�������� �������
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