#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <math.h>
#include <fstream>
#include <string>

namespace specf
{
	double fa�torial(double x);

	double SecondDerivative(std::function<double(double)>, double);

	//������ �������� - ����� ��������� ��� ������������ ��� ����������
	//����� ���� ��� ������������ ������� �� ����������� �������
	bool next_combination(std::vector<int>&, int);

	void Difference(double*, double*, double*, int);

	double NormInf(double*, int);
}
//�������� �������
namespace testf
{
	double fun0(double x);
	double fun1(double x);
	double fun2(double x);
	double fun3(double x);
	double fun4(double x);
	double fun5(double x);
	double fun6(double x);
	double fun7(double x);
	double fun8(double x);
	double fun9(double x);
}
