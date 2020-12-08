#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <string>
#include <chrono>

namespace specf
{
	double faсtorial(double x);

	//Второй аргемент - число элементов для перестановки без повторений
	//Число мест для перестановки берется из размерности вектора
	bool next_combination(std::vector<int>&, int);
}
//Тестовые функции
namespace testf
{
	double fun1(double x);
	double fun2(double x);
	double fun3(double x);
	double fun4(double x);
	double fun5(double x);
	double fun6(double x);

}
