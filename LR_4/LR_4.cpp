#include "function.h"

std::fstream InterpFunc::stream;

int main()
{
	//std::vector<std::function<double(double)>> vecFun = { testf::fun1,testf::fun2,testf::fun3,testf::fun4,testf::fun5 };
	//std::vector<std::function<double(double)>> vecFun = { testf::fun1 };

	//std::vector<int> numNod = { 4, 8, 16, 32, 64, 128 };
	std::vector<int> numNod = { 4, 8, 16, 32, 64, 128 };

	std::vector<InterpFunc> elements;
	//elements.emplace_back(testf::fun1, -1, 1, std::string("fun1"));
	//elements.emplace_back(testf::fun2, -1, 1, std::string("fun2"));
	//elements.emplace_back(testf::fun3, -3, 3, std::string("fun3"));
	//elements.emplace_back(testf::fun4, -1, 1, std::string("fun4"));
	//elements.emplace_back(testf::fun5, -1, 1, std::string("fun5"));
	//elements.emplace_back(testf::fun6, 0, 2, std::string("fun6"));
	elements.emplace_back(testf::fun6, -1, 1, std::string("fun7"));

	for (auto& e : elements)
	{
		e.StreamOpen("Results/");
		for (const auto& n : numNod)
		{
			e.NumNods(n);
			//e.OutToFile("Lagrangian polynomial: ");
			//e.LagrangePol('U');
			//e.LagrangePol('C');
			//e.OutToFile("Spline interpolation: ");
			e.SplineInterpolation('U');
			//e.SplineInterpolation('C');
		}
	}
}
