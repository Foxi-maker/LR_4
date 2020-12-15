#include "function.h"

std::fstream InterpFunc::stream;

int main()
{
	std::vector<int> numNod = { 4, 8, 16, 32, 64, 128 };
	//std::vector<int> numNod = { 4, 8, 16, 32 };

	std::vector<InterpFunc> elements;
	//elements.emplace_back(testf::fun0, -1, 1, std::string("fun0"));
	//elements.emplace_back(testf::fun1, -1, 1, std::string("fun1"));
	//elements.emplace_back(testf::fun2, -1, 1, std::string("fun2"));
	//elements.emplace_back(testf::fun3, -3, 3, std::string("fun3"));
	//elements.emplace_back(testf::fun4, -1, 1, std::string("fun4"));
	//elements.emplace_back(testf::fun5, -1, 1, std::string("fun5"));
	//elements.emplace_back(testf::fun6, 0, 2, std::string("fun6"));
	//elements.emplace_back(testf::fun7, -1, 1, std::string("fun7"));
	//elements.emplace_back(testf::fun8, -1, 1, std::string("fun8"));
	elements.emplace_back(testf::fun9, -1, 1, std::string("fun9"));

	for (auto& e : elements)
	{
		for (const auto& n : numNod)
		{
			for(int i=0;i<5;i++)
			e.NumNods(n);
			//e.OutToFile("Lagrangian polynomial: ");
			e.LagrangePol('U');
			//e.LagrangePol('C');
			//e.OutToFile("Spline interpolation: ");
			//e.SplineInterpolation('U');
			//e.SplineInterpolation('C');
		}
	}
}
