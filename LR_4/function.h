#pragma once
#include "Head.h"

class InterpFunc
{
	std::function<double(double)> fun;

	std::vector<double> unifPart;
	std::vector<double> chebPart;

	int left, right;

	std::string streamName;
	static std::fstream stream;
public:
	InterpFunc();
	//(”казатель на фукнцию, интервал а,b, название функции)
	InterpFunc(std::function<double(double)>,int,int,std::string);

	void NumNods(int n);

	//TODO: перенести в приват
	void UniformPartition(int);
	void ChebyshevPartition(int);

	void LagrangePol(const char ch);
	void SplineInterpolation(const char ch);

	void ErrorRate(double*,int,std::vector<double>&);

	void StreamOpen(std::string);
	void OutToFile(const std::string) const;

	void LRenderGridToFile(double*,int);
	void SRenderGridToFile(double**, std::vector<double>&);

	~InterpFunc();
};

