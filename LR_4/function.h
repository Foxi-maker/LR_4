#pragma once
#include "Head.h"

class InterpFunc
{
	std::function<double(double)> fun;

	std::vector<double> unifPart;
	std::vector<double> fuP;
	std::vector<double> chebPart;
	std::vector<double> fcP;

	int left, right;

	std::string streamName;
	static std::fstream stream;

	void UniformPartition(int,int);
	void ChebyshevPartition(int);

	double LagrangePoint(std::vector<double>&, std::vector<double>&, double);

	void LErrorRate(std::vector<double>&, std::vector<double>&);

	void LRenderGridToFile(std::vector<double>&, std::vector<double>&);
	void SRenderGridToFile(double**, std::vector<double>&);

	void LOrderOfConvergence(std::vector<double>&, std::vector<double>&);
	//void SOrderOfConvergence(double**, std::vector<double>&);
public:
	InterpFunc();
	//(”казатель на фукнцию, интервал а,b, название функции)
	InterpFunc(std::function<double(double)>,int,int,std::string);

	void NumNods(int n,int);

	void StreamOpen(std::string, std::string);
	void StreamAssociate(std::string);
	void OutToFile(const std::string) const;

	void LagrangePol(const char ch);
	void SplineInterpolation(const char ch);

	~InterpFunc();
};

