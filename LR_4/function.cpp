#include "function.h"

InterpFunc::InterpFunc()
{
	fun = nullptr;
}

InterpFunc::InterpFunc(std::function<double(double)> f, int a, int b, std::string name)
{
	fun = f;

	left = a;
	right = b;

	streamName = name;

	//UniformPartition(n);
	//ChebyshevPartition(n);
}

void InterpFunc::NumNods(int n)
{
	stream << "-----------------------\n"
		<< "Number of nodes: " << n << "\n";
	unifPart.clear();
	chebPart.clear();
	UniformPartition(n);
	ChebyshevPartition(n);
}

void InterpFunc::UniformPartition(int numNod)
{
	double h = (right - left) / (double)(numNod - 1);
	for (int i = 0; i < numNod; i++)
	{
		unifPart.push_back(left + i * h);
		//std::cout << unifPart[i] << "\n";
	}
}

void InterpFunc::ChebyshevPartition(int numNod)
{
	double h1 = (left + right) / 2;
	double h2 = (right - left) / 2;

	for (int i = 0; i < numNod; i++)
	{
		chebPart.push_back(h1 + h2 * cos(((2 * i + 1)*M_PI) / (2 * (numNod))));
		//std::cout << chebPart[i] << "\n";
	}
}

void InterpFunc::LagrangePol(const char ch)
{
	//---------------
	auto start = std::chrono::system_clock::now();
	//--------------

	std::vector<double> part;
	switch (ch)
	{
	case 'C':					//C- Chebyshev Partition
		part = chebPart;
		stream << "Chebyshev Partition\n";
		break;
	case 'U':					//U- Uniform Partition
		part = unifPart;
		stream << "Uniform Partition\n";
		break;
	default:
		std::cout << "Split selection error!\n";
		break;
	}

	int size = part.size();
	//В каждой ячейке массива l содержится коэф. перед x, начиная c x^0
	double* l = new double[size];
	for (int i = 0; i < size; i++)
		l[i] = 0.;

	std::vector<int> lIndex;
	for (int i = 0; i < (size - 1); i++)
		lIndex.push_back(i);

	double denominator;
	int indexDen = size;
	do
	{
		denominator = 1.;
		indexDen--;
		for (const auto& li : lIndex)
			denominator *= part[indexDen] - part[li];

		denominator = fun(part[indexDen]) / denominator;

		std::vector<double> tPart;
		for (const auto& li : lIndex)
			tPart.push_back(part[li]);

		for (int brackIndex = 0; brackIndex < size - 1; brackIndex++)
		{
			//Вектор индексов 
			std::vector<int> perIndex;
			for (int i = 0; i < (size - 1 - brackIndex); i++)
				perIndex.push_back(i);

			double sum = 0.;
			do
			{
				double temp = 1.;
				for (const auto& v : perIndex)
					temp *= tPart[v];

				sum -= temp;

			} while (specf::next_combination(perIndex, size - 1));
			l[brackIndex] += sum * denominator;
		}
		l[size - 1] += denominator;
	} while (specf::next_combination(lIndex, size));

	//for (int i = 0; i < size; i++)
	//	std::cout << l[i] << " ";
	//std::cout << "\n ";

	if (!(l[0] < 1.e-10))
		stream << l[0];
	for (int i = 1; i < size - 1; i++)
		if (l[i] < 1.e-10)
			continue;
		else
			if (l[i] < 0)
				stream << l[i] << "*x^(" << i << ")";
			else
				stream << "+" << l[i] << "*x^(" << i << ")";
	stream << "\n";

	ErrorRate(l,size,part);
	//----------------
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed = end - start;
	stream << "Elapsed time: " << elapsed.count() << "s\n";
	//----------------

	//RenderGridToFile(l,size);
	delete[] l;
}

void InterpFunc::SplineInterpolation(const char ch)
{
	//---------------
	auto start = std::chrono::system_clock::now();
	//--------------

	std::vector<double> part;
	switch (ch)
	{
	case 'C':					//C- Chebyshev Partition
		part = chebPart;
		stream << "Chebyshev Partition\n";
		break;
	case 'U':					//U- Uniform Partition
		part = unifPart;
		stream << "Uniform Partition\n";
		break;
	default:
		std::cout << "Split selection error!\n";
		break;
	}

	
	//----------------
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed = end - start;
	stream << "Elapsed time: " << elapsed.count() << "s\n";
	//----------------

}

void InterpFunc::ErrorRate(double* coef, int s, std::vector<double>& part)
{
	double max = 0., temp, dif;
	for (const auto& p : part)
	{
		temp = 0.;
		for (int j = 0; j < s; j++)
			temp += coef[j] * pow(p, j);
		dif = fabs(fun(p) - temp);
		if (dif > max)
			max = dif;
	}
	stream << "Error rate:" << max << "\n";
}

void InterpFunc::StreamOpen(std::string file)
{
	if (stream.is_open())
		stream.close();
	stream.open(file + streamName + ".txt", std::ios_base::out);
	if (!stream.is_open())
	{
		std::cout << "File " << streamName << " was not opened for writing!";
		exit(1);
	}
}

void InterpFunc::OutToFile(const std::string s) const
{
	stream << s << "\n";
}

void InterpFunc::RenderGridToFile(double* coef, int s)
{
	//Разбиение
	const int nods = 50;
	StreamOpen("Render_Grid/");

	double h = (right - left) / (double)(nods - 1);
	double temp;
	for (int i = 0; i < nods; i++)
	{
		temp = 0.;
		stream << left + i * h << " ";
		for (int j = 0; j < s; j++)
			temp += coef[j] * pow(left + i * h, j);
		stream << temp << "\n";
	}

}

InterpFunc::~InterpFunc()
{
	if (stream.is_open())
		stream.close();
}