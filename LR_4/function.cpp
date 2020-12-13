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
	//В каждой ячейке массива l содержится коэф. перед c, начиная c c^0
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
				stream << l[i] << "*c^(" << i << ")";
			else
				stream << "+" << l[i] << "*c^(" << i << ")";
	stream << "\n";

	ErrorRate(l, size, part);
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

	int size = part.size();

	//Массив кубических полиномов 
	double** splAr = new double*[size - 1];
	for (int i = 0; i < size - 1; i++)
		splAr[i] = new double[4];

	//a_i=y_i-1
	for (int i = 0; i < size - 1; i++)
		splAr[i][0] = fun(part[i]);

	std::cout << "part:\n";
	for (const auto& x : part)
		std::cout << x << " ";
	std::cout << "\n";

	std::cout << "a:\n";
	for (int i = 0; i < size - 1; i++)
		std::cout << splAr[i][0] << " ";
	std::cout << "\n";

	//c_i
	//Приведение к виду c = Cx + y
	//C=-(A-E)
	std::vector<double> topbotAr, midAr;
	std::vector<double> bAr;
	std::vector<double> h;

	double gPrev, gNext;

	h.push_back(part[1] - part[0]);
	gPrev = (fun(part[1]) - fun(part[0])) / h[0];

	splAr[0][1] = gPrev;		//b_i
	//std::cout << h.back() << "\n";
	//std::cout << gPrev << "\n";
	for (int i = 1; i < size - 1; i++)
	{
		h.push_back(part[i + 1] - part[i]);
		gNext = (fun(part[i + 1]) - fun(part[i])) / h.back();

		splAr[i][1] = gNext;			//b_i

		midAr.push_back(2 * (*(&h.back() - 1) + h.back()));
		topbotAr.push_back(h.back());

		bAr.push_back(3 * (gNext - gPrev));

		gPrev = gNext;
	}

	std::cout << "h:\n";
	for (const auto& x : h)
		std::cout << x << " ";
	std::cout << "\n";

	std::cout << "midAr:\n";
	for (const auto& x : midAr)
		std::cout << x << " ";
	std::cout << "\n";

	std::cout << "topbotAr:\n";
	for (const auto& x : topbotAr)
		std::cout << x << " ";
	std::cout << "\n";

	std::cout << "bAr:\n";
	for (const auto& x : bAr)
		std::cout << x << " ";
	std::cout << "\n";


	double* c = new double[size];

	//TODO: точные значения
	c[0] = 0.;
	c[size - 1] = 0.;

	double* c_K = new double[size];			//Начальное приближение

	for (int i = 0; i < size; i++)
		c_K[i] = 0.;

	double* dif = new double[size];

	do
	{
		c[1] = (-midAr[0] + 1.)*c_K[1] - topbotAr[0] * c_K[2];
		c[1] += bAr[0];

		for (int row = 2; row < size - 2; row++)
		{
			c[row] = -topbotAr[row - 2] * c_K[row - 1] + (-midAr[row - 1] + 1)*c_K[row] - topbotAr[row - 1] * c_K[row + 1];
			c[row] += bAr[row];
		}

		c[size - 2] = -topbotAr[size - 4] * c_K[size - 3] + (-midAr[size - 3] + 1.)*c_K[size - 2];
		c[size - 2] += bAr[size - 3];

		specf::Difference(c, c_K, dif, size);

		std::cout << "c:\n";
		for (int i = 0; i < size; i++)
			std::cout << c[i] << " ";
		std::cout << "\n";

		for (int i = 1; i < size - 1; i++)
			c_K[i] = c[i];

	} while (specf::NormInf(dif, size) > 1.e-5);


	for (int i = 0; i < size - 1; i++)
	{

		//b_i
		splAr[i][1] += (h[i] * (2 * c[i + 1] + c[i]) / 3);
		//c_i
		splAr[i][2] = c[i];
		//d_i
		splAr[i][3] = (c[i + 1] - c[i]) / (3 * h[i]);
	}

	for (int i = 0; i < size - 1; i++)
	{
		std::cout << "s_" << i << " = " << splAr[i][0] << " + "
			<< splAr[i][1] << "(x-x_" << i << ") + "
			<< splAr[i][2] << "(x-x_" << i << ")^2 + "
			<< splAr[i][3] << "(x-x_0)^3\n";
	}

	for (int i = 0; i < size - 1; i++)
	{
		stream << splAr[i][0] << " " << splAr[i][1] << " " << splAr[i][2] << " " << splAr[i][3] << "\n";
	}



	//----------------
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed = end - start;
	stream << "Elapsed time: " << elapsed.count() << "s\n";
	//----------------

	SRenderGridToFile(splAr, part);

	delete[] c;
	delete[] c_K;
	for (int i = 0; i < size - 1; i++)
	{
		delete[](splAr[i]);
	}
	delete[] splAr;
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
	stream << "Error rate: " << max << "\n";
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

void InterpFunc::LRenderGridToFile(double* coef, int s)
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

void InterpFunc::SRenderGridToFile(double** coef, std::vector<double>& part)
{
	//Разбиение

	StreamOpen("sRender_Grid/");

	for (int i = 0; i < part.size() - 1; i++)
	{
		const int nods = 4;
		double h = (part[i + 1] - part[i]) / (double)(nods - 1);
		double temp;

		for (int j = 0; j < nods; j++)
		{
			stream << part[i] + j * h << " ";
			temp = coef[i][0] + coef[i][1] * (j*h) + coef[i][2] * pow(j * h, 2) + coef[i][3] * pow(j * h, 3);
			stream << temp << "\n";

		}
	}
}

InterpFunc::~InterpFunc()
{
	if (stream.is_open())
		stream.close();
}