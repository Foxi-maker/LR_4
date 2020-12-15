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

	//очистка файла результатов
	if (stream.is_open())
		stream.close();
	stream.open("Results/" + streamName + ".txt", std::ios_base::out);
	if (!stream.is_open())
	{
		std::cout << "File " << streamName << " was not opened for writing!";
		exit(1);
	}
}

void InterpFunc::NumNods(int n, int wp = 0)
{
	StreamAssociate("Results/");
	stream << "-----------------------\n"
		<< "Number of nodes: " << n << "\n";
	unifPart.clear();
	fuP.clear();
	chebPart.clear();
	fcP.clear();
	UniformPartition(n,wp);
	ChebyshevPartition(n,wp);
}

void InterpFunc::UniformPartition(int numNod, int wp=0)
{
	double q = 0.5;
	if (!wp)
		stream << "q: " << q << "\n";
	double h = pow(q,wp)*(right - left) / (double)(numNod - 1);

	for (int i = 0; (left + i * h )< right; i++)
	{
		unifPart.push_back(left + i * h);
		//std::cout << unifPart[i] << "\n";
	}
	std::cout << "\n";
	for (const auto& u : unifPart)
	{
		fuP.push_back(fun(u));
		//std::cout << fuP.back() << "\n";
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
	for (const auto& c : chebPart)
		fcP.push_back(fun(c));
}

double InterpFunc::LagrangePoint(std::vector<double>& p, std::vector<double>& fp, double x)
{
	long double num, den;
	long double res = 0.;

	for (int index = 0; index < p.size(); index++)
	{
		num = 1.;
		den = 1.;

		for (int j = 0; j < index; j++)
		{
			num *= x - p[j];
			den *= p[index] - p[j];
		}
		for (int j = index + 1; j < p.size(); j++)
		{
			num *= x - p[j];
			den *= p[index] - p[j];
		}
		//std::cout << "num " << num << "\n";
		//std::cout << "den" << den << "\n";
		//TODO: переполнение
		res += fp[index] * num / den;
		//std::cout << "fp[index] " << fp[index] << "\n";
		//std::cout << "res " << res << "\n";
	}
	return res;
}

void InterpFunc::LagrangePol(const char ch)
{
	//StreamOpen("Results/");

	std::vector<double> part, fpart;
	switch (ch)
	{
	case 'C':					//C- Chebyshev Partition
		part = chebPart;
		fpart = fcP;
		stream << "Chebyshev Partition\n";
		break;
	case 'U':					//U- Uniform Partition
		part = unifPart;
		fpart = fuP;
		stream << "Uniform Partition\n";
		break;
	default:
		std::cout << "Split selection error!\n";
		break;
	}

	//LErrorRate(part, fpart);

	//LRenderGridToFile(part, fpart);
}

void InterpFunc::SplineInterpolation(const char ch)
{
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

	//c_i
	//Приведение к виду c = Cx + y
	//C=-(A-E)
	std::vector<double> midAr;
	std::vector<double> bAr;
	std::vector<double> h;

	double gPrev, gNext;

	h.push_back(part[1] - part[0]);
	gPrev = (fun(part[1]) - fun(part[0])) / h[0];

	splAr[0][1] = gPrev;		//b_i

	for (int i = 1; i < size - 1; i++)
	{
		h.push_back(part[i + 1] - part[i]);
		gNext = (fun(part[i + 1]) - fun(part[i])) / h.back();

		splAr[i][1] = gNext;			//b_i
		midAr.push_back(2 * (*(&h.back() - 1) + h.back()));
		bAr.push_back(3 * (gNext - gPrev));

		gPrev = gNext;
	}
	double* c = new double[size];
	double* c_K = new double[size];			//Начальное приближение

	c[0] = specf::SecondDerivative(fun, left);
	c[size - 1] = specf::SecondDerivative(fun, right);
	c_K[0] = c[0];
	c_K[size - 1] = c[size - 1];

	for (int i = 1; i < size - 1; i++)
		c_K[i] = bAr[i - 1];

	double* dif = new double[size];

	do
	{
		for (int row = 1; row < size - 1; row++)
			c[row] = (-h[row - 1] * c_K[row - 1] - h[row] * c_K[row + 1] + bAr[row - 1]) / midAr[row - 1];

		specf::Difference(c, c_K, dif, size);

		for (int i = 1; i < size - 1; i++)
			c_K[i] = c[i];

	} while (specf::NormInf(dif, size) > 1.e-5);


	for (int i = 0; i < size - 1; i++)
	{

		//b_i
		splAr[i][1] -= ((h[i] * (c[i + 1] + 2 * c[i])) / 3);
		//c_i
		splAr[i][2] = c[i];
		//d_i
		splAr[i][3] = (c[i + 1] - c[i]) / (3 * h[i]);
	}

	SRenderGridToFile(splAr, part);

	delete[] c;
	delete[] c_K;
	for (int i = 0; i < size - 1; i++)
	{
		delete[](splAr[i]);
	}
	delete[] splAr;
}

void InterpFunc::LErrorRate(std::vector<double>& p, std::vector<double>& fp)
{
	int size = p.size();

	//Разбиение
	const int nods = 256;

	double h = (right - left) / (double)(nods - 1); 

	double max = 0., point, dif;
	for (int i = 1; i < nods; i++)
	{
		point = left + i * h;
		dif = fabs(fun(point) - LagrangePoint(p, fp, point));
		if (dif > max)
			max = dif;
	}
	
	stream << "Error rate: " << max << "\n";
}

void InterpFunc::StreamOpen(std::string file, std::string fn = "")
{
	if (stream.is_open())
		stream.close();
	stream.open(file + streamName + "_" + fn + "nods.txt", std::ios_base::out);
	if (!stream.is_open())
	{
		std::cout << "File " << streamName << " was not opened for writing!";
		exit(1);
	}
}

void InterpFunc::StreamAssociate(std::string file)
{
	if (stream.is_open())
		stream.close();
	stream.open(file + streamName + ".txt", std::ios_base::app);
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

void InterpFunc::LRenderGridToFile(std::vector<double>& part, std::vector<double>& fpart)
{
	int size = part.size();
	//Разбиение
	const int nods = 256 / size;

	StreamOpen("Render_Grid/", std::to_string(size));

	double h;
	double point;
	for (int i = 0; i < size - 1; i++)
	{
		h = (part[i + 1] - part[i]) / (double)(nods - 1);
		for (int j = 0; j < nods; j++)
		{
			point = part[i] + j * h;

			stream << point << " ";
			stream << LagrangePoint(part, fpart, point) << "\n";
		}
	}
}

void InterpFunc::SRenderGridToFile(double** coef, std::vector<double>& part)
{
	StreamOpen("sRender_Grid/", std::to_string(part.size()));


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

void InterpFunc::LOrderOfConvergence(std::vector<double>& p, std::vector<double>& fp)
{
	int size = p.size();

	const int nods = 256;

	double h = (right - left) / (double)(nods - 1);

	double max = 0., point, dif;
	double prevErr;
	for (int i = 0; i < nods; i++)
	{
		point = left + i * h;
		dif = fabs(fun(point) - LagrangePoint(p, fp, point));
		if (dif > max)
			max = dif;
	}

	stream << "Error rate: " << max << "\n";
}
//
//void InterpFunc::SOrderOfConvergence(double **, std::vector<double>&)
//{
//	double q = 0.5;
//}

InterpFunc::~InterpFunc()
{
	if (stream.is_open())
		stream.close();
}