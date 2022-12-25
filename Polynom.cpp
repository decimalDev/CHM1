#include "Polynom.h"
#include<stdio.h>
/*
double expression(double x) {
	return exp(x) / (1 + x * x);
	//return atan(x)/(1 + x * x);
}
*/

double vect_to_polynom(std::vector<double>& P, double x) //преобразует последовательность коэффициентов в функцию-многочлен, возвращает значение многочлена в данной точке x
{
	double Px = 0;
	int N = P.size();
	for (int i = 0; i < N; i++)
	{
		Px += P[i] * pow(x, i);
	}
	return Px;
}

// У нас есть класс Polynom с полями коэффициентов (вектор vector<double> Coeffs) и степени (int Deg); у него есть дефолтный конструктор (нулевой вектор), конструктор по вектору коэффициентов (в порядке от младших к старших) и конструктор копирования; он умеет возвращать степень многочлена и возвращать его значение в произвольной точке (double value(double x))

//ОПРЕДЕЛЕНИЕ МЕТОДОВ КЛАССА Polynom

Polynom::Polynom() //дефолтный конструктор -- создаёт нулевой многочлен
{
	this->Coeffs = { 0 };
	this->Deg = 0;
	this->Grid = {};
}

Polynom::Polynom(vector<double>& X) //на вход подаётся вектор коэффициентов (от младших к старшим), по нему возвращается многочлен
{
	this->Coeffs = X;
	this->Deg = X.size() - 1;
	Grid = {};
}

Polynom::Polynom(vector<double>& X, vector<double>& G)
{
	this->Coeffs = X;
	this->Deg = G.size() - 1;
	this->Grid = G;
}

Polynom::Polynom(const Polynom& P) //копирование
{
	this->Coeffs = P.Coeffs;
	this->Deg = P.Deg;
	this->Grid = P.Grid;
}

Polynom& Polynom::operator=(const Polynom& P) //присваивание (не работает)
{
	this->Coeffs = P.Coeffs;
	this->Deg = P.Deg;
	this->Grid = P.Grid;
	return *this;
}

double Polynom::get_coeff(int i)
{
	return this->Coeffs[i];
}

double Polynom::value(double x)
{
	if (Grid.size() == 0) return(vect_to_polynom(Coeffs, x));
	else
	{
		double result = 0;
		for (int i = 0; i <= Deg; i++)
		{
			double product_i = 1;
			for (int j = 0; j <= Deg; j++)
			{
				if (i == j) continue;
				product_i *= (x - Grid[j]);
			}
			result += Coeffs[i] * product_i;
		}
		return result;
	}
}

int Polynom::get_Deg()
{
	return this->Deg;
}

bool Polynom::lagrange()
{
	if (Grid.size() == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

vector<double> Chebyshev_grid(int N, double a, double b) { //возвращает сетку нулей многочлена Чебышёва
	std::vector<double> X;
	for (int i = 0; i < N; i++) {
		X.push_back((b + a) / 2.0 + (b - a) / 2.0 * cos(3.141592653 * (2 * i + 1) / (2 * N)));
	}
	return X;
}

Polynom BUA_for_Polynom(Polynom& P, double a, double b) //принимает на вход многочлен P; узнаёт его степень и возвращает его наилучшее равномерное приближение степени N-1 на отрезке [a;b]
{
	cout.setf(ios::fixed);
	Polynom Q;
	int N = P.get_Deg();
	//if (P.lagrange() == false)
	//{
	//	N++;
	//}
	vector<double> G = Chebyshev_grid(N, a, b);
	double alpha = 1;
	vector<double> values;
	for (int i = 0; i < N; i++)
	{
		alpha = 1;
		for (int j = 0; j < N; j++)
		{
			if (i == j) continue;
			alpha *= 1.0 / (G[i] - G[j]);
		}
		values.push_back(P.value(G[i]) * alpha);
	}
	Q = Polynom(values, G);
	cout << Q.get_Deg() << " " << G.size() - 1 << endl;

	//double x, p, q, dev, max_dev;
	//max_dev = 0;
	//for (int i = 0; i <= 2000; i++)
	//{
	//	x = i * 0.001;
	//	p = P.value(x);
	//	q = Q.value(x);
	//	dev = abs(p - q);
	//	cout << setw(2) << setprecision(3) << x << " " << setprecision(15) << dev << " ";
	//	cout << endl;
	//}
	//cout << "\n" << endl;
	return Q;
}

double deviation(Polynom& Q)
{
	double a = 0.0;
	double b = 2.0;
	double x, f, q, dev, max_dev, min_dev, mean_dev;
	max_dev = 0;
	min_dev = 1.0;
	mean_dev = 0;
	for (int i = 0; i <= 2000; i++)
	{
		x = i * 0.001;
		f = expression(x);
		q = Q.value(x);
		dev = abs(f - q);
		if (max_dev < dev)
		{
			max_dev = dev;
		}
		if (min_dev > dev)
		{
			min_dev = dev;
		}
		mean_dev += dev;
		//cout << setw(2) << setprecision(3) << x << " " << setprecision(15) << dev << " " << endl;
	}
	mean_dev *= 0.002;
	cout << setprecision(10) << "minimum deviation: " << min_dev << "     mean deviation: " << mean_dev << "     maximum deviation: " << max_dev << endl;
	return max_dev;
}

Polynom BUA_for_Function()
{
	double a = 0.0; double b = 2.0;
	vector<double> X = { 0.9998529453, 1.003762849, -0.5449030727, -0.4998272697, -1.190897556, 7.538626711, -20.51685510, 46.17272397, -87.86076302, 134.3446405, -162.3206992, 155.4993829, -118.8785478, 72.73764734, -35.51448972, 13.69556335, -4.089527504, 0.9139089311, -0.1440991992, 0.01431562880, -0.0006746939898 };/*{ 0.349349163207821e-3, .994977456991976, 0.26776030251257e-1, -1.35412996071123, -.565523671624120, 5.60395583308993, -16.1192149183527, 42.5381062819594, -90.1442172348290, 142.232741966616, -168.521338373430, 153.150683261450, -108.305498239531, 59.8874369295559, -25.7608129773473, 8.48079451903690, -2.07035771939434, .353947182547754, -0.378910185027301e-1, 0.1914334529e-2 };*/
	Polynom P(X);
	Polynom Q = BUA_for_Polynom(P, a, b);
	double dev_new, dev_old; Polynom Buf;
	dev_new = deviation(Q);
	FILE* out;
	fopen_s(&out,"task5/difference.txt","w");
	while (dev_new <= 0.0005)
	//while (dev_new <= 0.001)
	{
		P = Q;
		dev_old = dev_new;
		Q = BUA_for_Polynom(P, a, b);
		dev_new = deviation(Q);
		fprintf(out,"%d %lf\n", P.get_Deg(), dev_new);
		//if (dev_old < dev_new) cout << "ALARM! " << endl;
	}
	cout << "THE PROCESS IS FINISHED. N = " << P.get_Deg() << ", dev = " << dev_old << endl;
	return P;
}

