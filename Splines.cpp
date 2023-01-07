#include "Splines.h"

/*
double expression(double x)
{
	return exp(x) / (1 + x * x);
}

std::vector<double> Uniform_grid(double a, double b, int n) { //строит равномерную n-сетку
	double h = (b - a) / (n);
	std::vector<double> X;
	for (int i = 0; i <= n; i++) {
		X.push_back(a + h * i);
	}
	return X;
}
*/
Spline::Spline(int n) //строит сплан функции expression на равномерной сетке размера n
{
	vector<double> Grid = Uniform_grid(0.0, 2.0, n);
	double h = 2.0 / (double)n;
	this->Grid = Grid;

	vector<double> Y;
	for (double x : Grid)
	{
		Y.push_back(expression(x));
	}

	this->ABCD.setcontent(n, 4);
	Matrix ABCD(n, 4);
	for (int i = 1; i <= n; i++)
	{
		ABCD.setelement(i, 1, Y.at(i-1));
	}

	//Приводим матрицу ABCD к трёхдиагональному виду
	Matrix ThreeDiag(n + 1, n + 1);
	ThreeDiag.setelement(1, 1, 1.0);
	ThreeDiag.setelement(n + 1, n + 1, 1.0);
	for (int i = 2; i <= n; i++)
	{
		ThreeDiag.setelement(i, i - 1, h);
		ThreeDiag.setelement(i, i, 2 * h);
		ThreeDiag.setelement(i, i + 1, h);
	}

	//Задаём столбец значений
	Vector B(n + 1);
	for (int i = 2; i <= n; i++)
	{
		B.setelement(i, (3 / h) * (Y[i] - 2 * Y.at(i-1) + Y.at(i-2)));
	}

	//Вычисляем коэффициенты c_i
	Vector C = Thomas(ThreeDiag, B);
	for (int i = 1; i <= n; i++)
	{
		cout << C.getelement(i) << "\t";
	}

	//Вычисляем коэффициенты d_i
	for (int i = 1; i < n; i++)
	{
		double d = (ABCD.getelement(i + 1, 3) - ABCD.getelement(i, 3)) / (3 * h);
		ABCD.setelement(i, 4, d);
	}
	ABCD.setelement(n, 4, -1.0 * ABCD.getelement(n, 3) / (3 * h));

	//Вычисляем коэффициенты b_i
	for (int i = 1; i < n; i++)
	{
		double b = (Y.at(i) - Y.at(i-1)) / h - (1.0 / 3.0) * h * (ABCD.getelement(i + 1, 3) + 2 * ABCD.getelement(i, 3));
		ABCD.setelement(i, 2, b);
	}
	double b = (Y.at(n) - Y.at(n-1)) / h - (2.0 / 3.0) * h * ABCD.getelement(n, 3);
	ABCD.setelement(n, 2, b);

	//Матрица ABCD вычислена
	this->ABCD = ABCD;
}

double Spline::value(double x)
{
	//определяем, какому отрезку сетки вида [x_{m - 1}; x_m] принадлежит точка x
	double h = 2.0 / Grid.size();
	double m;
	m = 0;
	do
	{
		m++;
		if (m == ABCD.getrows())
		{
			break;
		}
	} while ((double)m * h <= x);

	double res = ABCD.getelement(m, 1) + ABCD.getelement(m, 2) * (x - Grid[m - 1]) + ABCD.getelement(m, 3) * (x - Grid[m - 1]) * (x - Grid[m - 1]) + ABCD.getelement(m, 4) * (x - Grid[m - 1]) * (x - Grid[m - 1]) * (x - Grid[m - 1]);
	return res;
}

double Spline::dev_show(int gridcard)
{
	vector<double> X = Uniform_grid(0.0, 2.0, gridcard);
	double dev, max_dev, f, s;
	max_dev = 0;
	//FILE* out;
	//fopen_s(&out, "task6/difference.txt", "w");
	for (double x : X)
	{
		f = expression(x);
		s = value(x);
		dev = abs(f - s);
		if (max_dev < dev) max_dev = dev;
		//cout << x << "\t" << f << "\t" << s << "\t" << setprecision(15) << dev << endl;
		//fprintf(out,"%.9lf %.9lf\n",x,dev);
	}
	//cout << "FINISHED. Maximum deviation for grid " << gridcard << " equals " << setprecision(15) << max_dev << endl;

	return max_dev;
}

double Spline::dev(int gridcard)
{
	vector<double> X = Uniform_grid(0.0, 2.0, gridcard);
	double dev, max_dev, f, s;
	max_dev = 0;
	//FILE* out;
	//fopen_s(&out, "task6/difference.txt", "w");
	for(double x : X)
	{
		f = expression(x);
		s = value(x);
		dev = abs(f - s);
		if (max_dev < dev) max_dev = dev;
		//cout <<endl<< x << "\t" << f << "\t" << s << "\t" << setprecision(15) << dev << endl;
		//fprintf(out,"%.9lf %.9lf\n",x,dev);
	}
	cout << "FINISHED. Maximum deviation for grid " << gridcard << " equals " << setprecision(15) << max_dev << endl;

	return max_dev;
}
