#include "Matrix algebra.h"

// класс ћј“–»÷ј характеризуетс€ трем€ пол€ми: числом строк (rows), числом столбцов (columns) и двумерным массивом элементов (content)

int Matrix::getrows()
{
	return rows;
}

int Matrix::getcolumns()
{
	return columns;
}

void Matrix::setcontent(int m, int n) // по дефолту content -- нулева€ матрица соответствующего размера
{
	rows = m;
	columns = n;
	vector<vector<double>> Theta;
	Theta.resize(m);
	for (int i = 0; i < m; i++)
	{
		Theta[i].resize(n);
		fill(Theta[i].begin(), Theta[i].end(), 0.0);
	}
	content = Theta;
}

//метод, позвол€ющий заполнить выбранную €чейку нужным числом
void Matrix::setelement(int i, int j, double a)
{
	if (i < 1 || j < 1 || i > rows || j > columns)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		content[i-1][j-1] = a;
	}
}

//метод, извлекающий значение из нужной €чейки
double Matrix::getelement(int i, int j)
{
	if (i < 1 || j < 1 || i > rows || j > columns)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		return content[i-1][j-1];
	}
}

//дефолтный конструктор строит матрицу 3x3
Matrix::Matrix()
{
	rows = 3;
	columns = 3;
	setcontent(3,3);
}
	
//конструктор, стро€щий нулевую матрицу m x n
Matrix::Matrix(int m, int n)
{
	rows = m;
	columns = n;
	setcontent(m, n);
}

//конструктор, стро€щий матрицу m x n и заполн€ющий еЄ элементами данного массива C
Matrix::Matrix(vector<vector<double>>& C)
{
	int m = C.size();
	int n = C[0].size();
	setcontent(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			content[i][j] = C[i][j];
		}
	}
	cout << "Matrix created!" << endl;
}

//ѕоложить элемент в вектор
void Vector::setelement(int i, double a)
{
	if (i < 1 || i > rows)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		content[i-1][0] = a;
	}
}

//получить элемент вектора
double Vector::getelement(int i)
{
	if (i < 1 || i > rows)
	{
		throw "\n INCORRECT INPUT FOR MATRIX CELL. \n";
	}
	else
	{
		return content[i-1][0];
	}
}

//дефолтный конструктор вектора -- нулева€ матрица 3 x 1
Vector::Vector()
{
	rows = 3;
	columns = 1;
	setcontent(3, 1);
}

//нулевой вектор на m строках
Vector::Vector(int m)
{
	rows = m;
	columns = 1;
	setcontent(m, 1);
}

//вектор размера m, заполненный элементами массива C
Vector::Vector(vector<double>& C)
{
	int m = C.size();
	columns = 1;
	setcontent(m,1);
	for (int i = 0; i < m; i++)
	{
		content.at(i).at(0) = C.at(i);
	}
}

//оператор присваивани€ дл€ векторов
Vector& Vector::operator= (const Vector& v)
{
	if (rows != v.rows)
	{
		throw "Error.";
	}
	for (int i = 0; i < v.rows; i++)
	{
		content[i-1].at(0) = v.content[i-1].at(0);
	}
	return *this;
}

Vector Thomas(const Matrix& A, const Vector& B) //ѕ–ќ√ќЌ ј ƒЋя ѕ–ќ»«¬ќЋ№Ќќ… “–®’ƒ»ј√ќЌјЋ№Ќќ… ћј“–»÷џ
{
	Matrix A_mod = A;
	Vector D = B;
	int N = A_mod.getrows();
	int M = A_mod.getcolumns();
	for (int i = 1; i <= min(M, N); i++)
	{
		A_mod.setelement(i, i, -1.0 * A_mod.getelement(i, i)); //умножаем главную диагональ на минус, чтобы свести к каноническому виду
	}

	//получаем канонический вид дл€ уравнени€ трЄхдиагональной системы
	Matrix Coeffs(N, 3);
	Coeffs.setelement(1, 1, 0.0); //a_1 = 0;
	Coeffs.setelement(1, 2, A_mod.getelement(1, 1));
	Coeffs.setelement(1, 3, A_mod.getelement(1, 2));

	Coeffs.setelement(N, 1, A_mod.getelement(N, M - 1));
	Coeffs.setelement(N, 2, A_mod.getelement(N, M));
	Coeffs.setelement(N, 3, 0.0); //c_N = 0;

	for (int i = 2; i <= N-1; i++)
	{
		Coeffs.setelement(i, 1, A_mod.getelement(i, i - 1));
		Coeffs.setelement(i, 2, A_mod.getelement(i, i)); //столбец b заполн€етс€ элементами главной диагонали
		Coeffs.setelement(i, 3, A_mod.getelement(i, i + 1)); //полагаем, что на вход подаютс€ лишь "правильные" матрицы
	}

	Vector Xi(N + 1);
	Vector Eta(N + 1);
	Xi.setelement(1, 0.0);
	Eta.setelement(1, 0.0);
	double p;

	for (int i = 1; i <= N; i++)
	{
		p = (Coeffs.getelement(i, 2) - Coeffs.getelement(i, 1) * Xi.getelement(i));
		Xi.setelement(i + 1, Coeffs.getelement(i, 3) / p);
		Eta.setelement(i + 1, (Coeffs.getelement(i, 1) * Eta.getelement(i) - D.getelement(i)) / p);
	}

	Vector X(N + 1);
	X.setelement(N + 1, 0.0);
	for (int i = N; i >= 1; i--)
	{
		X.setelement(i, Xi.getelement(i + 1) * X.getelement(i + 1) + Eta.getelement(i + 1));
	}

	return X;
}