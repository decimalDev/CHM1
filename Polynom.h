#pragma once
#include<iostream>
#include <fstream>
#include<vector>
#include<stdio.h>
#include<cmath>
#include <iomanip>

using namespace std;

double expression(double x);

double vect_to_polynom(std::vector<double>& P, double x);

class Polynom {
private:
	vector<double> Coeffs; //вектор коэффициентов -- СОДЕРЖАНИЕ РАЗЛИЧНО ДЛЯ ОБЫЧНЫХ МНОГОЧЛЕНОВ И МНОГОЧЛЕНОВ ЛАГРАНЖА: в случае обычного многочлена (такой у нас в программе ровно один -- это многочлен Тейлора, плюс дефолтный нулевой многочлен) это стандартная значащая часть последовательности коэффициентов (от младших к старшим), определяющей многочлен (поле GRID = null); в случае многочлена Лагранжа это последовательность коэффициентов в записи многочлена Лагранжа (они не расписаны по степеням), в поле GRID хранится сетка, по которой строился многочлен
	int Deg;
	vector<double> Grid; //сетка в случае многочлена Лагранжа
public:
	Polynom(); //строит нулевой многочлен
	Polynom(const Polynom& P); //копирование
	Polynom(vector<double>& X); //строит многочлен по заданной последовательности коэффициентов (без сетки!)
	Polynom(vector<double>& X, vector<double>& G); //строит многочлен по заданной последовательности коэффициентов (без сетки!)
	Polynom& operator=(const Polynom& P);
	double get_coeff(int i);
	double value(double x);
	int get_Deg();
	bool lagrange();
};


vector<double> Chebyshev_grid(int N, double a, double b);

Polynom BUA_for_Polynom(Polynom& P, double a, double b);

double deviation(Polynom& Q);

Polynom BUA_for_Function();
