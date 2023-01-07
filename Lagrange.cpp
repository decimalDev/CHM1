//var 1 y=exp(x)/(1+x^2) [0,2] delta = 5*10^(-4)
//#pragma warning(disable : 6385)
//#pragma warning(disable : 26451)
#pragma once
#define _USE_MATH_DEFINES
#include<vector>
#include<stdio.h>
#include<cmath>
#include<iostream>
#include "Newton.hpp"
#define M_PI 3.14159265358979323846


double expression(double x) { //функция из варианта
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

void Lagrange_coeffs(double (*function)(double),std::vector<double> &X,std::vector<double> &Coeffisients){ //получает функцию, сетку, пустой вектор; заполняет вектор коэффициентами многочлена Лагранжа для данной функции по этой сетке
		double a = 1;
		size_t n = X.size();
		for(int i = 0;i<n;i++){
			a = 1;
			for(int j = 0;j<n;j++){
				if(i==j) continue;
				a *= 1.0/(X[i]-X[j]);
			}
			Coeffisients.push_back(function(X[i])*a);
		}	
}

double Lagrange_value(double parametr_x,std::vector<double> &X,std::vector<double> &Coeffisients){
	size_t n = Coeffisients.size(); //вычисление значения многочлена Лагранжа в заданной точке
	double result = 0;
	for(int i = 0;i<n;i++){
		double product_i = 1;
			for(int j = 0;j<n;j++){
				if(i==j) continue;
				product_i*=(parametr_x - X[j]);
			}
	result+=Coeffisients[i]*product_i;
	}
	return result;
}

double Lagrange_grid(double (*function)(double),double a,double b,double x,std::vector<double> &X){
	//строит на заданной сетке X для заданной функции лагранжеву интерполяцию (Lagrange_coeffs), затем вычисляет значение этой интерполяции в заданной точке (Lagrange_value); возвращает полученное значение
	std::vector<double> Coeffisients;
	Lagrange_coeffs(function,X,Coeffisients);
	double res = Lagrange_value(x,X,Coeffisients);
	return res;
}

double Lagrange(double (*function)(double), double a, double b, int n, double x) { //то же для равномерной сетки
	std::vector<double> X = Uniform_grid(a, b, n);
	return Lagrange_grid(function, a, b, x, X);
}

double Lagrange_max_deviation(double (*function)(double),double a,double b,int n){ //вычисляет наибольшее отклонение многочлена L_n (x) от начальной функции, перебирая значения на равномерной 10^5-сетке

	int grid = 100000;
	std::vector<double> X = Uniform_grid(a, b, grid);
	
	double max = 0;
	double value;
	for(int i = 0;i<grid;i++){
		value = function(X[i]) - Lagrange(function,a,b,n,X[i]);
		value = std::abs(value);
	
		if(value>max) max = value;
	}
	
	return max;
}

void task2(int n_max, int &n_opt, double &max_diff){ //для всех n от 1 до n_max вычисляет максимальное отклонение (Lagrange_max_deviation) многочлена Лагранжа от заданной функции; возвращает n_0; строит график отклонения f от n_0
	FILE* out;
	//fopen_s(&out,"Lagrange/max_difference.txt","w");
	fopen_s(&out, "Lagrange/max_difference.txt", "w");
	if (out == NULL) {
		printf("ERROR\n");
	}
	double *max_difference = new double[n_max];
	int n_0 = 1;
	//max_difference[0] = 0;
	
	
	for(int n = 1; n <= n_max; n++){
		int temp_n = n - 1;
		max_difference[temp_n] = Lagrange_max_deviation(expression,0,2,n);
		printf("n = %d\t max_dif = %lf\n",n,max_difference[temp_n]);
		if(out!=0)
		fprintf(out,"%d %.14lf\n", n, max_difference[temp_n]);

		if(max_difference[temp_n]<max_difference[n_0]) n_0 = n - 1;
	}
	printf("optimal n is %d\n",n_0+1);
	
	n_0++;
	n_opt = n_0;
	max_diff = max_difference[n_0 - 1];

	FILE* dev_0;
	fopen_s(&dev_0,"Lagrange/n0_difference.txt", "w");
	int grid = 100000;
	double h = (2.0) / grid;
	grid++;
	std::vector<double> X;
	for (int i = 0; i < grid; i++) X.push_back( i * h);
	for (int i = 0; i < grid; i++) {
		if(dev_0!=0)
		fprintf(dev_0, "%lf %.14lf\n", X[i], std::abs(expression(X[i]) - Lagrange(expression, 0, 2, n_0+1, X[i])));
	}
}

std::vector<double> Chebyshev_zeroes(int N, double a, double b) { //возвращает сетку нулей многочлена Чебышёва
	std::vector<double> X;
	for (int i = 0; i < N; i++) {
		X.push_back((b + a) / 2.0 + (b - a) / 2.0 * cos(M_PI * (2 * i + 1) / (2 * N)));
	}
	return X;
}

double Chebyshev(int N, double (*function)(double), double a, double b, double x) { //в заданной точке находит значение многочлена Лагранжа степени N, построенного на сетке из нулей Чебышёва
	std::vector<double> X = Chebyshev_zeroes(N, a, b);
	double res = Lagrange_grid(function, a, b, x, X);
	return res;
}

double Chebyshev_max_deviation(double (*function)(double), double a, double b,int n) {
	int grid = 10000;
	std::vector<double> X = Uniform_grid(a, b, grid);
	double max = 0;
	double value;
	for (int i = 0; i <= grid; i++) {
		value = function(X[i]) - Chebyshev(n, function, a, b, X[i]);
		value = std::abs(value);
		if (value > max) max = value;
	}
	return max;
}

double Chebyshev_compare(double (*function)(double), double a, double b, int n) {
	int grid = 100000;
	std::vector<double> X = Uniform_grid(a, b, grid);

	FILE* out;
	//fopen_s(&out,"Lagrange/max_difference.txt","w");
	fopen_s(&out, "Lagrange/Chebyshev_n0_deviation.txt", "w");
	if (out == NULL) {
		printf("ERROR\n");
	}

	double max = 0;
	double value;
	for (int i = 0; i < grid; i++) {
		value = function(X[i]) - Chebyshev(n, function, a, b, X[i]);
		
		value = std::abs(value);
		fprintf(out, "%.6lf %.14lf\n", X[i], value);
		if (value > max) max = value;
	}
	fclose(out);
	return max;
}

void compare_L_and_N(int n){
	double L = Lagrange_max_deviation(expression,0,2,n);
	double N = Newton_max_deviation(expression,0,2,n);
	printf("max deviation Lagrange:%.10lf and Newton:%.10lf\n",L,N);
}

int Lagrange_task(){
	int n = 10;
	int n_0; double max_div_uniform;
	task2(15, n_0, max_div_uniform);
	double max_div_chebyshev = Chebyshev_compare(expression, 0, 2, n_0);
	printf("max div uniform %.10lf\n", max_div_uniform);
	printf("max_div_chebyshev %.10lf\n", max_div_chebyshev);
	if (max_div_chebyshev < max_div_uniform)  printf("Chebyshev grid is better than the uniform one\n");
	else printf("Chebyshev grid is not better that the uniform one\n");
	compare_L_and_N(n_0);
	return 0;
}