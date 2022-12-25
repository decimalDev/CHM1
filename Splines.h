#pragma once
#include<iostream>
#include <fstream>
#include<vector>
#include<stdio.h>
#include<cmath>
#include <iomanip>
#include "Matrix algebra.h"

using namespace std;

double expression(double x);

std::vector<double> Uniform_grid(double a, double b, int n);

class Spline {
private:
	vector<double> Grid;
	Matrix ABCD;

public:
	Spline(int n);
	double value(double x);
	double dev(int gridcard);
};