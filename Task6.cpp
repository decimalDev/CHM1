#include<iostream>
#include <fstream>
#include<vector>
#include<stdio.h>
#include<cmath>
#include <iomanip>
#include "Splines.h"
#include "Matrix algebra.h"

#include "Splines.cpp"
#include "Matrix algebra.cpp"

using namespace std;

int task6()
{
	/*vector<double> R1 = { 2, 1, 0, 0 };
	vector<double> R2 = { 3, 2, 1, 0 };
	vector<double> R3 = { 0, 3, 2, 1 };
	vector<double> R4 = { 0, 0, 3, 2 };
	vector<vector<double>> content = { R1, R2, R3, R4 };
	Matrix A(content);

	vector<double> Y = { 1, 1, 1, 1 };
	Vector B(Y);
	Vector X = Thomas(A, B);

	for (int i = 1; i <= X.getrows(); i++)
	{
		cout << X.getelement(i) << endl;
	}*/

	Spline S(1000);
	cout.setf(ios::fixed);
	cout << S.value(1.0) << endl;
	S.dev(10000);
	return 0;
}
