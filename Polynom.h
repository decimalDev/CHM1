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
	vector<double> Coeffs; //������ ������������� -- ���������� �������� ��� ������� ����������� � ����������� ��������: � ������ �������� ���������� (����� � ��� � ��������� ����� ���� -- ��� ��������� �������, ���� ��������� ������� ���������) ��� ����������� �������� ����� ������������������ ������������� (�� ������� � �������), ������������ ��������� (���� GRID = null); � ������ ���������� �������� ��� ������������������ ������������� � ������ ���������� �������� (��� �� ��������� �� ��������), � ���� GRID �������� �����, �� ������� �������� ���������
	int Deg;
	vector<double> Grid; //����� � ������ ���������� ��������
public:
	Polynom(); //������ ������� ���������
	Polynom(const Polynom& P); //�����������
	Polynom(vector<double>& X); //������ ��������� �� �������� ������������������ ������������� (��� �����!)
	Polynom(vector<double>& X, vector<double>& G); //������ ��������� �� �������� ������������������ ������������� (��� �����!)
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
