#pragma once
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

#pragma once

class Vector;

// ����� ������� ��������������� ����� ������: ������ ����� (rows), ������ �������� (columns) � ��������� �������� ��������� (content)
class Matrix
{
protected:
	int rows;
	int columns;
	vector<vector<double>> content;

public:
	int getrows();

	int getcolumns();

	void setcontent(int m, int n); // �� ������� content -- ������� ������� ���������������� �������

	//�����, ����������� ��������� ��������� ������ ������ ������
	void setelement(int i, int j, double a);

	//�����, ����������� �������� �� ������ ������
	double getelement(int i, int j);

	//��������� ����������� ������ ������� 3x3
	Matrix();

	Matrix(int m, int n);

	Matrix(vector<vector<double>>& C);
};

class Vector : public Matrix
{
public:
	void setelement(int i, double a);

	double getelement(int i);

	Vector();

	Vector(int m);

	Vector(vector<double>& C);

	//friend Vector operator* (const RotationMatrix& A, const Vector& X);

	Vector& operator= (const Vector& v);
};

Vector Thomas(const Matrix& A, const Vector& B);