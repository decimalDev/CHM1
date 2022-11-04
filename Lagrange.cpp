//var 1 y=exp(x)/(1+x^2) [0,2] delta = 5*10^(-4)

#include<iostream>
#include<vector>
#include<stdio.h>
#include<cmath>

double expression(double x) { //функция из варианта
	return exp(x) / (1 + x * x);
} 



void Lagrange_coeffs(double (*function)(double),std::vector<double> &X,std::vector<double> &Coeffisients){ //получает функцию, сетку, пустой вектор; заполняет вектор коэффициентами многочлена Лагранжа для данной функции по этой сетке
		double a = 1;
		int n = X.size();
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
	int n = Coeffisients.size(); //вычисление значения многочлена Лагранжа в заданной точке
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

double Lagrange(double (*function)(double),double a,double b,int n,double x){
	//double a = 0,b = 2; //строит на отрезке [a;b] равномерную n-сетку, строит на этой сетке для заданной функции лагранжеву интерполяцию (Lagrange_coeffs), затем вычисляет значение этой интерполяции в заданной точке (Lagrange_value); возвращает полученное значение
	//int n = 10;
	double h = (b-a)/(n);
	n++;
	//double x = 0.5;
	std::vector<double> Coeffisients;
	std::vector<double> X;
	for(int i = 0;i<n;i++){
		X.push_back(a+h*i);
	}
	Lagrange_coeffs(function,X,Coeffisients);
	double res = Lagrange_value(x,X,Coeffisients);
	
	//std::cout<<"Josepf Lui Lagrange is tired. But he could complete the calculation and found number is "<<res<<std::endl;
	
	return res;
}

double Lagrange_max_deviation(double (*function)(double),double a,double b,int n){ //вычисляет наибольшее отклонение многочлена L_n (x) от начальной функции, перебирая значения на равномерной 10^5-сетке

	
	int grid = 100000;
	double h = (b-a)/grid;
	grid++;
	std::vector<double> X;
	for(int i = 0;i<grid;i++) X.push_back(a+i*h);
	
	double max = 0;
	double value;
	//printf("exp(x) = %lf\n",expression());
	for(int i = 0;i<grid;i++){
		value = function(X[i]) - Lagrange(function,a,b,n,X[i]);
		value = std::abs(value);
	
		if(value>max) max = value;
	}
	
	return max;
}

int task2(int n_max){ //для всех n от 1 до n_max вычисляет максимальное отклонение (Lagrange_max_deviation) многочлена Лагранжа от заданной функции; возвращает n_0; строит график отклонения f от n_0
	FILE* out = fopen("Lagrange/max_difference.txt","w");
	double max_difference[n_max];
	int n_0 = 1;
	//max_difference[0] = 0;
	
	
	for(int n = 1;n<=n_max;n++){
		max_difference[n-1] = Lagrange_max_deviation(expression,0,2,n);
		printf("n = %d\t max_dif = %lf\n",n,max_difference[n-1]);
		fprintf(out,"%d %lf\n",n,max_difference[n-1]);
		if(max_difference[n-1]<max_difference[n_0]) n_0 = n-1;
	}
	printf("optimal n is %d\n",n_0);

	FILE* dev_0 = fopen("Lagrange/n0_difference.txt", "w");
	int grid = 100000;
	double h = (2.0) / grid;
	grid++;
	std::vector<double> X;
	for (int i = 0; i < grid; i++) X.push_back( i * h);
	for (int i = 0; i < grid; i++) {
		fprintf(dev_0, "%lf %lf\n", X[i], expression(X[i]) - Lagrange(expression, 0, 2, n_0, X[i]));
	}
	return 0;
}

int main(){
	int n = 10;
	task2(15);
	return 0;
}