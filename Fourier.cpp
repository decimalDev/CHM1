#define _USE_MATH_DEFINES
#include<iostream>
#include<cmath>
#include<stdio.h>
#include<vector>


double alpha = 1/M_PI;
double betta = 0;

double expression(double x) {
	return exp(x) / (1 + x * x);
}

double expression2(double x) {
	return exp(x) / (1 + (alpha*x+betta)*(alpha*x+betta));
}

std::pair<double,double> Coeffisient(double (*function)(double), std::vector<double> X,int k){
	int n = X.size();//=2*N+1
	double sum_AK = 0;
	double sum_BK = 0;
	for(int j = 0;j<n;j++){
		sum_AK+=function(X[j])*cos(2*M_PI*k*j/n);
		sum_BK+=function(X[j])*sin(2*M_PI*k*j/n);
	}
	return std::pair<double,double>{2*sum_AK/n,2*sum_BK/n};
}

void Fourier_coeffs(double (*function)(double),std::vector<double> &X,std::vector<std::pair<double,double>> &Coeffisients){
	int n = X.size();
	for(int k = 0;k<n/2;k++){
		Coeffisients.push_back(Coeffisient(function,X,k));
	}
	Coeffisients[0].first = Coeffisients[0].first/2;
}

double Fourier_value(double parametr_x,std::vector<double> &X,std::vector<std::pair<double,double>> &Coeffisients){
	int n = Coeffisients.size();
	
	double result = 0;
	
	for(int k = 0;k<n;k++){
		result += Coeffisients[k].first*cos(k*parametr_x) + Coeffisients[k].second*sin(k*parametr_x);
		//printf("res=%lf\n",result);
	}

	return result;
}

void Fourier_show(double (*function)(double),int n){
	FILE* out;
	out = fopen("Fourier/Fourier.txt","w");
	
	int grid = 10000;
	double h = (M_PI)/grid;
	grid++;
	std::vector<double> X;
	for(int i = 0;i<grid;i++) X.push_back(i*h);

	double value;
	
	n = 2*n + 1;
	double h_interpolation = (2*M_PI)/n;
	std::vector<double> X_interpolation;
	for(int i = 1;i<=n;i++){
		X_interpolation.push_back(h_interpolation*(i-1));
	}
	std::vector<std::pair<double,double>> Coeffisients;
	Fourier_coeffs(function,X_interpolation,Coeffisients);
	
	for(int i = 0;i<grid;i++){

		value = Fourier_value(X[i],X_interpolation,Coeffisients);
		//std::cout<<"value="<<value<<" func="<<function(X[i])<<" intep="<<Newton_optimal(function,X_interpolation,Coeffisients,X[i])<<std::endl;
		
		fprintf(out,"%lf %lf\n",X[i],value);
	}
	
}

void Fourier_show2(double (*function)(double),int n){
	FILE* out;
	out = fopen("Fourier/Fourier.txt","w");
	
	int grid = 10000;
	double h = (M_PI)/grid;
	grid++;
	std::vector<double> X;
	for(int i = 0;i<grid;i++) X.push_back(i*h);

	double value;
	
	n = 2*n + 1;
	double h_interpolation = (2*M_PI)/n;
	std::vector<double> X_interpolation;
	for(int i = 1;i<=n;i++){
		X_interpolation.push_back(h_interpolation*(i-1));
	}
	std::vector<std::pair<double,double>> Coeffisients;
	Fourier_coeffs(function,X_interpolation,Coeffisients);
	
	for(int i = 0;i<grid;i++){

		value = Fourier_value(X[i],X_interpolation,Coeffisients);
		//std::cout<<"value="<<value<<" func="<<function(X[i])<<" intep="<<Newton_optimal(function,X_interpolation,Coeffisients,X[i])<<std::endl;
		
		fprintf(out,"%lf %lf\n",X[i],value);
	}
	
}

int main(){
	int n;
	std::cout<<"n = ";
	std::cin>>n;
	Fourier_show2(expression2,n);
	return 0;
}