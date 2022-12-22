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


std::pair<double,double> Coeffisient(double (*function)(double), std::vector<double> X,int k){
	size_t n = X.size();//=2*N+1
	double sum_AK = 0;
	double sum_BK = 0;
	for(int j = 0;j<n;j++){
		sum_AK+=function(X[j])*cos(2*M_PI*k*j/n);
		sum_BK+=function(X[j])*sin(2*M_PI*k*j/n);
	}
	return std::pair<double,double>{2*sum_AK/n,2*sum_BK/n};
}

void Fourier_coeffs(double (*function)(double),std::vector<double> &X,std::vector<std::pair<double,double>> &Coeffisients){
	size_t n = X.size();
	for(int k = 0;k<n/2;k++){
		Coeffisients.push_back(Coeffisient(function,X,k));
	}
	Coeffisients[0].first = Coeffisients[0].first/2;
}

double Fourier_value(double parametr_x,std::vector<double> &X,std::vector<std::pair<double,double>> &Coeffisients){
	size_t n = Coeffisients.size();
	
	double result = 0;
	
	for(int k = 0;k<n;k++){
		result += Coeffisients[k].first*cos(k*parametr_x) + Coeffisients[k].second*sin(k*parametr_x);
		//printf("res=%lf\n",result);
	}

	return result;
}




double expression2(double x) {
	return exp(alpha*x+betta) / (1 + (alpha*x+betta)*(alpha*x+betta));
}




std::vector<double> Uniform_grid(double a, double b, int n) { //строит равномерную n-сетку
	double h = (b - a) / (n);
	std::vector<double> X;
	for (int i = 0; i <= n; i++) {
		X.push_back(a + h * i);
	}
	return X;
}



double Fourier(double (*function)(double), double a, double b, int n, double x){	
	n = 2*n + 1;
	alpha = (b-a)/2/M_PI;
	betta = a;
	x = (x-betta)/alpha;
	double h_interpolation = (2*M_PI)/n;
	std::vector<double> X_interpolation;
	for(int i = 1;i<=n;i++){
		X_interpolation.push_back(h_interpolation*(i-1));
	}
	std::vector<std::pair<double,double>> Coeffisients;
	Fourier_coeffs(function,X_interpolation,Coeffisients);

	return Fourier_value(x,X_interpolation,Coeffisients);
	
}

double Fourier_optimal(double (*function)(double), double a, double b, int n, double x,std::vector<std::pair<double,double>> Coeffisients,std::vector<double> &X_interpolation){	
	//n = 2*n + 1;
	alpha = (b-a)/2/M_PI;
	betta = a;
	x = (x-betta)/alpha;

	return Fourier_value(x,X_interpolation,Coeffisients);
	
}

double Fourier_max_deviation(double (*function)(double),double a,double b,int n){ //вычисляет наибольшее отклонение многочлена L_n (x) от начальной функции, перебирая значения на равномерной 10^5-сетке
	
	
	int grid = 100000;
	double h = (b - a) / (grid);
	std::vector<double> X;
	for (int i = grid*5/8; i <= grid*6/8; i++) {
		X.push_back(a + h * i);
	}
	//std::cout<<"hi there"<<std::endl;
	double max = 0;
	double value;
	
	n = 2*n + 1;
	double h_interpolation = (2*M_PI)/n;
	std::vector<double> X_interpolation;
	for(int i = 1;i<=n;i++){
		X_interpolation.push_back(h_interpolation*(i-1));
	}
	std::vector<std::pair<double,double>> Coeffisients;
	Fourier_coeffs(function,X_interpolation,Coeffisients);
	
	for(int i = 0;i < grid/8+1;i++){
		value = function(X[i]) - Fourier_optimal(function,a,b,n,X[i],Coeffisients,X_interpolation);
		value = std::abs(value);
		//fprintf(out,"%lf %lf\n",i,value);
	
		if(value>max) max = value;
	}
	
	return max;
}

void Fourier_show(double (*function)(double),int n){//рисует для [0,2*pi] 
	FILE* out;
	fopen_s(out,"Fourier/Fourier.txt","w");
	
	int grid = 10000;
	double h = (2*M_PI)/grid;
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
	fopen_s(&out,"Fourier/Fourier.txt","w");
	
	int grid = 10000;
	double h = (2*M_PI)/grid;
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
		
		fprintf(out,"%lf %lf\n",X[i]*alpha+betta,value);
	}
	
}


void Fourier_show3(double (*function)(double),int n){
	FILE* out;
	fopen_s(out,"Fourier/Fourier.txt","w");
	
	int grid = 10000;
	double h = 2.0/grid;
	grid++;
	std::vector<double> X;
	for(int i = 0;i<grid;i++) X.push_back(i*h);

	double value;
	//int n_temp = n;
	
	n = 2*n + 1;
	double h_interpolation = (2*M_PI)/n;
	std::vector<double> X_interpolation;
	for(int i = 1;i<=n;i++){
		X_interpolation.push_back(h_interpolation*(i-1));
	}
	std::vector<std::pair<double,double>> Coeffisients;
	Fourier_coeffs(function,X_interpolation,Coeffisients);
	
	for(int i = 0;i<grid;i++){

		value = Fourier_optimal(function,0,2,n,X[i],Coeffisients,X_interpolation);
		//std::cout<<"value="<<value<<" func="<<function(X[i])<<" intep="<<Newton_optimal(function,X_interpolation,Coeffisients,X[i])<<std::endl;
		
		
	}
	
}

void task1(){
	std::cout<<"task1"<<std::endl;
	int n;
	std::cout<<"n = ";
	std::cin>>n;
	Fourier_show(expression,n);
	std::cout<<"write any letter"<<std::endl;
	char c;
	std::cin>>c;
}

void task2(){
	std::cout<<"task2"<<std::endl;
	int n;
	std::cout<<"n = ";
	std::cin>>n;
	Fourier_show2(expression2,n);
	std::cout<<"write any letter"<<std::endl;
	char c;
	std::cin>>c;
}

void task34(){
	std::cout<<"task3"<<std::endl;
	
	int n;
	double max = 1000;
	n = 2;
	
	FILE* out;
	fopen_s(&out,"Fourier/max_difference.txt","w");
	
	
	while(max>0.077){
		max = Fourier_max_deviation(expression2,0,2,n);
		std::cout<<"n = "<<n<<" max = "<<max<<std::endl;
		fprintf(out,"%d %lf\n",n,max);
		n++;
	}
	std::cout<<"optimal n is "<<n<<std::endl;
	std::cout<<"write any letter"<<std::endl;
	char c;
	std::cin>>c;
	fclose(out);
}


int NM_task4(){
	task1();
	task2();
	task34();
	return 0;
}