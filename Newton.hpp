//var 1 y=exp(x)/(1+x^2) [0,2] delta = 5*10^(-4)

#include<iostream>
#include<vector>
#include<stdio.h>
#include<cmath>


/*
double expression(double x) {
	return exp(x) / (1 + x * x);
}
*/

double divided_differences(double (*function)(double),std::vector<double> &X,int n,int i,int j){
	if(n==0){
		return function(X[i]);
	}
	return ( divided_differences(function,X,n-1,i,j-1) - divided_differences(function,X,n-1,i+1,j) )/(X[i]-X[j]);
}

void Newton_coeffs(double (*function)(double),std::vector<double> X,std::vector<double> &Coeffisients){
	int n = X.size();
	for(int i = 0;i<n;i++){
		Coeffisients.push_back(divided_differences(function,X,i,0,i));
	}
}

double Newton_value(double x,std::vector<double> X,std::vector<double> &Coeffisients){
	int n = Coeffisients.size();
	
	double result = Coeffisients[0];
	
	for(int i = 1;i<n;i++){
		double coeffisient_i = Coeffisients[i];
		for(int j = i-1;j>-1;j--)
			coeffisient_i*=(x-X[j]);
		result += coeffisient_i;
	}
	
	
	return result;
}

double Newton(double (*function)(double),double a,double b,int n,double x){
	double h = (b-a)/n;
	n++;
	std::vector<double> X;
	for(int i = 0;i<n;i++){
		X.push_back(a+h*i);
	}
	std::vector<double> Coeffisients;
	Newton_coeffs(function,X,Coeffisients);

	double res = Newton_value(x,X,Coeffisients);
	
	return res;
}




double Newton_optimal(double (*function)(double),std::vector<double> &X,std::vector<double> &Coeffisients,double x){
	double res = Newton_value(x,X,Coeffisients);
	
	return res;
}



double Newton_max_deviation(double (*function)(double),double a,double b,int n){

	int grid = 100000;
	double h = (b-a)/grid;
	grid++;
	std::vector<double> X;
	for(int i = 0;i<grid;i++) X.push_back(a+i*h);
	
	double max = 0;
	double value;
	
	double h_interpolation = (b-a)/n;
	n++;
	std::vector<double> X_interpolation;
	for(int i = 0;i<n;i++){
		X_interpolation.push_back(a+h_interpolation*i);
	}
	std::vector<double> Coeffisients;
	Newton_coeffs(function,X_interpolation,Coeffisients);
	
	for(int i = 0;i<grid;i++){

		value = function(X[i]) - Newton_optimal(function,X_interpolation,Coeffisients,X[i]);
		//std::cout<<"value="<<value<<" func="<<function(X[i])<<" intep="<<Newton_optimal(function,X_interpolation,Coeffisients,X[i])<<std::endl;
		value = std::abs(value);
		
		//fprintf(out,"%lf %lf\n",X[i],value);
	
		if(value>max) max = value;
	}
	
	return max;
}


double Newton_show(double (*function)(double),double a,double b,int n){
	FILE* out;
	//if(write_n==n)
	out = fopen("Newton/Newton.txt","w");
	
	
	int grid = 1000;
	double h = (b-a)/grid;
	grid++;
	std::vector<double> X;
	for(int i = 0;i<grid;i++) X.push_back(a+i*h);
	
	double max = 0;
	double value;

	
	
	for(int i = 0;i<grid;i++){

		value = Newton(function,a,b,n,X[i]);
		fprintf(out,"%lf %lf\n",X[i],value);
	}
	
	return max;
}



double Newton_show_optimal(double (*function)(double),double a,double b,int n){
	FILE* out;
	out = fopen("Newton/Newton.txt","w");
	
	
	int grid = 1000;
	double h = (b-a)/grid;
	grid++;
	std::vector<double> X;
	for(int i = 0;i<grid;i++) X.push_back(a+i*h);
	
	double max = 0;
	double value;
	
	
	
	double h_interpolation = (b-a)/n;
	n++;
	std::vector<double> X_interpolation;
	for(int i = 0;i<n;i++){
		X_interpolation.push_back(a+h_interpolation*i);
	}
	std::vector<double> Coeffisients;
	Newton_coeffs(function,X_interpolation,Coeffisients);
	
	for(int i = 0;i<grid;i++){

		value = Newton_optimal(function,X_interpolation,Coeffisients,X[i]);

		//std::cout<<"point x="<<X[i]<<std::endl;
		fprintf(out,"%lf %lf\n",X[i],value);
	
		//if(value>max) max = value;
	}
	
	return max;
}





/*
int task3(){
	//FILE* out = fopen("Newton.txt","w");
	int n0 = 2;
	double value;
	Newton_show_optimal(expression,0,2,n0);
	
	FILE* out = fopen("Newton/max_difference.txt","w");
	double max_difference[15];
	int min_i = 1;
	//max_difference[0] = 0;
	
	
	for(int n = 1;n<16;n++){
		max_difference[n-1] = Newton_max_deviation(expression,0,2,n);
		printf("n = %d\t max_dif = %lf\n",n,max_difference[n-1]);
		fprintf(out,"%d %lf\n",n,max_difference[n-1]);
		if(max_difference[n-1]<max_difference[min_i]) min_i = n-1;
	}
	printf("optimal n is %d\n",min_i);
	return 0;
	
}
int main(){
	//int n = 10;
	//Newton_show_optimal(expression,0,2,n);
	task3();
	return 0;
}
*/