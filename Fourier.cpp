#include<iostream>
#include<cmath>
#define _USE_MATH_DEFINES

double expression(double x) {
	return exp(x) / (1 + x * x);
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
	for(int i = 0;i<n/2;i++){
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

int main(){
	return 0;
}