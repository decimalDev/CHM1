//var 1 y=exp(x)/(1+x^2) [0,2] delta = 5*10^(-4)

#include<iostream>
#include<vector>

double expression(double x) {
	return exp(x) / (1 + x * x);
}

void Josepf_Lui_Lagrange(double (*function)(double),int a,int b,int n,std::vector<double> &A){
	double h = (b-a)/n;
	std::vector<double> X;
	for(int i = 0;i<n;i++) X.push_back(a+i*h);
	
	//std::vector<double> Coefs;
	//for(int i = 0;i<n;i++) Coefs.push_back( function(X[i]) );
	
		double a = 1;
		for(int i = 0;i<n;i++){
			a = 1;
			for(int j = 0;j<n;j++){
				if(i==j) continue;
				a *= 1/(X[i]-X[j]);
			}
			A.push_back(function(X[i])*a);
		}	
}
double L_function(std::vector<double> &A, double a,double b,std::vector<double> &W,std::vector<double> &Y){
	int n = A.size();
	
	double h_0 = (b-a)/n;
	
	std::vector<double> X_0;
	
	for(int i = 0;i<n;i++) X_0.push_back(a+i*h_0);
	
	
	for(int t = 0;t<W.size();t++){
		double result = 0;
		for(int i = 0;i<n;i++){
			double result_i = 1;
				for(int j = 0;j<n;j++){
				if(i==j) continue;
				result_i*=(W[t]-X_0[j]);
			}
		result+=A[i]*result_i;
		}
		Y.push_back(result);
	}
}

int main(){
	double h = (b-a)/(n+1);//вторая равномерная сетка для пункта 1
	std::vector<double> X;
	for(int i = 0;i<n;i++) X.push_back(a+i*h);
	
	h = (b-a)/100000;
	
	return 0;
}