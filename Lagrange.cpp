//var 1 y=exp(x)/(1+x^2) [0,2] delta = 5*10^(-4)

#include<iostream>
#include<vector>
#include<stdio.h>
#include<cmath>

double expression(double x) {
	return exp(x) / (1 + x * x);
}



void Josepf_Lui_Lagrange_is_thinking(double (*function)(double),std::vector<double> &X,std::vector<double> &Coeffisients){
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




double Josepf_Lui_Lagrange_is_calculating(double parametr_x,std::vector<double> &X,std::vector<double> &Coeffisients){
	int n = Coeffisients.size();
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

double Josepf_Lui_Lagrange_give_me_order(double (*function)(double),double a,double b,int n,double x){
	//double a = 0,b = 2;
	//int n = 10;
	double h = (b-a)/(n);
	n++;
	//double x = 0.5;
	std::vector<double> Coeffisients;
	std::vector<double> X;
	for(int i = 0;i<n;i++){
		X.push_back(a+h*i);
	}
	Josepf_Lui_Lagrange_is_thinking(function,X,Coeffisients);
	double res = Josepf_Lui_Lagrange_is_calculating(x,X,Coeffisients);
	
	//std::cout<<"Josepf Lui Lagrange is tired. But he could complete the calculation and found number is "<<res<<std::endl;
	
	return res;
}

int write_n;

double Josepf_Lui_Lagrange_max_difference(double (*function)(double),double a,double b,int n){

	
	int delta_n = 100000;
	double h = (b-a)/delta_n;
	delta_n++;
	std::vector<double> X;
	for(int i = 0;i<delta_n;i++) X.push_back(a+i*h);
	
	double max = 0;
	double value;
	//printf("exp(x) = %lf\n",expression());
	for(int i = 0;i<delta_n;i++){
		//printf("expr(x) = %lf and interpolation(x) = %lf\n",function(X[i]),Josepf_Lui_Lagrange_give_me_order(function,a,b,n,X[i]));
		value = function(X[i]) - Josepf_Lui_Lagrange_give_me_order(function,a,b,n,X[i]);
		value = std::abs(value);
	
		if(value>max) max = value;
	}
	
	return max;
}





int task1(){
	FILE* out = fopen("Lagrange/max_difference.txt","w");
	double max_difference[15];
	int min_i = 1;
	//max_difference[0] = 0;
	
	
	for(int n = 1;n<16;n++){
		max_difference[n-1] = Josepf_Lui_Lagrange_max_difference(expression,0,2,n);
		printf("n = %d\t max_dif = %lf\n",n,max_difference[n-1]);
		fprintf(out,"%d %lf\n",n,max_difference[n-1]);
		if(max_difference[n-1]<max_difference[min_i]) min_i = n-1;
	}
	printf("optimal n is %d\n",min_i);
	return 0;
}


double Josepf_Lui_Lagrange_show(double (*function)(double),double a,double b,int n){
	FILE* out;
	//if(write_n==n)
	out = fopen("Lagrange/Lagrange.txt","w");
	
	
	int delta_n = 1000;
	double h = (b-a)/delta_n;
	delta_n++;
	std::vector<double> X;
	for(int i = 0;i<delta_n;i++) X.push_back(a+i*h);
	
	double max = 0;
	double value;
	//printf("exp(x) = %lf\n",expression());
	for(int i = 0;i<delta_n;i++){
		//printf("expr(x) = %lf and interpolation(x) = %lf\n",function(X[i]),Josepf_Lui_Lagrange_give_me_order(function,a,b,n,X[i]));
		value = Josepf_Lui_Lagrange_give_me_order(function,a,b,n,X[i]);
	
		fprintf(out,"%lf %lf\n",X[i],value);
	
		//if(value>max) max = value;
	}
	
	return max;
}

int main(){
	int n = 10;
	Josepf_Lui_Lagrange_show(expression,0,2,n);
	task1();
	return 0;
}