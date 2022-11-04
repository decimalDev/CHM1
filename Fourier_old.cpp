//var 1 y=exp(x)/(1+x^2) [0,2] delta = 5*10^(-4)

#include<iostream>
#include<vector>
#include<math.h>

#define INTEGRATION_DIVIDE 300

double expression(double x) {
	return exp(x) / (1 + x * x);
}




double dot_product_cos(double (*function)(double),std::vector<double> &X,int k){
	int n = X.size();
	int n0 = 100;
	double a0 = X[0];
	double b0 = X[n-1];
	double h0 = (b0-a0)/n0;
	n0++;
	std::vector<double> X0;
	for(int i = 0;i<n0;i++){
		X0.push_back(a0+h0*i);
	}
	
	
	double result = 0;
	for(int i = 0;i<n0-1;i++){
		result+=(function(X0[i])*cos(k*X0[i]) + function(X0[i+1])*cos(k*X0[i+1]))*h0/2;
		//result+=function(X0[i])*cos(k*X0[i]);
	}
	return result/M_PI;
}

double dot_product_sin(double (*function)(double),std::vector<double> &X,int k){
	int n = X.size();
	int n0 = 100;
	double a0 = X[0];
	double b0 = X[n-1];
	double h0 = (b0-a0)/n0;
	n0++;
	std::vector<double> X0;
	for(int i = 0;i<n0;i++){
		X0.push_back(a0+h0*i);
	}
	
	
	double result = 0;
	for(int i = 0;i<n0-1;i++){
		result+=(function(X0[i])*sin(k*X0[i]) + function(X0[i+1])*sin(k*X0[i+1]))*h0/2;
		//result+=function(X0[i])*sin(k*X0[i]);
	}
	return result/M_PI;
}



///*

double dot_product_cos_with_displacement(double (*function)(double),std::vector<double> &X,int k,double displacement){
	int n = X.size();
	int n0 = INTEGRATION_DIVIDE;
	double a0 = X[0];
	double b0 = X[n-1];
	double h0 = (b0-a0)/n0;
	n0++;
	std::vector<double> X0;
	for(int i = 0;i<n0;i++){
		X0.push_back(a0+h0*i);
	}
	
	/*
	double result = 0;
	for(int i = 0;i<n0-1;i++){
		result+=(function(X0[i]+displacement)*cos(k*X0[i]*M_PI/(b0-a0)) + function(X0[i+1]+displacement)*cos(k*X0[i+1]*M_PI/(b0-a0)))*h0/2;
		//result+=function(X0[i])*cos(k*X0[i]);
	}
	return result/(b0-a0);
	*/
	
	double result = 0;
	for(int i = 1;i<n0-1;i++){
		result+=(function(X0[i]+displacement)*cos(k*X0[i]*M_PI/(b0-a0)))*h0;
		//result+=function(X0[i])*cos(k*X0[i]);
	}
	result+=(function(X0[0]+displacement)*cos(k*X0[0]*M_PI/(b0-a0)) + function(X0[n0-1]+displacement)*cos(k*X0[n0-1]*M_PI/(b0-a0)))*h0/2;
	return result/(b0-a0);
	
	
}



double dot_product_sin_with_displacement(double (*function)(double),std::vector<double> &X,int k,double displacement){
	int n = X.size();
	int n0 = INTEGRATION_DIVIDE;
	double a0 = X[0];
	double b0 = X[n-1];
	double h0 = (b0-a0)/n0;
	n0++;
	std::vector<double> X0;
	for(int i = 0;i<n0;i++){
		X0.push_back(a0+h0*i);
	}
	
	/*
	double result = 0;
	for(int i = 0;i<n0-1;i++){
		result+=(function(X0[i]+displacement)*sin(k*X0[i]*M_PI/(b0-a0)) + function(X0[i+1]+displacement)*sin(k*X0[i+1]*M_PI/(b0-a0)))*h0/2;
		//result+=function(X0[i])*cos(k*X0[i]);
	}
	return result/(b0-a0);
	*/
	
	double result = 0;
	for(int i = 1;i<n0-1;i++){
		result+=(function(X0[i]+displacement)*sin(k*X0[i]*M_PI/(b0-a0)))*h0;
		//result+=function(X0[i])*cos(k*X0[i]);
	}
	result+=(function(X0[0]+displacement)*sin(k*X0[0]*M_PI/(b0-a0)) + function(X0[n0-1]+displacement)*sin(k*X0[n0-1]*M_PI/(b0-a0)))*h0/2;
	return result/(b0-a0);
	
}

void Fourier_coeffs(double (*function)(double),std::vector<double> &X,std::vector<std::pair<double,double>> &Coeffisients){
	
	int n = X.size();
	Coeffisients.push_back(std::pair<double,double>{dot_product_cos(function,X,0)/2,dot_product_sin(function,X,0)/2});
	for(int i = 1;i<n;i++){
		Coeffisients.push_back(std::pair<double,double>{dot_product_cos(function,X,i),dot_product_sin(function,X,i)});
	}
	
}


void Fourier_coeffs_with_displacement(double (*function)(double),std::vector<double> &X,std::vector<std::pair<double,double>> &Coeffisients,double displacement){
	
	int n = X.size();
	Coeffisients.push_back(std::pair<double,double>{dot_product_cos_with_displacement(function,X,0,displacement)/2,dot_product_sin_with_displacement(function,X,0,displacement)/2});
	for(int i = 1;i<n;i++){
		Coeffisients.push_back(std::pair<double,double>{dot_product_cos_with_displacement(function,X,i,displacement),dot_product_sin_with_displacement(function,X,i,displacement)});
	}
	
}




double Fourier_value(double parametr_x,std::vector<double> &X,std::vector<std::pair<double,double>> &Coeffisients){
	int n = Coeffisients.size();
	
	double result = 0;
	
	for(int i = 0;i<n;i++){
		result += Coeffisients[i].first*cos(i*parametr_x) + Coeffisients[i].second*sin(i*parametr_x);
		//printf("res=%lf\n",result);
	}

	return result;
}


double Fourier_value_with_displacement(double parametr_x,std::vector<double> &X,std::vector<std::pair<double,double>> &Coeffisients,double displacement){
	int n = Coeffisients.size();
	
	double result = 0;
	
	for(int i = 0;i<n;i++){
		result += Coeffisients[i].first*cos(i*(parametr_x-displacement)*M_PI/(X[n-1]-X[0])) + Coeffisients[i].second*sin(i*(parametr_x-displacement)*M_PI/(X[n-1]-X[0]));
		//printf("res=%lf\n",result);
	}

	return result;
}




double Fourier(double (*function)(double),double a,double b,int n,double x){
	double h = (b-a)/n;
	n++;
	std::vector<double> X;
	for(int i = 0;i<n;i++){
		X.push_back(a+h*i);
	}
	std::vector<std::pair<double,double>> Coeffisients;
	Fourier_coeffs(function,X,Coeffisients);
	//std::cout<<"Isaac Newton is tired. And found number is "<<Isaac_Newton_is_calculating(x,Coeffisients,a,b,n)<<std::endl;
	double res = Fourier_value(x,X,Coeffisients);
	
	
	return res;
}



double Fourier_with_displacement(double (*function)(double),double a,double b,int n,double x){
	double h = (b-a)/n;
	double displacement = (a+b)/2;
	n++;
	std::vector<double> X;
	for(int i = 0;i<n;i++){
		X.push_back(a+h*i-displacement);
	}
	std::vector<std::pair<double,double>> Coeffisients;
	Fourier_coeffs_with_displacement(function,X,Coeffisients,displacement);
	//std::cout<<"Isaac Newton is tired. And found number is "<<Isaac_Newton_is_calculating(x,Coeffisients,a,b,n)<<std::endl;
	double res = Fourier_value_with_displacement(x,X,Coeffisients,displacement);
	
	
	return res;
}



double Fourier_with_displacement_optimal(double (*function)(double),std::vector<double> X,std::vector<std::pair<double,double>> &Coeffisients,double x,double displacement){
	

	//Fourier_coeffs_with_displacement(function,X,Coeffisients,displacement);
	//std::cout<<"Isaac Newton is tired. And found number is "<<Isaac_Newton_is_calculating(x,Coeffisients,a,b,n)<<std::endl;
	double res = Fourier_value_with_displacement(x,X,Coeffisients,displacement);
	
	
	return res;
}




double Fourier_show(double (*function)(double),double a,double b,int n){
	FILE* out;
	//if(write_n==n)
	out = fopen("Fourier/Fourier.txt","w");
	
	
	int grid = 1000;
	double h = (b-a)/grid;
	grid++;
	std::vector<double> X;
	for(int i = 0;i<grid;i++) X.push_back(a+i*h);
	
	double max = 0;
	double value;

	
	
	for(int i = 0;i<grid;i++){

		value = Fourier_with_displacement(function,a,b,n,X[i]);
		
		fprintf(out,"%lf %lf\n",X[i],value);
	}
	
	return max;
}


double Fourier_show_optimal(double (*function)(double),double a,double b,int n){
	FILE* out;
	//if(write_n==n)
	out = fopen("Fourier/Fourier.txt","w");
	

	int grid = 1000;
	double h = (b-a)/grid;
	grid++;
	std::vector<double> X;
	for(int i = 0;i<grid;i++) X.push_back(a+i*h);
	
	double max = 0;
	double value;

	double displacement = (a+b)/2;
	double h_interpolation = (b-a)/n;
	n++;
	std::vector<double> X_interpolation;
	for(int i = 0;i<n;i++){
		X_interpolation.push_back(a+h_interpolation*i-displacement);
	}
	std::vector<std::pair<double,double>> Coeffisients;
	Fourier_coeffs_with_displacement(function,X_interpolation,Coeffisients,displacement);
	
	for(int i = 0;i<grid;i++){
		
		
		value = Fourier_with_displacement_optimal(function,X_interpolation,Coeffisients,X[i],displacement);
		
		fprintf(out,"%lf %lf\n",X[i],value);
	}
	fclose(out);
	
	return max;
}



double Fourier_max_difference(double (*function)(double),double a,double b,int n){
	
	//FILE* out;
	//if(write_n==n)
	//out = fopen("Fourier/Fourier.txt","w");
	
	int grid = 100000;
	double h = (b-a)/grid;
	grid++;
	std::vector<double> X;
	for(int i = 0;i<grid;i++) X.push_back(a+i*h);
	
	double max = 0;
	double value;
	
	
	double displacement = (a+b)/2;
	double h_interpolation = (b-a)/n;
	n++;
	std::vector<double> X_interpolation;
	/*
	for(int i = -14;i<=n+13;i++){
		X_interpolation.push_back(a+h_interpolation*i-displacement);
	}
	*/
	
	for(int i = -2;i<=n+2;i++){
		X_interpolation.push_back(a+h_interpolation*i-displacement);
	}
	
	
	std::vector<std::pair<double,double>> Coeffisients;
	Fourier_coeffs_with_displacement(function,X_interpolation,Coeffisients,displacement);
	
	
	int i = 0;
	while(X[++i]<a);
	
	while(X[--grid]>b);
	
	for(;i<grid;i++){
		//printf("expr(x) = %lf and interpolation(x) = %lf\n",function(X[i]),Josepf_Lui_Lagrange_give_me_order(function,a,b,n,X[i]));
		value = Fourier_with_displacement_optimal(function,X_interpolation,Coeffisients,X[i],displacement);
		//fprintf(out,"%lf %lf\n",X[i],value);
		value -= function(X[i]);
		value = std::abs(value);
	
		if(value>max) max = value;
	}
	
	//fclose(out);
	
	return max;
}





int task4(){
	FILE* out = fopen("Fourier/max_difference.txt","w");
	double max_difference[99];
	int min_i = 1;
	//max_difference[0] = 0;
	double difference = 100;
	int n = INTEGRATION_DIVIDE/2;
	while(difference>0.0005&&n<500){
		difference = Fourier_max_difference(expression,0,2,n);
		printf("n = %d\t max_dif = %lf\n",n,difference);
		fprintf(out,"%d %lf\n",n,difference);
		n++;
	}
	printf("optimal is %d\n",n);
	fclose(out);
	return 0;
}


int main(){
	int n = 10;
	Fourier_show_optimal(expression,0,2,n);
	task4();
	return 0;
}