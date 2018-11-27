#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <Python.h>
#include "prog.h"
#include "t.h"

using namespace std;

const double PI  =3.141592653589793238463;


double u_x(double x)
{
	return sin(x);
}

double accurate_solve(double x, double t)
{
	return exp(-0.5*t)*sin(x);	
}
//запихнуть в заголовочник
void print_v(vector<double> &m)
{
	size_t n = m.size();
	for(size_t i = 0; i < n; i++) {
		cout<<m[i]<<endl;
	}
	cout<<endl;
}
void print_matr(vector<vector<double> > &m)
{
	size_t n = m[0].size();
	for(size_t i = 0; i < n; i++) {
		for(size_t j = 0; j < n; j++) {
			cout<<m[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
}
void printmatrix(int n, vector<vector<double> > &m, ofstream &in, int ind)
{
    for(int i = 0; i < n; i++) {
        in<<m[i][ind]<<endl;
    }
}
void printv(int n, vector<double> &m, ofstream &in)
{
    for(int i = 0; i < n; i++) {
        in<<m[i]<<endl;
    }
}

double f(double x, double t)
{
	return 0.5 * exp(-0.5 * t) * cos(x);
}

int ind = 5;
int n = 0;
int T = 0;
int appr = 0;

double calculate_error(vector<vector<double> > &m, vector<double> &x, vector<double> &t);

double yav()
{
	double h = PI / n;
	double tau = 0.05 / T;
	size_t count_of_t = (size_t) (0.05 / tau);
	vector<double> x(n + 1);
	vector<double> t(count_of_t + 1);
	cout<<x.size()<<endl;
	cout<<t.size()<<endl;
	for(size_t i = 0; i <= count_of_t; i++) {
		t[i] = tau * i;
	}
	for(size_t i = 0; i <= n; i++) {
		x[i] = h * i;
	}
	print_v(x);
	print_v(t);
	vector<vector<double> > u(n + 1);
	for(size_t i = 0; i < n + 1; i++) {
		u[i].resize(count_of_t + 1);
	}
	for(size_t i = 0; i < n + 1; i++) {
		u[i][0] = u_x(x[i]);
	}
	//print_matr(u);
	for(size_t j = 0; j < count_of_t; j++) {
		for(size_t i = 1; i < n; i++) {
			u[i][j + 1] = u[i][j] + (tau / (h * h))*(u[i + 1][j] - 2.0*u[i][j] + u[i - 1][j]) + (tau * f(x[i], t[j]));
		}
		if (appr == 0) {
			u[0][j + 1] = u[1][j + 1] - h*exp(-0.5 * t[j + 1]);
			u[n][j + 1] = u[n - 1][j + 1] - h*exp(-0.5 * t[j + 1]);
		}
		if (appr == 1) {
			u[0][j + 1] = (4 / 3.) * u[1][j + 1] - (1 / 3.) * u[2][j + 1] - (2 / 3.)*h*exp(-0.5*t[j + 1]);
			u[n][j + 1] = (2 / 3.)*h*(-exp(-0.5*t[j + 1])) - (1 / 3.)*u[n - 2][j + 1] + (4 / 3.)*u[n - 1][j + 1];
		}
		if (appr == 2) {
			double betta = (1 / h) + (h / (2 * tau));
			u[0][j + 1] = (1 / betta)*u[1][j + 1] + ((h*h) / (2*betta*tau))*u[0][j] + ((h*h) / (2*betta))*f(0, t[j + 1]) - (h / betta)*(exp(-0.5 * t[j + 1]));
			u[n][j + 1] = (1 / betta)*u[n - 1][j + 1] + ((h*h)/(2*tau*betta))*u[n][j] + ((h*h)/(2*betta))*f(x[n], t[j + 1]) + (h/betta)*(- h*exp(-0.5 * t[j + 1]));
		}
	}
	//print_matr(u);
	ofstream time;
	time.open("data/tau7.txt");
	time<<t[ind];
	time.close();
	setting1();
	ofstream in;
    ofstream in1;
    in.open("data/x7.txt");
    in1.open("data/y7.txt");
    printv(n + 1, x, in);
    printmatrix(n + 1, u, in1, ind);
    in.close();
    in1.close();
    double error = calculate_error(u, x, t);
    cout<<"Error "<<error<<endl;
	return error;
}


double neyav()
{
	double h = PI / n;
	double tau = 0.05 / T;
	size_t count_of_t = (size_t) (0.05 / tau);
	vector<double> x(n + 1);
	vector<double> t(count_of_t + 1);
	cout<<x.size()<<endl;
	cout<<t.size()<<endl;
	for(size_t i = 0; i <= count_of_t; i++) {
		t[i] = tau * i;
	}
	for(size_t i = 0; i <= n; i++) {
		x[i] = h * i;
	}
	print_v(x);
	print_v(t);
	vector<vector<double> > u(n + 1);
	for(size_t i = 0; i < n + 1; i++) {
		u[i].resize(count_of_t + 1);
	}
	for(size_t i = 0; i < n + 1; i++) {
		u[i][0] = u_x(x[i]);
	}
	vector<vector<double> > a(n + 1);
	for(size_t i = 0; i < n + 1; i++) {
		a[i].resize(n + 1);
	}
	double ksi = -(1 + ((2 * tau) / (h * h)));
	for(size_t j = 0; j < count_of_t; j++) {
		for(size_t i = 1; i < n; i++) {
			a[i][i - 1] = tau / (h * h);
			a[i][i] = ksi;
			a[i][i + 1] =  tau / (h * h);
		}
		vector<double> b(n + 1);
		if (appr == 0) {
			a[0][0] = -1;
			a[0][1] = 1;
			a[n][n - 1] = -1;
			a[n][n] = 1;
			b[0] = h*exp(-0.5*t[j + 1]);
			b[n] = -h*exp(-0.5*t[j + 1]);
		}
		if (appr == 1) {
			a[0][0] = -2;
			a[0][1] = 4 + ksi * ((h*h) / tau);
			a[n][n-1] = -4 - ksi * ((h*h) / tau);
			a[n][n] = 2;
			b[0] = 2*h*exp(-0.5*t[j + 1]) - ((h*h) / tau)*(u[1][j] + tau * f(x[1], t[j]));
			b[n] = -2*h*exp(-0.5*t[j + 1]) + ((h*h) / tau)*(u[n - 1][j] + tau * f(x[n - 1], t[j]));
		}
		if (appr == 2) {
			double betta = (1 / h) + (h / (2 * tau));
			a[0][0] = 1;
			a[0][1] = -1 / betta;
			b[0] = ((h*h)/(2*tau*betta))*u[0][j] + ((h*h)/ (2*betta))*f(0, t[j + 1]) - (h / betta)*exp(-0.5*t[j + 1]);
			a[n][n - 1] = -(1 / betta);
			a[n][n] = 1;
			b[n] =  ((h*h)/(2*tau*betta))*u[n][j] + ((h*h)/ (2*betta))*f(x[n], t[j + 1]) - (h / betta)*exp(-0.5*t[j + 1]);
		}
 		print_matr(a);
		for(size_t i = 1; i < n; i++) {
			b[i] = -(u[i][j] + tau * f(x[i], t[j]));
		}
		print_v(b);
		vector<double> solution = prog_solution(a, b);
		for(size_t i = 0; i < n + 1; i++) {
			u[i][j + 1] = solution[i];
		}
		for(size_t i = 0; i < n + 1; i++) {
			for(size_t z = 0; z < n + 1; z++) {
				a[i][z] = 0;
			}
		}
	}
	ofstream time;
	time.open("data/tau7.txt");
	time<<t[ind];
	time.close();
	setting1();
	ofstream in;
    ofstream in1;
    in.open("data/x17.txt");
    in1.open("data/y17.txt");
    if (!in.is_open()) {
    	cout<<"Not open!!!!!!!"<<endl;
    }
    printv(n + 1, x, in);
    printmatrix(n + 1, u, in1, ind);
    in.close();
    in1.close();
    double error = calculate_error(u, x, t);
    cout<<"Error "<<error<<endl;
	return error;

}


double krank_nikolson()
{
	double h = PI / n;
	double tau = 0.05 / T;
	size_t count_of_t = (size_t) (0.05 / tau);
	size_t count_of_x = n;
	vector<double> x(n + 1);
	vector<double> t(count_of_t + 1);
	cout<<x.size()<<endl;
	cout<<t.size()<<endl;
	for(size_t i = 0; i <= count_of_t; i++) {
		t[i] = tau * i;
	}
	for(size_t i = 0; i <= n; i++) {
		x[i] = h * i;
	}
	print_v(x);
	print_v(t);
	vector<vector<double> > u(n + 1);
	for(size_t i = 0; i < n + 1; i++) {
		u[i].resize(count_of_t + 1);
	}
	for(size_t i = 0; i < n + 1; i++) {
		u[i][0] = u_x(x[i]);
	}
	vector<vector<double> > a(n + 1);
	for(size_t i = 0; i < n + 1; i++) {
		a[i].resize(n + 1);
	}
	double r = tau / (2 * h * h);
	for(size_t j = 0; j < count_of_t; j++) {
		for(size_t i = 1; i < n; i++) {
			a[i][i - 1] = -r;
			a[i][i] = 1 + 2 * r;
			a[i][i + 1] =  -r;
		}
		vector<double> b(n + 1);
		if (appr == 0) {
			a[0][0] = -1;
			a[0][1] = 1;
			a[n][n - 1] = -1;
			a[n][n] = 1;
			b[0] = h*exp(-0.5*t[j + 1]);
			b[n] = -h*exp(-0.5*t[j + 1]);
		}
		if (appr == 1) {
			a[0][0] = -2;
			a[0][1] = 2 - (1 / r);
			a[n][n - 1] = -2 + (1 / r);
			a[n][n] = 2;
			b[0] = 2*h*exp(-0.5*t[j + 1]) - (1 / r) * (r*u[2][j] + (1 - 2*r)*u[1][j] +r*u[0][j] + tau * f(x[1], t[j] + (tau * 0.5)));
			b[n] = -2*h*exp(-0.5*t[j + 1]) + (1 / r) * (r*u[n][j] + (1 - 2*r)*u[n - 1][j] +r*u[n - 2][j] + tau * f(x[n - 1], t[j] + (tau * 0.5)));
		}
		if(appr == 2) {
			double betta = (1 / h) + (h / (2 * tau));
			a[0][0] = 1;
			a[0][1] = -1 / betta;
			b[0] = ((h*h)/(2*tau*betta))*u[0][j] + ((h*h)/ (2*betta))*f(0, t[j + 1]) - (h / betta)*exp(-0.5*t[j + 1]);
			a[n][n - 1] = -(1 / betta);
			a[n][n] = 1;
			b[n] =  ((h*h)/(2*tau*betta))*u[n][j] + ((h*h)/ (2*betta))*f(x[n], t[j + 1]) - (h / betta)*exp(-0.5*t[j + 1]);
		}
 		print_matr(a);
		for(size_t i = 1; i < n; i++) {
			b[i] = r*u[i + 1][j] + (1 - 2*r)*u[i][j] +r*u[i - 1][j] + tau * f(x[i], t[j] + (tau * 0.5)); 
		}
		vector<double> solution = prog_solution(a, b);
		for(size_t i = 0; i < n + 1; i++) {
			u[i][j + 1] = solution[i];
		}
		for(size_t i = 0; i < n + 1; i++) {
			for(size_t z = 0; z < n + 1; z++) {
				a[i][z] = 0;
			}
		}
	}
	ofstream time;
	time.open("data/tau7.txt");
	time<<t[ind];
	time.close();
	setting1();
	ofstream in;
    ofstream in1;
    in.open("data/x27.txt");
    in1.open("data/y27.txt");
    if (!in.is_open()) {
    	cout<<"Not open!!!!!!!"<<endl;
    }
    printv(n + 1, x, in);
    printmatrix(n + 1, u, in1, ind);
    in.close();
    in1.close();
    double error = calculate_error(u, x, t);
    cout<<"Error "<<error<<endl;
	return error;
}
/*
double calculate_error(vector<vector<double> > &m, vector<double> &x, vector<double> &t)
{
	vector<double> res(n + 1);
	for(size_t i = 0; i < n + 1; i++) {
		res[i] = accurate_solve(x[i], t[ind]) - m[i][ind];
	}
	double result = 0;
	for(size_t i = 0; i < n; i++) {
		result += res[i]*res[i];
	}
	return sqrt(result);
}
*/
double calculate_error(vector<vector<double> > &m, vector<double> &x, vector<double> &t)
{
	vector<double> res(n + 1);
	for(size_t i = 0; i < n + 1; i++) {
		res[i] = abs(accurate_solve(x[i], t[ind]) - m[i][ind]);
	}
	double result = 0;
	for(size_t i = 0; i < n; i++) {
		if (res[i] > result) {
			result = res[i];
		}
	}
	return result;
}





int main(int argc, char **argv)
{
	const char* script_path = "pp.py";
	wchar_t **changed_argv;
	changed_argv = (wchar_t**)malloc((argc)* sizeof(*changed_argv));
	for(int i = 0; i < argc; i++)
	{
    	changed_argv[i] =(wchar_t*) malloc(strlen(argv[i]) + 1);
    	mbstowcs(changed_argv[i], argv[i], strlen(argv[i]) + 1);
	}
	Py_Initialize();
	PySys_SetArgv(argc, changed_argv);
    PyRun_SimpleFile(fopen(script_path, "r"), script_path);
    for(int i = 0; i < argc; i++) {
    	free(changed_argv[i]);
    }
    free(changed_argv);
    ifstream data;
    data.open("data/pp.txt");
	int i = 0;
	data>>n;
	data>>T;
	data>>i;
	data>>appr;
	data>>ind;
	data.close();
	cout<<"i "<<i<<endl;
	cout<<n<<endl;
	cout<<T<<endl;
	cout<<endl;
	script_path = "47.py";
	if (i == 0) {
		yav();
	}
	if (i == 1) {
		neyav();
		script_path = "471.py";
	}
	if (i == 2) {
		krank_nikolson();
		script_path = "472.py";
	}
	cout<<script_path<<endl;
    PyRun_SimpleFile(fopen(script_path, "r"), script_path);
    Py_Finalize();
    n = 30;
    T = 30;
    ind = 6;
    appr = 0;
    double error11 = yav();
    double error12 = neyav();
    double error13 = krank_nikolson();
    appr = 1;
    double error21 = yav();
    double error22 = neyav();
    double error23 = krank_nikolson();
    appr = 2;
    double error31 = yav();
    double error32 = neyav();
    double error33 = krank_nikolson();
    cout<<error11<<" "<<error21<<" "<<error31<<endl;
    cout<<error12<<" "<<error22<<" "<<error32<<endl;
    cout<<error13<<" "<<error23<<" "<<error33<<endl;
   	return 0;
}