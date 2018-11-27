#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>


using namespace std;
double Pi  =3.141592653589793238463;

double tau = 0;

double func(double x)
{
	return exp(-4*Pi*Pi*tau) * sin(2.0*Pi*x);
}

double func1(double x)
{
    return exp(-0.5 * tau) * sin(x);
}


int setting(void)
{
    ifstream t;
    t.open("data/tau.txt");
    t>>tau;
    t.close();
	ofstream in;
    ofstream in1;
    in.open("data/t.txt");
    in1.open("data/t1.txt");
    for(double x = 0.0; x <= 1.0; x += 0.001) {
    	in<<x<<endl;
    	in1<<func(x)<<endl;
    }
    in.close();
    in1.close();
	return 0;
}

int setting1(void)
{
    ifstream t;
    t.open("data/tau7.txt");
    t>>tau;
    t.close();
    ofstream in;
    ofstream in1;
    in.open("data/t7.txt");
    in1.open("data/t17.txt");
    for(double x = 0.0; x <= 3.14; x += 0.001) {
        in<<x<<endl;
        in1<<func1(x)<<endl;
    }
    in.close();
    in1.close();
    return 0;
}