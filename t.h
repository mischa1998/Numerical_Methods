#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
//для построения графика точного решения создаются файлы с координатами (x,y) - которые получены из точного аналитического решения

using namespace std;
double Pi  =3.141592653589793238463;

double tau = 0;

double func1(double x)
{
    return exp(-0.5 * tau) * sin(x);
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
