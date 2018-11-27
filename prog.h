#include <iostream>
#include <vector>

using namespace std;


vector<double> prog(vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &d, int n)
{
    vector<double> p, q, x;
    x.resize(n);
    p.resize(n);
    q.resize(n);
    p[0] = - c[0] / b[0];
    q[0] = d[0] / b[0];
    for(int i = 1; i < n - 1; i++) {
        p[i] = - c[i] / (b[i] + a[i] * p[i - 1]);
        q[i] = (d[i] - a[i] * q[i- 1]) / (b[i] + a[i] * p[i - 1]);
    }
    p[n - 1] = 0;
    q[n - 1] = (d[n - 1] - a[n - 1] * q[n - 2]) / (b[n - 1] + a[n - 1] * p[n - 2]);
    x[n - 1] = q[n - 1];
    for(int i = n - 2; i >= 0; i--) {
        x[i] = p[i] * x[i + 1] + q[i];
    }
    return x;
}


vector<double> prog_solution(vector<vector<double> > &m, vector<double> &d)
{
    size_t size = m[0].size();
    if (size == 0) {
        return m[0];
    }
    vector<double> a, b, c;
    a.resize(size);
    b.resize(size);
    c.resize(size);
    b[0] = m[0][0];
    c[0] = m[0][1];
    for(size_t i = 1; i < size - 1; i++) {
        a[i] = m[i][i - 1];
        b[i] = m[i][i];
        c[i] = m[i][i + 1];
    }
    a[size - 1] = m[size - 1][size - 2];
    b[size - 1] = m[size - 1][size - 1];
    return prog(a, b, c, d, size);
}