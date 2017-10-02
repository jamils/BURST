#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <fstream>
#include <algorithm>

#define _USE_MATH_DEFINES

using namespace std;

int main()
{
    ifstream infile("CVC.txt");

    double a, b;
    while (infile >> a >> b)
    {
        cout << a << ", " << b << endl;
    }

return 0;
}