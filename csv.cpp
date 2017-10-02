#include <fstream>

int main()
{
    std::ifstream infile("C:\Users\Jamil\Documents\BURST\CVC.txt");

    double a, b;
    while (infile >> a >> b)
    {
        cout << a << ", " << b << endl;
    }

return 0;
}