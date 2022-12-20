#pragma once
#include<vector>
#include<string>
using std::vector;
using std::string;
void extrapolateZeroGrad(vector<vector<double>>& p, double IFb, string boundary) {
    if (boundary == "north")
    {
        int j = p[1].size() - 1;
        for (int i = 1; i < p.size() - 1; i++) {
            p[i][j] = p[i][j - 1] + (1-IFb) * (p[i][j - 1] - p[i][j - 2]);
            //p[i][j] = p[i][j - 1];
        }
    }
    if (boundary == "south")
    {
        int j = 0;
        for (int i = 1; i < p.size() - 1; i++) {
            p[i][j] = p[i][j + 1] +IFb * (p[i][j + 1] - p[i][j + 2]);

        }
    }
    if (boundary == "east")
    {
        int i = p.size() - 1;
        for (int j = 1; j < p.size() - 1; j++) {
            p[i][j] = p[i - 1][j] +IFb * (p[i - 1][j] - p[i - 2][j]);
        }
    }
    if (boundary == "west")
    {
        int i = 0;
        for (int j = 1; j < p.size() - 1; j++) {
            p[i][j] = p[i + 1][j] +IFb * (p[i + 1][j] - p[i + 2][j]);
        }
    }
}




 

