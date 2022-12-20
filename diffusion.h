#pragma once
#include "forLoop.h"
#include <vector>
using std::vector;
void diffusion(double mu, double dy, double dx, vector<vector<double>>& AE, vector<vector<double>>& AW,
    vector<vector<double>>& AN, vector<vector<double>>& AS, vector<vector<double>>& x, vector<vector<double>>& y) {
    forAllInternalUCVs(x) {
        AE[i][j]     += -mu * dy / (x[i + 1][j] - x[i][j]);
        AW[i + 1][j] += AE[i][j];
    }

    //diffusion U
    //north
    forAllInternalVCVs(y) {
        AN[i][j]     += -mu * dx / (y[i][j + 1] - y[i][j]);
        AS[i][j + 1] += AN[i][j];
    }
}
