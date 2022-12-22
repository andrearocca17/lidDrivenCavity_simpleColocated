#pragma once
#include "forLoop.h"
#include <vector>
using std::vector;
void diffusion(double mu, double dy, double dx, vector<vector<double>>& AE, vector<vector<double>>& AW,
    vector<vector<double>>& AN, vector<vector<double>>& AS, vector<double>& xc, vector<double>& yc) {
    forAllInternalUCVs(AE) {
        AE[i][j]     += -mu * dy / (xc[i + 1] - xc[i]);
        AW[i + 1][j] += AE[i][j];
    }

    //north
    forAllInternalVCVs(AN) {
        AN[i][j]     += -mu * dx / (yc[j + 1] - yc[j]);
        AS[i][j + 1] += AN[i][j];
    }
}
