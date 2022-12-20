#pragma once


#include "forLoop.h"
#include <vector>
#include <cmath>
using std::vector;
void convection(double IFi, vector<vector<double>>& AE, vector<vector<double>>& AW,
    vector<vector<double>>& AN, vector<vector<double>>& AS, vector<vector<double>>& F1, vector<vector<double>>& F2, vector<vector<double>>& velocityComponent, vector<vector<double>>& B) {

    double blendFactor = 1 ;
    //east
    forAllInternalVCVs(velocityComponent) {
        AE[i][j]     += std::min(F1[i][j], 0.0);
        AW[i + 1][j] += -std::max(F1[i][j], 0.0);
        double UDS    = std::max(F1[i][j], 0.0) * velocityComponent[i][j] + std::min(F1[i][j], 0.0) * velocityComponent[i + 1][j];
        double CDS    = F1[i][j] * (velocityComponent[i + 1][j] * IFi + velocityComponent[i][j] * (1- IFi));
        B[i][j]      += blendFactor * (UDS - CDS);
        B[i + 1][j]  += -blendFactor * (UDS - CDS);
    }

    //north
    forAllInternalUCVs(velocityComponent) {
        AN[i][j]     += std::min(F2[i][j], 0.0);
        AS[i][j + 1] += -std::max(F2[i][j], 0.0);
        double UDS    = std::max(F2[i][j], 0.0) * velocityComponent[i][j] + std::min(F2[i][j], 0.0) * velocityComponent[i][j + 1];
        double CDS    = F2[i][j] * (velocityComponent[i][j + 1] * IFi + velocityComponent[i][j] * (1-IFi));
        B[i][j]      += blendFactor * (UDS - CDS);
        B[i][j+1]    += -blendFactor * (UDS - CDS);
    }
}