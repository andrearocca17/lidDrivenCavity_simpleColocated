#pragma once
#include "forLoop.h"
#include <vector>
#include <cmath>
using std::vector;
void convection(vector<double> IFx, vector<double> IFy, vector<vector<double>>& AE, vector<vector<double>>& AW,
    vector<vector<double>>& AN, vector<vector<double>>& AS, vector<vector<double>>& F1, vector<vector<double>>& F2, vector<vector<double>>& velocityComponent, vector<vector<double>>& B) {

    double blendFactor = 1 ;
    //east
    forAllInternalUCVs(velocityComponent) {
        AE[i][j]      +=  std::min(F1[i][j], 0.0);
        AW[i + 1][j]  += -std::max(F1[i][j], 0.0);
        double UDS    =  std::max(F1[i][j], 0.0) * velocityComponent[i][j] + std::min(F1[i][j], 0.0) * velocityComponent[i + 1][j];
        double CDS    =  F1[i][j] * (velocityComponent[i + 1][j] * IFx[i] + velocityComponent[i][j] * (1 - IFx[i]));
        B[i][j]       +=  blendFactor * (UDS - CDS);
        B[i + 1][j]   += -blendFactor * (UDS - CDS);
    }

    //north
    forAllInternalVCVs(velocityComponent) {
        AN[i][j]     += std::min(F2[i][j], 0.0);
        AS[i][j + 1] += -std::max(F2[i][j], 0.0);
        double UDS    = std::max(F2[i][j], 0.0) * velocityComponent[i][j] + std::min(F2[i][j], 0.0) * velocityComponent[i][j + 1];
        double CDS    = F2[i][j] * (velocityComponent[i][j + 1] * IFy[j] + velocityComponent[i][j] * (1-IFy[j]));
        B[i][j]      += blendFactor * (UDS - CDS);
        B[i][j+1]    += -blendFactor * (UDS - CDS);
    }
}