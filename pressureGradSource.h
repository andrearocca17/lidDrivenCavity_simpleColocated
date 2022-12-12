#pragma once
#include "forLoop.h"
#include <vector>
#include <string>
using std::string;
using std::vector;
void pressureGradSource(double IFi, vector<vector<double>>& p, vector<vector<double>> Bu, vector<vector<double>>& Bv, double dx, double dy) {
    forAllInternal(p) {

        double pressureEastFace = (p[i + 1][j] * IFi) + (p[i][j] * (1.0 - IFi));
        double pressureWestFace = (p[i][j] * IFi) + (p[i - 1][j] * (1.0 - IFi));
        double pressureNorthFace = (p[i][j + 1] * IFi) + (p[i][j] * (1.0 - IFi));
        double pressureSouthFace = (p[i][j] * IFi) + (p[i][j - 1] * (1.0 - IFi));


        double pressureEastGrad = (pressureEastFace - pressureWestFace) / dx;
        double pressureNorthGrad = (pressureNorthFace - pressureSouthFace) / dy;

        // il gradiente di pressione lavora sulle source.

            Bu[i][j] += -pressureEastGrad * dx * dy;
        

            Bv[i][j] += -pressureNorthGrad * dx * dy;
        

    }   // end forAllInternal
}