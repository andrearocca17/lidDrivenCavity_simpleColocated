#pragma once
#include "forLoop.h"
#include <vector>
#include <string>
using std::string;
using std::vector;
void pressureGradSource(vector<double> IFx, vector<double> IFy, vector<vector<double>>& p, vector<vector<double>>& Bu, vector<vector<double>>& Bv, double dx, double dy) {
    forAllInternal(p) {

        double pressureEastFace = (p[i + 1][j] * IFx[i]) + (p[i][j] * (1.0 - IFx[i]));
        double pressureWestFace = (p[i][j] * IFx[i-1]) + (p[i - 1][j] * (1.0 - IFx[i-1]));
        double pressureNorthFace = (p[i][j + 1] * IFy[j]) + (p[i][j] * (1.0 - IFy[j]));
        double pressureSouthFace = (p[i][j] * IFy[j-1]) + (p[i][j - 1] * (1.0 - IFy[j - 1]));


        double pressureEastGrad = (pressureEastFace - pressureWestFace) / dx;
        double pressureNorthGrad = (pressureNorthFace - pressureSouthFace) / dy;
     //   cout << pressureEastGrad << " ";

        // il gradiente di pressione lavora sulle source.

            Bu[i][j] += -pressureEastGrad * dx * dy;
        

            Bv[i][j] += -pressureNorthGrad * dx * dy;

    }   // end forAllInternal
}