#pragma once
#include "forLoop.h"
#include <vector>
#include <string>
using std::string;
using std::vector;
void noWallShearY(double Nx, double Ny, double dx, double dy, double mu, double V_wall, vector<double>& xc, vector<double>& yc,
    vector<vector<double>>& APbv, vector<vector<double>>& Bv) {
    double dx1 = dx;
    double dy1 = yc[1] - yc[0];
    double dx2 = xc[1] - xc[0];
    double dy2 = dy;
    //cout << dx1 << dx2 << dy1 << dy2 << endl;
    // east boundary. no wall shear in Y direction. Solo V
    for (int j = 1; j < Ny - 1; j++) {
        APbv[1][j] += mu * dy2 / dx2;     // west boundary
        APbv[Nx - 2][j] += mu * dy2 / dx2;     // east boundary
        Bv[Nx - 2][j] += mu * dy2 / dx2 * V_wall;
        Bv[1][j] += mu * dy2 / dx2 * V_wall;
    }
}