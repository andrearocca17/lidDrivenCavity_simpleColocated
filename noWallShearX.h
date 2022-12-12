#pragma once
#include "forLoop.h"
#include <vector>
#include <string>
using std::string;
using std::vector;
void noWallShearX(double Nx, double Ny, double dx, double dy, double mu, double U_north, double U_wall, vector<vector<double>>& x, vector<vector<double>>&y, 
    vector<vector<double>>& APbu, vector<vector<double>>& Bu) {
    double dx1 = dx;
    double dy1 = y[0][1] - y[0][0];
    double dx2 = x[1][0] - x[0][0];
    double dy2 = dy;
    // North boundary. no wall shear in X direction. Solo U
    for (int i = 1; i < Nx - 1; i++) {
        APbu[i][Ny - 2] += mu * dx1 / dy1;     // north boundary
        APbu[i][1] += mu * dx1 / dy1;     // South boundary
        Bu[i][Ny - 2] += mu * dx1 / dy1 * U_north;
        Bu[i][1] += mu * dx1 / dy1 * U_wall;
    }
}