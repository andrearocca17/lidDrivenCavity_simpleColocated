#pragma once
#include "forLoop.h"
#include <vector>
#include <string>
using std::string;
using std::vector;

void assembleMatrixUV(vector<vector<double>>&  AP, vector<vector<double>>& APb, vector<vector<double>>& AE,
    vector<vector<double>>& AW, vector<vector<double>>& AS, vector<vector<double>>& AN) {


    // ap final
    forAllInternal(AP) {
        AP[i][j] += APb[i][j] - AE[i][j] - AW[i][j] - AS[i][j] - AN[i][j];
    }
}