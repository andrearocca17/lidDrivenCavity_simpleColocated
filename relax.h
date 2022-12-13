#pragma once


void relaxEq(double URF, vector<vector<double>>&  AP, vector<vector<double>>& BRelaxed, vector<vector<double>>& B, vector<vector<double>>& velComp) {
    forAllInternal(AP){
            AP[i][j] = (1.0 / URF) * AP[i][j];
            BRelaxed[i][j] = B[i][j] + (1.0 - URF) * (AP[i][j] * velComp[i][j]);
        }
    
}