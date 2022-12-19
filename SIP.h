#pragma once


vector<vector<double>>  SIP(vector<vector<double>>& phi, double alpha, vector<vector<double>>& AS, vector<vector<double>>& AN, vector<vector<double>>& AW, vector<vector<double>>& AE, vector<vector<double>>& AP, vector<vector<double>>& Brelaxed, int& iterations, string EqnName)
{
    vector<vector<double>> phitemp(phi.size(), vector<double>(phi[0].size()));
    vector<vector<double>> LW(phi.size(), vector<double>(phi[0].size()));
    vector<vector<double>> LS(phi.size(), vector<double>(phi[0].size()));
    vector<vector<double>> UN(phi.size(), vector<double>(phi[0].size()));
    vector<vector<double>> UE(phi.size(), vector<double>(phi[0].size()));
    vector<vector<double>> LPR(phi.size(), vector<double>(phi[0].size()));
    vector<vector<double>> RES(phi.size(), vector<double>(phi[0].size()));

    double Residual = 0.0;
    double value = 0.0;
    double RSM = 0.0;
    double RESOR = 0;
    double rURF;
    
    double SOR = 0.2;


    //..COEFFICIENTS OF UPPER AND LOWER TRIANGULAR MATRICES
    for (unsigned int i = 1; i < phi.size() - 1; i++)
    {
        for (unsigned int j = 1; j < phi[i].size() - 1; j++)
        {
            LW[i][j] = AW[i][j] / (1.0 + (alpha * UN[i - 1][j]));
            LS[i][j] = AS[i][j] / (1.0 + (alpha * UE[i][j - 1]));
            double P1 = alpha * LW[i][j] * UN[i - 1][j];
            double P2 = alpha * LS[i][j] * UE[i][j - 1];
            LPR[i][j] = 1.0 / (AP[i][j] + P1 + P2 - LW[i][j] * UE[i - 1][j] - LS[i][j] * UN[i][j - 1]);
            UN[i][j] = (AN[i][j] - P1) * LPR[i][j];
            UE[i][j] = (AE[i][j] - P2) * LPR[i][j];
        }
    }


    int Liter = iterations;
    //.....CALCULATE RESIDUAL AND OVERWRITE IT BY INTERMEDIATE VECTOR
    for (int L = 0; L < Liter; L++)
    {
        Residual = 0.0;

        for (unsigned int i = 1; i < phi.size() - 1; i++)
        {
            for (unsigned int j = 1; j < phi[i].size() - 1; j++)
            {
                RES[i][j] = Brelaxed[i][j] - (AN[i][j] * phi[i][j + 1]) - (AS[i][j] * phi[i][j - 1]) -
                    (AE[i][j] * phi[i + 1][j]) - (AW[i][j] * phi[i - 1][j]) - (AP[i][j] * phi[i][j]);
                Residual += abs(RES[i][j]);
                RES[i][j] = (RES[i][j] - (LS[i][j] * RES[i][j - 1]) - (LW[i][j] * RES[i - 1][j])) * LPR[i][j];
                //  cout << std::setprecision(3)  << RES[i][j] << ' ';

            }
            //     cout << endl;
        }

        double small = 1e-20;
        if (L == 0)
        {
            RESOR = Residual;
        }
        RSM = Residual / (RESOR + small);


        cout << EqnName << " Inner It: " << L+1 << " and Residual --> " << Residual << " RSM " << RSM << endl;

        // Back Subsitution and Correction
        for (unsigned int i = phi.size() - 2; i >= 1; --i)
        {
            for (unsigned int j = phi[i].size() - 2; j >= 1; --j)
            {
                RES[i][j] = RES[i][j] - (UN[i][j] * RES[i][j + 1]) - (UE[i][j] * RES[i + 1][j]);
                phi[i][j] = phi[i][j] + RES[i][j];
            }
        }

    }

    for (int i = 0; i < phi.size(); i++) {
        for (int j = 0; j < phi[1].size(); j++) {
            phitemp[i][j] = phi[i][j];
        }
    }
    return phitemp;
}
