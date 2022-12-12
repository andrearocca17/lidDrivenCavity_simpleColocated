#pragma once
#pragma once
// 15/11/2022 inizio di un nuovo progetto lid driven cavity. Obiettivo arrivare a lid driven cavity cuda
// finite volume method colocated



#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>

using std::vector;
using std::cout;
using std::endl;
using std::string;

void extaraPolateZeroGrad(vector<vector<double>>& p, double IFb, string boundary) {

    if (boundary == "north")
    {
        int j = p[1].size() - 1;
        for (int i = 1; i < p.size() - 1; i++) {
            p[i][j] = p[i][j - 1] + IFb * (p[i][j - 1] - p[i][j - 2]);
        }
    }

    if (boundary == "south")
    {
        int j = 0;
        for (int i = 1; i < p.size() - 1; i++) {
            p[i][j] = p[i][j + 1] + IFb * (p[i][j + 1] - p[i][j + 2]);
        }
    }


    if (boundary == "east")
    {
        int i = p.size() - 1;
        for (int j = 1; j < p.size() - 1; j++) {
            p[i][j] = p[i - 1][j] + IFb * (p[i - 1][j] - p[i - 2][j]);
        }
    }

    if (boundary == "west")
    {
        int i = 0;
        for (int j = 1; j < p.size() - 1; j++) {
            p[i][j] = p[i + 1][j] + IFb * (p[i + 1][j] - p[i + 2][j]);
        }
    }


}


void  printMat(vector<vector<double>>& mat) {
    int row = mat[0].size();
    int col = mat[1].size();
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            cout << mat[i][j] << " ";

        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------" << endl;
}

void  printMat(vector<vector<double>>& mat, string nameVar) {
    cout << "------------------------" << nameVar << "------------------------" << endl;

    int row = mat[0].size();
    int col = mat[1].size();
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            cout << mat[i][j] << " ";

        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------" << endl;
}

void  printVec(vector<double>& vec) {

    int row = vec.size();
    for (int i = 0; i < row; i++) {
        cout << vec[i] << endl;
    }
    cout << "---------------------------------------------------------------" << endl;
}

void luDecomposition(vector<vector<double>>& a, vector<vector<double>>& u, vector<vector<double>>& l, int n) {

    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (a[j][i] == 0)
                l[j][i] == 0;
            else if (j < i)
                l[j][i] = 0;
            else {
                l[j][i] = a[j][i];
                for (k = 0; k < i; k++) {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
        for (j = 0; j < n; j++) {

            if (a[i][j] == 0)
                u[i][j] = 0;

            else if (j < i)
                u[i][j] = 0;
            else if (j == i)
                u[i][j] = 1;
            else {
                u[i][j] = a[i][j] / l[i][i];
                for (k = 0; k < i; k++) {
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
                }
            }
        }
    }

}

vector<vector<double>>  SIP(vector<vector<double>>& phi, double alpha, vector<vector<double>>& AS, vector<vector<double>>& AN, vector<vector<double>>& AW, vector<vector<double>>& AE, vector<vector<double>>& AP, vector<vector<double>>& Brelaxed, int& iterations)
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
    string EqnName = "U-Momentum";
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


        cout << EqnName << " Inner It: " << L << " and Residual --> " << Residual << " RSM " << RSM << endl;

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


#define forAll(tempvector) \
    for(unsigned int i=0; i<tempvector.size(); i++) \
        for(unsigned int j=0; j<tempvector[i].size(); j++)
#define forAllInternal(tempvector) \
    for(unsigned int i=1; i<tempvector.size()-1; i++) \
        for(unsigned int j=1; j<tempvector[i].size()-1; j++)
#define forAllBackwardInternal(tempvector) \
    for(unsigned int i=tempvector.size()-2; i >=1; --i) \
        for(unsigned int j=tempvector[i].size()-2; j >=1; --j)
#define forAllInternalUCVs(tempvector) \
       for(unsigned int i=1; i<tempvector.size()-2; i++) \
        for(unsigned int j=1; j<tempvector[i].size()-1; j++)
#define forAllInternalVCVs(tempvector) \
       for(unsigned int j=1; j<tempvector.size()-2; j++) \
        for(unsigned int i=1; i<tempvector[j].size()-1; i++)
#define forAllInternalX(tempvector) \
    for(unsigned int i=1; i<tempvector.size()-1; i++) \
        for(unsigned int j=1; j<tempvector[i].size(); j++)
#define forAllInternalY(tempvector) \
    for(unsigned int i=1; i<tempvector.size(); i++) \
        for(unsigned int j=1; j<tempvector[i].size()-1; j++)
#define forSouthBoundary(tempvector) \
    for(unsigned int i=1; i<tempvector.size()-1; i++) \
        for(unsigned int j=0; j<1; j++)
#define forNorthBoundary(tempvector) \
    for(unsigned int i=1; i<tempvector.size()-1; i++) \
        for(unsigned int j=tempvector[i].size()-1; j<tempvector[i].size(); j++)
#define forWestBoundary(tempvector) \
  for(unsigned int j=1; j<tempvector[0].size(); j++) \
    for(unsigned int i=1; i<2; i++)
#define forEastBoundary(tempvector) \
  for(unsigned int j=1; j<tempvector[0].size(); j++) \
        for(unsigned int i=tempvector[j].size()-2; i<tempvector[j].size()-1; i++)
#define forAllOppInternal(tempvector) \
    for(unsigned int j=1; j<tempvector.size()-1; j++) \
        for(unsigned int i=1; i<tempvector[i].size()-1; i++)



int main() {
    double L = 1.0; //lunghezza del lato della cavità
    double mu = 0.001;
    double density = 1.0;
    double dt = 0.001;
    int Nx = 7, Ny = 7; // numero di nodi interni su x e su y 
    int nb = 2; //numero di boundaries
    int Nxb = Nx + 2;
    int Nyb = Ny + 2;
    double dx = L / Nx; //dimensione-x del volume di controllo
    double dy = L / Ny;
    int Nb = (Nx + 2) * (Ny + 2); //numero di nodi totali
    int N = Nx * Ny; //numero di incognite
    int Nfx = Nx + 1;
    int Nfy = Ny + 1;
    double IFi = 0.5; //interpolation factor interno
    double IFb = 1.0; //interpolation factor al bordo

    // 0 ----------------------- 
    // |          //west      |
    // |                      |
    // |                      |
    // |                      |north
    // |                      |
    // |                      |
    // |                      |
    // |                      |
    // |                      |
    // -----------------------
                //east

    // ^------> j
    // |
    // |       
    // |
    // \/ i   

    // oriento nel verso di stampa di un ciclo for. 

    vector<double> xe(N); // coordinata x centri delle facce est
    vector<double> ye(N); // coordinata y centri delle facce est
    vector<double> xn(N); // coordinata x centri delle facce nord
    vector<double> yn(N); // coordinata y centri delle facce nord
    vector<double> xw(N); // coordinata x centri delle facce ovest
    vector<double> yw(N); // coordinata y centri delle facce ovest
    vector<double> xs(N); // coordinata x centri delle facce sud
    vector<double> ys(N); // coordinata y centri delle facce sud

    vector<vector<double>> x(Nxb, vector<double>(Nyb));
    vector<vector<double>> y(Nxb, vector<double>(Nyb));

    vector<vector<double>> F1(Nx, vector<double>(Ny)); //flussi x
    vector<vector<double>> F2(Nx, vector<double>(Ny)); //flussi y

    vector<vector<double>> U(Nx, vector<double>(Ny));
    vector<vector<double>> p(Nxb, vector<double>(Nxb));


    vector<vector<double>> AE(Nx, vector<double>(Ny));
    vector<vector<double>> AW(Nx, vector<double>(Ny));
    vector<vector<double>> AS(Nx, vector<double>(Ny));
    vector<vector<double>> AN(Nx, vector<double>(Ny));
    vector<vector<double>> APtime(Nx, vector<double>(Ny));
    vector<vector<double>> APb(Nx, vector<double>(Ny));
    vector<vector<double>> AP(Nx, vector<double>(Ny));
    vector<vector<double>> B(Nx, vector<double>(Ny));
    vector<vector<double>> Brelaxed(Nx, vector<double>(Ny));
    vector<vector<double>> rAP(Nx, vector<double>(Ny));
    // source


    for (int i = 0; i < Nxb; i++) {
        for (int j = 0; j < Nyb; j++) {
            if (i == 0) { x[i][j] = 0; }
            else if (i == Nxb - 1) { x[i][j] = 1; }
            else x[i][j] = dx / 2 + (i - 1) * dx;
            //  y[n] = dy / 2 + j * dy;
        }
    }
    for (int i = 0; i < Nxb; i++) {
        for (int j = 0; j < Nyb; j++) {
            //  cout << x[i][j] << " ";

        }
        cout << endl;
    }
    cout << "---------------------------------------------------------------" << endl;

    for (int i = 0; i < Nxb; i++) {
        for (int j = 0; j < Nyb; j++) {
            if (j == 0) { y[i][j] = 0; }
            else if (j == Nyb - 1) { y[i][j] = 1; }
            else y[i][j] = dy / 2 + (j - 1) * dy;
            //  y[n] = dy / 2 + j * dy;
        }
    }
    for (int i = 0; i < Nxb; i++) {
        for (int j = 0; j < Nyb; j++) {
            //   cout << y[i][j] << " ";

        }
        // cout << endl;
    }
    cout << "---------------------------------------------------------------" << endl;

    //diffusion U
    //east
    for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 2; i++) {
            //    n = i + (Nx)*j;
            AE[i][j] += -mu * dy / (x[i + 1][j] - x[i][j]);
            AW[i + 1][j] += AE[i][j];
            // AN[i][j] += -mu * dx / (yn[n + 1] - yn[n]);
            // AS[i][j + 1] += AN[i][j];
        }
    }
    //diffusion U
    //north
    for (int j = 1; j < Ny - 2; j++) {
        for (int i = 1; i < Nx - 1; i++) {
            AN[i][j] += -mu * dx / (y[i][j + 1] - y[i][j]);
            AS[i][j + 1] += AN[i][j];
        }
    }
    // printMat(AE);
    // printVec(U);


     //convection U
    double blendFactor = 0.8;
    //east
    for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 2; i++) {
            int n = j * (Ny - 1) + i;
            AE[i][j] += std::min(F1[i][j], 0.0);
            AW[i + 1][j] += -std::max(F1[i][j], 0.0);
            double UDS = std::max(F1[i][j], 0.0) * U[i][j] + std::min(F1[i][j], 0.0) * U[i + 1][j];
            double CDS = F1[i][j] * (U[i + 1][j] * IFi + U[i][j] * IFi);
            B[i][j] = blendFactor * (UDS - CDS);
            B[i + 1][j] = -blendFactor * (UDS - CDS);
        }
    }

    //north
    for (int j = 1; j < Ny - 2; j++) {
        for (int i = 1; i < Nx - 1; i++) {

            AN[i][j] += std::min(F2[i][j], 0.0);
            AS[i][j + 1] += -std::max(F2[i][j], 0.0);
            double UDS = std::max(F2[i][j], 0.0) * U[i][j] + std::min(F2[i][j], 0.0) * U[i][j + 1];
            double CDS = F2[i][j] * (U[i][j + 1] * IFi + U[i][j] * IFi);
            B[i][j] += blendFactor * (UDS - CDS);
            B[i + 1][j] += -blendFactor * (UDS - CDS);
        }
    }

    //pressure grad source term

  // 1. ASSEMBLE AND SOLVE MOMENTUM EQUATION 



    extaraPolateZeroGrad(p, IFb, "north");
    extaraPolateZeroGrad(p, IFb, "south");

    extaraPolateZeroGrad(p, IFb, "east");
    extaraPolateZeroGrad(p, IFb, "west");

    printMat(p, "pressure");

    forAllInternal(p) {

        double pressureEastFace = (p[i + 1][j] * IFi) + (p[i][j] * (1.0 - IFi));
        double pressureWestFace = (p[i][j] * IFi) + (p[i - 1][j] * (1.0 - IFi));
        double pressureNorthFace = (p[i][j + 1] * IFi) + (p[i][j] * (1.0 - IFi));
        double pressureSouthFace = (p[i][j] * IFi) + (p[i][j - 1] * (1.0 - IFi));


        double pressureEastGrad = (pressureEastFace - pressureWestFace) / dx;
        double pressureNorthGrad = (pressureNorthFace - pressureSouthFace) / dy;

        // il gradiente di pressione lavora sulle source.

       // B[i][j] = - pressureEastGrad * dx * dy;



          //  B[i][j] = -pressureNorthGrad * dx*dy;


    }   // end forAllInternal







    //boundayy conditions 


    double U_north = 1;
    double U_wall = 0;
    double dx1 = dx;
    double dy1 = y[0][1] - y[0][0];
    double dx2 = x[1][0] - x[0][0];
    double dy2 = dy;
    // North boundary. no wall shear in X direction. Solo U
    for (int i = 0; i < Nx; i++) {
        APb[i][Ny - 1] += mu * dx1 / dy1;     // north boundary
        APb[i][0] += mu * dx1 / dy1;     // South boundary
        B[i][Ny - 1] += mu * dx1 / dy1 * U_north;
        B[i][0] += mu * dx1 / dy1 * U_wall;
    }
    // ap final
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            AP[i][j] += APb[i][j] - AE[i][j] - AW[i][j] - AS[i][j] - AN[i][j];
        }
    }

    //printMat(APb);

    /*int  n = 5;
    vector<vector<double>> a(n, vector<double>(n));
    vector<vector<double>> u(n, vector<double>(n));
    vector<vector<double>> l(n, vector<double>(n));

    a= {         {1, 1, 3, 0, 2}, {2, -3, 3, 0, -6}, {1, 7, 1, 9, 0}, {1, -2, -3, 2, 4}, {0, -1, 2, 1, 0}       };
   // printMat(a);
 */
 //luDecomposition(a, u, l, n);

    double URF = 0.8;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            AP[i][j] = (1.0 / URF) * AP[i][j];
            Brelaxed[i][j] = B[i][j] + (1.0 - URF) * (AP[i][j] * U[i][j]);
            rAP[i][j] = 1.0 / AP[i][j];
        }
    }

    /*printMat(AS, "AS");
    printMat(AN, "AN");
    printMat(AE, "AE");
    printMat(AW, "AW");
    printMat(AP, "AP");*/

    printMat(Brelaxed, "B");

    int iterations = 5;
    U = SIP(U, 0.8, AS, AN, AW, AE, AP, Brelaxed, iterations);
    // printMat(U, "U");
    return 0;
}

