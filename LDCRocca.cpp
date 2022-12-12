// 15/11/2022 inizio di un nuovo progetto lid driven cavity. Obiettivo arrivare a lid driven cavity cuda
// finite volume method colocated



#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include "forLoop.h"
#include "printTools.h"
#include "extrapolatePressure.h"
#include "SIP.h"
#include "diffusion.h"
#include "convection.h"
#include "pressureGradSource.h"
#include "noWallShearX.h"
#include "noWallShearY.h"
#include "assembleMatrixUV.h"
#include "relax.h"
using std::vector;
using std::cout;
using std::endl;
using std::string;

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

    vector<vector<double>> x(Nxb, vector<double>(Nyb));
    vector<vector<double>> y(Nxb, vector<double>(Nyb));
    vector<vector<double>> U(Nxb, vector<double>(Nyb));
    vector<vector<double>> UU(Nxb, vector<double>(Nyb));
    vector<vector<double>> V(Nxb, vector<double>(Nyb));
    vector<vector<double>> VV(Nxb, vector<double>(Nyb));
    vector<vector<double>> p(Nxb, vector<double>(Nxb));
    vector<vector<double>> pp(Nxb, vector<double>(Nxb));
    vector<vector<double>> Vn(Nxb, vector<double>(Nyb));
    vector<vector<double>> Ue(Nxb, vector<double>(Nyb));
    vector<vector<double>> Fe(Nxb, vector<double>(Nyb)); //flussi x
    vector<vector<double>> Fn(Nxb, vector<double>(Nyb)); //flussi y
    vector<vector<double>> Fw(Nxb, vector<double>(Nyb)); //flussi x
    vector<vector<double>> Fs(Nxb, vector<double>(Nyb)); //flussi y
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
            if (j == 0) { y[i][j] = 0; }
            else if (j == Nyb - 1) { y[i][j] = 1; }
            else y[i][j] = dy / 2 + (j - 1) * dy;
            //  y[n] = dy / 2 + j * dy;
        }
    }
    for (int outerIt = 1; outerIt < 15; outerIt++) {
        cout << "ITERATION NUMBER " << outerIt << endl;
        cout << endl;

        vector<vector<double>> AEu(Nxb, vector<double>(Nyb));
        vector<vector<double>> AWu(Nxb, vector<double>(Nyb));
        vector<vector<double>> ASu(Nxb, vector<double>(Nyb));
        vector<vector<double>> ANu(Nxb, vector<double>(Nyb));
        vector<vector<double>> APbu(Nxb, vector<double>(Nyb));
        vector<vector<double>> APu(Nxb, vector<double>(Nyb));
        vector<vector<double>> Bu(Nxb, vector<double>(Nyb));
        vector<vector<double>> BuRelaxed(Nxb, vector<double>(Nyb));
        vector<vector<double>> rAPu(Nxb, vector<double>(Nyb));
        vector<vector<double>> AEv(Nxb, vector<double>(Nyb));
        vector<vector<double>> AWv(Nxb, vector<double>(Nyb));
        vector<vector<double>> ASv(Nxb, vector<double>(Nyb));
        vector<vector<double>> ANv(Nxb, vector<double>(Nyb));
        vector<vector<double>> APbv(Nxb, vector<double>(Nyb));
        vector<vector<double>> APv(Nxb, vector<double>(Nyb));
        vector<vector<double>> Bv(Nxb, vector<double>(Nyb));
        vector<vector<double>> BvRelaxed(Nxb, vector<double>(Nyb));
        vector<vector<double>> rAPv(Nxb, vector<double>(Nyb));
        vector<vector<double>> AEp(Nxb, vector<double>(Nyb));
        vector<vector<double>> AWp(Nxb, vector<double>(Nyb));
        vector<vector<double>> ASp(Nxb, vector<double>(Nyb));
        vector<vector<double>> ANp(Nxb, vector<double>(Nyb));
        vector<vector<double>> APbp(Nxb, vector<double>(Nyb));
        vector<vector<double>> APp(Nxb, vector<double>(Nyb));
        vector<vector<double>> Bp(Nxb, vector<double>(Nyb));

 

        extaraPolateZeroGrad(p, IFb, "north");
        extaraPolateZeroGrad(p, IFb, "south");
        extaraPolateZeroGrad(p, IFb, "east");
        extaraPolateZeroGrad(p, IFb, "west");
        pressureGradSource(IFi, p, Bu, Bv, dx, dy);
        diffusion(mu, dy, dx, AEu, AWu, ANu, ASu, x, y);    //diffusion U
        convection(IFi, AEu, AWu, ANu, ASu, Fe, Fn, U, Bu);

        double U_north = 1;
        double U_wall = 0.0;
        double V_wall = 0.0;

        noWallShearX(Nx, Ny, dx, dy, mu, U_north, U_wall, x, y, APbu, Bu);


        //boundayy conditions 


        assembleMatrixUV(APu, APbu, AEu, AWu, ASu, ANu);
        double URF = 0.8;

        relaxEq(URF, APu, BuRelaxed, Bu, rAPu, U);



        //U=SIP(U, 0.92, AS, AN, AW, AE, AP, Brelaxed, iterations);
        int iteration = 7;
        UU = SIP(U, 0.7, ASu, ANu, AWu, AEu, APu, BuRelaxed, iteration, "U-Momentum");


        diffusion(mu, dy, dx, AEv, AWv, ANv, ASv, x, y);    //diffusion U
        convection(IFi, AEv, AWv, ANv, ASv, Fe, Fn, V, Bv);
        noWallShearY(Nx, Ny, dx, dy, mu, V_wall, x, y, APbv, Bv);

        //boundayy conditions 
        assembleMatrixUV(APv, APbv, AEv, AWv, ASv, ANv);
        relaxEq(URF, APv, BvRelaxed, Bv, rAPv, V);

        //  printMat(VV, "VV");
        VV = SIP(V, 0.7, ASv, ANv, AWv, AEv, APv, BvRelaxed, iteration, "V-Momentum");

        double linearTerm;
        double c1;
        double c2;
        double c3;

        forAllInternalUCVs(Ue) {
            linearTerm = IFi * UU[i + 1][j] + (1 - IFi) * UU[i][j];
            c1 = -(URF * dy * (p[i + 1][j] - p[i][j]))*(IFi / APu[i + 1][j] + (1 - IFi) / APu[i][j]);
            c2 = (IFi * URF * dy / APu[i + 1][j]) * ((IFi * p[i + 2][j] + (1 - IFi) * p[i+1][j]) - (IFi * p[i + 1][j] + (1 - IFi) * p[i][j]));
            c3 = ((1 - IFi) * URF * dy / APu[i][j] )* ((IFi * p[i + 1][j] + (1 - IFi) * p[i][j]) - (IFi * p[i-1][j] + (1 - IFi) * p[i][j])) ;
            Ue[i][j] = linearTerm + c1 + c2 + c3;
                
            Fe[i][j] = density * dy * Ue[i][j];
            Fw[i][j] = Fe[i - 1][j];
        }

        
        forAllInternalVCVs(Vn) {
            linearTerm = IFi * VV[i][j+1] + (1 - IFi) * VV[i][j];
            c1 = -(URF * dx * (p[i][j+1] - p[i][j])) * (IFi / APv[i][j+1] + (1 - IFi) / APv[i][j]);
            c2 = (IFi * URF * dx / APv[i][j+1]) * ((IFi * p[i][j+2] + (1 - IFi) * p[i][j+1]) - (IFi * p[i ][j+1] + (1 - IFi) * p[i][j]));
            c3 = ((1 - IFi) * URF * dy / APv[i][j]) * ((IFi * p[i][j+1] + (1 - IFi) * p[i][j]) - (IFi * p[i][j-1] + (1 - IFi) * p[i][j]));
            Vn[i][j] = linearTerm + c1 + c2 + c3;

            Fn[i][j] = density * dy * Vn[i][j];
            Fs[i][j] = Fn[i][j-1];
        }

     //   printMat(Vn, "Vn");

        // pressure corretion equation 
        forAllInternalUCVs(p) {
            AEp[i][j]     += -density * dy * dy / (APu[i + 1][j] * IFi + APu[i][j] * (1 - IFi));
            AWp[i + 1][j] += -density * dy * dy / (APu[i + 1][j] * IFi + APu[i][j] * (1 - IFi));
        }
        // pressure corretion equation 
        forAllInternalVCVs(p) {
            ANp[i][j]     += -density * dx * dx / (APv[i][j + 1] * IFi + APv[i][j] * (1 - IFi));
            ASp[i][j + 1] += -density * dx * dx / (APv[i][j + 1] * IFi + APv[i][j] * (1 - IFi));
        }
        forAllInternal(p) {
            APp[i][j] += - (ANp[i][j] + ASp[i][j] + AEp[i][j] + AWp[i][j]);
        }
        forAllInternal(Bp) {
            Bp[i][j] += Fw[i][j] - Fe[i][j] + Fs[i][j] - Fn[i][j];
        }

        iteration = 7;
        pp = SIP(p, 0.9, ASp, ANp, AWp, AEp, APp, Bp, iteration, "Pressure-correction");
        //   printMat(pp, "pp");
        extaraPolateZeroGrad(pp, IFb, "north");
        extaraPolateZeroGrad(pp, IFb, "south");
        extaraPolateZeroGrad(pp, IFb, "east");
        extaraPolateZeroGrad(pp, IFb, "west");
        //   printMat(pp, "pp");

        double URFp = 0.3;
        //correct p

        forAll(p) {
            p[i][j] = p[i][j] + URFp * pp[i][j];
        }


        //correct u e v

        forAllInternal(U) {
            UU[i][j] = U[i][j] + (URF * dy / APu[i][j]) * ( (pp[i - 1][j] * IFi + (1 - IFi) * pp[i][j]) - (pp[i + 1][j] * IFi + (1 - IFi) * pp[i][j]));
            VV[i][j] = V[i][j] + (URF * dx / APv[i][j]) * ( (pp[i][j - 1] * IFi + (1 - IFi) * pp[i][j]) - (pp[i][j + 1] * IFi + (1 - IFi) * pp[i][j]));
        }

        //correct face velocities
        forAllInternalUCVs(Ue) {
            Ue[i][j] = Ue[i][j] + (URF * dy) * (IFi/ APu[i + 1][j] + (1 - IFi) / APu[i][j]) * (pp[i][j] - pp[i + 1][j]);
        } 
        forAllInternalVCVs(Vn) {
            Vn[i][j] = Vn[i][j] + (URF * dx) * (IFi / APv[i][j + 1] + (1 - IFi)/ APv[i][j]) * (pp[i][j] - pp[i][j + 1]);
        } 

        forAllInternal(U) {
            U[i][j] = UU[i][j];
            V[i][j] = VV[i][j];
        }

        forAllInternalUCVs(Ue) {
            Fe[i][j] = density * dy * Ue[i][j];
            Fw[i][j] = Fe[i - 1][j];
        }

        forAllInternalVCVs(Vn) {
            Fn[i][j] = density * dy * Vn[i][j];
            Fs[i][j] = Fn[i][j - 1];
        }


    }
    return 0;
}

