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
    int Nx = 5, Ny = 5; // numero di nodi interni su x e su y 
    int nb = 2; //numero di boundaries
    int Nxb = Nx + 2;
    int Nyb = Ny + 2;
    double dx = L / Nx; //dimensione-x del volume di controllo
    double dy = L / Ny;
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
    for (int outerIt = 0; outerIt < 1; outerIt++) {
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

 

        extrapolateZeroGrad(p, IFb, "north");
        extrapolateZeroGrad(p, IFb, "south");
        extrapolateZeroGrad(p, IFb, "east");
        extrapolateZeroGrad(p, IFb, "west");
        printMat(p, "p");
        pressureGradSource(IFi, p, Bu, Bv, dx, dy);
        diffusion(mu, dy, dx, AEu, AWu, ANu, ASu, x, y);    //diffusion U
        convection(IFi, AEu, AWu, ANu, ASu, Fe, Fn, U, Bu);

        double U_north = 1;
        double U_wall = 0.0;
        double V_wall = 0.0;

        noWallShearX(Nxb, Nyb, dx, dy, mu, U_north, U_wall, x, y, APbu, Bu);
        assembleMatrixUV(APu, APbu, AEu, AWu, ASu, ANu);
        printMat(APu, "APu");
        printMat(AEu, "AEu");
        printMat(ANu, "ANu");
        printMat(AWu, "AWu");
        printMat(ASu, "ASu");

        double URF = 0.8;

        relaxEq(URF, APu, BuRelaxed, Bu, rAPu, U);



        //U=SIP(U, 0.92, AS, AN, AW, AE, AP, Brelaxed, iterations);
        int iteration = 1;
        UU = SIP(U, 0.92, ASu, ANu, AWu, AEu, APu, BuRelaxed, iteration, "U-Momentum");
        printMat(UU, "UU");

        diffusion(mu, dy, dx, AEv, AWv, ANv, ASv, x, y);    //diffusion U
        convection(IFi, AEv, AWv, ANv, ASv, Fe, Fn, V, Bv);
        noWallShearY(Nxb, Nyb, dx, dy, mu, V_wall, x, y, APbv, Bv);

        //boundayy conditions 
        assembleMatrixUV(APv, APbv, AEv, AWv, ASv, ANv);
        relaxEq(URF, APv, BvRelaxed, Bv, rAPv, V);

        //  printMat(VV, "VV");
        VV = SIP(V, 0.92, ASv, ANv, AWv, AEv, APv, BvRelaxed, iteration, "V-Momentum");
        printMat(VV, "VV");

        double linearTerm;
        double c1;
        double c2;

        forAllInternalUCVs(Ue) {
            linearTerm = IFi * UU[i + 1][j] + (1 - IFi) * UU[i][j];
            c1         = - dy * (IFi / APu[i + 1][j] + (1 - IFi) / APu[i][j]);
            c2         = ( p[i+1][j] - p[i][j] ) - 0.5 * ( (IFi * (p[i + 2][j] - p[i][j]) + (1 - IFi) * ( p[i+1][j]-p[i][j] )));
            Ue[i][j]   = linearTerm + c1 * c2;

            Fe[i][j] = density * dy * Ue[i][j];
            Fw[i][j] = Fe[i - 1][j];
        }

        forAllInternalVCVs(Vn) {
            linearTerm = IFi * VV[i][j+1] + (1 - IFi) * VV[i][j];
            c1         = - dx * ( IFi / APv[i][j+1] + (1 - IFi) / APv[i][j]);
            c2         = ( p[i][j+1] - p[i][j] ) - 0.5 * ( (IFi * (p[i][j+2] - p[i][j]) + (1 - IFi) * ( p[i][j+1] - p[i][j] )));

            Vn[i][j] = linearTerm + c1 * c2;

            Fn[i][j] = density * dy * Vn[i][j];
            Fs[i][j] = Fn[i][j-1];
        }

        double sum = 0;
        // pressure corretion equation 
        forAllInternalUCVs(p) {
            AEp[i][j]     += - density * dy * dy * ( IFi/ APu[i + 1][j]  + (1 - IFi) /APu[i][j]);
            AWp[i + 1][j] +=   AEp[i][j];
        }
        // pressure corretion equation 
        forAllInternalVCVs(p) {
            ANp[i][j]     += -density * dx * dx * (IFi / APv[i][j + 1]  + (1 - IFi) / APv[i][j]);
            ASp[i][j + 1] +=  ANp[i][j];
        }
        forAllInternal(p) {
            APp[i][j] += - (ANp[i][j] + ASp[i][j] + AEp[i][j] + AWp[i][j]);
        }
        forAllInternal(Bp) {
            Bp[i][j] += Fe[i-1][j] - Fe[i][j] + Fn[i][j-1] - Fn[i][j];
            sum += Bp[i][j];
        }
        cout << "SUM" << sum << endl;
        printMat(APp);
        iteration = 6;
        pp = SIP(p, 0.92, ASp, ANp, AWp, AEp, APp, Bp, iteration, "Pressure-correction");
        printMat(pp, "pp");
        extrapolateZeroGrad(pp, IFb, "north");
        extrapolateZeroGrad(pp, IFb, "south");
        extrapolateZeroGrad(pp, IFb, "east");
        extrapolateZeroGrad(pp, IFb, "west");
        printMat(pp, "pp-extrapolate");

        double URFp = 0.2;
        //correct p

        forAllInternal(p) {
            p[i][j] = p[i][j] + URFp * pp[i][j];
        }


        //correct u e v

        forAllInternal(U) {
            UU[i][j] = UU[i][j] + ( dy / APu[i][j]) * ( ( pp[i - 1][j] * IFi + (1 - IFi) * pp[i][j] ) - ( pp[i + 1][j] * IFi + (1 - IFi) * pp[i][j]));
            VV[i][j] = VV[i][j] + ( dx / APv[i][j]) * ( ( pp[i][j - 1] * IFi + (1 - IFi) * pp[i][j] ) - ( pp[i][j + 1] * IFi + (1 - IFi) * pp[i][j]));
        }

        //correct face velocities
        forAllInternalUCVs(Ue) {
            Ue[i][j] = Ue[i][j] + dy * (IFi / APu[i + 1][j] + (1 - IFi) / APu[i][j]) * (pp[i][j] - pp[i + 1][j]);
        } 
        forAllInternalVCVs(Vn) {
            Vn[i][j] = Vn[i][j] + dx * (IFi / APv[i][j + 1] + (1 - IFi)/ APv[i][j]) * (pp[i][j] - pp[i][j + 1]);
        } 

        forAllInternal(U) {
            U[i][j] = UU[i][j];
            V[i][j] = VV[i][j];
        }

        forAllInternalUCVs(Ue) {
            Fe[i][j] = density * dy * Ue[i][j];
            //Fw[i][j] = Fe[i - 1][j];
        }

        forAllInternalVCVs(Vn) {
            Fn[i][j] = density * dy * Vn[i][j];
            //Fs[i][j] = Fn[i][j - 1];
        }

        printMat(U, "Unew");
        printMat(V, "Vnew");
        

    }
    return 0;
}

