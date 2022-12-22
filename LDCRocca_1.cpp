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
    vector<double> IFx(Nxb); // interpolation factor x
    vector<double> IFy(Nyb); // interpolation factor y
    vector<double> X(Nxb); // interpolation factor y
    vector<double> XC(Nxb); // interpolation factor y
    vector<double> Y(Nyb); // interpolation factor y
    vector<double> YC(Nyb); // interpolation factor y



    for (int i = 0; i < Nxb; i++) {
        for (int j = 0; j < Nyb; j++) {
            if (i == 0) { x[i][j] = 0; }
            else if (i == Nxb - 1) { x[i][j] = 1; }
            else x[i][j] = dx / 2 + (i - 1) * dx;
            //  y[n] = dy / 2 + j * dy;
        }
    }


    for (int i = 1; i < X.size(); i++)
    {
        X[i] = X[i - 1] + dx;
    }
    X[X.size() - 1] = X[X.size() - 2];


    for (int i = 1; i < XC.size(); i++)
    {
        XC[i] = 0.5 * (X[i] + X[i - 1]);
    }
    XC[0] = X[0];
    XC[X.size() - 1] = X[X.size() - 2];



    for (int i = 1; i < Y.size(); i++)
    {
        Y[i] = Y[i - 1] + dy;
    }
    Y[Y.size() - 1] = Y[Y.size() - 2];


    for (int i = 1; i < Y.size(); i++)
    {
        YC[i] = 0.5 * (Y[i] + Y[i - 1]);
    }
    YC[0] = Y[0];
    YC[YC.size() - 1] = Y[Y.size() - 2];

    for (int i = 0; i < X.size() - 1; i++)
    {
        IFx[i] = (X[i] - XC[i]) / (XC[i + 1] - XC[i]);
    }

    for (int i = 0; i < Y.size() - 1; i++)
    {
        IFy[i] = (Y[i] - YC[i]) / (YC[i + 1] - YC[i]);
    }


    printVec(X, "X");
    printVec(XC, "XC");
    printVec(IFx, "IFx");



    for (int i = 0; i < Nxb; i++) {
        for (int j = 0; j < Nyb; j++) {
            if (j == 0) { y[i][j] = 0; }
            else if (j == Nyb - 1) { y[i][j] = 1; }
            else y[i][j] = dy / 2 + (j - 1) * dy;
            //  y[n] = dy / 2 + j * dy;
        }
    }
    for (int outerIt = 1; outerIt < 30; outerIt++) {
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



        extrapolateZeroGrad(p, IFy, "north");
        extrapolateZeroGrad(p, IFy, "south");
        extrapolateZeroGrad(p, IFx, "east");
        extrapolateZeroGrad(p, IFx, "west");
        printMat(p, "p");
        diffusion(mu, dy, dx, AEu, AWu, ANu, ASu, XC, YC);    //diffusion U

        double U_north = 1;
        double U_wall = 0.0;
        double V_wall = 0.0;
        noWallShearX(Nxb, Nyb, dx, dy, mu, U_north, U_wall, XC, YC, APbu, Bu);
        printMat(APbu, "APBu");

        convection(IFx, IFy, AEu, AWu, ANu, ASu, Fe, Fn, U, Bu);
        pressureGradSource(IFx, IFy, p, Bu, Bv, dx, dy);
       // printMat(Bu, "PRESSURE GRAD SOURCE");
        printMat(Bu, "TOTAL SOURCE");

        assembleMatrixUV(APu, APbu, AEu, AWu, ASu, ANu);
        printMat(APu, "APu");
        printMat(AEu, "AEu");
        printMat(ANu, "ANu");
        printMat(AWu, "AWu");
        printMat(ASu, "ASu");
        printMat(Bu, "Bu");

        double URF = 0.8;

        relaxEq(URF, APu, BuRelaxed, Bu, U);
        int iteration = 1;
        UU = SIP(U, 0.92, ASu, ANu, AWu, AEu, APu, BuRelaxed, iteration, "U-Momentum");
        printMat(UU, "UU");

        diffusion(mu, dy, dx, AEv, AWv, ANv, ASv, XC, YC);    //diffusion U
        convection(IFx, IFy, AEv, AWv, ANv, ASv, Fe, Fn, V, Bv);
        noWallShearY(Nxb, Nyb, dx, dy, mu, V_wall, XC, YC, APbv, Bv);

        //boundayy conditions 
        assembleMatrixUV(APv, APbv, AEv, AWv, ASv, ANv);

        printMat(APv, "APv");
        printMat(AEv, "AEv");
        printMat(ANv, "ANv");
        printMat(AWv, "AWv");
        printMat(ASv, "ASv");
        printMat(Bv, "Bv");

        relaxEq(URF, APv, BvRelaxed, Bv,  V);

        //  printMat(VV, "VV");
        VV = SIP(V, 0.92, ASv, ANv, AWv, AEv, APv, BvRelaxed, iteration, "V-Momentum");
        printMat(VV, "VV");

        double linearTerm;
        double c1;
        double c2;

        forAllInternalUCVs(Ue) {
            linearTerm = IFx[i] * UU[i + 1][j] + (1 - IFx[i]) * UU[i][j];
            c1 = -dy * (IFx[i] / APu[i + 1][j] + (1 - IFx[i]) / APu[i][j]);
            c2 = (p[i + 1][j] - p[i][j]) - 0.5 * ((IFx[i] * (p[i + 2][j] - p[i][j]) + (1 - IFx[i]) * (p[i + 1][j] - p[i][j])));
            Ue[i][j] = linearTerm + c1 * c2;

            Fe[i][j] = density * dy * Ue[i][j];
            Fw[i][j] = Fe[i - 1][j];
        }

        forAllInternalVCVs(Vn) {
            linearTerm = IFy[j] * VV[i][j + 1] + (1 - IFy[j]) * VV[i][j];
            c1 = -dx * (IFy[j] / APv[i][j + 1] + (1 - IFy[j]) / APv[i][j]);
            c2 = (p[i][j + 1] - p[i][j]) - 0.5 * ((IFy[j] * (p[i][j + 2] - p[i][j]) + (1 - IFy[j]) * (p[i][j + 1] - p[i][j])));

            Vn[i][j] = linearTerm + c1 * c2;

            Fn[i][j] = density * dy * Vn[i][j];
            Fs[i][j] = Fn[i][j - 1];
        }

        double sum = 0;
        // pressure corretion equation 
        forAllInternalUCVs(p) {
            AEp[i][j] += -density * dy * dy * (IFx[i] / APu[i + 1][j] + (1 - IFx[i]) / APu[i][j]);
            AWp[i + 1][j] += AEp[i][j];
        }
        // pressure corretion equation 
        forAllInternalVCVs(p) {
            ANp[i][j] += -density * dx * dx * (IFy[j] / APv[i][j + 1] + (1 - IFy[j]) / APv[i][j]);
            ASp[i][j + 1] += ANp[i][j];
        }
        forAllInternal(p) {
            APp[i][j] += -(ANp[i][j] + ASp[i][j] + AEp[i][j] + AWp[i][j]);
        }
        forAllInternal(Bp) {
            Bp[i][j] += Fe[i - 1][j] - Fe[i][j] + Fn[i][j - 1] - Fn[i][j];
            sum += Bp[i][j];
        }
        cout << "SUM" << sum << endl;
        printMat(APp);
        iteration = 6;
        pp = SIP(pp, 0.92, ASp, ANp, AWp, AEp, APp, Bp, iteration, "Pressure-correction");
        printMat(pp, "pp");
        extrapolateZeroGrad(pp, IFy, "north");
        extrapolateZeroGrad(pp, IFy, "south");
        extrapolateZeroGrad(pp, IFx, "east");
        extrapolateZeroGrad(pp, IFx, "west");
        printMat(pp, "pp-extrapolate");

        double URFp = 0.2;
        //correct p

     
        printMat(p, "p before correction");

        //correct u e v

        forAllInternal(U) {

            dx = X[i] - X[i - 1];
            dy = Y[j] - Y[j - 1];

            double PPE = (pp[i + 1][j] * IFx[i]) + (pp[i][j] * (1.0 - IFx[i]));
            double PPW = (pp[i][j] * IFx[i - 1]) + (pp[i - 1][j] * (1.0 - IFx[i-1]));
            double PPN = (pp[i][j + 1] * IFy[j]) + (pp[i][j] * (1.0 - IFy[j]));
            double PPS = (pp[i][j] * IFy[j-1])   + (pp[i][j - 1] * (1.0 - IFy[j-1]));

            U[i][j] = UU[i][j] - (PPE - PPW) * dy * (1 / APu[i][j]);
            V[i][j] = VV[i][j] - (PPN - PPS) * dx * (1 / APv[i][j]);
            p[i][j] = p[i][j] + URFp * (pp[i][j]);// -pp[4][4]);

            // UU[i][j] = UU[i][j] + ( dy / APu[i][j]) * ( ( pp[i - 1][j] * IFi + (1 - IFi) * pp[i][j] ) - ( pp[i + 1][j] * IFi + (1 - IFi) * pp[i][j]));
            // VV[i][j] = VV[i][j] + ( dx / APv[i][j]) * ( ( pp[i][j - 1] * IFi + (1 - IFi) * pp[i][j] ) - ( pp[i][j + 1] * IFi + (1 - IFi) * pp[i][j]));

        }

        //correct face velocities
        forAllInternalUCVs(Ue) {
      //      Ue[i][j] = Ue[i][j] + dy * (IFx[i] / APu[i + 1][j] + (1 - IFx[i]) / APu[i][j]) * (pp[i][j] - pp[i + 1][j]);
            //10.1 East
            Fe[i][j] = Fe[i][j] + (AEp[i][j] * (pp[i + 1][j] - pp[i][j]));
   
//            Fe[i][j] = density * dy * Ue[i][j];

        }
        forAllInternalVCVs(Vn) {
          //  Vn[i][j] = Vn[i][j] + dx * (IFy[j] / APv[i][j + 1] + (1 - IFy[j]) / APv[i][j]) * (pp[i][j] - pp[i][j + 1]);
           // Fn[i][j] = density * dy * Vn[i][j];
            Fn[i][j] = Fn[i][j] + (ANp[i][j] * (pp[i][j + 1] - pp[i][j]));
        }

        printMat(U, "Unew");
        printMat(V, "Vnew");
        printMat(p, "Pnew");

        printMat(Fe, "Fe NEW");
        printMat(Fn, "Fn NEW");

    }
    return 0;
}

