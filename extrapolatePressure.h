#pragma once
#include<vector>
#include<string>
#include<iostream>
using std::vector;
using std::string;
//void extrapolateZeroGrad(vector<vector<double>>& p,  vector<double> IFb, string boundary)
/* {
    if (boundary == "north")
    {
        int j = p[1].size() - 1;
        for (int i = 1; i < p.size() - 1; i++) {
            p[i][j] = p[i][j - 1] + (1.0 - IFb[j - 1]) * (p[i][j - 1] - p[i][j - 2]);
            //p[i][j] = p[i][j - 1];
        }
    }
    if (boundary == "south")
    {
        int j = 0;
        for (int i = 1; i < p.size() - 1; i++) {
            p[i][j] = p[i][j + 1] + IFb[2] * (p[i][j + 1] - p[i][j + 2]);

        }
    }
    if (boundary == "east")
    {
        int i = p.size() - 1;
        for (int j = 1; j < p.size() - 1; j++) {
            p[i][j] = p[i - 1][j] + IFb[i - 1] * (p[i - 1][j] - p[i - 2][j]);
        }
    }
    if (boundary == "west")
    {
        int i = 0;
        for (int j = 1; j < p.size() - 1; j++) {
            p[i][j] = p[i + 1][j] + IFb[2] * (p[i + 1][j] - p[i + 2][j]);
        }
    }
}*/
void extrapolateZeroGrad(vector<vector<double>>& vec, vector<double> FXvec, string direction)

{

    // EAST
    if (direction == "east")
    {
        for (int j = 1; j < vec[0].size() - 1; j++)
        {
            int i = vec[0].size() - 1;
            vec[i][j] = vec[i - 1][j] + (vec[i - 1][j] - vec[i - 2][j]) * FXvec[i - 1];
        }
    }

    else if (direction == "west")
    {
        for (int j = 1; j < vec[0].size() - 1; j++)
        {
            int i = 0;
            vec[i][j] = vec[i + 1][j] + (vec[i + 1][j] - vec[i + 2][j]) * FXvec[2];
        }
    }

    else if (direction == "south")
    {
        for (int i = 1; i < vec.size() - 1; i++)
        {
            int j = 0;
            vec[i][j] = vec[i][j + 1] + (vec[i][j + 1] - vec[i][j + 2]) * FXvec[2];
        }
    }

    else  if (direction == "north")
    {
        for (int i = 1; i < vec.size() - 1; i++)
        {
            int j = vec.size() - 1;
            vec[i][j] = vec[i][j - 1] +(vec[i][j - 1] - vec[i][j - 2]) * (1.0 - FXvec[j - 1]);
         //   std::cout << "vec[i][j]=" << vec[i][j] << std::endl;
          //  std::cout << "vec[i][j-1]=" << vec[i][j] << std::endl;

        }
    }

    else
    {
        std::cout << " direction is wrong " << std::endl;
    }
}

