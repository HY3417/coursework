/**
* Import libraries
*/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

using namespace std;

#include <cblas.h>

int main(){

/**
* Defining Materials
*/
double kx = 250.0; // Thermal conductivity [W/mK]
double ky = 250.0; // Thermal conductivity [W/mK]
double kxy = 0.0;  // Thermal conductivity [W/mK]

/**
* Defining Section
*/
double th = 0.2;   // Thickness [m]

/**
* Integration scheme
*/
const int gaussorder = 2; // number of Gauss points

/**
* Defining Geometry
*
* h = a * x**2 + b * x + h1: Beam height as a function of x
*/
double a = 0;                   // constant in the polynomial describing the beam height
double h1 = 1.;                    //[m] Height at beam left edge
double h2 = h1 * 1.;              // [m] Height at beam right edge
double L = 2. * h1;                // [m] Beam length
double b = -a * L + (h2 - h1) / L; // constant in the polynomial describing the beam height

/**
* Meshing Geometry
*/
int nelem_y = 5;      // Number of element in y-direction
int nelem_x = 10;     // Number of element in x-direction
int nnode_elem = 4;   // Number of nodes in each element
int vNodeDof[4] = {1,1,1,1}; // Number of DoF per node (1 = Temperature only)
int neDof = 4;        // Total number of DoF per element

/**
* Calculations
*/
int nelem = nelem_x * nelem_y;             // Total number of elements
int nnode = (nelem_x + 1) * (nelem_y + 1); // Total number of nodes

/**
*
* Calculation of Nodal coordinate matrix
* h(x) = a * x**2 + b * x + h1 Beam height as a function of x
*/
double vX[nelem_x + 1];
vX[0] = 0.;
for (int i = 1; i < nelem_x + 1; i++){
    vX[i] = vX[i-1] + L / nelem_x;
}

double vH[nelem_x + 1];
for(int i = 0; i < nelem_x + 1; i++){
    vH[i] = a * vX[i] * vX[i] + b * vX[i] + h1;
}

double vY[nelem_y + 1][nelem_x + 1];
for(int i = 0; i < nelem_y + 1; i++){
    for(int j = 0; j < nelem_x + 1; j++){
    vY[i][j] = -0.5 * vH[j] + ((double) i) / ((double) nelem_y) * vH[j];
    }
}

double vCoord[nnode][2];  // Coordinate of the nodes [x, y]
for(int i = 0; i < nelem_y + 1; i++){
    for(int j = 0; j < nelem_x + 1; j++){
       vCoord[j*(nelem_y + 1) + i][0] = vX[j];
    }
}
for(int j = 0; j < nelem_x + 1; j++){
    for(int i = 0; i < nelem_y + 1; i++){
       vCoord[j*(nelem_y + 1) + i][1] = vY[i][j];
    }
}

/**
* Calculation of topology matrix NodeTopo
*/
int vNodeTopo[nelem_y + 1][nelem_x+1];
for(int i = 0; i < nelem_y + 1; i++){
    for(int j = 0; j < nelem_x + 1; j++){
       vNodeTopo[i][j] = j * (nelem_y + 1) + i;
    }
}

/**
* Calculation of topology matrix Elemnode
*/
int vElemNode[nelem][nnode_elem+1];
int elemnr = 0;
for(int colnr = 0; colnr < nelem_x; colnr++){
    for(int rownr = 0; rownr < nelem_y; rownr++){
    vElemNode[elemnr][0] = elemnr;
    vElemNode[elemnr][1] = vNodeTopo[rownr][colnr];
    vElemNode[elemnr][2] = vNodeTopo[rownr][colnr+1];
    vElemNode[elemnr][3] = vNodeTopo[rownr+1][colnr+1];
    vElemNode[elemnr][4] = vNodeTopo[rownr+1][colnr];
    elemnr += 1;
    }
}

double vElemX[nelem][nnode_elem];  // Element x nodal coordinates
for(int i = 0; i < nelem; i++){
    for(int j =1; j < nnode_elem + 1; j++){
    vElemX[i][j-1] = vCoord[vElemNode[i][j]][0];
    }
}

double vElemY[nelem][nnode_elem];  // Element y nodal coordinates
for(int i = 0; i < nelem; i++){
    for(int j =1; j < nnode_elem + 1; j++){
    vElemY[i][j-1] = vCoord[vElemNode[i][j]][1];
    }
}

int vGlobDof[nnode][2];  // nDof/node, Dof number
for(int i = 0; i < nnode; i++){
    for(int j = 0; j < nnode_elem; j++){
        vGlobDof[i][j] = 0;
    }
}

for(int i = 0; i < nelem; i++){
    for(int j = 0; j < nnode_elem; j++){  // loop over element nodes creats the first column
        int nNode = vElemNode[i][j+1];

        if(vGlobDof[nNode][0] < vNodeDof[j]){ // if the already existing ndof of the present node is less than
            vGlobDof[nNode][0] = vNodeDof[j]; // the present elements ndof then replace the ndof for that node
        }
    }
}

int nDof = 0; // counting the global dofs and inserting in vGlobDof
for(int i = 0; i < nnode; i++){
    int eDof = vGlobDof[i][0];
    for(int j = 0; j < eDof; j++){
        vGlobDof[i][j+1] = nDof;
        nDof += 1;
    }
}

/**
*Assumbly of global stiffness matrix K
*/
int gauss = gaussorder;
double vGP[2] = {-1 / sqrt(3), 1 / sqrt(3)}; // Points
int vW[2] = {1, 1}; // Weights
double vD[gauss*gauss] = {kx, kxy, kxy, ky}; // Conductivity matrix D
double vK[nDof*nDof] = {}; // Initiation of global stiffness matrix K

for(int i = 0; i < nelem; i++){
    int vNodes[4];  // Element nodes
    double eCoord[nnode_elem*gauss] = {}; // node coordinates

    int iCoord = 0;
    for(int index = 1; index < nnode_elem + 1; index++){

        vNodes[index-1] = vElemNode[i][index];

        eCoord[iCoord] = vCoord[vNodes[index-1]][0];
        eCoord[iCoord+1] = vCoord[vNodes[index-1]][1];

        iCoord = iCoord + 2;
    }

    int gDof[nnode_elem];
    for(int j = 0; j < nnode_elem; j++){
        gDof[j] = vNodes[j]; // global dof for node j

    }

    double vKe[nnode_elem*nnode_elem] = {}; // Local stiffness matrix
    double DetJ; // For storing the determinants of J
    double vInvJ[gauss*gauss] = {}; // Inverse Jacobian matrix
    double vGN[gauss*nnode_elem] = {}; // Derivative (gradient) of the shape functions
    double vB[gauss*nnode_elem] = {}; // Gradient of the shape function


    for(int k = 0; k < gauss; k++){
        for(int j = 0; j < gauss; j++){
            double eta = vGP[k];
            double xi = vGP[j];


            double vN[4] = {0.25 * (1 - xi) * (1 - eta),
                            0.25 * (1 + xi) * (1 - eta),
                            0.25 * (1 + xi) * (1 + eta),
                            0.25 * (1 - xi) * (1 + eta)};
            double vGN[gauss*nnode_elem] = {0.25 * -(1 - eta), 0.25 * (1 - eta),
                                            0.25 * (1 + eta), 0.25 * -(1 + eta),
                                            0.25 * -(1 - xi), 0.25 * -(1 + xi),
                                            0.25 * (1 + xi), 0.25 * (1 - xi)};

            double vJ[gauss*gauss] = {}; //Jacobian matrix
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, gauss, gauss, nnode_elem, 1.0, vGN, nnode_elem, eCoord, gauss, 0.0, vJ, gauss);


            double DetJ = vJ[0] * vJ[3] - vJ[1] * vJ[2]; // Jacobian determinant

            vInvJ[0] = vJ[3] / (vJ[0] * vJ[3] - vJ[1] * vJ[2]);
            vInvJ[1] = -vJ[1] / (vJ[0] * vJ[3] - vJ[1] * vJ[2]);
            vInvJ[2] = -vJ[2] / (vJ[0] * vJ[3] - vJ[1] * vJ[2]);
            vInvJ[3] = vJ[0] / (vJ[0] * vJ[3] - vJ[1] * vJ[2]);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, gauss, nnode_elem, gauss, 1.0, vInvJ, gauss, vGN, nnode_elem, 0.0, vB, nnode_elem);

            double vBT[nnode_elem*gauss] = {}; // the transposition matrix of vB
            vBT[0] = vB[0];
            vBT[1] = vB[4];
            vBT[2] = vB[1];
            vBT[3] = vB[5];
            vBT[4] = vB[2];
            vBT[5] = vB[6];
            vBT[6] = vB[3];
            vBT[7] = vB[7];

            double vBD[nnode_elem*gauss] = {}; // the matrix generated by timing transposition of matrix vB with matrix vD
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nnode_elem, gauss, gauss, 1.0, vBT, gauss, vD, gauss, 0.0, vBD, gauss);

            double vBDB[nnode_elem*nnode_elem] = {}; // the matrix generated by timing matrix vBD and MATRIX vB
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nnode_elem, nnode_elem, gauss, 1.0, vBD, gauss, vB, nnode_elem, 0.0, vBDB, nnode_elem);

            for(int n = 0; n < nnode_elem*nnode_elem; n++){
                vKe[n] = vKe[n] + vBDB[n] * th * DetJ * vW[k] * vW[j];
            }
        }
    }

    for(int n = 0; n < nnode_elem; n++){
        for(int k = 0; k < nnode_elem; k++){
            vK[gDof[n]*nDof+gDof[k]] = vK[gDof[n]*nDof+gDof[k]] + vKe[n*nnode_elem+k];
            }
       }
}

return 0;
}
