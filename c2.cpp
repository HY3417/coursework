/**
*
* Import libraries
*
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <fstream>
using namespace std;

#include <cblas.h>
#define F77NAME(x) x##_

/**
*
*LAPACK routine for solving systems of linear equations
*
*/
extern "C"{
    void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                        const int& lda, int * ipiv, double * B,
                        const int& ldb, int& info);
}

int main(int argc, char** argv){
double a = atof(argv[1]);     // constant in the polynomial describing the beam height
double h1 = atof(argv[2]);    // [m] Height at beam left edge
double h2 = atof(argv[3]);    // [m] Height at beam right edge
double L = atof(argv[4]);     // [m] Beam length
double tp = atof(argv[5]);    // Thickness [m]
double kx = atof(argv[6]);    // Thermal conductivity [W/mK]
double ky = atof(argv[7]);    // Thermal conductivity [W/mK]
double kxy = atof(argv[8]);   // Thermal conductivity [W/mK]
int nelem_x = atoi(argv[9]);  // Number of element in x-direction
int nelem_y = atoi(argv[10]); // Number of element in y-direction

/**
* Integration scheme
*/
const int gaussorder = 2; // number of Gauss points

/**
* Defining Geometry
*
* h = a * x**2 + b * x + h1: Beam height as a function of x
*/
double b = -a * L + (h2 - h1) / L; // constant in the polynomial describing the beam height

/**
* Meshing Geometry
*/
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
                vKe[n] = vKe[n] + vBDB[n] * tp * DetJ * vW[k] * vW[j];
            }
        }
    }

    for(int n = 0; n < nnode_elem; n++){
        for(int k = 0; k < nnode_elem; k++){
            vK[gDof[n]*nDof+gDof[k]] = vK[gDof[n]*nDof+gDof[k]] + vKe[n*nnode_elem+k];
            }
       }
}

/**
* Compute nodal boundary flux vector
* Natural boundary condition
*/

/**
*Defined on edges
*/
int vFluxNodes[nelem_x + 1] = {};
int nFluxNodes = 0;
for(int i = 0; i < nelem_x + 1; i++){
    vFluxNodes[i] = vNodeTopo[nelem_y][i];  // Nodes at the top edge of the beam
    nFluxNodes += 1;  // Number of nodes on the top edge of the beam
}

/**
*Defining load
*/
int q = 2500; // Constant flux at top edge of the beam
double vN_bc[4][nFluxNodes-1];
for(int i = 0; i < nFluxNodes - 1; i++){
    vN_bc[0][i] = vFluxNodes[i];   // Node 1
    vN_bc[1][i] = vFluxNodes[i+1]; // Node 2
    vN_bc[2][i] = q; // Flux value at node 1
    vN_bc[3][i] = q; // FLux value at node 2
}

double vF[nDof] = {}; // Initialize nodal flux vector
for(int i = 0; i < nFluxNodes - 1; i++){
    double vFq[2] = {}; // Initialize the nodal source vector

    int node1; // First node
    int node2; // Second node
    node1 = vN_bc[0][i];
    node2 = vN_bc[1][i];

    double vN_bce[2] = {};
    for(int m = 0; m < 2; m++){
        vN_bce[m] = vN_bc[m+2][i]; // FLux value at an edge
    }

    double x1; // x coord of the first node
    double x2; // x coord of the second node
    double y1; // y coord of the first node
    double y2; // y coord of the second node
    x1 = vCoord[node1][0];
    y1 = vCoord[node1][1];
    x2 = vCoord[node2][0];
    y2 = vCoord[node2][1];

    double leng; // Edge length
    double detJ; // 1D Jacobian
    leng = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    detJ = leng / 2;

    for(int j = 0; j < gauss; j++){ // Integrate in x direction (1D integration)
        double x;
        double vNN[2] = {}; // 1D shape functions in parent domain
        double flux;
        x = vGP[j];
        vNN[0] = 0.5 * (1 - x);
        vNN[1] = 0.5 * (1 + x);
        flux = cblas_ddot(2, vNN, 1, vN_bce, 1);
        vFq[0] = vFq[0] + vW[j] * vNN[0] * flux * detJ * tp; // Nodal flux
        vFq[1] = vFq[1] + vW[j] * vNN[1] * flux * detJ * tp;
    }

    for(int j = 0; j < gauss; j++){
        vFq[j] = -vFq[j]; // Define flux as negative integrals
    }
    vF[node1] += vFq[0];
    vF[node2] += vFq[1];

}

/**
*
* Apply boundary conditions
* Essential boundary conditions
*
*/
int vTempNodes[nelem_x + 1] = {};
for(int i = 0; i < nelem_x+1; i++){
    vTempNodes[i] = vNodeTopo[0][i]; // Nodes at the bottom edge of the beam
}

int nTempNodes = nelem_x + 1; // NUmber of nodes with temp boundary condition

double vBC[nTempNodes*2] = {}; // Initialize the nodal temperature vector
int T0 = 10; // Temperature at boundary
for(int i = 0; i < nTempNodes; i++){
    vBC[i * 2] = vTempNodes[i];
    vBC[i * 2 + 1] = T0;
}

/**
*
* Assembling global "Force" vector
*
*/
int vOrgDof[nDof] = {}; // Original Dof number
double vT[nDof] = {}; // Initialize nodal temperature vector
int rDof = nDof; // Reduced number of DOF

int vInd[nTempNodes] = {};
for(int i = 0; i < nTempNodes; i++){
    vInd[i] = vBC[i * 2];
}
for(int i = 0; i < nTempNodes; i++){
    vOrgDof[vInd[i]] = -1;
    vT[vInd[i]] = vBC[i * 2 + 1];
}

rDof = rDof - nTempNodes;

int vRedDof[rDof] = {};
int counter1 = 0;
for(int j = 0; j < nDof; j++){
    if(vOrgDof[j] == 0){
        vOrgDof[j] = counter1;
        vRedDof[counter1] = j;
        counter1 += 1;
    }
}

/**
*
*  Partition matrices
*
*/
int lenT = sizeof(vT) / sizeof(vT[0]);
bool *mask_E = new bool[lenT];
int ee = 0;
for(int i = 0; i < lenT; i++){
    mask_E[i] = 0;
    for(int j = 0; j < nTempNodes; j++){
        if(i == vTempNodes[j]){
            mask_E[i] = 1;  // Known temperature Dof
            ee++;
            break;
        }
    }
}

double vT_E[ee] = {};
double vF_F[nDof - ee] = {};
for(int i = 0, j = 0, k = 0; i < nDof; i++){
    if(mask_E[i] == 1)
        vT_E[j++] = vT[i];
    else
        vF_F[k++] = vF[i];
}

double vK_EE[ee * ee] = {};
double vK_FF[(nDof - ee) * (nDof - ee)] = {};
double vK_EF[ee * (nDof - ee)] = {};
for(int i = 0, k = 0, m = 0, n = 0; i < nDof; i++){
    if(mask_E[i] == 1)
        for(int j = 0; j < nDof; j++){
            if(mask_E[j] == 1)
                vK_EE[k++] = vK[i * nDof + j];
            else
                vK_EF[m++] = vK[i * nDof + j];
        }
    else
        for(int l = 0; l < nDof; l++){
            if(mask_E[l] == 0){
                vK_FF[n++] = vK[i * nDof + l];
            }
        }
}

/**
*
* Solve for d_F
*
*/
double vRhs[nDof - ee] = {};
double vA[nDof - ee] = {};
double vK_EFT[(nDof - ee) * ee] = {};
for(int i = 0; i < nDof - ee; i++){
    for(int j = 0; j < ee; j++){
        vK_EFT[i * ee + j] = vK_EF[j * (nDof - ee) + i];
    }
}
cblas_dgemv(CblasRowMajor, CblasNoTrans, nDof - ee, ee, 1.0, vK_EFT, ee, vT_E, 1, 0.0, vA, 1);
for(int i = 0; i < nDof - ee; i++){
    vRhs[i] = vF_F[i] - vA[i];
}

double vT_F[nDof - ee] = {};
int vIpiv[nDof - ee] = {};
int info = 0;
F77NAME(dgesv)(nDof - ee, 1, vK_FF, nDof - ee, vIpiv, vRhs, nDof - ee, info);
for(int i = 0; i < nDof - ee; i++){
    vT_F[i] = vRhs[i];
}

/**
*
*Reconstruct the global displacement d
*
*/
for(int i = 0, m = 0, n = 0; i < nDof; i++){
    if(mask_E[i] == 1)
        vT[i] = vT_E[m++];
    else
        vT[i] = vT_F[n++];
}

/**
*
*Compute the reaction f_E
*
*/
double vB[ee];
double vC[ee];
double vF_E[ee];
cblas_dgemv(CblasRowMajor, CblasNoTrans, ee, ee,1.0, vK_EE, ee, vT_E, 1, 0.0, vB, 1);
cblas_dgemv(CblasRowMajor, CblasNoTrans, ee, nDof - ee, 1.0, vK_EF, nDof - ee, vT_F, 1, 0.0, vC, 1);
for(int i = 0; i < ee; i++){
    vF_E[i] = vB[i] + vC[i];
}

/**
*
*Reconstruct the global reactions f
*
*/
for(int i = 0, m = 0, n = 0; i < nDof; i++){
    if(mask_E[i] == 1)
        vF[i] = vF_E[m++];
    else
        vF[i] = vF_F[n++];

}
ofstream myfile2 ("disp2.vtk");
if (myfile2.is_open()){
    myfile2 << "# vtk DataFile Version 4.0\n";
    myfile2 << "vtk output\n";
    myfile2 << "ASCII\n";
    myfile2 << "DATASET UNSTRUCTURED_GRID\n";
    myfile2 << "POINTS " << nnode << " double\n";
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < nDof; j++){
            myfile1 << fixed <<setprecision(3)<< vCoord[i][j] << " ";
        }
        myfile2  << "0.0 ";
    }
    myfile2 << "\n";
    myfile2 << "CELLS " << nelem << " "<< (nnode_elem+1)*nelem << "\n";


    for (int i = 0; i < nelem; i++) {
        for (int j = 0; j < nnode_elem + 1; j++) {
            myfile2 << ElemNode1[i][j] << " ";
        }
        myfile2 << "\n";
    }
    myfile2 << "CELL_TYPES " << nelem <<"\n";
    for (int i = 0; i < nelem; i++) {
        myfile2 << "9\n";
    }

    myfile2 << "POINT_DATA " << nDof <<"\n";
    myfile2 << "FIELD FieldData 2" << "\n";
    myfile2 << "disp 2 " << nDof << " double\n";
    for (int i = 0; i < nDof; i++) {
        myfile2 << fixed << setprecision(5) << T[i] <<" ";
    }
    }
    myfile2.close();
return 0;
}
