/**
*
* Import libraries
*
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

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

int c1{
/**
* Compute nodal boundary flux vector
* Natural boundary condition
*/

/**
*Defined on edges
*/
int vFluxNodes[nelem_y + 1] = {};
int nFluxNodes = 0;
for(int i = 0; i < nelem_y + 1; i++){
    vFluxNodes[i] = vNodeTopo[0][i]; // Nodes at the bottom edge of the beam
    nFluxNodes += 1;  // Number of nodes on the bottom edge of the beam
}

/**
*Defining load
*/
int q = 2500; // Constant flux at bottom edge of the beam
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
        vFq[0] = vFq[0] + vW[j] * vNN[0] * flux * detJ * th; // Nodal flux
        vFq[1] = vFq[1] + vW[j] * vNN[1] * flux * detJ * th;
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
int vTempNodes[nelem_y+1] = {};
for(int i = 0; i < nelem_y+1; i++){
    vTempNodes[i] = vNodeTopo[i][0]; // Nodes at the left edge of the beam
}

int nTempNodes = nelem_y + 1; // NUmber of nodes with temp boundary condition

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


return 0;
}
