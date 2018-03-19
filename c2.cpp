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

int c2{
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
    vFluxNodes[i] = vNodeTopo[0][i];  // Nodes at the top edge of the beam
    nFluxNodes += 1;  // Number of nodes on the top edge of the beam
}

/**
*Defining load
*/
int q = 2500; // Constant flux at right edge of the beam
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
int vTempNodes[nelem_x + 1] = {};
for(int i = 0; i < nelem_x+1; i++){
    vTempNodes[i] = vNodeTopo[nelem_y][i]; // Nodes at the bottom edge of the beam
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
int mask_E = nTempNodes; // Known temperature Dof
int mmask_E = nDof - nTempNodes;
double vT_E[mask_E] = {};
double vF_F[mmask_E] = {};
double vK_EE[mask_E * mask_E] = {};
double vK_FF[mmask_E * mmask_E] = {};
double vK_EF[mask_E * mmask_E] = {};
double vK_EFT[mmask_E * mask_E] = {};

for(int i = 0; i < mask_E; i++){
    vT_E[i] = vT[i];
}

for(int i = 0; i < mmask_E; i++){
    vF_F[i] = vF[i + mask_E];
}

for(int i = 0; i < mask_E; i++){
    for(int j = 0; j < mask_E; j++){
        vK_EE[i * mask_E + j] = vK[i * nDof + j];
    }
}

for(int i = 0; i < mmask_E; i++){
    for(int j = 0; j < mmask_E; j++){
        vK_FF[i * mmask_E + j] = vK[(i + mask_E) * nDof + (j + mask_E)];
    }
}

for(int i = 0; i < mask_E; i++){
    for(int j = 0; j < mmask_E; j++){
        vK_EF[i * mmask_E + j] = vK[i * nDof + (j + mask_E)];
    }
}

for(int i = 0; i < mmask_E; i++){
    for(int j = 0; j < mask_E; j++){
        vK_EFT[i * mask_E + j] = vK_EF[j * mmask_E + i];
    }
}

/**
*
* Solve for d_F
*
*/
double vRhs[mmask_E] = {};
double vA[mmask_E] = {};
cblas_dgemv(CblasRowMajor, CblasNoTrans, mmask_E, mask_E, 1.0, vK_EFT, mask_E, vT_E, 1, 0.0, vA, 1);
for(int i = 0; i < mmask_E; i++){
    vRhs[i] = vF_F[i] - vA[i];
}

double vT_F[mmask_E] = {};
int vIpiv[mmask_E] = {};
int info = 0;
F77NAME(dgesv)(mmask_E, 1, vK_FF, mmask_E, vIpiv, vRhs, mmask_E, info);
for(int i = 0; i < mmask_E; i++){
    vT_F[i] = vRhs[i];
}

/**
*
*Reconstruct the global displacement d
*
*/
for(int i = 0; i < mask_E; i++){
    vT[i] = vT_E[i];
}
for(int i = 0; i < mmask_E; i++){
    vT[i + mask_E] = vT_F[i];
}

/**
*
*Compute the reaction f_E
*
*/
double vB[mask_E];
double vC[mask_E];
double vF_E[mask_E];
cblas_dgemv(CblasRowMajor, CblasNoTrans, mask_E, mask_E,1.0, vK_EE, mask_E, vT_E, 1, 0.0, vB, 1);
cblas_dgemv(CblasRowMajor, CblasNoTrans, mask_E, mmask_E, 1.0, vK_EF, mmask_E, vT_F, 1, 0.0, vC, 1);
for(int i = 0; i < mask_E; i++){
    vF_E[i] = vB[i] + vC[i];
}

/**
*
*Reconstruct the global reactions f
*
*/
for(int i = 0; i < mask_E; i++){
    vF[i] = vF_E[i];
}
for(int i = 0; i < mmask_E; i++){
    vF[i + mask_E] = vF_F[i];
}

return 0;
}
