/**
* Import libraries
*/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

using namespace std;

#include <cblas.h>

/**
 * @brief Prints a square matrix supplied in column-major full storage.
 *
 * @param   n     Matrix dimension
 * @param   H     Pointer to matrix storage of size n*n
 */
void PrintMatrix(int n, double *H){
    cout.precision(4);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << setw(6) << H[j*n+i] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

/**
 * @brief Prints a vector
 *
 * @param   n       Vector dimension
 * @param   u       Pointer to vector storage of length n
 * @param   dx      Grid spacing
 */
void PrintVector(int n, double *u, double dx){
    cout.precision(4);
    for (int i = 0; i < n; ++i) {
        cout << setw(6) << i*dx << "  " << u[i] << endl;
    }
    cout << endl;
}

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
double a = 0.25;                   // constant in the polynomial describing the beam height
double h1 = 1.;                    //[m] Height at beam left edge
double h2 = h1 * 1.3;              // [m] Height at beam right edge
double L = 3. * h1;                // [m] Beam length
double b = -a * L + (h2 - h1) / L; // constant in the polynomial describing the beam height

/**
* Meshing Geometry
*/
int nelem_y = 8;      // Number of element in y-direction
int nelem_x = 15;     // Number of element in x-direction
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
    vY[i][j] = -0.5 * vH[j] + ((double) i) / ((double) 8) * vH[j];
    }
}

double vCoord[144][2];  // Coordinate of the nodes [x, y]
for(int i = 0; i < 9; i++){
    for(int j = 0; j < 16; j++){
       vCoord[j*9 + i][0] = vX[j];
    }
}
for(int j = 0; j < 16; j++){
    for(int i = 0; i < 9; i++){
       vCoord[j*9 + i][1] = vY[i][j];
    }
}

/**
* Calculation of topology matrix NodeTopo
*/
int vNodeTopo[nelem_y + 1][nelem_x+1];
for(int i = 0; i < nelem_y + 1; i++){
    for(int j = 0; j < nelem_x + 1; j++){
       vNodeTopo[i][j] = j * 9 + i;
    }
}

/**
* Calculation of topology matrix Elemnode
*/
int vElemNode[nelem][5];
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
    for(int j =1; j < 5; j++){
    vElemX[i][j-1] = vCoord[vElemNode[i][j]][0];
    }
}

double vElemY[nelem][nnode_elem];  // Element y nodal coordinates
for(int i = 0; i < nelem; i++){
    for(int j =1; j < 5; j++){
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
double vGP[2] = {-1 / sqrt(3), 1 / sqrt(3)}; // Points
int vW[2] = {1, 1}; // Weights
double vD[4] = {kx, kxy, kxy, ky}; // Conductivity matrix D
double vK[nDof*nDof] = {}; // Initiation of global stiffness matrix K
int gauss = gaussorder;

for(int i = 0; i < nelem; i++){
    int vNodes[4];  // Element nodes
    double eCoord[nnode_elem*gauss] = {}; // node coordinates

    int iCoord = 0;
    for(int index = 1; index < 5; index++){

        vNodes[index-1] = vElemNode[i][index];

        eCoord[iCoord] = vCoord[vNodes[index-1]][0];
        eCoord[iCoord+1] = vCoord[vNodes[index-1]][1];

        iCoord = iCoord + 2;
    }
    //cout << eCoord[0] << " " << eCoord[1] << " " << eCoord[2] << " " << eCoord[3] << " " << eCoord[4] << " " << eCoord[5] << " " << eCoord[6] << " " << eCoord[7] << endl;
    int gDof[nnode_elem];
    for(int j = 0; j < nnode_elem; j++){
        gDof[j] = vNodes[j]; // global dof for node j
    }

    double vKe[nnode_elem*nnode_elem] = {}; // Local stiffness matrix
    double vDetJ[gauss*gauss] = {}; // For storing the determinants of J
    double vGN[gauss*nnode_elem] = {};


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
            //cout << vGN[0] << " " << vGN[1] << " " << vGN[2] << " " << vGN[3] << " " << vGN[4] << " " << vGN[5] << " " << vGN[6] << " " << vGN[7] << endl;
            double J[gauss*gauss] = {}; //Jacobian matrix
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, gauss, gauss, nnode_elem, 1.0, vGN, nnode_elem, eCoord, gauss, 0.0, J, gauss);

            cout << (int)(J[0] * 10000.0) / 10000.0 << " " << (int)(J[1] * 10000.0) / 10000.0 << endl;
            cout << J[2] << " " << J[3] << endl;
            cout << endl;
        }
    }



}



return 0;
}
