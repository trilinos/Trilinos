// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert Kirby (robert.c.kirby@ttu.edu)
//
// ************************************************************************
// @HEADER

/** \file
\brief  Convergence tests for high-order bases using various integration rules.
\author Created by D. Ridzal and M. Keegan.
*/
#include <iostream>

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_Flops.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_SerialComm.h"
#include "Epetra_CrsSingletonFilter.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_Reindex_LinearProblem.h"

// Teuchos includes
#include "Teuchos_RefCountPtr.hpp"

// Amesos includes
#include "Amesos.h"
#include "Amesos_Klu.h"
#include "Amesos_Umfpack.h"

// Intrepid includes
#include "Intrepid_Types.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CubatureTensor.hpp"
#include "Intrepid_CubatureSparse.hpp"
#include "Intrepid_CubatureGenSparse.hpp"
#include "Intrepid_F0_HEX_DD.hpp"
#include "Intrepid_LocalForm0.hpp"

using namespace std;
using namespace Intrepid;

// forward declarations
bool FEVector2MATLAB( const Epetra_FEVector &, ostream &);
bool CrsMatrix2MATLAB(const Epetra_CrsMatrix &, ostream &);
// used to select one of several RHS functionals
enum RHSType {
  SINSINSIN,
  SINSINSINEXP,
  POLY_Q,
  POLY_P
};
double rhsFunc(Point<double>, RHSType, int);
double dirichletFunc(RHSType, int, double, Point<double>, int);
double u_exact(Point<double>, RHSType, int);
void   u_grad(double*, Point<double>, RHSType, int);

// implements RHS functionals
double rhsFunc(Point<double> point, RHSType rhstype, int deg) {
  double val = 0.0; 
  switch (rhstype) {
    case SINSINSIN:
      val = 3.0*M_PI*M_PI*sin(M_PI*point[0])*sin(M_PI*point[1])*sin(M_PI*point[2]);
      break;
    case SINSINSINEXP: {
      double exyz = exp(point[0]+point[1]+point[2]);
      double s0 = sin(M_PI*point[0]), s1 = sin(M_PI*point[1]), s2 = sin(M_PI*point[2]);
      val =  3.0*M_PI*M_PI*s0*s1*s2*exyz
            -2.0*M_PI*cos(M_PI*point[0])*s1*s2*exyz
            -2.0*M_PI*s0*cos(M_PI*point[1])*s2*exyz
            -2.0*M_PI*s0*s1*cos(M_PI*point[2])*exyz
            -3.0*s0*s1*s2*exyz;
      break;
    }
    case POLY_Q: {
      double px=0, py=0, pz=0, pxdd=0, pydd=0, pzdd=0;
      for (int i=0; i<=deg; i++) {
        px += std::pow(point[0], i);
        py += std::pow(point[1], i);
        pz += std::pow(point[2], i);
      }      
      for (int i=0; i<deg-1; i++) {
        pxdd += (i+1)*(i+2)*std::pow(point[0], i);
        pydd += (i+1)*(i+2)*std::pow(point[1], i);
        pzdd += (i+1)*(i+2)*std::pow(point[2], i);
      }
      val = - pxdd*py*pz
            - pydd*px*pz
            - pzdd*px*py;
      break;
    }
    case POLY_P: {
      for (int i=0; i<=deg; i++) {
        for (int j=0; j<=deg-i; j++) {
          for (int k=0; k<=deg-i-j; k++) {
            if (i>1)
              val -= i*(i-1)*std::pow(point[0], i-2)*std::pow(point[1], j)*std::pow(point[2], k);
            if (j>1)
              val -= j*(j-1)*std::pow(point[0], i)*std::pow(point[1], j-2)*std::pow(point[2], k);
            if (k>1)
              val -= k*(k-1)*std::pow(point[0], i)*std::pow(point[1], j)*std::pow(point[2], k-2);
          }
        }
      }      
      break;
    }
    default:
       TEST_FOR_EXCEPTION( ( 1 ),
                           std::invalid_argument,
                           ">>> ERROR (rhsFunc): Invalid type of RHS.");
  }
  return val;
}


// implements Dirichlet bdry values
double dirichletFunc(RHSType rhstype, int deg, double sideval, Point<double> point, int side) {
  double val = 0.0;
  switch (rhstype) {
    case SINSINSIN:
    case SINSINSINEXP:
      break;
    case POLY_Q: {
      double px=0, py=0, pz=0, mult=0;
      for (int i=0; i<=deg; i++) {
        px   += std::pow(point[0], i);
        py   += std::pow(point[1], i);
        pz   += std::pow(point[2], i);
        mult += std::pow(sideval, i);
      }
      switch (side) {
        case 0:
        case 2:
          val = mult*px*pz;
          break;
        case 1:
        case 3:
          val = mult*py*pz;
          break;
        case 4:
        case 5:
          val = mult*px*py;
          break;
      }
      break;
    }
    case POLY_P: {
      val = u_exact(point, rhstype, deg);
      break;
    }
    default:
       TEST_FOR_EXCEPTION( ( 1 ),
                           std::invalid_argument,
                           ">>> ERROR (rhsFunc): Invalid type of RHS.");
  }
  return val;
}


// implements exact solution of PDE
double u_exact(Point<double> point, RHSType rhstype, int deg) {
  double val = 0.0; 
  switch (rhstype) {
    case SINSINSIN:
      val = sin(M_PI*point[0])*sin(M_PI*point[1])*sin(M_PI*point[2]);
      break;
    case SINSINSINEXP:
      val = sin(M_PI*point[0])*sin(M_PI*point[1])*sin(M_PI*point[2])*exp(point[0]+point[1]+point[2]);
      break;
    case POLY_Q: {
      double px=0, py=0, pz=0;
      for (int i=0; i<=deg; i++) {
        px += std::pow(point[0], i);
        py += std::pow(point[1], i);
        pz += std::pow(point[2], i);
      }      
      val = px * py * pz;
      break;
    }
    case POLY_P: {
      for (int i=0; i<=deg; i++) {
        for (int j=0; j<=deg-i; j++) {
          for (int k=0; k<=deg-i-j; k++) {
            val += std::pow(point[0], i)*std::pow(point[1], j)*std::pow(point[2], k);
          }
        }
      }      
      break;
    }
    default:
      TEST_FOR_EXCEPTION( ( 1 ),
                          std::invalid_argument,
                          ">>> ERROR (u_exact): Invalid type of RHS.");
  }
  return val;
}


// implements exact gradients of PDE solution
void u_grad(double* grad, Point<double> point, RHSType rhstype, int deg) {
  switch (rhstype) {
    case SINSINSIN:
      //val = sin(M_PI*point[0])*sin(M_PI*point[1])*sin(M_PI*point[2]);
      grad[0] = M_PI*cos(M_PI*point[0])*sin(M_PI*point[1])*sin(M_PI*point[2]);
      grad[1] = M_PI*cos(M_PI*point[1])*sin(M_PI*point[0])*sin(M_PI*point[2]);
      grad[2] = M_PI*cos(M_PI*point[2])*sin(M_PI*point[0])*sin(M_PI*point[1]);
      break;
    case SINSINSINEXP: {
      double val = sin(M_PI*point[0])*sin(M_PI*point[1])*sin(M_PI*point[2])*exp(point[0]+point[1]+point[2]);
      grad[0] = M_PI*cos(M_PI*point[0])*sin(M_PI*point[1])*sin(M_PI*point[2])*exp(point[0]+point[1]+point[2]) + val;
      grad[1] = M_PI*cos(M_PI*point[1])*sin(M_PI*point[0])*sin(M_PI*point[2])*exp(point[0]+point[1]+point[2]) + val;
      grad[2] = M_PI*cos(M_PI*point[2])*sin(M_PI*point[0])*sin(M_PI*point[1])*exp(point[0]+point[1]+point[2]) + val;
      break;
    }
    case POLY_Q: {
      double px=0, py=0, pz=0, pxd=0, pyd=0, pzd=0;
      for (int i=0; i<=deg; i++) {
        px += std::pow(point[0], i);
        py += std::pow(point[1], i);
        pz += std::pow(point[2], i);
      }
      for (int i=0; i<=deg-1; i++) {
        pxd += (i+1)*std::pow(point[0], i);
        pyd += (i+1)*std::pow(point[1], i);
        pzd += (i+1)*std::pow(point[2], i);
      }
      grad[0] = pxd*py*pz;
      grad[1] = pyd*px*pz;
      grad[2] = pzd*px*py;
      break;
    }
    case POLY_P: {
      grad[0] = 0.0;
      grad[1] = 0.0;
      grad[2] = 0.0;
      for (int i=0; i<=deg; i++) {
        for (int j=0; j<=deg-i; j++) {
          for (int k=0; k<=deg-i-j; k++) {
            if (i>0)
              grad[0] += i*std::pow(point[0], i-1)*std::pow(point[1], j)*std::pow(point[2], k);
            if (j>0)
              grad[1] += j*std::pow(point[0], i)*std::pow(point[1], j-1)*std::pow(point[2], k);
            if (k>0)
              grad[2] += k*std::pow(point[0], i)*std::pow(point[1], j)*std::pow(point[2], k-1);
          }
        }
      }
      break;
    }
    default:
      TEST_FOR_EXCEPTION( ( 1 ),
                          std::invalid_argument,
                          ">>> ERROR (u_exact): Invalid type of RHS.");
  }
}


int main(int argc, char *argv[]) {
  // Check number of arguments.
  TEST_FOR_EXCEPTION( ( argc < 8 ),
                      std::invalid_argument,
                      ">>> ERROR (example_01): Invalid number of arguments. See source for proper calling sequence.");

  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 7)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  *outStream \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|              Example: Convergence study for high-order H1 spaces            |\n" \
  << "|                                                                             |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
  << "|                      Robert Kirby (robert.c.kirby@ttu.edu).                 |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n";

  Epetra_SerialComm Comm;

  try {
    // Gather command-line parameters.
    int cubType  = atoi(argv[1]); // cubature type
    int basisDeg = atoi(argv[2]); // basis degree
    int cubDeg   = atoi(argv[3]); // cubature degree
    int dataDeg  = atoi(argv[4]); // polynomial degree of PDE solution
    int xint     = atoi(argv[5]); // num intervals in x direction (assumed box domain)
    int yint     = atoi(argv[6]); // num intervals in y direction (assumed box domain)
    int zint     = atoi(argv[7]); // num intervals in z direction (assumed box domain)

    std::cout << "Problem parameters:\n";
    std::cout << "   CubType" << "   BasDeg" << "   CubDeg" << "   PolyDeg" << "    NX" << "   NY" << "   NZ\n";
    std::cout << std::setw(7) << cubType <<
                 std::setw(9) << basisDeg <<
                 std::setw(10) << cubDeg <<
                 std::setw(9) << dataDeg <<
                 std::setw(9) << xint <<
                 std::setw(5) << yint <<
                 std::setw(5) << zint << "\n\n";
	
    /*** BEGIN SELECTION OF PROBLEM PARAMETERS ***/

    // Choose basis:
    Basis_F0_HEX_DD<double> myBasis( basisDeg );
    myBasis.initialize();
    Teuchos::RCP<Basis<double> > hexBasis = Teuchos::rcp( &myBasis, false );

    // Choose integration rule:
    Teuchos::Array<Teuchos::Array< Teuchos::RCP<Cubature<double> > > > hexCubature;
    Teuchos::RCP<Cubature<double> > cellCub;
    switch (cubType) {
      case 1:
        cellCub = Teuchos::rcp(new CubatureTensor<double>(CELL_HEX, cubDeg) );
        break;
      case 2:
        cellCub = Teuchos::rcp(new CubatureSparse<double,3>(cubDeg) );
        break;
      case 3:
        cellCub = Teuchos::rcp(new CubatureGenSparse<double,3>(cubDeg) );
        break;
      default:
        TEST_FOR_EXCEPTION( ( 1 ),
                            std::invalid_argument,
                            ">>> ERROR (cubType): Invalid type of cubature rule.");
    }
    hexCubature.resize(1);
    hexCubature[0].resize(1);
    hexCubature[0][0] = cellCub;

    // Choose computational engine:
    ECompEngine compEng = COMP_BLAS;

    // Initialize local 0-form:
    LocalForm0<double> hexH1form(hexBasis, hexCubature, compEng);

    // Set domain geometry:
    double leftX = 0.0, rightX = 1.0;
    double leftY = 0.0, rightY = 1.0;
    double leftZ = 0.0, rightZ = 1.0;
    int NX = xint, NY = yint, NZ = zint;
    // Dirichlet boundary sides
    bool side0  = true; // y==leftY
    bool side1  = true; // x==rightX
    bool side2  = true; // y==rightY
    bool side3  = true; // x==leftX
    bool side4  = true; // z=leftZ, bottom
    bool side5  = true; // z=rightZ, top

    // RHS function specification:
    RHSType rhstype = SINSINSINEXP; 

    /*** END SELECTION OF PROBLEM PARAMETERS ***/



    /*** Begin data structure preparation ***/
    //
    //
      // extract integration points for later use (e.g. in evaluating RHS data)
      int                             numCubPoints;
      Teuchos::Array< Point<double> > cubPoints;
      Teuchos::Array<double>          cubWeights;
      cellCub->getCubature(numCubPoints, cubPoints, cubWeights);

      // compute mesh spacings
      double hx = (rightX-leftX)/((double)NX);
      double hy = (rightY-leftY)/((double)NY);
      double hz = (rightZ-leftZ)/((double)NZ);
      
      // DOFs
      int localDOFs   = hexBasis->getNumLocalDof();
      int localDOFs1D = (int)std::floor(std::pow((double)localDOFs, (1.0/3.0))+0.1);
      int globalDOFsX = (localDOFs1D-1)*NX+1;
      int globalDOFsY = (localDOFs1D-1)*NY+1;
      int globalDOFsZ = (localDOFs1D-1)*NZ+1;
      int globalDOFs  = globalDOFsX*globalDOFsY*globalDOFsZ;

      // other variables
      double hexNodes[24];
      FieldContainer<double> localGradGrad;
      FieldContainer<double> localRHS;
      FieldContainer<double> trialData(1, numCubPoints);
      Point<double>          transPoint(3);
      FieldContainer<int>    insertList(localDOFs);
      int xd = globalDOFsX, yd = globalDOFsY, zd = globalDOFsZ, ld = localDOFs1D;
      std::vector<int> dirichletDOFs;
      std::vector<int> freeDOFs;

      // operators and functionals
      Epetra_Map allDOFmap(globalDOFs, 0, Comm);
      Teuchos::RCP<Epetra_FECrsMatrix> A = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, allDOFmap, (ld+1)*(ld+1)*(ld+1) ));
      Teuchos::RCP<Epetra_FEVector>    b = Teuchos::rcp(new Epetra_FEVector(allDOFmap));
      Teuchos::RCP<Epetra_FEVector>   uh = Teuchos::rcp(new Epetra_FEVector(allDOFmap));
      Teuchos::RCP<Epetra_FEVector>   Au = Teuchos::rcp(new Epetra_FEVector(allDOFmap));
    //
    // 
    /*** End data structure preparation ***/




    // Loop over the mesh and assemble local operators:
    // reverse lexicographical order
    std::cout << "\nAssembling ...\n";
    for (int i=0; i<NX; i++) {
      for (int j=0; j<NY; j++) {
        for (int k=0; k<NZ; k++) {
          double ih, jh, kh, ihh, jhh, khh;
          ih = (double)i*hx; jh = (double)j*hy; kh = (double)k*hz;
          ihh = ih+hx;       jhh = jh+hy;       khh = kh+hz;
          hexNodes[0]  = ih;  hexNodes[1]  = jh;  hexNodes[2]  = kh;
          hexNodes[3]  = ihh; hexNodes[4]  = jh;  hexNodes[5]  = kh;
          hexNodes[6]  = ihh; hexNodes[7]  = jhh; hexNodes[8]  = kh;
          hexNodes[9]  = ih;  hexNodes[10] = jhh; hexNodes[11] = kh;
          hexNodes[12] = ih;  hexNodes[13] = jh;  hexNodes[14] = khh;
          hexNodes[15] = ihh; hexNodes[16] = jh;  hexNodes[17] = khh;
          hexNodes[18] = ihh; hexNodes[19] = jhh; hexNodes[20] = khh;
          hexNodes[21] = ih;  hexNodes[22] = jhh; hexNodes[23] = khh;

          MultiCell<double> mCell(1, CELL_HEX, hexNodes);
          
          // compute local stiffness contribution
          hexH1form.getOperator(localGradGrad, OPERATOR_GRAD, OPERATOR_GRAD, mCell);

          // compute local RHS contribution
          for (int qp=0; qp<numCubPoints; qp++) {
            // first, transform integration points into physical space
            transPoint = mCell.mapToPhysicalCell(0, cubPoints[qp]);
            // second, evaluate RHS functional at transformed integration points
            trialData(0, qp) = rhsFunc(transPoint, rhstype, dataDeg);
          }
          hexH1form.getFunctional(localRHS, trialData, OPERATOR_VALUE, mCell);
 
          // insert into global operator
          int count = 0;
          int aid   = i*(yd*zd*(ld-1)) + j*(zd*(ld-1)) + k*(ld-1);
          for (int li=0; li<ld; li++) {
            for (int lj=0; lj<ld; lj++) {
              for (int lk=0; lk<ld; lk++) {

                int gid = aid + li*(yd*zd) + lj*zd + lk;
                insertList[count] = gid;
                count++;

                Point<double> dirPoint(hexNodes[0] + (double)li*(hx/(ld-1)),
                                       hexNodes[1] + (double)lj*(hy/(ld-1)),
                                       hexNodes[2] + (double)lk*(hz/(ld-1)));

                // Dirichlet boundary conditions
                if (side0 && (j==0) && (lj==0)) {       // y==leftY
                  dirichletDOFs.push_back(gid);
                  double dVal = dirichletFunc(rhstype, dataDeg, leftY, dirPoint, 0);
                  uh->ReplaceGlobalValues(1, &gid, &dVal);
                }
                if (side1 && (i==NX-1) && (li==ld-1)) { // x==rightX
                  dirichletDOFs.push_back(gid);
                  double dVal = dirichletFunc(rhstype, dataDeg, rightX, dirPoint, 1);
                  uh->ReplaceGlobalValues(1, &gid, &dVal);
                }
                if (side2 && (j==NY-1) && (lj==ld-1)) { // y==rightY
                  dirichletDOFs.push_back(gid);
                  double dVal = dirichletFunc(rhstype, dataDeg, rightY, dirPoint, 2);
                  uh->ReplaceGlobalValues(1, &gid, &dVal);
                }
                if (side3 && (i==0) && (li==0)) {       // x==leftX
                  dirichletDOFs.push_back(gid);
                  double dVal = dirichletFunc(rhstype, dataDeg, leftX, dirPoint, 3);
                  uh->ReplaceGlobalValues(1, &gid, &dVal);
                }
                if (side4 && (k==0) && (lk==0)) {       // z=leftZ, bottom
                  dirichletDOFs.push_back(gid);
                  double dVal = dirichletFunc(rhstype, dataDeg, leftZ, dirPoint, 4);
                  uh->ReplaceGlobalValues(1, &gid, &dVal);
                }
                if (side5 && (k==NZ-1) && (lk==ld-1)) { // z=rightZ, top
                  dirichletDOFs.push_back(gid);
                  double dVal = dirichletFunc(rhstype, dataDeg, rightZ, dirPoint, 5);
                  uh->ReplaceGlobalValues(1, &gid, &dVal);
                }

              }
            }
          }
          A->InsertGlobalValues(localDOFs, &insertList[0], &localGradGrad[0], Epetra_FECrsMatrix::ROW_MAJOR);
          b->SumIntoGlobalValues(localDOFs, &insertList[0], &localRHS[0]);

          *outStream << mCell;
          *outStream << localGradGrad;
        }
      }
    }

    A->GlobalAssemble();
    b->GlobalAssemble();
    uh->GlobalAssemble();
    Au->GlobalAssemble();

    // compute 'free' DOFs
    std::sort(dirichletDOFs.begin(), dirichletDOFs.end());
    std::vector<int>::iterator p = std::unique(dirichletDOFs.begin(), dirichletDOFs.end());
    dirichletDOFs.erase(p, dirichletDOFs.end());
    std::vector<int> allDOFs(globalDOFs);
    freeDOFs.resize(globalDOFs-(int)dirichletDOFs.size());
    for (unsigned ii=0; ii<globalDOFs; ii++) {
      allDOFs[ii] = ii;
    }
    std::set_difference(allDOFs.begin(), allDOFs.end(), dirichletDOFs.begin(), dirichletDOFs.end(), freeDOFs.begin());

    // Modify b to account for Dirichlet conditions
    A->Multiply(false, *uh, *Au);
    b->Update(-1.0, *Au, 1.0);

    // Prepare reduced data structures (Dirichlet DOFs are excluded)
    std::cout << "\nEliminating Dirichlet data ...\n";
    Epetra_Map reducedMap(freeDOFs.size(), freeDOFs.size(), &freeDOFs[0], 0, Comm);
    // it is important to specify the column map for the reduced matrix
    Teuchos::RCP<Epetra_FECrsMatrix> A_red = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, reducedMap, reducedMap, localDOFs));
    Teuchos::RCP<Epetra_FEVector>    b_red = Teuchos::rcp(new Epetra_FEVector(reducedMap));
    Teuchos::RCP<Epetra_FEVector>    uh_red = Teuchos::rcp(new Epetra_FEVector(reducedMap));
    Epetra_Import importer(reducedMap, allDOFmap);
    A_red->Import(*A, importer, Insert);
    b_red->Import(*b, importer, Insert);

    A_red->GlobalAssemble();
    b_red->GlobalAssemble();
    uh_red->GlobalAssemble();

    CrsMatrix2MATLAB(*A, *outStream);
    CrsMatrix2MATLAB(*A_red, *outStream);

    // Set up reduced linear problem
    std::cout << "\nSolving ...\n";
    Epetra_LinearProblem Problem;
    Problem.SetOperator(&(*A_red));
    Problem.SetRHS(&(*b_red));
    Problem.SetLHS(&(*uh_red));

    // Reindex linear problem (Amesos_Umfpack doesn't understand non-contiguous indices)
    EpetraExt::LinearProblem_Reindex reindex(NULL);
    Epetra_LinearProblem newProblem = reindex(Problem);

    // Pick solver and solve reindexed reduced problem
    Amesos_Umfpack solver(newProblem);
    AMESOS_CHK_ERR(solver.SymbolicFactorization());
    solver.PrintStatus();
    AMESOS_CHK_ERR(solver.NumericFactorization());
    AMESOS_CHK_ERR(solver.Solve());
    solver.PrintTiming();

    // Assemble full solution vector
    uh->ReplaceGlobalValues(freeDOFs.size(), &freeDOFs[0], uh_red->Values() );




    /*** Begin error computation ***/

      // Choose integration rule, possibly different from the above:
      Teuchos::RCP<Cubature<double> > errCub = Teuchos::rcp(new CubatureTensor<double>(CELL_HEX, basisDeg*2) );

      // Choose computational engine:
      compEng = COMP_BLAS;

      // extract integration points
      errCub->getCubature(numCubPoints, cubPoints, cubWeights);

      // set up basis functions and gradients at integration points
      FieldContainer<double> u_refvals;
      FieldContainer<double> u_refgrads;
      hexBasis->getValues(u_refvals, cubPoints, OPERATOR_VALUE);
      hexBasis->getValues(u_refgrads, cubPoints, OPERATOR_GRAD);

      FieldContainer<double> exactSoln(numCubPoints);
      FieldContainer<double> approxSoln(numCubPoints);
      FieldContainer<double> exactGrads(numCubPoints, 3);
      FieldContainer<double> approxGrads(numCubPoints, 3);
      FieldContainer<double> approxGradsTmp(numCubPoints, 3);

      double errL2 = 0.0;
      double errH1 = 0.0;
      double errLinf = 0.0;
      
      std::cout << "\nComputing L2, Linf, H1 error ...\n";
      for (int i=0; i<NX; i++) {
        for (int j=0; j<NY; j++) {
          for (int k=0; k<NZ; k++) {
            double ih, jh, kh, ihh, jhh, khh;
            ih = (double)i*hx; jh = (double)j*hy; kh = (double)k*hz;
            ihh = ih+hx;       jhh = jh+hy;       khh = kh+hz;
            hexNodes[0]  = ih;  hexNodes[1]  = jh;  hexNodes[2]  = kh;
            hexNodes[3]  = ihh; hexNodes[4]  = jh;  hexNodes[5]  = kh;
            hexNodes[6]  = ihh; hexNodes[7]  = jhh; hexNodes[8]  = kh;
            hexNodes[9]  = ih;  hexNodes[10] = jhh; hexNodes[11] = kh;
            hexNodes[12] = ih;  hexNodes[13] = jh;  hexNodes[14] = khh;
            hexNodes[15] = ihh; hexNodes[16] = jh;  hexNodes[17] = khh;
            hexNodes[18] = ihh; hexNodes[19] = jhh; hexNodes[20] = khh;
            hexNodes[21] = ih;  hexNodes[22] = jhh; hexNodes[23] = khh;

            MultiCell<double> mCell(1, CELL_HEX, hexNodes);
            mCell.setAtlas();
            mCell.initializeMeasures(3, 0, cubPoints, cubWeights);
            Teuchos::Array<double> weightedDetJacs = mCell.getWeightedMeasure(0, 3, 0);
            Teuchos::Array<Matrix<double> > JacInvTransp = mCell.getJacobian(0, 3, 0);

            double errL2_loc = 0.0;
            double errH1_loc = 0.0;
            for (int qp=0; qp<numCubPoints; qp++) {
              // first, transform integration points into physical space
              transPoint = mCell.mapToPhysicalCell(0, cubPoints[qp]);
              // second, evaluate exact solution and gradients at transformed integration points
              exactSoln(qp) = u_exact(transPoint, rhstype, dataDeg);
              u_grad(&exactGrads(qp,0), transPoint, rhstype, dataDeg);
              // third, transform Jacobian
              JacInvTransp[qp].transpose();
              JacInvTransp[qp].invert();

              // get approximate solution data
              approxSoln(qp) = 0.0;
              approxGradsTmp(qp,0) = 0.0;
              approxGradsTmp(qp,1) = 0.0;
              approxGradsTmp(qp,2) = 0.0;
              int count = 0;
              int aid   = i*(yd*zd*(ld-1)) + j*(zd*(ld-1)) + k*(ld-1);
              for (int li=0; li<ld; li++) {
                for (int lj=0; lj<ld; lj++) {
                  for (int lk=0; lk<ld; lk++) {
                    int gid = aid + li*(yd*zd) + lj*zd + lk;
                    approxSoln(qp)  += uh->Values()[gid] * u_refvals(qp, count);
                    approxGradsTmp(qp,0) += uh->Values()[gid] * u_refgrads(qp, count, 0);
                    approxGradsTmp(qp,1) += uh->Values()[gid] * u_refgrads(qp, count, 1);
                    approxGradsTmp(qp,2) += uh->Values()[gid] * u_refgrads(qp, count, 2);
                    JacInvTransp[qp].multiplyLeft(&approxGrads(qp,0), &approxGradsTmp(qp,0));
                    count++;
                  }
                }
              }

              errL2_loc += (exactSoln(qp) - approxSoln(qp))*(exactSoln(qp) - approxSoln(qp))*weightedDetJacs[qp];
              errH1_loc += ( (exactGrads(qp,0) - approxGrads(qp,0))*(exactGrads(qp,0) - approxGrads(qp,0)) +
                             (exactGrads(qp,1) - approxGrads(qp,1))*(exactGrads(qp,1) - approxGrads(qp,1)) +
                             (exactGrads(qp,2) - approxGrads(qp,2))*(exactGrads(qp,2) - approxGrads(qp,2)) )*weightedDetJacs[qp];
              errLinf    = max(errLinf, abs(exactSoln(qp) - approxSoln(qp)));

            }

            errL2 += errL2_loc;
            errH1 += errH1_loc;

          }
        }
      }
      errH1 = sqrt(errH1+errL2);
      errL2 = sqrt(errL2);

      cout << "\nerrL2   = " << errL2 << "\n";
      cout << "\nerrH1   = " << errH1 << "\n";
      cout << "\nerrLinf = " << errLinf << "\n\n\n";

    /*** End error computation ***/





    // Output in Matlab format
    CrsMatrix2MATLAB(*A, *outStream);
    FEVector2MATLAB(*b, *outStream);
    FEVector2MATLAB(*uh, *outStream);
    CrsMatrix2MATLAB(*A_red, *outStream);
    FEVector2MATLAB(*b_red, *outStream);
    FEVector2MATLAB(*uh_red, *outStream);

    // Output in native Epetra format
    *outStream << *A;
    *outStream << *b;
    *outStream << *uh;
    *outStream << *A_red;
    *outStream << *b_red;
    *outStream << *uh_red;
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
  };

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return 0;
}



// Output functions (Matlab format)

/* ======== ================ *
 * function CrsMatrix2MATLAB *
 * ======== ================ *
 *
 * Print out a CrsMatrix in a MATLAB format. Each processor prints out
 * its part, starting from proc 0 to proc NumProc-1. The first line of
 * each processor's output states the number of local rows and of
 * local nonzero elements.
 *
 *
 * Return code:        true if matrix has been printed out
 * -----------         false otherwise
 *
 * Parameters:
 * ----------
 *
 * - Epetra_CrsMatrix  reference to the distributed CrsMatrix to
 *                     print out
 * - ostream &         reference to output stream
 */
bool CrsMatrix2MATLAB(const Epetra_CrsMatrix & A, ostream & outfile) {

  int MyPID = A.Comm().MyPID();
  int NumProc = A.Comm().NumProc();

  // work only on transformed matrices;
  if( A.IndicesAreLocal() == false ) {
    if( MyPID == 0 ) {
      cerr << "ERROR in "<< __FILE__ << ", line " << __LINE__ << endl;
      cerr << "Function CrsMatrix2MATLAB accepts\n";
      cerr << "transformed matrices ONLY. Please call A.TransformToLoca()\n";
      cerr << "on your matrix A to that purpose.\n";
      cerr << "Now returning...\n";
    }
    return false;
  }

  int NumMyRows = A.NumMyRows(); // number of rows on this process
  int NumNzRow;   // number of nonzero elements for each row
  int NumEntries; // number of extracted elements for each row
  int NumGlobalRows; // global dimensio of the problem
  int GlobalRow;  // row in global ordering
  int NumGlobalNonzeros; // global number of nonzero elements

  NumGlobalRows = A.NumGlobalRows();
  NumGlobalNonzeros = A.NumGlobalNonzeros();

  // print out on cout if no filename is provided

  int IndexBase = A.IndexBase(); // MATLAB starts from 0
  if( IndexBase == 0 )
    IndexBase = 1;
  else if ( IndexBase == 1)
    IndexBase = 0;

  // write on file the dimension of the matrix

  if( MyPID==0 ) {
    outfile << "A = spalloc(";
    outfile << NumGlobalRows << ',' << NumGlobalRows;
    outfile << ',' << NumGlobalNonzeros << ");\n";
  }

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {
    A.Comm().Barrier();
    if( MyPID == Proc ) {

      outfile << "\n\n% On proc " << Proc << ": ";
      outfile << NumMyRows << " rows and ";
      outfile << A.NumMyNonzeros() << " nonzeros\n";

      // cycle over all local rows to find out nonzero elements
      for( int MyRow=0 ; MyRow<NumMyRows ; ++MyRow ) {

        GlobalRow = A.GRID(MyRow);

        NumNzRow = A.NumMyEntries(MyRow);
        double *Values = new double[NumNzRow];
        int *Indices = new int[NumNzRow];

        A.ExtractMyRowCopy(MyRow, NumNzRow,
                           NumEntries, Values, Indices);
        // print out the elements with MATLAB syntax
        for( int j=0 ; j<NumEntries ; ++j ) {
          outfile << "A(" << GlobalRow  + IndexBase
                  << "," << A.GCID(Indices[j]) + IndexBase
                  << ") = " << Values[j] << ";\n";
        }

        delete Values;
        delete Indices;
      }

    }
    A.Comm().Barrier();
  }

  return true;

}



/* ======== =============== *
 * function FEVector2MATLAB *
 * ======== =============== *
 *
 * Print out a Epetra_Vector in a MATLAB format. Each processor prints out
 * its part, starting from proc 0 to proc NumProc-1. The first line of
 * each processor's output states the number of local rows and of
 * local nonzero elements.
 *
 * Return code:        true if vector has been printed out
 * -----------         false otherwise
 *
 * Parameters:
 * ----------
 *
 * - Epetra_FEVector   reference to FE vector
 * - ostream &         reference to output stream
 */
bool FEVector2MATLAB( const Epetra_FEVector & v, ostream & outfile) {

  int MyPID = v.Comm().MyPID();
  int NumProc = v.Comm().NumProc();
  int MyLength = v.MyLength();
  int GlobalLength = v.GlobalLength();

  // print out on cout if no filename is provided

  // write on file the dimension of the matrix

  if( MyPID == 0 ) outfile << "v = zeros(" << GlobalLength << ",1)\n";

  int NumMyElements = v.Map().NumMyElements();
  // get update list
  int * MyGlobalElements = v.Map().MyGlobalElements( );

  int Row;

  int IndexBase = v.Map().IndexBase(); // MATLAB starts from 0
  if( IndexBase == 0 )
    IndexBase = 1;
  else if ( IndexBase == 1)
    IndexBase = 0;

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {
    v.Comm().Barrier();
    if( MyPID == Proc ) {

      outfile << "% On proc " << Proc << ": ";
      outfile << MyLength << " rows of ";
      outfile << GlobalLength << " elements\n";

      for( Row=0 ; Row<MyLength ; ++Row ) {
        outfile << "v(" << MyGlobalElements[Row] + IndexBase
             << ") = " << v[0][Row] << ";\n";
      }

    }

    v.Comm().Barrier();
  }

  return true;

} /* FEVector2MATLAB */
