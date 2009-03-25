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
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file test_02.cpp
\brief  Unit tests for the Intrepid::LocalForm0  class,
        with Basis_F0_TRI_C1_FEM_DEFAULT basis.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_LocalForm0.hpp"
#include "Intrepid_DefaultBasisFactory.hpp"
#include "Intrepid_CubatureDirect.hpp"
#include "Intrepid_MultiCell.hpp"
#include "Intrepid_Utils.hpp"
#include "Intrepid_DefaultFieldFactory.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  *outStream \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                 Unit Test (Basis_F0_TRI_C1_FEM_DEFAULT)                     |\n" \
  << "|                                                                             |\n" \
  << "|     1) Correctness of Operators                                             |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 1: Correctness of mass and stiffness matrices: no Jacobian reuse       |\n"\
  << "===============================================================================\n";

  int errorFlag = 0;

  try {
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > basis = \
      BFactory.create(FIELD_FORM_0, 
                      CELL_TRI, 
                      RECONSTRUCTION_SPACE_COMPLETE, 
                      1, 
                      BASIS_FEM_DEFAULT, 
                      COORDINATES_CARTESIAN);

    Teuchos::Array<Teuchos::Array< Teuchos::RCP<Cubature<double> > > > allCubs;
    Teuchos::RCP<Cubature<double> > cellCub = Teuchos::rcp(new CubatureDirect<double>(CELL_TRI,2) );
    allCubs.resize(1); 
    allCubs[0].resize(1);
    allCubs[0][0] = cellCub;
    /*************** use if you need face cubatures
    Teuchos::RCP<Cubature<double> > faceCub = Teuchos::rcp(new CubatureDirect<double>(CELL_EDGE,2) );
    allCubs.resize(2); 
    allCubs[0].resize(1);
    allCubs[0][0] = cellCub;
    allCubs[1].resize(3);
    allCubs[1][0] = faceCub; 
    allCubs[1][1] = faceCub; 
    allCubs[1][2] = faceCub;
    ****************/

    string basedir = "./data";

    int nCells = 4;

    double triNodes[] = {
      // tri 0
      0.0, 0.0,
      1.0, 0.0,
      0.0, 1.0,
      // tri 1
      4.0, 5.0,
      -6.0, 2.0,
      4.0, -3.0,
      // tri 2
      -6.0, -3.0,
      9.0, 2.0,
      8.9, 2.1,
      // tri 3
      -6.0, -3.0,
      12.0, 3.0,
      2.9, 0.1
    };

    // This MCell and container are to test matrix assembly without Jacobian reuse
    MultiCell<double> mCell(nCells,          // number of cells (triangles) in the multicell instance
                            CELL_TRI,        // generating cell type
                            triNodes);       // array with interleaved node coordinates
    FieldContainer<double> massMatrices;
    FieldContainer<double> stiffMatrices;

    // This MCell and container are for to test matrix assembly with Jacobian reuse.
    MultiCell<double> mCellReuse(nCells,          
                            CELL_TRI,        
                            triNodes);       
    FieldContainer<double> massMatricesReuse;
    FieldContainer<double> stiffMatricesReuse;
    
    // This MCell and container are for to test matrix assembly with auxiliary "right" LocalForm0
    MultiCell<double> mCellWithRight(nCells,          
                                     CELL_TRI,        
                                     triNodes);       
    FieldContainer<double> massMatricesWithRight;
    FieldContainer<double> stiffMatricesWithRight;

    for (ECompEngine compEng = COMP_CPP; compEng < COMP_ENGINE_MAX; compEng++) {
      for (int cubDeg=2; cubDeg<=20; cubDeg++) {
        
        // Delete precomputed measures before changing cubature rule!!
        mCellReuse.deleteSubcellMeasures(2,0);
        mCellWithRight.deleteSubcellMeasures(2,0);

        // set cubature
        cellCub = Teuchos::rcp(new CubatureDirect<double>(CELL_TRI,cubDeg) );
        allCubs[0][0] = cellCub;

        // create local form
        LocalForm0<double> form0(basis, allCubs, compEng);
        
        // Create "right" local form for testing purposes using the same basis and cubature
        LocalForm0<double> rightForm0(basis, allCubs, compEng);
        
        // compute mass and stiffness matrices without Jacobian reuse
        form0.getOperator(massMatrices, OPERATOR_VALUE, OPERATOR_VALUE, mCell);
        form0.getOperator(stiffMatrices, OPERATOR_GRAD, OPERATOR_GRAD, mCell);
        
        // compute mass and stiffness matrices with Jacobian reuse
        form0.getOperator(massMatricesReuse, OPERATOR_VALUE, OPERATOR_VALUE, mCellReuse, true);
        form0.getOperator(stiffMatricesReuse, OPERATOR_GRAD, OPERATOR_GRAD, mCellReuse, true);
        
        // compute mass and stiffness  matrices using an "auxiliary" right local field:
        form0.getOperator(massMatricesWithRight, OPERATOR_VALUE, OPERATOR_VALUE, rightForm0, mCellWithRight, true);
        form0.getOperator(stiffMatricesWithRight, OPERATOR_GRAD, OPERATOR_GRAD, rightForm0, mCellWithRight, true);
         
        *outStream << "\nComputational engine: " << ECompEngineToString(compEng) << "\n";
        *outStream << "Cubature degree:      " << cubDeg << "\n";

        // compare mass matrices to analytic
        for (int cellId = 0; cellId < nCells; cellId++) {
          *outStream << "\n Cell Id = " << cellId << " -----------\n\n";
          stringstream namestream;
          string filename;
          namestream <<  basedir << "/mass_TRI_FEM_P1" << "_" << "0" << cellId+1 << ".dat";
          namestream >> filename;

          // Temp arrays to load mass matrix entries from FieldContainer file
          Teuchos::Array<Teuchos::Array<double> > cellMass;
          Teuchos::Array<Teuchos::Array<double> > cellMassReuse;
          Teuchos::Array<Teuchos::Array<double> > cellMassWithRight;
          
          // fill mass matrix for this cell
          int numLbf = massMatrices.dimension(1);
          int numRbf = massMatrices.dimension(2);
          cellMass.resize(numLbf);
          cellMassReuse.resize(numLbf);
          cellMassWithRight.resize(numLbf);
          
          for (int i=0; i<numLbf; i++) {
            cellMass[i].resize(numRbf);
            cellMassReuse[i].resize(numRbf);
            cellMassWithRight[i].resize(numRbf);
            for (int j=0; j<numRbf; j++) {
              cellMass[i][j]          = massMatrices(cellId, i, j);
              cellMassReuse[i][j]     = massMatricesReuse(cellId, i, j);
              cellMassWithRight[i][j] = massMatricesWithRight(cellId, i, j);
            }
          }
          
          // Compare mass matrix without reuse and entries from data file
          ifstream massfile(&filename[0]);
          if (massfile.is_open()) {
            *outStream << " Mass matrix without Jacobian reuse:  \n";
            if (compareToAnalytic<double>(cellMass, massfile, 1e-10, iprint) > 0) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << std::setw(70) << " Wrong mass matrix when not reusing Jacobians" << "\n";
            }
            massfile.close();
            massfile.clear();
          }
          else {
            errorFlag = -999;
          }
          
          // Compare mass matrix with reuse and entries from data file. Need to open file again
          massfile.open(&filename[0]);
          if (massfile.is_open()) {              
            *outStream << " Mass matrix with Jacobian reuse: \n";
            if (compareToAnalytic<double>(cellMassReuse, massfile, 1e-10, iprint) > 0) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << std::setw(70) << " Wrong mass matrix when reusing Jacobians" << "\n";
            }
            massfile.close();
            massfile.clear();
          }
          else {
            errorFlag = -999;
          }
          
          // Compare mass matrix with auxiliary right field and entries from data file. Need to open file again
          massfile.open(&filename[0]);
          if (massfile.is_open()) {              
            *outStream << " Mass matrix with auxiliary right op: \n";
            if (compareToAnalytic<double>(cellMassWithRight, massfile, 1e-10, iprint) > 0) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << std::setw(70) << " Wrong mass matrix when using auxiliary LocalForm0 right field" << "\n";
            }
            massfile.close();
            massfile.clear();
          }
          else {
            errorFlag = -999;
          }
          
        } // for(cellId) mass matrices

//#ifdef QWERTYUIO
        // compare stiffness matrices to analytic
        for (int cellId = 0; cellId < nCells; cellId++) {
          *outStream << "\n Cell Id = " << cellId << " -----------\n\n";
          stringstream namestream;
          string filename;
          namestream <<  basedir << "/stiff_TRI_FEM_P1" << "_" << "0" << cellId+1 << ".dat";
          namestream >> filename;

          // Temp arrays to load stiffness matrix entries from FieldContainer file
          Teuchos::Array<Teuchos::Array<double> > cellStiff;
          Teuchos::Array<Teuchos::Array<double> > cellStiffReuse;
          Teuchos::Array<Teuchos::Array<double> > cellStiffWithRight;
          
          // fill stiffness matrix for this cell
          int numLbf = stiffMatrices.dimension(1);
          int numRbf = stiffMatrices.dimension(2);
          cellStiff.resize(numLbf);
          cellStiffReuse.resize(numLbf);
          cellStiffWithRight.resize(numLbf);
          
          for (int i=0; i<numLbf; i++) {
            cellStiff[i].resize(numRbf);
            cellStiffReuse[i].resize(numRbf);
            cellStiffWithRight[i].resize(numRbf);
            for (int j=0; j<numRbf; j++) {
              Teuchos::Array<int> mIndex(3);
              cellStiff[i][j]          = stiffMatrices(cellId, i, j);
              cellStiffReuse[i][j]     = stiffMatricesReuse(cellId, i, j);
              cellStiffWithRight[i][j] = stiffMatricesWithRight(cellId, i, j);
            }
          }
          
          // Compare stiffness matrix without reuse and entries from data file
          ifstream stifffile(&filename[0]);
          if (stifffile.is_open()) {
            *outStream << " Stiffness matrix without Jacobian reuse:  \n";
            if (compareToAnalytic<double>(cellStiff, stifffile, 1e-10, iprint) > 0) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << std::setw(70) << " Wrong stiffness matrix when not reusing Jacobians" << "\n";
            }
            stifffile.close();
            stifffile.clear();
          }
          else {
            errorFlag = -999;
          }
          
          // Compare stiffness matrix with reuse and entries from data file. Need to open file again
          stifffile.open(&filename[0]);
          if (stifffile.is_open()) {              
            *outStream << " Stiffness matrix with Jacobian reuse: \n";
            if (compareToAnalytic<double>(cellStiffReuse, stifffile, 1e-10, iprint) > 0) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << std::setw(70) << " Wrong stiffness matrix when reusing Jacobians" << "\n";
            }
            stifffile.close();
            stifffile.clear();
          }
          else {
            errorFlag = -999;
          }
          
          
          // Compare stiffness matrix with auxiliary right field and entries from data file. Need to open file again
          stifffile.open(&filename[0]);
          if (stifffile.is_open()) {              
            *outStream << " Stiffness matrix with auxiliary right op: \n";
            if (compareToAnalytic<double>(cellStiffWithRight, stifffile, 1e-10, iprint) > 0) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << std::setw(70) << " Wrong stiffness matrix when using auxiliary LocalForm0 right field" << "\n";
            }
            stifffile.close();
            stifffile.clear();
          }
          else {
            errorFlag = -999;
          }
          
        } // for(cellId) stiffness matrices
//#endif

      } // for(cubDeg)
    }// for(compEngine)
  }
  catch (std::logic_error err) {
      *outStream << err.what() << "\n\n";
          errorFlag = -999;
  };

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
