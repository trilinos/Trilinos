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
        with Basis_F0_TRI_C2_FEM_DEFAULT basis.
\author Created by P. Bochev and D. Ridzal.
*/

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
  << "|                 Unit Test (Basis_F0_TRI_C2_FEM_DEFAULT)                     |\n" \
  << "|                                                                             |\n" \
  << "|     1) Basis creation, computation of basis function values                 |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 1: Correctness of mass and stiffness matrices                          |\n"\
  << "===============================================================================\n";

  int errorFlag = 0;

  try {
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

    MultiCell<double> mCell(nCells,          // number of cells (triangles) in the multicell instance
                            CELL_TRI,        // generating cell type
                            triNodes);       // array with interleaved node coordinates

    FieldContainer<double> massMatrices;
    DefaultFieldFactory<double> FFactory;

    for (ECompEngine compEng = COMP_CPP; compEng < COMP_ENGINE_MAX; compEng++) {
      for (int cubDeg=4; cubDeg<=20; cubDeg++) {
        Teuchos::RCP<LocalField<double> > form0 = FFactory.create(FIELD_FORM_0,
                                                                  CELL_TRI,
                                                                  RECONSTRUCTION_SPACE_COMPLETE,
                                                                  2,
                                                                  BASIS_FEM_DEFAULT,
                                                                  COORDINATES_CARTESIAN,
                                                                  cubDeg);
        // compute mass matrices
        form0->getOperator(massMatrices, OPERATOR_VALUE, OPERATOR_VALUE, mCell);

        *outStream << "\nComputational engine: " << ECompEngineToString(compEng) << "\n";
        *outStream << "Cubature degree:      " << cubDeg << "\n";

        for (int cell_id = 0; cell_id < nCells; cell_id++) {
          stringstream namestream;
          string filename;
          namestream <<  basedir << "/mass_TRI_FEM_P2" << "_" << "0" << cell_id+1 << ".dat";
          namestream >> filename;

          ifstream massfile(&filename[0]);
          if (massfile.is_open()) {
            Teuchos::Array<Teuchos::Array<double> > cellMass;
            // fill mass matrix for this cell
            int numLbf = massMatrices.getDimension(1);
            int numRbf = massMatrices.getDimension(2);
            cellMass.resize(numLbf);
            for (int i=0; i<numLbf; i++) {
              cellMass[i].resize(numRbf);
              for (int j=0; j<numRbf; j++) {
                cellMass[i][j] = massMatrices(cell_id, i, j);
              }
            }
            if (compareToAnalytic<double>(cellMass, massfile, 1e-10, iprint) > 0) {
              errorFlag++;
            }
            massfile.close();
          }
          else {
            errorFlag = -999;
          }
        }
      }
    }
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
