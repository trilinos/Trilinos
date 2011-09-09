// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   example_05.cpp
    \brief  Demonstrate diagonalized mass matrices for H(grad) elements in 1d
            using Gauss-Legendre quadrature

    \author Created by P. Bochev, R. Kirby D. Ridzal and K. Peterson.

    
     \remark Usage
     \verbatim

     ./Intrepid_example_Drivers_Example_03.exe max_deg verbose

        int min_deg         - beginning polynomial degree to check 
	int max_deg         - maximum polynomial degree to check
        verbose (optional)  - any character, indicates verbose output

     \endverbatim

    \remark Sample command line: checks mass matrix of degree 1,2,3
    \code   ./Intrepid_example_Drivers_Example_03.exe 1 3 \endcode
*/

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HGRAD_LINE_Cn_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid_HGRAD_HEX_Cn_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_CubaturePolylib.hpp"
#include "Intrepid_CubatureTensor.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialComm.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

  //Check number of arguments
    if (argc < 2) {
       std::cout <<"\n>>> ERROR: Invalid number of arguments.\n\n";
       std::cout <<"Usage:\n\n";
       std::cout <<"  ./Intrepid_example_Drivers_Example_05.exe min_deg max_deg verbose\n\n";
       std::cout <<" where \n";
       std::cout <<"   int min_deg         - beginning polynomial degree to check \n";
       std::cout <<"   int max_deg         - final polynomial degree to check \n";
       std::cout <<"   verbose (optional)  - any character, indicates verbose output \n\n";
       exit(1);
    }
    
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 1)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|  Example: Check diagonalization of reference mass matrix                    |\n" \
    << "|           on line, quad, and hex                                            |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";


// ************************************ GET INPUTS **************************************
  int min_degree = atoi(argv[1]);
  int max_degree = atoi(argv[2]);
  

// *********************************** CELL TOPOLOGY **********************************
  
  // Get cell topology for base line
  typedef shards::CellTopology    CellTopology;
  CellTopology line_2(shards::getCellTopologyData<shards::Line<2> >() );
  CellTopology quad_4(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
  CellTopology hex_8(shards::getCellTopologyData<shards::Hexahedron<8> >() );
  
  std::vector<CellTopology> cts(3);

  // loop over degrees
  for (int deg=min_degree;deg<=max_degree;deg++) {
    std::vector<Teuchos::RCP<Basis<double,FieldContainer<double> > > > bases;
    bases.push_back( Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<double, FieldContainer<double> >(deg,POINTTYPE_SPECTRAL)) );
    bases.push_back( Teuchos::rcp(new Basis_HGRAD_QUAD_Cn_FEM<double, FieldContainer<double> >(deg,POINTTYPE_SPECTRAL)) );
    bases.push_back( Teuchos::rcp(new Basis_HGRAD_HEX_Cn_FEM<double, FieldContainer<double> >(deg,POINTTYPE_SPECTRAL)) );

    // ************************************ CUBATURE **************************************
    Teuchos::RCP<Cubature<double,FieldContainer<double>,FieldContainer<double> > > glcub
      = Teuchos::rcp(new CubaturePolylib<double,FieldContainer<double>,FieldContainer<double> >(2*deg-1,PL_GAUSS_LOBATTO) );
    
    std::vector<Teuchos::RCP<Cubature<double,FieldContainer<double>,FieldContainer<double> > > >
      cub_to_tensor;

    // now we loop over spatial dimensions
    for (int sd=1;sd<=3;sd++) {
      // ************************************ CUBATURE **************************************
      cub_to_tensor.push_back( glcub );
      CubatureTensor<double,FieldContainer<double>,FieldContainer<double> > cubcur( cub_to_tensor );

      int cubDim       = cubcur.getDimension();
      int numCubPoints = cubcur.getNumPoints();
    
      FieldContainer<double> cubPoints(numCubPoints, cubDim);
      FieldContainer<double> cubWeights(numCubPoints);
    
      cubcur.getCubature(cubPoints, cubWeights);

      // ************************************ BASIS AT QP**************************************
      Teuchos::RCP<Basis<double,FieldContainer<double> > > basis_cur = bases[sd-1];
      const int numFields = basis_cur->getCardinality();
      FieldContainer<double> bf_at_cub_pts( numFields, numCubPoints );
      FieldContainer<double> trans_bf_at_cub_pts( 1 , numFields,numCubPoints );
      FieldContainer<double> w_trans_bf_at_cub_pts( 1, numFields , numCubPoints );
      basis_cur->getValues( bf_at_cub_pts , cubPoints , OPERATOR_VALUE );
      
      // *********************************** FORM MASS MATRIX *****************************
       FunctionSpaceTools::HGRADtransformVALUE<double>( trans_bf_at_cub_pts ,
							bf_at_cub_pts );
       cubWeights.resize(1,numCubPoints);
       FunctionSpaceTools::multiplyMeasure<double>( w_trans_bf_at_cub_pts ,
 						   cubWeights ,
 						   trans_bf_at_cub_pts );
       cubWeights.resize(numCubPoints);

      FieldContainer<double> mass_matrix( 1 , numFields, numFields );
      FunctionSpaceTools::integrate<double>( mass_matrix ,
					     trans_bf_at_cub_pts ,
					     w_trans_bf_at_cub_pts ,
					     COMP_BLAS );

      // now we loop over the mass matrix and compare the nondiagonal entries to zero
      double max_offdiag = 0.0;
      for (int i=0;i<numFields;i++) {
	for (int j=0;j<numFields;j++) {
	  if (i != j) {
	    if ( abs(mass_matrix(0,i,j)) >= max_offdiag) {
	      max_offdiag = abs(mass_matrix(0,i,j));
	    }
	  }
	}
      }
      double min_diag = mass_matrix(0,0,0);
      for (int i=0;i<numFields;i++) {
	if ( mass_matrix(0,i,i) <= min_diag ) {
	  min_diag = mass_matrix(0,i,i);
	}
      }
      *outStream << "Degree = " << deg << " and dimension = " << sd << "; Max offdiagonal"
		 << " element is " << max_offdiag << " and min diagonal element is " << min_diag
		 << ".\n";
	  
    }

  }    
    
   std::cout << "End Result: TEST PASSED\n";

  std::cout.copyfmt(oldFormatState);

return 0;
}

