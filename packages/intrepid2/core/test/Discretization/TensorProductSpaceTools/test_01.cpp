/*
// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Intrepid2_TensorProductSpaceTools.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"
#include "Intrepid2_CubaturePolylib.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Types.hpp"


using Teuchos::Array;
using Intrepid2::FieldContainer;
using Intrepid2::Basis;
using Intrepid2::TensorBasis;

#define INTREPID_TEST_COMMAND( S )                                                                                  \
{                                                                                                                   \
  try {                                                                                                             \
    S ;                                                                                                             \
  }                                                                                                                 \
  catch (std::logic_error err) {                                                                                    \
      *outStream << "Expected Error ----------------------------------------------------------------\n";            \
      *outStream << err.what() << '\n';                                                                             \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";    \
  };                                                                                                                \
}


int main( int argc , char **argv )
{  
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
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
    << "|                      Unit Test TensorProductSpace Tools                     |\n" \
    << "|                                                                             |\n" \
    << "|     Tests sum-factored polynomial evaluation and integration                |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
    << "|                      Robert Kirby (robert.c.kirby@ttu.edu)                  |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";



  int errorFlag = 0;

  Array<RCP<TensorBasis<double,FieldContainer<double> > > > basesByDim(4);
  Array<RCP<FieldContainer<double> > > cubPtsByDim(4);

  Intrepid2::CubaturePolylib<double> cpl(2,Intrepid2::PL_GAUSS_LOBATTO);

  FieldContainer<double> cubPoints( cpl.getNumPoints() ,1 );
  FieldContainer<double> cubWeights( cpl.getNumPoints() );

  cpl.getCubature( cubPoints, cubWeights );

  basesByDim[2] = Teuchos::rcp( new Intrepid2::Basis_HGRAD_QUAD_Cn_FEM<double,FieldContainer<double> >( 2 , Intrepid2::POINTTYPE_SPECTRAL ) );
  basesByDim[3] = Teuchos::rcp( new Intrepid2::Basis_HGRAD_HEX_Cn_FEM<double,FieldContainer<double> >( 2 , Intrepid2::POINTTYPE_SPECTRAL ) );


  // get points
  FieldContainer<double> quadPts( cpl.getNumPoints() * cpl.getNumPoints() , 2 );
  for (int j=0;j<cpl.getNumPoints();j++)
    {
      for (int i=0;i<cpl.getNumPoints();i++)
	{
	  int index = j*cpl.getNumPoints() + i;
	  quadPts(index,0) = cubPoints(i,0);
	  quadPts(index,1) = cubPoints(j,0);
	}
    }

  FieldContainer<double> cubPts( cpl.getNumPoints() * cpl.getNumPoints() * cpl.getNumPoints() , 3 );
  for (int k=0;k<cpl.getNumPoints();k++)
    {
      for (int j=0;j<cpl.getNumPoints();j++)
	{
	  for (int i=0;i<cpl.getNumPoints();i++)
	    {
	      int index = k* cpl.getNumPoints() * cpl.getNumPoints() + j*cpl.getNumPoints() + i;
	      cubPts(index,0) = cubPoints(i,0);
	      cubPts(index,1) = cubPoints(j,0);
	      cubPts(index,2) = cubPoints(k,0);
	    }
	}
    }

  cubPtsByDim[2] = Teuchos::rcp( &quadPts , false );
  cubPtsByDim[3] = Teuchos::rcp( &cubPts , false );

  int space_dim = 2;

  Array<Array<RCP<Basis<double,FieldContainer<double> > > > > &bases = basesByDim[space_dim]->getBases();

  FieldContainer<double> coeff(1,1,basesByDim[space_dim]->getCardinality());


  
  Array<RCP<FieldContainer<double> > > pts( space_dim );
  pts[0] = Teuchos::rcp( &cubPoints, false );
  for (int i=1;i<space_dim;i++)
    {
      pts[i] = pts[0];
    }

  Array<RCP<FieldContainer<double> > > wts(space_dim);
  wts[0] = Teuchos::rcp( &cubWeights , false );
  for (int i=1;i<space_dim;i++)
    {
      wts[i] = wts[0];
    }

  FieldContainer<double> Phix(bases[0][0]->getCardinality(),
			      cpl.getNumPoints() );
  FieldContainer<double> Phiy(bases[0][1]->getCardinality(),
			      cpl.getNumPoints() );
  FieldContainer<double> DPhix(bases[0][0]->getCardinality(),
			       cpl.getNumPoints(), 1 );
  FieldContainer<double> DPhiy(bases[0][1]->getCardinality(),
			       cpl.getNumPoints(), 1 );

  bases[0][0]->getValues( Phix , cubPoints, Intrepid2::OPERATOR_VALUE );
  bases[0][1]->getValues( Phiy , cubPoints, Intrepid2::OPERATOR_VALUE );
  bases[0][0]->getValues( DPhix , cubPoints, Intrepid2::OPERATOR_D1 );
  bases[0][1]->getValues( DPhiy , cubPoints, Intrepid2::OPERATOR_D1 );

  Array<RCP<FieldContainer<double> > > basisVals(2);
  basisVals[0] = Teuchos::rcp( &Phix , false );
  basisVals[1] = Teuchos::rcp( &Phiy , false );

  Array<RCP<FieldContainer<double> > > basisDVals(2);
  basisDVals[0] = Teuchos::rcp( &DPhix , false );
  basisDVals[1] = Teuchos::rcp( &DPhiy , false );

  FieldContainer<double> vals(1,1,pts[0]->size() * pts[1]->size() );

  // first basis function is the polynomial.
  coeff(0,0,0) = 1.0;

  Intrepid2::TensorProductSpaceTools::evaluate<double,FieldContainer<double>,FieldContainer<double>,FieldContainer<double> >( vals , coeff , basisVals );

  FieldContainer<double> grads( 1 , 1, pts[0]->size() * pts[1]->size() , 2 );

  Intrepid2::TensorProductSpaceTools::evaluateGradient<double,FieldContainer<double>,FieldContainer<double>,FieldContainer<double> >( grads , 
																     coeff , 
																     basisVals ,
																     basisDVals );

  // confirm by comparing to actual gradients
  FieldContainer<double> fullVals(basesByDim[space_dim]->getCardinality(),
				  basesByDim[space_dim]->getCardinality());
  FieldContainer<double> fullGrads(basesByDim[space_dim]->getCardinality(),
				   basesByDim[space_dim]->getCardinality(),
				   space_dim );


  basesByDim[space_dim]->getValues( fullVals ,
		    quadPts ,
		    Intrepid2::OPERATOR_VALUE );
  basesByDim[space_dim]->getValues( fullGrads ,
		    quadPts ,
		    Intrepid2::OPERATOR_GRAD );

  for (int i=0;i<fullVals.dimension(1);i++)
    {
      if (std::abs( fullVals(0,i) - vals(0,0,i) ) > Intrepid2::INTREPID_TOL ) 
	{
	  errorFlag++;
	  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          
	  // Output the multi-index of the value where the error is:
	  *outStream << " Evaluating first bf At multi-index { ";
	  *outStream << i;
	  *outStream << "}  brute force value: " << fullVals(0,i)
		     << " but tensor-product  value: " << vals(0,i) << "\n";
	  *outStream << "Difference: " << std::abs( fullVals(0,i) - vals(0,i) ) << "\n";
          }
    }

  for (int i=0;i<fullGrads.dimension(1);i++)
    {
      for (int j=0;j<fullGrads.dimension(2);j++)
	{
	  if (std::abs( fullGrads(0,i,j) - grads(0,0,i,j) ) > Intrepid2::INTREPID_TOL ) 
	    {
	      errorFlag++;
	      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
	      
	      // Output the multi-index of the value where the error is:
	      *outStream << " Evaluating first bf At multi-index { ";
	      *outStream << i << " " << j;
	      *outStream << "}  brute force value: " << fullGrads(0,i,j)
			 << " but tensor-product  value: " << grads(0,0,i,j) << "\n";
	      *outStream << "Difference: " << std::abs( fullGrads(0,i,j) - grads(0,0,i,j) ) << "\n";
	    }
	}
    }

  
  // now test moments.  
  // I've already evaluated the first basis function at the quadrature points.
  // why not use it?

  FieldContainer<double> momentsNaive(1,basesByDim[2]->getCardinality());
  for (int i=0;i<basesByDim[2]->getCardinality();i++)
    {
      momentsNaive(0,i) = 0.0;
      for (int qpty=0;qpty<cubPoints.dimension(0);qpty++)
	{
	  for (int qptx=0;qptx<cubPoints.dimension(0);qptx++)
	    {
	      momentsNaive(0,i) += cubWeights(qpty) * cubWeights(qptx) *
		vals( 0, 0, qpty*cubPoints.dimension(0)+qptx )
		* fullVals(i,qpty*cubPoints.dimension(0)+qptx);
	    }
	}
    }

  FieldContainer<double> momentsClever(1,1,basesByDim[space_dim]->getCardinality());
  Intrepid2::TensorProductSpaceTools::moments<double,FieldContainer<double>,FieldContainer<double>,FieldContainer<double>,FieldContainer<double> >( momentsClever , 
																		   vals ,
																		   basisVals ,
																		   wts );
  for (int j=0;j<momentsClever.dimension(0);j++)
    {
      if (std::abs( momentsClever(0,0,j) - momentsNaive(0,j) ) > Intrepid2::INTREPID_TOL ) 
	{
	  errorFlag++;
	  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
	  
	  // Output the multi-index of the value where the error is:
	  *outStream << " At multi-index { ";
	  *outStream << " " << j;
	  *outStream << "}  brute force value: " << momentsNaive(0,j)
		     << " but sum-factored value: " << momentsClever(0,0,j) << "\n";
	  *outStream << "Difference: " << std::abs( momentsNaive(0,j) - momentsClever(0,0,j) ) << "\n";
	}
    }

  if (errorFlag != 0)
    {
      std::cout << "End Result: TEST FAILED\n";
    }
  else
    {
      std::cout << "End Result: TEST PASSED\n";
    }
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  Kokkos::finalize();
  return errorFlag;

}
