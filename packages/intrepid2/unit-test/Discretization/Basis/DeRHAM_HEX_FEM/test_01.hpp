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


/** \file
    \brief  Test for checking the De Rham complex for FEM on HEX
            Testing that HGRAD_Cn --(grad)--> HCURL_In --(curl)--> HDIV_In --(div)--> HVOL_C(n-1)
            And that curl grad f = 0, and div curl f =0.
    \author Created by M. Perego
 */


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"
#include "Intrepid2_HDIV_HEX_In_FEM.hpp"
#include "Intrepid2_HVOL_HEX_Cn_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"


namespace Intrepid2 {

namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }

//Warning: This is not supposed to be build with CUDA.


template<typename ValueType, typename DeviceSpaceType>
int DeRHAM_HEX_FEM_Test01(const bool verbose) {

  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  if (verbose)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs,       false);

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  *outStream
  <<"\n"
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|                 Unit Test checking properties of DeRHAM complex             |\n"
  << "|                                                                             |\n"
  << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
  << "|                      Robert Kirby  (robert.c.kirby@ttu.edu),                |\n"
  << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
  << "|                      Kara Peterson (kjpeter@sandia.gov),                    |\n"
  << "|                      Kyungjoo Kim  (kyukim@sandia.gov),                     |\n"
  << "|                      Mauro Perego  (mperego@sandia.gov).                    |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";

  typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

  const ValueType tol = tolerence();
  int errorFlag = 0;

  // for virtual function, value and point types are declared in the class
  typedef ValueType outputValueType;
  typedef ValueType pointValueType;

  typedef Basis_HGRAD_HEX_Cn_FEM<DeviceSpaceType,outputValueType,pointValueType> HexHGRADBasisType;
  typedef Basis_HCURL_HEX_In_FEM<DeviceSpaceType,outputValueType,pointValueType> HexHCURLBasisType;
  typedef Basis_HDIV_HEX_In_FEM<DeviceSpaceType,outputValueType,pointValueType> HexHDIVBasisType;
  typedef Basis_HVOL_HEX_Cn_FEM<DeviceSpaceType,outputValueType,pointValueType> HexHVOLBasisType;

  constexpr ordinal_type dim = 3;
  constexpr ordinal_type maxOrder = Intrepid2::Parameters::MaxOrder < 5 ? Intrepid2::Parameters::MaxOrder : 5;

  // In order to obtain nontrivial FE functions, we interpolate the following (vector) function on the space spanned by the fE basis considered.
  struct Fun {
    pointValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const pointValueType& x, const pointValueType& y, const pointValueType& z, const int comp=0) {
      switch (comp) {
      case 0:
        return std::sin(x)*std::sin(y)*std::sin(z)+2.0*std::cos(x)*std::cos(y)*std::cos(z);
        break;
      case 1:
        return std::exp(x)*std::exp(y)*std::exp(z);
        break;
      case 2:
        return std::log(x+2)*std::log(y+2)*std::log(z+2);
        break;
      default:
        return 0;
      }
    }
  };



  try {

    *outStream
    << "\n"
    << "=======================================================================================================\n"
    << "| TEST 1: Testing that grad of HGRAD_Cn function is in the kernel the curl operator in HCURL_In space |\n"
    << "=======================================================================================================\n";

    /*
     * For every order n<=maxOrder, we check that the gradient of a function in HGRAD_Cn space is in the HCURL_In space.
     * Then we take its curl, in the HCURL_In space, and verify is zero.
     */

    for (ordinal_type order=1; order <= maxOrder; ++order ) {

      //Getting ref coordinates (HGRAD dof coordinates):
      //Point evaluations at ref coordinates uniquely identify a polynomial of degree n.
      HexHGRADBasisType hexHGradBasisSpec(order, POINTTYPE_WARPBLEND);
      const ordinal_type refCardinality = hexHGradBasisSpec.getCardinality();
      DynRankView ConstructWithLabel(refCoords, refCardinality, dim);
      hexHGradBasisSpec.getDofCoords(refCoords);

      //Creating HGRAD and HCURL basis
      HexHGRADBasisType hexHGradBasis(order, POINTTYPE_EQUISPACED);  //any point type is fine
      HexHCURLBasisType hexHCurlBasis(order, POINTTYPE_EQUISPACED); //any point type is fine

      const ordinal_type hgradCardinality = hexHGradBasis.getCardinality();
      const ordinal_type hcurlCardinality = hexHCurlBasis.getCardinality();

      //Getting DOF coordinates for HGRAD and HCURL elements
      DynRankView ConstructWithLabel(hgradDofCoords, hgradCardinality , dim);
      hexHGradBasis.getDofCoords(hgradDofCoords);

      DynRankView ConstructWithLabel(hcurlDofCoords, hcurlCardinality , dim);
      hexHCurlBasis.getDofCoords(hcurlDofCoords);

      //Getting DOF coefficients for HGRAD and HCURL elements
      DynRankView ConstructWithLabel(hgradDofCoeffs, hgradCardinality);
      hexHGradBasis.getDofCoeffs(hgradDofCoeffs);

      DynRankView ConstructWithLabel(hcurlDofCoeffs, hcurlCardinality, dim);
      hexHCurlBasis.getDofCoeffs(hcurlDofCoeffs);

      //Evaluating the function at HGRAD dof coordinates
      DynRankView ConstructWithLabel(funAtHGradDofCoords, hgradCardinality);

      Fun fun;
      for(int i=0;i<hgradCardinality;i++)
        funAtHGradDofCoords(i) = fun(hgradDofCoords(i,0), hgradDofCoords(i,1), hgradDofCoords(i,2));

      //Interpolating the function in the HGRAD space by computing the degrees of freedom for the function
      DynRankView ConstructWithLabel(funHGradCoeffs, hgradCardinality);
      for(int i=0;i<hgradCardinality;i++)
        funHGradCoeffs(i) = funAtHGradDofCoords(i)*hgradDofCoeffs(i); //not really needed for HGRAD, hgradDofCoeffs = 1

      //Computing the gradient of the (hgrad) interpolated function (ifun) at HGRAD dof coordinates
      DynRankView ConstructWithLabel(gradOfHGradBasisAtRefCoords, hgradCardinality , refCardinality, dim);
      hexHGradBasis.getValues(gradOfHGradBasisAtRefCoords, refCoords, OPERATOR_GRAD);
      DynRankView ConstructWithLabel(ifunGradAtRefCoords, refCardinality, dim);
      for(int i=0;i<refCardinality;i++)
        for(int j=0;j<dim;j++)
          for(int k=0;k<hgradCardinality;k++)
            ifunGradAtRefCoords(i,j) += funHGradCoeffs(k)*gradOfHGradBasisAtRefCoords(k,i,j);


      //Computing the gradient of the (hgrad) interpolated function (ifun) at HCURL dof coordinates
      DynRankView ConstructWithLabel(gradOfHGradBasisAtHCurlDofCoords, hgradCardinality , hcurlCardinality, dim);
      hexHGradBasis.getValues(gradOfHGradBasisAtHCurlDofCoords, hcurlDofCoords, OPERATOR_GRAD);
      DynRankView ConstructWithLabel(ifunGradAtHCurlDofCoords, hcurlCardinality,dim);
      for(int i=0;i<hcurlCardinality;i++)
        for(int j=0;j<dim;j++)
          for(int k=0;k<hgradCardinality;k++)
            ifunGradAtHCurlDofCoords(i,j) += funHGradCoeffs(k)*gradOfHGradBasisAtHCurlDofCoords(k,i,j);


      //Interpolating the gradient of ifun in the HCURL space by computing the degrees of freedom for the function
      DynRankView ConstructWithLabel(ifunGradHCurlCoeffs, hcurlCardinality);
      for(int i=0;i<hcurlCardinality;i++)
        for(int j=0;j<dim;j++)
          ifunGradHCurlCoeffs(i) += ifunGradAtHCurlDofCoords(i,j)*hcurlDofCoeffs(i,j);

      //Evaluating the gradient of ifun in the HCURL space at HGRAD dof coordinates and compare it with the gradient of ifun at the same coords
      //We note if the two representations of the gradient of ifun are equal at the HGRAD dof coordinates, they are equal everywhere.
      DynRankView ConstructWithLabel(hcurlBasisAtRefCoords, hcurlCardinality , refCardinality, dim);
      hexHCurlBasis.getValues(hcurlBasisAtRefCoords, refCoords, OPERATOR_VALUE);
      DynRankView ConstructWithLabel(ifunGradInHCurlSpaceAtRefCoords, refCardinality,dim);
      pointValueType diffErr(0);
      for(int i=0;i<refCardinality;i++)
        for(int j=0;j<dim;j++) {
          for(int k=0;k<hcurlCardinality;k++)
            ifunGradInHCurlSpaceAtRefCoords(i,j) += ifunGradHCurlCoeffs(k)*hcurlBasisAtRefCoords(k,i,j);
          diffErr = std::max(diffErr, std::abs(ifunGradInHCurlSpaceAtRefCoords(i,j) - ifunGradAtRefCoords(i,j)));
        }

      //Compute the curl of the gradient of ifun in HCURL space, and evaluate it at HGRAD dof coordinates.
      //We note if the curl of the gradient of ifun is zero at the HGRAD dof coordinates, is identically zero.
      DynRankView ConstructWithLabel(curlOfHCurlBasisAtRefCoords, hcurlCardinality , refCardinality, dim);
      hexHCurlBasis.getValues(curlOfHCurlBasisAtRefCoords, hgradDofCoords, OPERATOR_CURL);
      DynRankView ConstructWithLabel(ifunCurlGradAtRefCoords, refCardinality,dim);
      pointValueType maxNorm(0);

      for(int i=0;i<refCardinality;i++)
        for(int j=0;j<dim;j++) {
          for(int k=0;k<hcurlCardinality;k++)
            ifunCurlGradAtRefCoords(i,j) += ifunGradHCurlCoeffs(k)*curlOfHCurlBasisAtRefCoords(k,i,j);
          maxNorm = std::max(maxNorm, std::abs(ifunCurlGradAtRefCoords(i,j)));
        }

      //Check that the two representations of the gradient of ifun are consistent
      if(diffErr > pow(7, order-1)*tol) { //heuristic relation on how round-off error depends on order
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << "Grad of HGRAD_C" << order << " function does not belong to HCURL_I" << order << "."<<
            "\nIn fact the HCURL interpolation of the function gradient is different from the function gradient."<<
            "\nThe max norm of the difference at HGRAD DOF coordinates is: " <<  diffErr << std::endl;
      }

      //Check that the curl of the grad of ifun is zero
      if(maxNorm > pow(7, order)*tol) { //heuristic relation on how round-off error depends on order
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << "The Curl Grad of a function is 0, but taking the gradient of the function in the HGRAD_C" << order << " space,"
            "\ninterpolating it on the HCURL_I" << order << " space and taking the curl gives a nonzero function."
            "\nIts max norm at HGRAD DOF coordinates is: " << maxNorm << std::endl;;
      }
    }
  } catch (std::exception err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };


  try {

    *outStream
    << "\n"
    << "=======================================================================================================\n"
    << "| TEST 2: Testing that curl of HCURL_In function is in the kernel the div operator in HDIV_In space   |\n"
    << "=======================================================================================================\n";

    /*
     * For every order n<=maxOrder, we check that the curl of a function in HCURL_In space is in the HDIV_In space.
     * Then we take its div, in the HDIV_In space, and verify is zero.
     */

    for (ordinal_type order=1; order <= maxOrder; ++order ) {

      //Getting ref coordinates (HGRAD dof coordinates):
      //Point evaluations at ref coordinates uniquely identify a polynomial of degree n.
      HexHGRADBasisType hexHGradBasis(order, POINTTYPE_WARPBLEND);
      const ordinal_type refCardinality = hexHGradBasis.getCardinality();
      DynRankView ConstructWithLabel(refCoords, refCardinality, dim);
      hexHGradBasis.getDofCoords(refCoords);

      //Creating HCURL and HDIV basis
      HexHCURLBasisType hexHCurlBasis(order, POINTTYPE_EQUISPACED);  //any point type is fine
      HexHDIVBasisType hexHDivBasis(order, POINTTYPE_EQUISPACED);   //any point type is fine

      const ordinal_type hcurlCardinality = hexHCurlBasis.getCardinality();
      const ordinal_type hdivCardinality = hexHDivBasis.getCardinality();


      //Getting DOF coordinates for HCURL and HDIV elements
      DynRankView ConstructWithLabel(hcurlDofCoords, hcurlCardinality, dim);
      hexHCurlBasis.getDofCoords(hcurlDofCoords);

      DynRankView ConstructWithLabel(hdivDofCoords, hdivCardinality, dim);
      hexHDivBasis.getDofCoords(hdivDofCoords);

      //Getting DOF coefficients for HCURL and HDIV elements
      DynRankView ConstructWithLabel(hcurlDofCoeffs, hcurlCardinality, dim);
      hexHCurlBasis.getDofCoeffs(hcurlDofCoeffs);

      DynRankView ConstructWithLabel(hdivDofCoeffs, hdivCardinality, dim);
      hexHDivBasis.getDofCoeffs(hdivDofCoeffs);

      //Evaluating the function at HCURL dof coordinates
      DynRankView ConstructWithLabel(funAtHCurlDofCoords, hcurlCardinality, dim);

      Fun fun;
      for(int i=0;i<hcurlCardinality;i++)
        for(int j=0;j<dim;j++)
          funAtHCurlDofCoords(i,j) = fun(hcurlDofCoords(i,0), hcurlDofCoords(i,1), hcurlDofCoords(i,2), j);

      //Interpolating the function in the HCURL space by computing the degrees of freedom for the function
      DynRankView ConstructWithLabel(funHCurlCoeffs, hcurlCardinality);
      for(int i=0;i<hcurlCardinality;i++)
        for(int j=0;j<dim;j++)
          funHCurlCoeffs(i) += funAtHCurlDofCoords(i,j)*hcurlDofCoeffs(i,j);

      //Computing the curl of the (hcurl) interpolated function (ifun) at ref coordinates
      DynRankView ConstructWithLabel(curlOfHCurlBasisAtRefCoords, hcurlCardinality , refCardinality, dim);
      hexHCurlBasis.getValues(curlOfHCurlBasisAtRefCoords, refCoords, OPERATOR_CURL);
      DynRankView ConstructWithLabel(ifunCurlAtRefCoords, refCardinality, dim);
      for(int i=0;i<refCardinality;i++)
        for(int j=0;j<dim;j++)
          for(int k=0;k<hcurlCardinality;k++)
            ifunCurlAtRefCoords(i,j) += funHCurlCoeffs(k)*curlOfHCurlBasisAtRefCoords(k,i,j);


      //Computing the curl of the (hcurl) interpolated function (ifun) at HDIV dof coordinates
      DynRankView ConstructWithLabel(curlOfHCurlBasisAtHDivDofCoords, hcurlCardinality , hdivCardinality, dim);
      hexHCurlBasis.getValues(curlOfHCurlBasisAtHDivDofCoords, hdivDofCoords, OPERATOR_CURL);
      DynRankView ConstructWithLabel(ifunCurlAtHDivDofCoords, hdivCardinality,dim);
      for(int i=0;i<hdivCardinality;i++)
        for(int j=0;j<dim;j++)
          for(int k=0;k<hcurlCardinality;k++)
            ifunCurlAtHDivDofCoords(i,j) += funHCurlCoeffs(k)*curlOfHCurlBasisAtHDivDofCoords(k,i,j);


      //Interpolating the curl of ifun in the HDIV space by computing the degrees of freedom for the function
      DynRankView ConstructWithLabel(ifunCurlHDivCoeffs, hdivCardinality);
      for(int i=0;i<hdivCardinality;i++)
        for(int j=0;j<dim;j++)
          ifunCurlHDivCoeffs(i) += ifunCurlAtHDivDofCoords(i,j)*hdivDofCoeffs(i,j);

      //Evaluating the curl of ifun in the HDIV space at HCURL dof coordinates and compare it with the curl of ifun at the same coords
      //We note if the two representations of the curl of ifun are equal at the HCURL dof coordinates, they are equal everywhere.
      DynRankView ConstructWithLabel(hdivBasisAtRefCoords, hdivCardinality , refCardinality, dim);
      hexHDivBasis.getValues(hdivBasisAtRefCoords, refCoords, OPERATOR_VALUE);
      DynRankView ConstructWithLabel(ifunCurlInHDivSpaceAtRefCoords, refCardinality,dim);
      pointValueType diffErr(0);
      for(int i=0;i<refCardinality;i++)
        for(int j=0;j<dim;j++) {
          for(int k=0;k<hdivCardinality;k++)
            ifunCurlInHDivSpaceAtRefCoords(i,j) += ifunCurlHDivCoeffs(k)*hdivBasisAtRefCoords(k,i,j);
          diffErr = std::max(diffErr, std::abs(ifunCurlInHDivSpaceAtRefCoords(i,j) - ifunCurlAtRefCoords(i,j)));
        }

      //Compute the div of the curl of ifun in HDIV space, and evaluate it at HCURL dof coordinates.
      //We note if the div of the curl of ifun is zero at the HCURL dof coordinates, is identically zero.
      DynRankView ConstructWithLabel(divOfHDivBasisAtRefCoords, hdivCardinality , refCardinality);
      hexHDivBasis.getValues(divOfHDivBasisAtRefCoords, refCoords, OPERATOR_DIV);
      DynRankView ConstructWithLabel(ifunDivCurlAtRefCoords, refCardinality);
      pointValueType maxNorm(0);

      for(int i=0;i<refCardinality;i++) {
        for(int k=0;k<hdivCardinality;k++)
          ifunDivCurlAtRefCoords(i) += ifunCurlHDivCoeffs(k)*divOfHDivBasisAtRefCoords(k,i);
        maxNorm = std::max(maxNorm, std::abs(ifunDivCurlAtRefCoords(i)));
      }

      //Check that the two representations of the curl of ifun are consistent
      if(diffErr > pow(7, order-1)*tol) { //heuristic relation on how round-off error depends on order
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << "Curl of HCURL_I" << order << " function does not belong to HDIV_I" << order << "."<<
            "\nIn fact the HDIV interpolation of the function curl is different from the function curl."<<
            "\nThe max norm of the difference at HCURL DOF coordinates is: " <<  diffErr << std::endl;
      }

      //Check that the div of the curl of ifun is zero
      if(maxNorm > pow(7, order)*tol) { //heuristic relation on how round-off error depends on order
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << "The Div Curl of a function is 0, but taking the curl of the function in the HCURL_I" << order << " space,"
            "\ninterpolating it on the HDIV_I" << order << " space and taking the div gives a nonzero function."
            "\nIts max norm at HCURL DOF coordinates is: " << maxNorm << std::endl;;
      }
    }
  } catch (std::exception err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };


  try {

    *outStream
    << "\n"
    << "=======================================================================================================\n"
    << "| TEST 3: Testing that the divergence of HDIV_In functions belong to HVOL_C(n-1) space                  |\n"
    << "=======================================================================================================\n";

    /*
     * For every order n<=maxOrder, we check that the div of a function in HDIV_In space is in the HVOL_C(n-1) space.
     */

    for (ordinal_type order=1; order <= maxOrder; ++order ) {

      //Getting ref coordinates (HGRAD dof coordinates):
      //Point evaluations at ref coordinates uniquely identify a polynomial of degree n.
      HexHGRADBasisType hexHGradBasis(order, POINTTYPE_WARPBLEND);
      const ordinal_type refCardinality = hexHGradBasis.getCardinality();
      DynRankView ConstructWithLabel(refCoords, refCardinality, dim);
      hexHGradBasis.getDofCoords(refCoords);

      //Creating HDIV and HVOL basis
      HexHDIVBasisType hexHDivBasis(order, POINTTYPE_EQUISPACED);  //any point type is fine
      HexHVOLBasisType hexHVOLBasis(order-1, POINTTYPE_EQUISPACED);   //any point type is fine

      const ordinal_type hdivCardinality = hexHDivBasis.getCardinality();
      const ordinal_type hvolCardinality = hexHVOLBasis.getCardinality();

      //Getting DOF coordinates for HDIV and HVOL elements
      DynRankView ConstructWithLabel(hdivDofCoords, hdivCardinality, dim);
      hexHDivBasis.getDofCoords(hdivDofCoords);

      DynRankView ConstructWithLabel(hvolDofCoords, hvolCardinality, dim);
      hexHVOLBasis.getDofCoords(hvolDofCoords);

      //Getting DOF coefficients for HDIV and HVOL elements
      DynRankView ConstructWithLabel(hdivDofCoeffs, hdivCardinality, dim);
      hexHDivBasis.getDofCoeffs(hdivDofCoeffs);

      DynRankView ConstructWithLabel(hvolDofCoeffs, hvolCardinality);
      hexHVOLBasis.getDofCoeffs(hvolDofCoeffs);

      //Evaluating the function at HDIV dof coordinates
      DynRankView ConstructWithLabel(funAtHDivDofCoords, hdivCardinality, dim);

      Fun fun;
      for(int i=0;i<hdivCardinality;i++)
        for(int j=0;j<dim;j++)
          funAtHDivDofCoords(i,j) = fun(hdivDofCoords(i,0), hdivDofCoords(i,1), hdivDofCoords(i,2), j);

      //Interpolating the function in the HDIV space by computing the degrees of freedom for the function
      DynRankView ConstructWithLabel(funHDivCoeffs, hdivCardinality);
      for(int i=0;i<hdivCardinality;i++)
        for(int j=0;j<dim;j++)
          funHDivCoeffs(i) += funAtHDivDofCoords(i,j)*hdivDofCoeffs(i,j);

      //Computing the div of the (hdiv) interpolated function (ifun) at HDIV dof coordinates
      DynRankView ConstructWithLabel(divOfHDivBasisAtRefCoords, hdivCardinality , refCardinality);
      hexHDivBasis.getValues(divOfHDivBasisAtRefCoords, refCoords, OPERATOR_DIV);
      DynRankView ConstructWithLabel(ifunDivAtRefCoords, refCardinality);
      for(int i=0;i<refCardinality;i++)
        for(int k=0;k<hdivCardinality;k++)
          ifunDivAtRefCoords(i) += funHDivCoeffs(k)*divOfHDivBasisAtRefCoords(k,i);


      //Computing the div of the (hdiv) interpolated function (ifun) at HVOL dof coordinates
      DynRankView ConstructWithLabel(divOfHDivBasisAtHVOLDofCoords, hdivCardinality , hvolCardinality);
      hexHDivBasis.getValues(divOfHDivBasisAtHVOLDofCoords, hvolDofCoords, OPERATOR_DIV);
      DynRankView ConstructWithLabel(ifunDivAtHVOLDofCoords, hvolCardinality);
      for(int i=0;i<hvolCardinality;i++)
        for(int k=0;k<hdivCardinality;k++)
          ifunDivAtHVOLDofCoords(i) += funHDivCoeffs(k)*divOfHDivBasisAtHVOLDofCoords(k,i);


      //Interpolating the div of ifun in the HVOL space by computing the degrees of freedom for the function
      DynRankView ConstructWithLabel(ifunDivHVOLCoeffs, hvolCardinality);
      for(int i=0;i<hvolCardinality;i++)
        ifunDivHVOLCoeffs(i) += ifunDivAtHVOLDofCoords(i)*hvolDofCoeffs(i); //not really needed for HVOL, hvolDofCoeffs = 1

      //Evaluating the div of ifun in the HVOL space at HDIV dof coordinates and compare it with the div of ifun at the same coords
      //We note if the two representations of the div of ifun are equal at the HDIV dof coordinates, they are equal everywhere.
      DynRankView ConstructWithLabel(hvolBasisAtRefCoords, hvolCardinality , refCardinality);
      hexHVOLBasis.getValues(hvolBasisAtRefCoords, refCoords, OPERATOR_VALUE);
      DynRankView ConstructWithLabel(ifunDivInHVOLSpaceAtRefCoords, refCardinality);
      pointValueType diffErr(0);
      for(int i=0;i<refCardinality;i++) {
        for(int k=0;k<hvolCardinality;k++)
          ifunDivInHVOLSpaceAtRefCoords(i) += ifunDivHVOLCoeffs(k)*hvolBasisAtRefCoords(k,i);
        diffErr = std::max(diffErr, std::abs(ifunDivInHVOLSpaceAtRefCoords(i) - ifunDivAtRefCoords(i)));
      }

      //Check that the two representations of the div of ifun are consistent
      if(diffErr > pow(7, order-1)*tol) { //heuristic relation on how round-off error depends on order
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << "Div of HDIV_I" << order << " function does not belong to HVOL_C" << order-1 << "."<<
            "\nIn fact the HVOL interpolation of the function div is different from the function div."<<
            "\nThe max norm of the difference at HDIV DOF coordinates is: " <<  diffErr << std::endl;
      }

    }
  } catch (std::exception err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
}
}

