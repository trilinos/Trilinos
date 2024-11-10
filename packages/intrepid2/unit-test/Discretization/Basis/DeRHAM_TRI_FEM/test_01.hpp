// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Test for checking the De Rham complex for FEM on TRI
            Testing that HGRAD_Cn --(grad)--> HCURL_In --(curl)--> HDIV_In --(div)--
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
#include "Intrepid2_HGRAD_TRI_Cn_FEM.hpp"
#include "Intrepid2_HCURL_TRI_In_FEM.hpp"
#include "Intrepid2_HDIV_TRI_In_FEM.hpp"
#include "Intrepid2_HVOL_TRI_Cn_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "packages/intrepid2/unit-test/Discretization/Basis/Macros.hpp"

namespace Intrepid2 {

namespace Test {

//Warning: This is not supposed to be build with CUDA.


template<typename ValueType, typename DeviceSpaceType>
int DeRHAM_TRI_FEM_Test01(const bool verbose) {

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

  const ValueType tol = tolerence();
  int errorFlag = 0;

  // for virtual function, value and point types are declared in the class
  typedef ValueType outputValueType;
  typedef ValueType pointValueType;

  typedef Basis_HGRAD_TRI_Cn_FEM<DeviceSpaceType,outputValueType,pointValueType> TriHGRADBasisType;
  typedef Basis_HCURL_TRI_In_FEM<DeviceSpaceType,outputValueType,pointValueType> TriHCURLBasisType;
  typedef Basis_HDIV_TRI_In_FEM<DeviceSpaceType,outputValueType,pointValueType> TriHDIVBasisType;
  typedef Basis_HVOL_TRI_Cn_FEM<DeviceSpaceType,outputValueType,pointValueType> TriHVOLBasisType;

  constexpr ordinal_type dim = 2;
  constexpr ordinal_type maxOrder = Intrepid2::Parameters::MaxOrder < 5 ? Intrepid2::Parameters::MaxOrder : 5;

  // In order to obtain nontrivial FE functions, we interpolate the following (vector) function on the space spanned by the fE basis considered.
  struct Fun {
    pointValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const pointValueType& x, const pointValueType& y, const int comp=0) {
      switch (comp) {
      case 0:
        return std::sin(x)*std::sin(y)+2.0*std::cos(x)*std::cos(y);
      case 1:
        return std::exp(x)*std::exp(y);
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
      TriHGRADBasisType triHGradBasisWarp(order, POINTTYPE_WARPBLEND);
      const ordinal_type refCardinality = triHGradBasisWarp.getCardinality();
      DynRankView ConstructWithLabel(refCoords, refCardinality, dim);
      triHGradBasisWarp.getDofCoords(refCoords);

      //Creating HGRAD and HCURL basis
      TriHGRADBasisType triHGradBasis(order, POINTTYPE_EQUISPACED);  //any point type is fine
      TriHCURLBasisType triHCurlBasis(order, POINTTYPE_WARPBLEND); //any point type is fine

      const ordinal_type hgradCardinality = triHGradBasis.getCardinality();
      const ordinal_type hcurlCardinality = triHCurlBasis.getCardinality();

      //Getting DOF coordinates for HGRAD and HCURL elements
      DynRankView ConstructWithLabel(hgradDofCoords, hgradCardinality , dim);
      triHGradBasis.getDofCoords(hgradDofCoords);

      DynRankView ConstructWithLabel(hcurlDofCoords, hcurlCardinality , dim);
      triHCurlBasis.getDofCoords(hcurlDofCoords);

      //Getting DOF coefficients for HGRAD and HCURL elements
      DynRankView ConstructWithLabel(hgradDofCoeffs, hgradCardinality);
      triHGradBasis.getDofCoeffs(hgradDofCoeffs);

      DynRankView ConstructWithLabel(hcurlDofCoeffs, hcurlCardinality, dim);
      triHCurlBasis.getDofCoeffs(hcurlDofCoeffs);

      //Evaluating the function at HGRAD dof coordinates
      DynRankView ConstructWithLabel(funAtHGradDofCoords, hgradCardinality);

      Fun fun;
      for(int i=0;i<hgradCardinality;i++)
        funAtHGradDofCoords(i) = fun(hgradDofCoords(i,0), hgradDofCoords(i,1));

      //Interpolating the function in the HGRAD space by computing the degrees of freedom for the function
      DynRankView ConstructWithLabel(funHGradCoeffs, hgradCardinality);
      for(int i=0;i<hgradCardinality;i++)
        funHGradCoeffs(i) = funAtHGradDofCoords(i)*hgradDofCoeffs(i); //not really needed for HGRAD, hgradDofCoeffs = 1

      //Computing the gradient of the (hgrad) interpolated function (ifun) at HGRAD dof coordinates
      DynRankView ConstructWithLabel(gradOfHGradBasisAtRefCoords, hgradCardinality , refCardinality, dim);
      triHGradBasis.getValues(gradOfHGradBasisAtRefCoords, refCoords, OPERATOR_GRAD);
      DynRankView ConstructWithLabel(ifunGradAtRefCoords, refCardinality, dim);
      for(int i=0;i<refCardinality;i++)
        for(int j=0;j<dim;j++)
          for(int k=0;k<hgradCardinality;k++)
            ifunGradAtRefCoords(i,j) += funHGradCoeffs(k)*gradOfHGradBasisAtRefCoords(k,i,j);


      //Computing the gradient of the (hgrad) interpolated function (ifun) at HCURL dof coordinates
      DynRankView ConstructWithLabel(gradOfHGradBasisAtHCurlDofCoords, hgradCardinality , hcurlCardinality, dim);
      triHGradBasis.getValues(gradOfHGradBasisAtHCurlDofCoords, hcurlDofCoords, OPERATOR_GRAD);
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
      triHCurlBasis.getValues(hcurlBasisAtRefCoords, refCoords, OPERATOR_VALUE);
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
      DynRankView ConstructWithLabel(curlOfHCurlBasisAtRefCoords, hcurlCardinality , refCardinality);
      triHCurlBasis.getValues(curlOfHCurlBasisAtRefCoords, hgradDofCoords, OPERATOR_CURL);
      DynRankView ConstructWithLabel(ifunCurlGradAtRefCoords, refCardinality);
      pointValueType maxNorm(0);

      for(int i=0;i<refCardinality;i++) {
         for(int k=0;k<hcurlCardinality;k++)
            ifunCurlGradAtRefCoords(i) += ifunGradHCurlCoeffs(k)*curlOfHCurlBasisAtRefCoords(k,i);
          maxNorm = std::max(maxNorm, std::abs(ifunCurlGradAtRefCoords(i)));
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
  } catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };


  try {

    *outStream
    << "\n"
    << "=======================================================================================================\n"
    << "| TEST 2: Testing that curl of HGRAD_Cn function is in the kernel the div operator in HDIV_In space   |\n"
    << "=======================================================================================================\n";

    /*
     * For every order n<=maxOrder, we check that the curl (grad perp) of a function in HGRAD_In space is in the 
     * HDIV_In space. Then we take its div, in the HDIV_In space, and verify is zero.
     * Note that for 2-D scalar fields, curl u := (-du/dy, du/dx)
     */

    for (ordinal_type order=1; order <= maxOrder; ++order ) {

      //Getting ref coordinates (HGRAD dof coordinates):
      //Point evaluations at ref coordinates uniquely identify a polynomial of degree n.
      TriHGRADBasisType triHGradBasisWarp(order, POINTTYPE_WARPBLEND);
      const ordinal_type refCardinality = triHGradBasisWarp.getCardinality();
      DynRankView ConstructWithLabel(refCoords, refCardinality, dim);
      triHGradBasisWarp.getDofCoords(refCoords);

      //Creating HGRAD and HDIV basis
      TriHGRADBasisType triHGradBasis(order, POINTTYPE_EQUISPACED);  //any point type is fine
      TriHDIVBasisType triHDivBasis(order, POINTTYPE_WARPBLEND);   //any point type is fine

      const ordinal_type hgradCardinality = triHGradBasis.getCardinality();
      const ordinal_type hdivCardinality = triHDivBasis.getCardinality();


      //Getting DOF coordinates for HGRAD and HDIV elements
      DynRankView ConstructWithLabel(hgradDofCoords, hgradCardinality, dim);
      triHGradBasis.getDofCoords(hgradDofCoords);

      DynRankView ConstructWithLabel(hdivDofCoords, hdivCardinality, dim);
      triHDivBasis.getDofCoords(hdivDofCoords);

      //Getting DOF coefficients for HGRAD and HDIV elements
      DynRankView ConstructWithLabel(hgradDofCoeffs, hgradCardinality);
      triHGradBasis.getDofCoeffs(hgradDofCoeffs);

      DynRankView ConstructWithLabel(hdivDofCoeffs, hdivCardinality, dim);
      triHDivBasis.getDofCoeffs(hdivDofCoeffs);

      //Evaluating the function at HGRAD dof coordinates
      DynRankView ConstructWithLabel(funAtHGradDofCoords, hgradCardinality);

      Fun fun;
      for(int i=0;i<hgradCardinality;i++)
          funAtHGradDofCoords(i) = fun(hgradDofCoords(i,0), hgradDofCoords(i,1));

      //Interpolating the function in the HGRAD space by computing the degrees of freedom for the function
      DynRankView ConstructWithLabel(funHGradCoeffs, hgradCardinality);
      for(int i=0;i<hgradCardinality;i++)
        funHGradCoeffs(i) = funAtHGradDofCoords(i)*hgradDofCoeffs(i); //not really needed for HGRAD, hgradDofCoeffs = 1

      //Computing the curl of the (hgrad) interpolated function (ifun) at ref coordinates
      DynRankView ConstructWithLabel(curlOfHGradBasisAtRefCoords, hgradCardinality , refCardinality, dim);
      triHGradBasis.getValues(curlOfHGradBasisAtRefCoords, refCoords, OPERATOR_CURL);
      DynRankView ConstructWithLabel(ifunCurlAtRefCoords, refCardinality, dim);
      for(int i=0;i<refCardinality;i++) {
        for(int k=0;k<hgradCardinality;k++)
          for(int j=0;j<dim;j++)
            ifunCurlAtRefCoords(i,j) += funHGradCoeffs(k)*curlOfHGradBasisAtRefCoords(k,i,j);
      }

      //Computing the curl of the (hgrad) interpolated function (ifun) at HDIV dof coordinates
      DynRankView ConstructWithLabel(curlOfHGradBasisAtHDivDofCoords, hgradCardinality , hdivCardinality, dim);
      triHGradBasis.getValues(curlOfHGradBasisAtHDivDofCoords, hdivDofCoords, OPERATOR_CURL);
      DynRankView ConstructWithLabel(ifunCurlAtHDivDofCoords, hdivCardinality,dim);
      for(int i=0;i<hdivCardinality;i++) {
        for(int k=0;k<hgradCardinality;k++)
          for(int j=0;j<dim;j++)
            ifunCurlAtHDivDofCoords(i,j) += funHGradCoeffs(k)*curlOfHGradBasisAtHDivDofCoords(k,i,j);
      }

      //Interpolating the curl of ifun in the HDIV space by computing the degrees of freedom for the function
      DynRankView ConstructWithLabel(ifunCurlHDivCoeffs, hdivCardinality);
      for(int i=0;i<hdivCardinality;i++)
        for(int j=0;j<dim;j++)
          ifunCurlHDivCoeffs(i) += ifunCurlAtHDivDofCoords(i,j)*hdivDofCoeffs(i,j);

      //Evaluating the curl of ifun in the HDIV space at HGRAD dof coordinates and compare it with the curl of ifun at the same coords
      //We note if the two representations of the curl of ifun are equal at the HGRAD dof coordinates, they are equal everywhere.
      DynRankView ConstructWithLabel(hdivBasisAtRefCoords, hdivCardinality , refCardinality, dim);
      triHDivBasis.getValues(hdivBasisAtRefCoords, refCoords, OPERATOR_VALUE);
      DynRankView ConstructWithLabel(ifunCurlInHDivSpaceAtRefCoords, refCardinality,dim);
      pointValueType diffErr(0);
      for(int i=0;i<refCardinality;i++)
        for(int j=0;j<dim;j++) {
          for(int k=0;k<hdivCardinality;k++)
            ifunCurlInHDivSpaceAtRefCoords(i,j) += ifunCurlHDivCoeffs(k)*hdivBasisAtRefCoords(k,i,j);
          diffErr = std::max(diffErr, std::abs(ifunCurlInHDivSpaceAtRefCoords(i,j) - ifunCurlAtRefCoords(i,j)));
        }

      //Compute the div of the curl of ifun in HDIV space, and evaluate it at HGRAD dof coordinates.
      //We note if the div of the curl of ifun is zero at the HGRAD dof coordinates, is identically zero.
      DynRankView ConstructWithLabel(divOfHDivBasisAtRefCoords, hdivCardinality , refCardinality);
      triHDivBasis.getValues(divOfHDivBasisAtRefCoords, refCoords, OPERATOR_DIV);
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
        *outStream << "Curl of HGRAD_C" << order << " function does not belong to HDIV_I" << order << "."<<
            "\nIn fact the HDIV interpolation of the function curl is different from the function curl."<<
            "\nThe max norm of the difference at HGRAD DOF coordinates is: " <<  diffErr << std::endl;
      }

      //Check that the div of the curl of ifun is zero
      if(maxNorm > pow(7, order)*tol) { //heuristic relation on how round-off error depends on order
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << "The Div Curl of a function is 0, but taking the curl of the function in the HGRAD_C" << order << " space,"
            "\ninterpolating it on the HDIV_I" << order << " space and taking the div gives a nonzero function."
            "\nIts max norm at HGRAD DOF coordinates is: " << maxNorm << std::endl;;
      }
    }
  } catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };


  try {

    *outStream
    << "\n"
    << "=======================================================================================================\n"
    << "| TEST 3: Testing that the divergence of HDIV_In functions and the curl of HCURL_In functions         |\n"
    << "|         belong to HVOL_C(n-1) space                                                                   |\n"
    << "=======================================================================================================\n";

    /*
     * For every order n<=maxOrder, we check that the div of a function in HDIV_In space and the curl of a
     * function in HCURL_In space is in the HVOL_C(n-1) space.
     */

    for (ordinal_type order=1; order <= maxOrder; ++order ) {

      //Getting ref coordinates (HGRAD dof coordinates):
      //Point evaluations at ref coordinates uniquely identify a polynomial of degree n.
      TriHGRADBasisType triHGradBasis(order, POINTTYPE_WARPBLEND);
      const ordinal_type refCardinality = triHGradBasis.getCardinality();
      DynRankView ConstructWithLabel(refCoords, refCardinality, dim);
      triHGradBasis.getDofCoords(refCoords);

      //Creating HDIV, HCURL and HVOL basis
      TriHDIVBasisType triHDivBasis(order, POINTTYPE_EQUISPACED);  //any point type is fine
      TriHCURLBasisType triHCurlBasis(order, POINTTYPE_EQUISPACED);  //any point type is fine
      TriHVOLBasisType triHVOLBasis(order-1, POINTTYPE_WARPBLEND);   //any point type is fine

      const ordinal_type hdivCardinality = triHDivBasis.getCardinality();
      const ordinal_type hcurlCardinality = triHCurlBasis.getCardinality();
      const ordinal_type hvolCardinality = triHVOLBasis.getCardinality();

      //Getting DOF coordinates for HDIV, HCURL and HVOL elements
      DynRankView ConstructWithLabel(hdivDofCoords, hdivCardinality, dim);
      triHDivBasis.getDofCoords(hdivDofCoords);

      DynRankView ConstructWithLabel(hcurlDofCoords, hcurlCardinality, dim);
      triHCurlBasis.getDofCoords(hcurlDofCoords);

      DynRankView ConstructWithLabel(hvolDofCoords, hvolCardinality, dim);
      triHVOLBasis.getDofCoords(hvolDofCoords);

      //Getting DOF coefficients for HDIV, HCURL and HVOL elements
      DynRankView ConstructWithLabel(hdivDofCoeffs, hdivCardinality, dim);
      triHDivBasis.getDofCoeffs(hdivDofCoeffs);

      DynRankView ConstructWithLabel(hcurlDofCoeffs, hcurlCardinality, dim);
      triHCurlBasis.getDofCoeffs(hcurlDofCoeffs);

      DynRankView ConstructWithLabel(hvolDofCoeffs, hvolCardinality);
      triHVOLBasis.getDofCoeffs(hvolDofCoeffs);

      //Evaluating the function at HDIV dof coordinates
      DynRankView ConstructWithLabel(funAtHDivDofCoords, hdivCardinality, dim);

      Fun fun;
      for(int i=0;i<hdivCardinality;i++)
        for(int j=0;j<dim;j++)
          funAtHDivDofCoords(i,j) = fun(hdivDofCoords(i,0), hdivDofCoords(i,1), j);

      //Interpolating the function in the HDIV space by computing the degrees of freedom for the function
      DynRankView ConstructWithLabel(funHDivCoeffs, hdivCardinality);
      for(int i=0;i<hdivCardinality;i++)
        for(int j=0;j<dim;j++)
          funHDivCoeffs(i) += funAtHDivDofCoords(i,j)*hdivDofCoeffs(i,j);

      //Computing the div of the (hdiv) interpolated function (ifun) at HDIV dof coordinates
      DynRankView ConstructWithLabel(divOfHDivBasisAtRefCoords, hdivCardinality , refCardinality);
      triHDivBasis.getValues(divOfHDivBasisAtRefCoords, refCoords, OPERATOR_DIV);
      DynRankView ConstructWithLabel(ifunDivAtRefCoords, refCardinality);
      for(int i=0;i<refCardinality;i++)
        for(int k=0;k<hdivCardinality;k++)
          ifunDivAtRefCoords(i) += funHDivCoeffs(k)*divOfHDivBasisAtRefCoords(k,i);


      //Computing the div of the (hdiv) interpolated function (ifun) at HVOL dof coordinates
      DynRankView ConstructWithLabel(divOfHDivBasisAtHVOLDofCoords, hdivCardinality , hvolCardinality);
      triHDivBasis.getValues(divOfHDivBasisAtHVOLDofCoords, hvolDofCoords, OPERATOR_DIV);
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
      triHVOLBasis.getValues(hvolBasisAtRefCoords, refCoords, OPERATOR_VALUE);
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

      //Evaluating the function at HCURL dof coordinates
      DynRankView ConstructWithLabel(funAtHCurlDofCoords, hcurlCardinality, dim);

      for(int i=0;i<hcurlCardinality;i++)
        for(int j=0;j<dim;j++)
           funAtHCurlDofCoords(i,j) = fun(hcurlDofCoords(i,0), hcurlDofCoords(i,1), j);

      //Interpolating the function in the HCURL space by computing the degrees of freedom for the function
      DynRankView ConstructWithLabel(funHCurlCoeffs, hcurlCardinality);
      for(int i=0;i<hcurlCardinality;i++)
        for(int j=0;j<dim;j++)
          funHCurlCoeffs(i) += funAtHCurlDofCoords(i,j)*hcurlDofCoeffs(i,j);
     
      //Computing the curl of the (hcurl) interpolated function (ifun) at HCURL dof coordinates
      DynRankView ConstructWithLabel(curlOfHCurlBasisAtRefCoords, hcurlCardinality , refCardinality);
      triHCurlBasis.getValues(curlOfHCurlBasisAtRefCoords, refCoords, OPERATOR_CURL);
      DynRankView ConstructWithLabel(ifunCurlAtRefCoords, refCardinality);
      for(int i=0;i<refCardinality;i++)
        for(int k=0;k<hcurlCardinality;k++)
          ifunCurlAtRefCoords(i) += funHCurlCoeffs(k)*curlOfHCurlBasisAtRefCoords(k,i);
    
      //Computing the curl of the (hcurl) interpolated function (ifun) at HVOL dof coordinates
      DynRankView ConstructWithLabel(curlOfHCurlBasisAtHVOLDofCoords, hcurlCardinality , hvolCardinality);
      triHCurlBasis.getValues(curlOfHCurlBasisAtHVOLDofCoords, hvolDofCoords, OPERATOR_CURL);
      DynRankView ConstructWithLabel(ifunCurlAtHVOLDofCoords, hvolCardinality);
      for(int i=0;i<hvolCardinality;i++)
        for(int k=0;k<hcurlCardinality;k++)
          ifunCurlAtHVOLDofCoords(i) += funHCurlCoeffs(k)*curlOfHCurlBasisAtHVOLDofCoords(k,i);
      

      //Interpolating the curl of ifun in the HVOL space by computing the degrees of freedom for the function
      DynRankView ConstructWithLabel(ifunCurlHVOLCoeffs, hvolCardinality);
      for(int i=0;i<hvolCardinality;i++)
        ifunCurlHVOLCoeffs(i) += ifunCurlAtHVOLDofCoords(i)*hvolDofCoeffs(i); //not really needed for HVOL, hvolDofCoeffs = 1
    
      //Evaluating the curl of ifun in the HVOL space at HCURL dof coordinates and compare it with the curl of ifun at the same coords
      //We note if the two representations of the curl of ifun are equal at the HCURL dof coordinates, they are equal everywhere.
      DynRankView ConstructWithLabel(ifunCurlInHVOLSpaceAtRefCoords, refCardinality);
      pointValueType diffErrCurl(0);
      for(int i=0;i<refCardinality;i++) {
        for(int k=0;k<hvolCardinality;k++)
          ifunCurlInHVOLSpaceAtRefCoords(i) += ifunCurlHVOLCoeffs(k)*hvolBasisAtRefCoords(k,i);
        diffErrCurl = std::max(diffErrCurl, std::abs(ifunCurlInHVOLSpaceAtRefCoords(i) - ifunCurlAtRefCoords(i)));
      } 

      //Check that the two representations of the div of ifun are consistent
      if(diffErrCurl > pow(7, order-1)*tol) { //heuristic relation on how round-off error depends on order
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << "Curl of HCURL_I" << order << " function does not belong to HVOL_C" << order-1 << "."<<
            "\nIn fact the HVOL interpolation of the function curl is different from the function curl."<<
            "\nThe max norm of the difference at HCURL DOF coordinates is: " <<  diffErrCurl << std::endl;
      }
      
    }
  } catch (std::exception &err) {
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

