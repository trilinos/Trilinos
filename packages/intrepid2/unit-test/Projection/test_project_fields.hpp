// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Test for checking the correctness of function ProjectionTools::projectField.  

    The test considers a uniform and structured mesh of the hypercube [-1,1]^dim, with dim=2 or 3.
    The hypercube is meshed with Hex, Tets, Wedges (in 3d) and Quad and Tri (in 2d).

    A finite element field is defined on the mesh and it is initialized with an analytic function.
    The field is the projected to an higher-order finite element space and then it is projected back to the original space.
    We then check that the original field is the same as the one projected back up to a tight tolerance.

    The tests consider different finite element spaces in HGRAD, HCURL, HDIV and HVOL. 

    \author Created by Mauro Perego
 */

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_ProjectionTools.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid2_HGRAD_WEDGE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_WEDGE_C2_FEM.hpp"
#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid2_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TET_C2_FEM.hpp"
#include "Intrepid2_HCURL_TET_I1_FEM.hpp"
#include "Intrepid2_HDIV_TET_I1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"
#include "Intrepid2_HDIV_QUAD_I1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid2_HCURL_TRI_I1_FEM.hpp"
#include "Intrepid2_HDIV_TRI_I1_FEM.hpp"
#include "Intrepid2_HVOL_C0_FEM.hpp"
#include "struct_mesh_utils.hpp"


#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include <array>
#include <set>
#include <random>
#include <algorithm>

namespace Intrepid2 {

namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception &err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }

template<typename ValueType, typename DeviceType>
int ProjectFields(const bool verbose) {

  using ExecSpaceType = typename DeviceType::execution_space;

  typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
  typedef Kokkos::DynRankView<ordinal_type,DeviceType> DynRankViewInt;

  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  if (verbose)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs,       false);

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  using HostSpaceType = Kokkos::DefaultHostExecutionSpace;

  *outStream << "DeviceSpace::  ";   ExecSpaceType().print_configuration(*outStream, false);
  *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);
  *outStream << "\n";

  ordinal_type errorFlag = 0;

  struct Fun {
    KOKKOS_INLINE_FUNCTION
    ValueType
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const ordinal_type comp=0) {
      ValueType f0 = sin(x*2)*sin(y*2)*sin(z*2)+sin(x*y*z*8);
      ValueType f1 = cos(x*2)*cos(y*2)*cos(z*2);
      ValueType f2 = cos(x*y*z*8);
      //fun = f + a \times x
      switch (comp) {
      case 0:
        return f0;
      case 1:
        return f1;
      case 2:
        return f2;
      default:
        return 0;
      }
    }
  };

  using ots = OrientationTools<DeviceType>;
  using pts = ProjectionTools<DeviceType>;
  using ProjStruct = ProjectionStruct<DeviceType,ValueType>;
  using basisPtrType = BasisPtr<DeviceType,ValueType,ValueType>;
  using CG_NBasis = DerivedNodalBasisFamily<DeviceType,ValueType,ValueType>;

  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|           Projecting fields to higher-dregree FE spaces and back            |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";

  try {

    // ************* Selecte basis for projections ***********

    //vector of basis pairs. we will project from the first basis of each the pair to its second basis, and back
    std::vector<std::pair<basisPtrType,basisPtrType> > basis_pair_set;

    //HEX HGRAD
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_HEX_C1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_HEX(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_HEX_C1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_HEX(2))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_HEX_C2_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_HEX(2))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_HEX_C2_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_HEX(3))));

    //HEX HCURL
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HCURL_HEX_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HCURL_HEX(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HCURL_HEX_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HCURL_HEX(2))));
    
    //HEX HDIV
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HDIV_HEX_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HDIV_HEX(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HDIV_HEX_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HDIV_HEX(2))));

    //HEX HVOL
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HVOL_C0_FEM<DeviceType,ValueType,ValueType>(shards::getCellTopologyData<shards::Hexahedron<8>>())), Teuchos::rcp(new typename  CG_NBasis::HVOL_HEX(1))));



    //TET HGRAD
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_TET_C1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_TET(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_TET_C1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_TET(2))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_TET_C2_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_TET(2))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_TET_C2_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_TET(3))));

    //TET HCURL
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HCURL_TET_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HCURL_TET(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HCURL_TET_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HCURL_TET(2))));

    //TET HDIV
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HDIV_TET_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HDIV_TET(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HDIV_TET_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HDIV_TET(2))));

    //TET HVOL
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HVOL_C0_FEM<DeviceType,ValueType,ValueType>(shards::getCellTopologyData<shards::Tetrahedron<4>>())), Teuchos::rcp(new typename  CG_NBasis::HVOL_TET(1))));


    //WEDGE HGRAD
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_WEDGE_C1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_WEDGE(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_WEDGE_C1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_WEDGE(2))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_WEDGE_C2_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_WEDGE(2))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_WEDGE_C2_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_WEDGE(3))));

    //WEDGE HCURL
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HCURL_WEDGE_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HCURL_WEDGE(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HCURL_WEDGE_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HCURL_WEDGE(2))));

    //WEDGE HDIV
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HDIV_WEDGE_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HDIV_WEDGE(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HDIV_WEDGE_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HDIV_WEDGE(2))));

    //WEDGE HVOL
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HVOL_C0_FEM<DeviceType,ValueType,ValueType>(shards::getCellTopologyData<shards::Wedge<6>>())), Teuchos::rcp(new typename  CG_NBasis::HVOL_WEDGE(1))));
      

    //QUAD HGRAD
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_QUAD_C1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_QUAD(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_QUAD_C1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_QUAD(2))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_QUAD_C2_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_QUAD(2))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_QUAD_C2_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_QUAD(3))));

    //QUAD HCURL
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HCURL_QUAD_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HCURL_QUAD(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HCURL_QUAD_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HCURL_QUAD(2))));

    //QUAD HDIV
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HDIV_QUAD_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HDIV_QUAD(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HDIV_QUAD_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HDIV_QUAD(2))));

    //QUAD HVOL
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HVOL_C0_FEM<DeviceType,ValueType,ValueType>(shards::getCellTopologyData<shards::Quadrilateral<4>>())), Teuchos::rcp(new typename  CG_NBasis::HVOL_QUAD(1))));

    

    //TRI HGRAD
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_TRI_C1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_TRI(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_TRI_C1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_TRI(2))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_TRI_C2_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_TRI(2))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HGRAD_TRI_C2_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HGRAD_TRI(3))));

    //TRI HCURL
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HCURL_TRI_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HCURL_TRI(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HCURL_TRI_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HCURL_TRI(2))));

    //TRI HDIV
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HDIV_TRI_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HDIV_TRI(1))));

    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HDIV_TRI_I1_FEM<DeviceType,ValueType,ValueType>()), Teuchos::rcp(new typename  CG_NBasis::HDIV_TRI(2))));

    //TRI HVOL
    basis_pair_set.push_back(std::make_pair(
      Teuchos::rcp(new Basis_HVOL_C0_FEM<DeviceType,ValueType,ValueType>(shards::getCellTopologyData<shards::Triangle<3>>())), Teuchos::rcp(new typename  CG_NBasis::HVOL_TRI(1))));


    DynRankView nodeCoords;
    DynRankViewInt elemNodes;
    Kokkos::DynRankView<Orientation,DeviceType> elemOrts;
    
    // ******** looping through basis pairs **************
    unsigned meshTopoKey = 0;
    for (auto basis_pair:basis_pair_set) {
      auto srcBasisPtr = basis_pair.first;
      auto dstBasisPtr = basis_pair.second;

      //output the projection we are preforming
      *outStream << "\n" << srcBasisPtr->getName() << " -> " << dstBasisPtr->getName() << " -> " << srcBasisPtr->getName() <<std::endl;
      ordinal_type srcBasisCardinality = srcBasisPtr->getCardinality();

      //Create mesh and orientations if not already computed;
      const auto& cellTopo = srcBasisPtr->getBaseCellTopology();
      if(cellTopo.getKey() != meshTopoKey) {
        ordinal_type NX(2), NY(2), NZ(2);
        createStructMesh(nodeCoords, elemNodes, cellTopo, NX, NY, NZ, false, *outStream);
        elemOrts = Kokkos::DynRankView<Orientation,DeviceType>("elemOrts", elemNodes.extent(0));
        ots::getOrientation(elemOrts, elemNodes, cellTopo);
        meshTopoKey = srcBasisPtr->getBaseCellTopology().getKey();    
      }
      ordinal_type numElems = elemNodes.extent(0);
      ordinal_type numNodesPerElem = elemNodes.extent(1);
      const ordinal_type dim = nodeCoords.extent(1);
      ordinal_type fieldDim = (srcBasisPtr->getFunctionSpace() == Intrepid2::FUNCTION_SPACE_HCURL || srcBasisPtr->getFunctionSpace() == Intrepid2::FUNCTION_SPACE_HDIV) ? dim : 1;


      // ******** compute projection-based interpolation of the function into the src FE space *********************
      DynRankView srcBasisCoeffs("srcBasisCoeffs", numElems, srcBasisPtr->getCardinality());
      {
        ordinal_type srcCubDegree(srcBasisPtr->getDegree());

        ProjStruct projStruct;
        projStruct.createL2ProjectionStruct(srcBasisPtr.get(), srcCubDegree);
        
        auto evaluationPoints = projStruct.getAllEvalPoints();
        ordinal_type numPoints = evaluationPoints.extent(0);

        DynRankView functionAtEvalPoints; 
        
        if(fieldDim == 1)
          functionAtEvalPoints = DynRankView("functionAtEvalPoints", numElems, numPoints);
        else 
          functionAtEvalPoints = DynRankView("functionAtEvalPoints", numElems, numPoints, dim);

        //(x,y,z) coords also in 2D. In 2D, z=0.
        DynRankView physEvalPoints("physEvalPoints", numElems, numPoints, 3); 
        {
          DynRankView linearBasisValuesAtEvalPoint("linearBasisValuesAtEvalPoint", numElems, numNodesPerElem);

          if(cellTopo.getKey() == shards::Hexahedron<8>::key) {
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
            KOKKOS_LAMBDA (const ordinal_type &i) {
              auto basisValuesAtEvalPoint = Kokkos::subview(linearBasisValuesAtEvalPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j){
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_HEX_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalPoints(i,j,d) += nodeCoords(elemNodes(i,k),d)*basisValuesAtEvalPoint(k);
              }
            });
          } else if(cellTopo.getKey() == shards::Tetrahedron<4>::key) {
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
            KOKKOS_LAMBDA (const ordinal_type &i) {
              auto basisValuesAtEvalPoint = Kokkos::subview(linearBasisValuesAtEvalPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j){
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TET_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalPoints(i,j,d) += nodeCoords(elemNodes(i,k),d)*basisValuesAtEvalPoint(k);
              }
            });
          } else if(cellTopo.getKey() == shards::Wedge<6>::key) {
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
            KOKKOS_LAMBDA (const ordinal_type &i) {
              auto basisValuesAtEvalPoint = Kokkos::subview(linearBasisValuesAtEvalPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j){
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_WEDGE_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalPoints(i,j,d) += nodeCoords(elemNodes(i,k),d)*basisValuesAtEvalPoint(k);
              }
            });
          } else if(cellTopo.getKey() == shards::Quadrilateral<4>::key) {
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
            KOKKOS_LAMBDA (const ordinal_type &i) {
              auto basisValuesAtEvalPoint = Kokkos::subview(linearBasisValuesAtEvalPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j){
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_QUAD_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalPoints(i,j,d) += nodeCoords(elemNodes(i,k),d)*basisValuesAtEvalPoint(k);
              }
            });
          } else if(cellTopo.getKey() == shards::Triangle<3>::key) {
            Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
            KOKKOS_LAMBDA (const ordinal_type &i) {
              auto basisValuesAtEvalPoint = Kokkos::subview(linearBasisValuesAtEvalPoint,i,Kokkos::ALL());
              for(ordinal_type j=0; j<numPoints; ++j){
                auto evalPoint = Kokkos::subview(evaluationPoints,j,Kokkos::ALL());
                Impl::Basis_HGRAD_TRI_C1_FEM::template Serial<OPERATOR_VALUE>::getValues(basisValuesAtEvalPoint, evalPoint);
                for(ordinal_type k=0; k<numNodesPerElem; ++k)
                  for(ordinal_type d=0; d<dim; ++d)
                    physEvalPoints(i,j,d) += nodeCoords(elemNodes(i,k),d)*basisValuesAtEvalPoint(k);
              }
            });
          } else {
            INTREPID2_TEST_FOR_EXCEPTION
              (true, std::runtime_error, "Intrepid2::Test::ProjectFields: Topology not supported");
          }
          ExecSpaceType().fence();
        }

        //transform the target function to the reference element (inverse of pullback operator)
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
        KOKKOS_LAMBDA (const ordinal_type &ic) {
          Fun fun;
          for(ordinal_type i=0;i<numPoints;i++)
            for(ordinal_type d=0;d<fieldDim;d++)
              functionAtEvalPoints.access(ic,i,d) = fun(physEvalPoints(ic,i,0), physEvalPoints(ic,i,1), physEvalPoints(ic,i,2), d);
        });
        ExecSpaceType().fence();

        pts::getL2BasisCoeffs(srcBasisCoeffs,
              functionAtEvalPoints,
              elemOrts,
              srcBasisPtr.get(),
              &projStruct);
      }


      // ******** project fields into higher-order finite element space and back **************

      //project from source to destination basis
      DynRankView dstBasisCoeffs("dstBasisCoeffs", numElems, dstBasisPtr->getCardinality());

      //Note, no need to create subviews. This is only to test that projectField works with subviews.
      auto range = std::make_pair(0,numElems);
      pts::projectField(Kokkos::subview(dstBasisCoeffs, range, Kokkos::ALL()), dstBasisPtr.get(), Kokkos::subview(srcBasisCoeffs, range, Kokkos::ALL()), srcBasisPtr.get(), Kokkos::subview(elemOrts, range));

      //project back
      DynRankView srcBasisCoeffs2("srcBasisCoeffs2", numElems, srcBasisPtr->getCardinality());
      pts::projectField(srcBasisCoeffs2, srcBasisPtr.get(), dstBasisCoeffs, dstBasisPtr.get(), elemOrts);


      // ******** compare results **************

      //compute error in l1 norm of original coefficients with those after projecting to higher FE space and back
      ValueType norm1(0);
      Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpaceType>(0,numElems),
      KOKKOS_LAMBDA (const ordinal_type &i, double &norm1Update) {
        for(ordinal_type j=0; j<srcBasisCardinality; ++j)
          norm1Update += std::abs(srcBasisCoeffs2(i,j) - srcBasisCoeffs(i,j));
      }, norm1);

      ExecSpaceType().fence();

      // this should be close to machine precision. Error out if the error is too large
      if(norm1 > 1000.0 * epsilon()){
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << "Error too large: " << norm1;
        *outStream << std::endl;
      }
    }
  } catch (std::exception &err) {
    std::cout << " Exception\n";
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED = " << errorFlag << "\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  return errorFlag;
}
}
}

