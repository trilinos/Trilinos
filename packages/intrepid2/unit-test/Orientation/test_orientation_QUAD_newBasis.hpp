// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file test_orientation_QUAD_newBasis.hpp
    \brief Test for checking orientation tools for the Nodal and Hierarchical Derived basis on Quadrilaterals

    The test considers two quadrilaterals in the physical space sharing a common edge. 
    In order to test significant configurations, we consider 4 mappings of the reference quadrilateral
    to the first (physical) quadrilateral, so that the common edge is mapped from all the 4 edges
    of the reference quadrilateral.
    Then, for each of the mappings, the global ids of the vertices of the common edge are permuted.
    This gives a total of 8 combinations

    The test considers HGRAD, HCURL and HDIV Lagrangian and Hierarchical basis functions, and for each of them:
    1. Computes the oriented basis.
    2. computes the basis coefficients, separately for each quad, for functions belonging to the H-space spanned by the basis.
    3. checks that the dofs shared between the two quads are equivalent (this ensures that the orientation works correctly)
    4. checks that the functions are indeed exactly reproduced

    \author Created by Mauro Perego
 */

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"

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

template<typename edgeT>
bool areEdgesEqual(const edgeT& edge0, const edgeT& edge1)
{
  return ((edge0[0] == edge1[0]) && (edge0[1] == edge1[1])) || ((edge0[0] == edge1[1]) && (edge0[1] == edge1[0]));
}

template<typename ValueType, typename DeviceType>
int OrientationQuadNewBasis(const bool verbose) {

  typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
  typedef Kokkos::DynRankView<ordinal_type,DeviceType> DynRankViewInt;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

  static Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  if (verbose)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs,       false);

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  int errorFlag = 0;
  const ValueType tol = tolerence();

  struct Fun {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y) {
      return x*(x+1)*(y-2);//*(z+3)*(x + 2*y +5*z+ 1.0/3.0);
    }
  };

  struct FunDiv {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const int comp=0) {
      ValueType a = 2*x*y+x*x;
      ValueType f0 = 5+y+x*x;
      ValueType f1 = -7-2*x+y*y;
      //fun = f + a x
      switch (comp) {
      case 0:
        return f0 + a*x;
      case 1:
        return f1 + a*y;
      default:
        return 0;
      }
    }
  };

  struct FunCurl {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const int comp=0) {
      ValueType a2 = y-7+y*y;
      ValueType f0 = 2+x+x*y;
      ValueType f1 = 3-3*x;
      //fun = f + a \times x
      switch (comp) {
      case 0:
        return f0 - a2*y;
      case 1:
        return f1 + a2*x;
      default:
        return 0;
      }
    }
  };


  class BasisFunctionsSystem{
  public:
    BasisFunctionsSystem(const ordinal_type basisCardinality_, const ordinal_type numRefCoords_, const ordinal_type  dim_) :
      basisCardinality(basisCardinality_),
      numRefCoords(numRefCoords_),
      dim(dim_),
      work("lapack_work", basisCardinality+dim*numRefCoords, 1) ,
      cellMassMat("basisMat", dim*numRefCoords,basisCardinality),
      cellRhsMat("rhsMat",dim*numRefCoords, 1) {
    };

    std::vector<int> computeBasisCoeffs(DynRankView basisCoeffs, ordinal_type& errorFlag, const DynRankView transformedBasisValuesAtRefCoordsOriented, const DynRankView funAtPhysRefCoords) {

      ordinal_type numCells = basisCoeffs.extent(0);
      std::vector<int> info(numCells);
      for(ordinal_type ic=0; ic<numCells; ++ic) {
        for(ordinal_type i=0; i<numRefCoords; ++i) {
          if (dim==1) {
            for(ordinal_type j=0; j<basisCardinality; ++j) {
              cellMassMat(i,j) = transformedBasisValuesAtRefCoordsOriented(ic,j,i);
            }
            cellRhsMat(i,0) = funAtPhysRefCoords(ic,i);
          } else
            for (ordinal_type id=0; id<dim; ++id) {
              for(ordinal_type j=0; j<basisCardinality; ++j) {
                cellMassMat(i+id*numRefCoords,j) = transformedBasisValuesAtRefCoordsOriented(ic,j,i,id);
              }
              cellRhsMat(i+id*numRefCoords,0) = funAtPhysRefCoords(ic,i,id);
            }
        }

        lapack.GELS('N', dim*numRefCoords, basisCardinality,1,
            cellMassMat.data(),
            cellMassMat.stride_1(),
            cellRhsMat.data(),
            cellRhsMat.stride_1(),
            work.data(),
            basisCardinality+dim*numRefCoords,
            &info[ic]);

        for(ordinal_type i=0; i<basisCardinality; ++i){
          basisCoeffs(ic,i) = cellRhsMat(i,0);
        }
      }
      return info;
    }
  private:
    Teuchos::LAPACK<ordinal_type,ValueType> lapack;
    ordinal_type basisCardinality, numRefCoords, dim;
    Kokkos::View<ValueType**,Kokkos::LayoutLeft,Kokkos::HostSpace> work;
    Kokkos::View<ValueType**,Kokkos::LayoutLeft,Kokkos::HostSpace> cellMassMat;
    Kokkos::View<ValueType**,Kokkos::LayoutLeft,Kokkos::HostSpace> cellRhsMat;
  };


  typedef std::array<ordinal_type,2> edgeType;
  typedef CellTools<DeviceType> ct;
  typedef OrientationTools<DeviceType> ots;
  typedef RealSpaceTools<DeviceType> rst;
  typedef FunctionSpaceTools<DeviceType> fst;

  using  basisType = Basis<DeviceType,ValueType,ValueType>;

  constexpr ordinal_type dim = 2;
  constexpr ordinal_type numCells = 2;
  constexpr ordinal_type numElemVertexes = 4;
  constexpr ordinal_type numTotalVertexes = 6;
  constexpr ordinal_type numSharedVertexes = 2;
  constexpr ordinal_type numSharedEdges = 1;

  ValueType  vertices_orig[numTotalVertexes][dim] = {{-1,-1},{1,-1},{1,1},{-1,1}, {-1,2},{1,2}};
  ordinal_type quads_orig[numCells][numElemVertexes] = {{0,1,2,3},{3,2,5,4}};
  static edgeType common_edge = {{2,3}};
  ordinal_type quads_rotated[numCells][numElemVertexes];

  static ordinal_type shared_vertexes[numCells][numSharedVertexes];
  static ordinal_type edgeIndexes[numCells][numSharedEdges];
  static ordinal_type faceIndex[numCells];


  class TestResults
  {
  private:
    const DynRankView basisCoeffs, transformedBasisValuesAtRefCoordsOriented, funAtPhysRefCoords;
    const basisType* basis;


  public:

    TestResults(const DynRankView basisCoeffs_,
        const DynRankView transformedBasisValuesAtRefCoordsOriented_,
        const DynRankView funAtPhysRefCoords_,
        const basisType* basis_) :
          basisCoeffs(basisCoeffs_),
          transformedBasisValuesAtRefCoordsOriented(transformedBasisValuesAtRefCoordsOriented_),
          funAtPhysRefCoords(funAtPhysRefCoords_),
          basis(basis_){}

    //check that fun values are consistent at the common vertexes
    void test(ordinal_type& errorFlag, ValueType tol){

      auto numVertexDOFs = basis->getDofCount(0,0);
      if(numVertexDOFs >0) {
        bool areDifferent(false);


        for(ordinal_type j=0;j<numSharedVertexes && !areDifferent;j++)
          areDifferent = std::abs(basisCoeffs(0,basis->getDofOrdinal(0,shared_vertexes[0][j],0))
              - basisCoeffs(1,basis->getDofOrdinal(0,shared_vertexes[1][j],0))) > 10*tol;

        if(areDifferent) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "Function  DOFs on shared vertexes computed using Cell 0 basis functions are not consistent with those computed using Cell 1 bssis functions\n";
          *outStream << "Function DOFs for Cell 0 are:";
          for(ordinal_type j=0;j<numSharedVertexes;j++)
            *outStream << " " << basisCoeffs(0,basis->getDofOrdinal(0,shared_vertexes[0][j],0));
          *outStream << "\nFunction DOFs for Cell 1 are:";
          for(ordinal_type j=0;j<numSharedVertexes;j++)
            *outStream << " " << basisCoeffs(1,basis->getDofOrdinal(0,shared_vertexes[1][j],0));
          *outStream << std::endl;
        }

      }


      //check that fun values are consistent on shared edges dofs
      auto numEdgeDOFs = basis->getDofCount(1,0);
      if(numEdgeDOFs>0)
      {

        bool areDifferent(false);
        for(std::size_t iEdge=0;iEdge<numSharedEdges;iEdge++) {
        for(ordinal_type j=0;j<numEdgeDOFs && !areDifferent;j++) {
          areDifferent = std::abs(basisCoeffs(0,basis->getDofOrdinal(1,edgeIndexes[0][iEdge],j))
              - basisCoeffs(1,basis->getDofOrdinal(1,edgeIndexes[1][iEdge],j))) > 10*tol;
        }
        if(areDifferent) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "Function DOFs on shared edge " << iEdge << " computed using Cell 0 basis functions are not consistent with those computed using Cell 1 basis functions\n";
          *outStream << "Function DOFs for Cell 0 are:";
          for(ordinal_type j=0;j<numEdgeDOFs;j++)
            *outStream << " " << basisCoeffs(0,basis->getDofOrdinal(1,edgeIndexes[0][iEdge],j));
          *outStream << "\nFunction DOFs for Cell 1 are:";
          for(ordinal_type j=0;j<numEdgeDOFs;j++)
            *outStream << " " << basisCoeffs(1,basis->getDofOrdinal(1,edgeIndexes[1][iEdge],j));
          *outStream << std::endl;
        }
        }
      }


      //check that fun values are consistent on common face dofs
      auto numFaceDOFs = basis->getDofCount(2,0);
      if(numFaceDOFs > 0 && dim>2)
      {
        bool areDifferent(false);
        for(ordinal_type j=0;j<numFaceDOFs && !areDifferent;j++) {
          areDifferent = std::abs(basisCoeffs(0,basis->getDofOrdinal(2,faceIndex[0],j))
              - basisCoeffs(1,basis->getDofOrdinal(2,faceIndex[1],j))) > 10*tol;
        }

        if(areDifferent) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "Function DOFs on common face computed using Hex 0 basis functions are not consistent with those computed using Hex 1\n";
          *outStream << "Function DOFs for Hex 0 are:";
          for(ordinal_type j=0;j<numFaceDOFs;j++)
            *outStream << " " << basisCoeffs(0,basis->getDofOrdinal(2,faceIndex[0],j));
          *outStream << "\nFunction DOFs for Hex 1 are:";
          for(ordinal_type j=0;j<numFaceDOFs;j++)
            *outStream << " " << basisCoeffs(1,basis->getDofOrdinal(2,faceIndex[1],j));
          *outStream << std::endl;
        }
      }

      ordinal_type numRefCoords = funAtPhysRefCoords.extent(1);
      ordinal_type basisCardinality = basisCoeffs.extent(1);
      ordinal_type basis_dim = (transformedBasisValuesAtRefCoordsOriented.rank()==3) ? 1 : dim;

      DynRankView ConstructWithLabel(funAtRefCoordsOriented, numCells, numRefCoords, basis_dim);
      for(ordinal_type ic=0; ic<numCells; ++ic) {
        ValueType error=0;
        for(ordinal_type j=0; j<numRefCoords; ++j) {
          if (basis_dim==1) {
            for(ordinal_type k=0; k<basisCardinality; ++k)
              funAtRefCoordsOriented(ic,j,0) += basisCoeffs(ic,k)*transformedBasisValuesAtRefCoordsOriented(ic,k,j);
            error = std::max(std::abs( funAtPhysRefCoords(ic,j) - funAtRefCoordsOriented(ic,j,0)), error);
          } else {
            for(ordinal_type d=0; d<basis_dim; ++d) {
              for(ordinal_type k=0; k<basisCardinality; ++k)
                funAtRefCoordsOriented(ic,j,d) += basisCoeffs(ic,k)*transformedBasisValuesAtRefCoordsOriented(ic,k,j,d);
              error = std::max(std::abs( funAtPhysRefCoords(ic,j,d) - funAtRefCoordsOriented(ic,j,d)), error);
            }
          }
        }

        if(error>100*tol) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "Function values at reference points differ from those computed using basis functions on Cell " << ic << "\n";
          *outStream << "Function values at reference points are:\n";
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            if(basis_dim==1)
              *outStream << " (" << funAtPhysRefCoords(ic,j);
            else
              *outStream << " (" << funAtPhysRefCoords(ic,j,0);
            for(ordinal_type d=1; d<basis_dim; d++)
              *outStream <<  ", " << funAtPhysRefCoords(ic,j,d);
            *outStream  << ")";
          }
          *outStream << "\nFunction values at reference points computed using basis functions are\n";
          for(ordinal_type j=0; j<numRefCoords; ++j) {
            *outStream << " (" << funAtRefCoordsOriented(ic,j,0);
            for(ordinal_type d=1; d<basis_dim; d++)
              *outStream <<  ", " << funAtRefCoordsOriented(ic,j,d);
            *outStream  << ")";
          }
          *outStream << std::endl;
        }
      }
    }
  };


  try {

    const ordinal_type order = 3;
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5};

    do {
      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type quads[numCells][numElemVertexes];
      std::copy(&quads_orig[0][0], &quads_orig[0][0]+numCells*numElemVertexes, &quads_rotated[0][0]);
      for (ordinal_type shift=0; shift<4; ++shift) {
        std::rotate_copy(&quads_orig[0][0], &quads_orig[0][0]+shift, &quads_orig[0][0]+4, &quads_rotated[0][0]);
        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            quads[i][j] = reorder[quads_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Quad 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << quads[0][j] << " ";
        *outStream << "] and Quad 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << quads[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, quad.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<quad.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[quads[i][j]][k];

        //computing common and edge
        {
              ordinal_type vertex;
              edgeType edge{};
              for(ordinal_type i=0; i<numCells; ++i) {
                for (std::size_t iv=0; iv<quad.getNodeCount(); ++iv) {
                  vertex = quads_rotated[i][quad.getNodeMap(0,iv,0)];
                  for (std::size_t isv=0; isv<common_edge.size(); ++isv)
                    if(common_edge[isv] == vertex)
                      shared_vertexes[i][isv] = iv;
                }
                //compute common edge
                for (std::size_t ie=0; ie<quad.getSideCount(); ++ie) {
                  for (std::size_t k=0; k<quad.getNodeCount(1,ie); ++k)
                    edge[k]= quads_rotated[i][quad.getNodeMap(1,ie,k)];

                  if(areEdgesEqual(edge,common_edge)) edgeIndexes[i][0]=ie;
                }
              }
        }

        using CG_NBasis = NodalBasisFamily<DeviceType,ValueType,ValueType>;
        using CG_DNBasis = DerivedNodalBasisFamily<DeviceType,ValueType,ValueType>;
        using CG_HBasis = HierarchicalBasisFamily<DeviceType,ValueType,ValueType>;
        std::vector<basisType*> basis_set;


        //compute reference points
        typename CG_NBasis::HGRAD_QUAD warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&quads[0][0], numCells, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, quad);

        DynRankView ConstructWithLabel(physRefCoords, numCells, numRefCoords, dim);
        {
          Basis_HGRAD_QUAD_C1_FEM<DeviceType,ValueType,ValueType> quadLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(quadLinearBasisValuesAtRefCoords, quad.getNodeCount(), numRefCoords);
          quadLinearBasis.getValues(quadLinearBasisValuesAtRefCoords, refPoints);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type d=0; d<dim; ++d) {
              for(ordinal_type j=0; j<numRefCoords; ++j)
                for(std::size_t k=0; k<quad.getNodeCount(); ++k)
                  physRefCoords(i,j,d) += vertices[quads[i][k]][d]*quadLinearBasisValuesAtRefCoords(k,j);
            }
        }

        //HGRAD BASIS
        {
          Fun fun;
          DynRankView ConstructWithLabel(funAtPhysRefCoords, numCells, numRefCoords);
          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j)
              funAtPhysRefCoords(i,j) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1));
          }

          basis_set.push_back(new typename  CG_NBasis::HGRAD_QUAD(order));
          basis_set.push_back(new typename  CG_DNBasis::HGRAD_QUAD(order));
          basis_set.push_back(new typename  CG_HBasis::HGRAD_QUAD(order));

          for (auto basisPtr:basis_set) {
            auto& basis = *basisPtr;
            auto name = basis.getName();
            *outStream << " " << name << std::endl;
            ordinal_type basisCardinality = basis.getCardinality();

            //check that fun values at reference points coincide with those computed using basis functions
            DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords);
            DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords);

            DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords);
            basis.getValues(basisValuesAtRefCoords, refPoints);

            // modify basis values to account for orientations
            ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
                basisValuesAtRefCoords,
                elemOrts,
                &basis);

            // transform basis values
            deep_copy(transformedBasisValuesAtRefCoordsOriented,
                basisValuesAtRefCoordsOriented);

            DynRankView ConstructWithLabel(basisCoeffs, numCells, basisCardinality);

            BasisFunctionsSystem  basisFunctionsSystem(basisCardinality, numRefCoords, 1);
            auto info = basisFunctionsSystem.computeBasisCoeffs(basisCoeffs, errorFlag, transformedBasisValuesAtRefCoordsOriented, funAtPhysRefCoords);

            for (int ic =0; ic < numCells; ++ic) {
              if(info[ic] != 0) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "LAPACK error flag for cell " << ic << " is: " << info[ic] << std::endl;
              }
            }

            TestResults testResults(basisCoeffs, transformedBasisValuesAtRefCoordsOriented, funAtPhysRefCoords, basisPtr);
            testResults.test(errorFlag,tol);
            delete basisPtr;
          }
        }


        //HCURL Case
        {
          FunCurl fun;
          DynRankView ConstructWithLabel(funAtPhysRefCoords, numCells, numRefCoords, dim);
          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<dim; ++k)
                funAtPhysRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), k);
            }
          }

          basis_set.clear();
          basis_set.push_back(new typename  CG_NBasis::HCURL_QUAD(order));
          basis_set.push_back(new typename  CG_DNBasis::HCURL_QUAD(order));
          basis_set.push_back(new typename  CG_HBasis::HCURL_QUAD(order));

          for (auto basisPtr:basis_set) {
            auto& basis = *basisPtr;
            auto name = basis.getName();
            *outStream << " " << name << std::endl;

            ordinal_type basisCardinality = basis.getCardinality();
            DynRankView ConstructWithLabel(basisCoeffs, numCells, basisCardinality);

            //check that fun values at reference points coincide with those computed using basis functions
            DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
            DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
            DynRankView basisValuesAtRefCoordsCells("inValues", numCells, basisCardinality, numRefCoords, dim);


            DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords, dim);
            basis.getValues(basisValuesAtRefCoords, refPoints);
            rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

            // modify basis values to account for orientations
            ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
                basisValuesAtRefCoordsCells,
                elemOrts,
                &basis);

            // transform basis values
            DynRankView ConstructWithLabel(jacobianAtRefCoords, numCells, numRefCoords, dim, dim);
            DynRankView ConstructWithLabel(jacobianAtRefCoords_inv, numCells, numRefCoords, dim, dim);
            ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, quad);
            ct::setJacobianInv (jacobianAtRefCoords_inv, jacobianAtRefCoords);
            fst::HCURLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
                jacobianAtRefCoords_inv,
                basisValuesAtRefCoordsOriented);

            BasisFunctionsSystem  basisFunctionsSystem(basisCardinality, numRefCoords, dim);
            auto info = basisFunctionsSystem.computeBasisCoeffs(basisCoeffs, errorFlag, transformedBasisValuesAtRefCoordsOriented, funAtPhysRefCoords);

            for (int ic =0; ic < numCells; ++ic) {
              if(info[ic] != 0) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "LAPACK error flag for cell " << ic << " is: " << info[ic] << std::endl;
              }
            }

            TestResults testResults(basisCoeffs, transformedBasisValuesAtRefCoordsOriented, funAtPhysRefCoords, basisPtr);
            testResults.test(errorFlag,tol);
            delete basisPtr;
          }
        }


        //HDIV Case
        {
          FunDiv fun;
          DynRankView ConstructWithLabel(funAtPhysRefCoords, numCells, numRefCoords, dim);
          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j) {
              for(ordinal_type k=0; k<dim; ++k)
                funAtPhysRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), k);
            }
          }
          basis_set.clear();
          basis_set.push_back(new typename  CG_NBasis::HDIV_QUAD(order));
          basis_set.push_back(new typename  CG_DNBasis::HDIV_QUAD(order));
          basis_set.push_back(new typename  CG_HBasis::HDIV_QUAD(order));

          for (auto basisPtr:basis_set) {
            auto& basis = *basisPtr;
            auto name = basis.getName();
            *outStream << " " << name << std::endl;
            ordinal_type basisCardinality = basis.getCardinality();

            //check that fun values at reference points coincide with those computed using basis functions
            DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
            DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
            DynRankView basisValuesAtRefCoordsCells("inValues", numCells, basisCardinality, numRefCoords, dim);

            DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords, dim);
            basis.getValues(basisValuesAtRefCoords, refPoints);
            rst::clone(basisValuesAtRefCoordsCells,basisValuesAtRefCoords);

            // modify basis values to account for orientations
            ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
                basisValuesAtRefCoordsCells,
                elemOrts,
                &basis);

            // transform basis values
            DynRankView ConstructWithLabel(jacobianAtRefCoords, numCells, numRefCoords, dim, dim);
            DynRankView ConstructWithLabel(jacobianAtRefCoords_det, numCells, numRefCoords);
            ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, quad);
            ct::setJacobianDet (jacobianAtRefCoords_det, jacobianAtRefCoords);
            fst::HDIVtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
                jacobianAtRefCoords,
                jacobianAtRefCoords_det,
                basisValuesAtRefCoordsOriented);

            DynRankView ConstructWithLabel(basisCoeffs, numCells, basisCardinality);


            BasisFunctionsSystem  basisFunctionsSystem(basisCardinality, numRefCoords, dim);
            auto info = basisFunctionsSystem.computeBasisCoeffs(basisCoeffs, errorFlag, transformedBasisValuesAtRefCoordsOriented, funAtPhysRefCoords);

            for (int ic =0; ic < numCells; ++ic) {
              if(info[ic] != 0) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "LAPACK error flag for cell " << ic << " is: " << info[ic] << std::endl;
              }
            }

            TestResults testResults(basisCoeffs, transformedBasisValuesAtRefCoordsOriented, funAtPhysRefCoords, basisPtr);
            testResults.test(errorFlag,tol);
            delete basisPtr;
          }
        }
      } //rotation of first cell vertices
    } while(std::next_permutation(&reorder[0]+2, &reorder[0]+4)); //reorder vertices of common edge

  } catch (std::exception &err) {
    std::cout << " Exeption\n";
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

