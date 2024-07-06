// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file test_orientation_HEX_newBasis.hpp
    \brief  Test for checking orientation tools for the Nodal and Hierarchical Derived basis on Hexahedra

    The test considers two hexahedra in the physical space sharing a common face. 
    In order to test significant configurations, we consider 6 mappings of the reference hexahedron
    to the first (physical) hexahedron, so that the common face is mapped from all the 6 faces
    of the reference hexahedron.
    Then, for each of the mappings, the global ids of the vertices of the common face are permuted.
    This gives a total of 144 combinations

    The test considers HGRAD, HCURL and HDIV Lagrangian and Hierarchical basis functions, and for each of them:
    1. Computes the oriented basis.
    2. computes the basis coefficients, separately for each hexa, for functions belonging to the H-space spanned by the basis.
    3. checks that the dofs shared between the two hexas are equivalent (this ensures that the orientation works correctly)
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

template<typename ValueType, typename DeviceType>
int OrientationHexNewBasis(const bool verbose) {

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
    operator()(const ValueType& x, const ValueType& y, const ValueType& z) {
      return (x+1)*(y-2);//*(z+3)*(x + 2*y +5*z+ 1.0/3.0);
    }
  };

  struct FunDiv {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      ValueType a = 2*x*y+x*x;
      ValueType f0 = 5+y+x*x+z*z;
      ValueType f1 = -7-2*z+x+y*y+z*z;
      ValueType f2 = 0.5+z*z+x*x;
      //fun = f + a x
      switch (comp) {
      case 0:
        return f0 + a*x;
      case 1:
        return f1 + a*y;
      case 2:
        return f2 + a*z;
      default:
        return 0;
      }
    }
  };

  struct FunCurl {
    ValueType
    KOKKOS_INLINE_FUNCTION
    operator()(const ValueType& x, const ValueType& y, const ValueType& z, const int comp=0) {
      ValueType a0 = y-7+z*z;
      ValueType a1 = 2*z-1+z*x;
      ValueType a2 = z-2+x*x;
      ValueType f0 = 2+x+z+x*y;
      ValueType f1 = 3-3*z;
      ValueType f2 = -5+x;
      //fun = f + a \times x
      switch (comp) {
      case 0:
        return f0 + (a1*z-a2*y);//2*x+y-z + (x+2*(y+z);
      case 1:
        return f1 + (a2*x-a0*z);//y+2*(z+x);
      case 2:
        return f2 + (a0*y-a1*x);//z+2*(x+y);
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
  typedef std::array<ordinal_type,4> faceType;
  typedef CellTools<DeviceType> ct;
  typedef OrientationTools<DeviceType> ots;
  typedef RealSpaceTools<DeviceType> rst;
  typedef FunctionSpaceTools<DeviceType> fst;

  using  basisType = Basis<DeviceType,ValueType,ValueType>;

  constexpr ordinal_type dim = 3;
  constexpr ordinal_type numCells = 2;
  constexpr ordinal_type numElemVertexes = 8;
  constexpr ordinal_type numTotalVertexes = 12;
  constexpr ordinal_type numSharedVertexes = 4;
  constexpr ordinal_type numSharedEdges = 4;


  ValueType  vertices_orig[numTotalVertexes][dim] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}, {-1,-1,2},{1,-1,2},{1,1,2},{-1,1,2}};
  ordinal_type hexas_orig[numCells][numElemVertexes] = {{0,1,2,3,4,5,6,7},{4,5,6,7,8,9,10,11}};  //common face is {4,5,6,7}
  faceType common_face = {{4,5,6,7}};
  faceType faceLeft = {{0, 3, 7, 4}};
  faceType faceRight = {{1, 2, 6, 5}};
  faceType faceFront = {{0, 4, 5, 1}};
  faceType faceBack = {{2, 3, 7, 6}};
  ordinal_type hexas_rotated[numCells][numElemVertexes];
  faceType faceLeftOriented, faceRightOriented, faceBackOriented, faceFrontOriented;

  static std::set<edgeType> common_edges;
  common_edges.insert(edgeType({{4,5}})); common_edges.insert(edgeType({{5,6}})); common_edges.insert(edgeType({{6,7}})); common_edges.insert(edgeType({{4,7}}));

  static ordinal_type shared_vertexes[numCells][numSharedVertexes];
  static ordinal_type edgeIndexes[numCells][numSharedEdges];
  static ordinal_type faceIndex[numCells];


  class TestResults
  {
  private:
    const DynRankView basisCoeffs, transformedBasisValuesAtRefCoords, transformedBasisValuesAtRefCoordsOriented, funAtPhysRefCoords;
    const Kokkos::DynRankView<Orientation,DeviceType> elemOrts;
    const basisType* basis;


  public:

    TestResults(const DynRankView basisCoeffs_,
                const DynRankView transformedBasisValuesAtRefCoords_,
                const DynRankView transformedBasisValuesAtRefCoordsOriented_,
                const DynRankView funAtPhysRefCoords_,
                const Kokkos::DynRankView<Orientation,DeviceType> elemOrts_,
                const basisType* basis_) :
          basisCoeffs(basisCoeffs_),
          transformedBasisValuesAtRefCoords(transformedBasisValuesAtRefCoords_),
          transformedBasisValuesAtRefCoordsOriented(transformedBasisValuesAtRefCoordsOriented_),
          funAtPhysRefCoords(funAtPhysRefCoords_),
          elemOrts(elemOrts_),
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
              - basisCoeffs(1,basis->getDofOrdinal(2,faceIndex[1],j))) > 100*tol;
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
      
      // check that global (oriented) basis coefficients contracted with the transformed oriented basis values agrees with the transpose-oriented basis coefficients contracted with the unoriented transformed basis values
      DynRankView ConstructWithLabel(localBasisCoeffs, numCells, basisCardinality);
      OrientationTools<DeviceType>::modifyBasisByOrientationTranspose(localBasisCoeffs,
                                                                      basisCoeffs,
                                                                      elemOrts,
                                                                      basis);
      
      DynRankView ConstructWithLabel(basisCoeffsModified, numCells, basisCardinality);
      OrientationTools<DeviceType>::modifyBasisByOrientation(basisCoeffsModified,
                                                             basisCoeffs,
                                                             elemOrts,
                                                             basis);
      
      
      for(ordinal_type ic=0; ic<numCells; ++ic) {
        ValueType error=0;
        for(ordinal_type j=0; j<numRefCoords; ++j) {
          for (int d=0; d<transformedBasisValuesAtRefCoords.extent_int(3); d++)
          {
            ValueType globalValue_d = 0; // oriented/global
            ValueType localValue_d  = 0; // unoriented/local
            for(ordinal_type k=0; k<basisCardinality; ++k)
            {
              auto orientedBasisValue = transformedBasisValuesAtRefCoordsOriented(ic,k,j,d);
              auto basisValue         = transformedBasisValuesAtRefCoords        (ic,k,j,d);
              
              globalValue_d += basisCoeffs(ic,k)      * orientedBasisValue;
              localValue_d  += localBasisCoeffs(ic,k) * basisValue;
            }
            error = std::max(std::abs( globalValue_d - localValue_d), error);
          }
        }

        if(error>100*tol) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "global values do not agree with local on cell " << ic << "\n";
        }
      }
    }
  };


  try {

    const ordinal_type order = 3;
    ordinal_type reorder[numTotalVertexes] = {0,1,2,3,4,5,6,7,8,9,10,11};

    do {
      ordinal_type orderback[numTotalVertexes];
      for(ordinal_type i=0;i<numTotalVertexes;++i) {
        orderback[reorder[i]]=i;
      }
      ValueType vertices[numTotalVertexes][dim];
      ordinal_type hexas[numCells][numElemVertexes];
      std::copy(&hexas_orig[0][0], &hexas_orig[0][0]+numCells*numElemVertexes, &hexas_rotated[0][0]);
      for (ordinal_type shift=0; shift<6; ++shift) {
        if(shift <4){
          std::rotate_copy(faceLeft.begin(), faceLeft.begin()+shift, faceLeft.end(), faceLeftOriented.begin());
          std::rotate_copy(faceRight.begin(), faceRight.begin()+shift, faceRight.end(), faceRightOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            hexas_rotated[0][faceLeft[ii]] = hexas_orig[0][faceLeftOriented[ii]];
            hexas_rotated[0][faceRight[ii]] = hexas_orig[0][faceRightOriented[ii]];
          }
        } else {
          ordinal_type iirot = (shift==4) ? 1 : 3;
          std::rotate_copy(faceFront.begin(), faceFront.begin()+iirot, faceFront.end(), faceFrontOriented.begin());
          std::rotate_copy(faceBack.begin(), faceBack.begin()+iirot, faceBack.end(), faceBackOriented.begin());
          for(ordinal_type ii=0; ii<4; ii++) {
            hexas_rotated[0][faceFront[ii]] = hexas_orig[0][faceFrontOriented[ii]];
            hexas_rotated[0][faceBack[ii]] = hexas_orig[0][faceBackOriented[ii]];
          }
        }

        for(ordinal_type i=0; i<numCells;++i)
          for(ordinal_type j=0; j<numElemVertexes;++j)
            hexas[i][j] = reorder[hexas_rotated[i][j]];

        for(ordinal_type i=0; i<numTotalVertexes;++i)
          for(ordinal_type d=0; d<dim;++d)
            vertices[i][d] = vertices_orig[orderback[i]][d];

        *outStream <<  "Considering Hex 0: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << hexas[0][j] << " ";
        *outStream << "] and Hex 1: [ ";
        for(ordinal_type j=0; j<numElemVertexes;++j)
          *outStream << hexas[1][j] << " ";
        *outStream << "]\n";

        shards::CellTopology hexa(shards::getCellTopologyData<shards::Hexahedron<8> >());
        shards::CellTopology quad(shards::getCellTopologyData<shards::Quadrilateral<4> >());
        shards::CellTopology line(shards::getCellTopologyData<shards::Line<2> >());

        //computing vertices coords
        DynRankView ConstructWithLabel(physVertexes, numCells, hexa.getNodeCount(), dim);
        for(ordinal_type i=0; i<numCells; ++i)
          for(std::size_t j=0; j<hexa.getNodeCount(); ++j)
            for(ordinal_type k=0; k<dim; ++k)
              physVertexes(i,j,k) = vertices[hexas[i][j]][k];



        //computing common face and edges


        {
          faceType face={};
          edgeType edge={};
          //bool faceOrientation[numCells][4];
          for(ordinal_type i=0; i<numCells; ++i) {
            for (std::size_t iv=0; iv<hexa.getNodeCount(); ++iv) {
              auto vertex = hexas_rotated[i][hexa.getNodeMap(0,iv,0)];
              for (std::size_t isv=0; isv<common_face.size(); ++isv)
                if(common_face[isv] == vertex)
                  shared_vertexes[i][isv] = iv;
            }
            //compute faces' tangents
            for (std::size_t is=0; is<hexa.getSideCount(); ++is) {
              for (std::size_t k=0; k<hexa.getNodeCount(2,is); ++k)
                face[k]= hexas_rotated[i][hexa.getNodeMap(2,is,k)];

              //rotate and flip
              auto minElPtr= std::min_element(face.begin(), face.end());
              std::rotate(face.begin(),minElPtr,face.end());
              if(face[3]<face[1]) {auto tmp=face[1]; face[1]=face[3]; face[3]=tmp;}

              if(face == common_face) faceIndex[i]=is;
            }
            //compute edges' tangents
            for (std::size_t ie=0; ie<hexa.getEdgeCount(); ++ie) {
              for (std::size_t k=0; k<hexa.getNodeCount(1,ie); ++k)
                edge[k]= hexas_rotated[i][hexa.getNodeMap(1,ie,k)];
              std::sort(edge.begin(),edge.end());
              auto it=common_edges.find(edge);
              if(it !=common_edges.end()){
                auto edge_lid = std::distance(common_edges.begin(),it);
                edgeIndexes[i][edge_lid]=ie;
              }
            }
          }
        }

        using CG_NBasis = NodalBasisFamily<DeviceType,ValueType,ValueType>;
        using CG_DNBasis = DerivedNodalBasisFamily<DeviceType,ValueType,ValueType>;
        using CG_HBasis = HierarchicalBasisFamily<DeviceType,ValueType,ValueType>;
        std::vector<basisType*> basis_set;

        //compute reference points
        typename CG_NBasis::HGRAD_HEX warpBasis(order,POINTTYPE_WARPBLEND); //used only for computing reference points
        ordinal_type numRefCoords = warpBasis.getCardinality();
        DynRankView ConstructWithLabel(refPoints, numRefCoords, dim);
        warpBasis.getDofCoords(refPoints);

        // compute orientations for cells (one time computation)
        DynRankViewInt elemNodes(&hexas[0][0], numCells, numElemVertexes);
        Kokkos::DynRankView<Orientation,DeviceType> elemOrts("elemOrts", numCells);
        ots::getOrientation(elemOrts, elemNodes, hexa);

        //Compute physical Dof Coordinates and Reference coordinates
        DynRankView ConstructWithLabel(physRefCoords, numCells, numRefCoords, dim);
        {
          Basis_HGRAD_HEX_C1_FEM<DeviceType,ValueType,ValueType> hexaLinearBasis; //used for computing physical coordinates
          DynRankView ConstructWithLabel(hexaLinearBasisValuesAtRefCoords, hexa.getNodeCount(), numRefCoords);
          hexaLinearBasis.getValues(hexaLinearBasisValuesAtRefCoords, refPoints);
          for(ordinal_type i=0; i<numCells; ++i)
            for(ordinal_type d=0; d<dim; ++d) {
              for(ordinal_type j=0; j<numRefCoords; ++j)
                for(std::size_t k=0; k<hexa.getNodeCount(); ++k)
                  physRefCoords(i,j,d) += vertices[hexas[i][k]][d]*hexaLinearBasisValuesAtRefCoords(k,j);
            }
        }


        //HGRAD BASIS
        {
          Fun fun;
          DynRankView ConstructWithLabel(funAtPhysRefCoords, numCells, numRefCoords);
          for(ordinal_type i=0; i<numCells; ++i) {
            for(ordinal_type j=0; j<numRefCoords; ++j)
              funAtPhysRefCoords(i,j) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2));
          }

          basis_set.push_back(new typename  CG_NBasis::HGRAD_HEX(order));
          basis_set.push_back(new typename  CG_DNBasis::HGRAD_HEX(order));
          basis_set.push_back(new typename  CG_HBasis::HGRAD_HEX(order));

          for (auto basisPtr:basis_set) {
            auto& basis = *basisPtr;
            auto name = basis.getName();
            *outStream << " " << name << std::endl;
            ordinal_type basisCardinality = basis.getCardinality();

            //check that fun values at reference points coincide with those computed using basis functions
            DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords);
            DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords);
            DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoords, numCells, basisCardinality, numRefCoords);

            DynRankView ConstructWithLabel(basisValuesAtRefCoords, basisCardinality, numRefCoords);
            basis.getValues(basisValuesAtRefCoords, refPoints);

            // modify basis values to account for orientations
            ots::modifyBasisByOrientation(basisValuesAtRefCoordsOriented,
                basisValuesAtRefCoords,
                elemOrts,
                &basis);
            
            // transform basis values
            fst::HGRADtransformVALUE(transformedBasisValuesAtRefCoords, basisValuesAtRefCoords);
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

            TestResults testResults(basisCoeffs, transformedBasisValuesAtRefCoords, transformedBasisValuesAtRefCoordsOriented, funAtPhysRefCoords, elemOrts, basisPtr);
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
                funAtPhysRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2), k);
            }
          }

          basis_set.clear();
          basis_set.push_back(new typename  CG_NBasis::HCURL_HEX(order));
          basis_set.push_back(new typename  CG_DNBasis::HCURL_HEX(order));
          basis_set.push_back(new typename  CG_HBasis::HCURL_HEX(order));

          for (auto basisPtr:basis_set) {
            auto& basis = *basisPtr;
            auto name = basis.getName();
            *outStream << " " << name << std::endl;

            ordinal_type basisCardinality = basis.getCardinality();
            DynRankView ConstructWithLabel(basisCoeffs, numCells, basisCardinality);

            //check that fun values at reference points coincide with those computed using basis functions
            DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
            DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
            DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoords, numCells, basisCardinality, numRefCoords, dim);
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
            ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, hexa);
            ct::setJacobianInv (jacobianAtRefCoords_inv, jacobianAtRefCoords);
            fst::HCURLtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
                jacobianAtRefCoords_inv,
                basisValuesAtRefCoordsOriented);
            fst::HCURLtransformVALUE(transformedBasisValuesAtRefCoords,
                jacobianAtRefCoords_inv,
                basisValuesAtRefCoords);

            BasisFunctionsSystem  basisFunctionsSystem(basisCardinality, numRefCoords, dim);
            auto info = basisFunctionsSystem.computeBasisCoeffs(basisCoeffs, errorFlag, transformedBasisValuesAtRefCoordsOriented, funAtPhysRefCoords);

            for (int ic =0; ic < numCells; ++ic) {
              if(info[ic] != 0) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << "LAPACK error flag for cell " << ic << " is: " << info[ic] << std::endl;
              }
            }

            TestResults testResults(basisCoeffs, transformedBasisValuesAtRefCoords, transformedBasisValuesAtRefCoordsOriented, funAtPhysRefCoords, elemOrts, basisPtr);
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
                funAtPhysRefCoords(i,j,k) = fun(physRefCoords(i,j,0), physRefCoords(i,j,1), physRefCoords(i,j,2), k);
            }
          }
          basis_set.clear();
          basis_set.push_back(new typename  CG_NBasis::HDIV_HEX(order));
          basis_set.push_back(new typename  CG_DNBasis::HDIV_HEX(order));
          basis_set.push_back(new typename  CG_HBasis::HDIV_HEX(order));

          for (auto basisPtr:basis_set) {
            auto& basis = *basisPtr;
            auto name = basis.getName();
            *outStream << " " << name << std::endl;
            ordinal_type basisCardinality = basis.getCardinality();

            //check that fun values at reference points coincide with those computed using basis functions
            DynRankView ConstructWithLabel(basisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
            DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoordsOriented, numCells, basisCardinality, numRefCoords, dim);
            DynRankView ConstructWithLabel(transformedBasisValuesAtRefCoords, numCells, basisCardinality, numRefCoords, dim);
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
            ct::setJacobian(jacobianAtRefCoords, refPoints, physVertexes, hexa);
            ct::setJacobianDet (jacobianAtRefCoords_det, jacobianAtRefCoords);
            fst::HDIVtransformVALUE(transformedBasisValuesAtRefCoordsOriented,
                jacobianAtRefCoords,
                jacobianAtRefCoords_det,
                basisValuesAtRefCoordsOriented);
            fst::HDIVtransformVALUE(transformedBasisValuesAtRefCoords,
                jacobianAtRefCoords,
                jacobianAtRefCoords_det,
                basisValuesAtRefCoords);

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

            TestResults testResults(basisCoeffs, transformedBasisValuesAtRefCoords, transformedBasisValuesAtRefCoordsOriented, funAtPhysRefCoords, elemOrts, basisPtr);
            testResults.test(errorFlag,tol);
            delete basisPtr;
          }
        }
      } //rotation of first cell vertices
    } while(std::next_permutation(&reorder[0]+4, &reorder[0]+8)); //reorder vertices of common face

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

