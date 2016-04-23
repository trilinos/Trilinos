1;2c// @HEADER
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


/** \file   Intrepid_CellToolsDef.hpp
    \brief  Definition file for the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_HPP__
#define __INTREPID2_CELLTOOLS_DEF_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  
  //============================================================================================//          
  //                                                                                            //          
  //                      Reference-to-physical frame mapping and its inverse                   //          
  //                                                                                            //          
  //============================================================================================//   
  template<class Scalar>
  template<class ArrayPhysPoint, class ArrayRefPoint, class ArrayCell>
  void CellTools<Scalar>::mapToPhysicalFrame(ArrayPhysPoint      &        physPoints,
                                             const ArrayRefPoint &        refPoints,
                                             const ArrayCell     &        cellWorkset,
                                             const shards::CellTopology & cellTopo,
                                             const int &                  whichCell)
  {
    INTREPID2_VALIDATE(validateArguments_mapToPhysicalFrame( physPoints, refPoints, cellWorkset, cellTopo, whichCell) );

    ArrayWrapper<Scalar,ArrayPhysPoint, Rank<ArrayPhysPoint >::value, false>physPointsWrap(physPoints);
    ArrayWrapper<Scalar,ArrayRefPoint, Rank<ArrayRefPoint >::value, true>refPointsWrap(refPoints);
    ArrayWrapper<Scalar,ArrayCell, Rank<ArrayCell >::value,true>cellWorksetWrap(cellWorkset);

    index_type spaceDim  = (index_type)cellTopo.getDimension();
    index_type numCells  = static_cast<index_type>(cellWorkset.dimension(0));
    //points can be rank-2 (P,D), or rank-3 (C,P,D)
    index_type numPoints = (getrank(refPoints) == 2) ? static_cast<index_type>(refPoints.dimension(0)) : static_cast<index_type>(refPoints.dimension(1));

    // Mapping is computed using an appropriate H(grad) basis function: define RCP to the base class
    Teuchos::RCP<Basis<Scalar, FieldContainer<Scalar> > > HGRAD_Basis;
    // Choose the H(grad) basis depending on the cell topology. \todo define maps for shells and beams
    switch( cellTopo.getKey() ){

      // Standard Base topologies (number of cellWorkset = number of vertices)
    case shards::Line<2>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_LINE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;

    case shards::Triangle<3>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TRI_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
    case shards::Quadrilateral<4>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_QUAD_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
    case shards::Tetrahedron<4>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
    case shards::Hexahedron<8>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
    case shards::Wedge<6>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
    case shards::Pyramid<5>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_PYR_C1_FEM<Scalar, FieldContainer<Scalar> >() );
      break;

      // Standard Extended topologies
    case shards::Triangle<6>::key:    
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TRI_C2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
    case shards::Quadrilateral<9>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_QUAD_C2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
    case shards::Tetrahedron<10>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_C2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;

    case shards::Tetrahedron<11>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_TET_COMP12_FEM<Scalar, FieldContainer<Scalar> >() );
      break;

    case shards::Hexahedron<20>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_I2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;

    case shards::Hexahedron<27>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_HEX_C2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;

    case shards::Wedge<15>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_I2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
    case shards::Wedge<18>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_WEDGE_C2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;

    case shards::Pyramid<13>::key:
      HGRAD_Basis = Teuchos::rcp( new Basis_HGRAD_PYR_I2_FEM<Scalar, FieldContainer<Scalar> >() );
      break;
      
      // These extended topologies are not used for mapping purposes
    case shards::Quadrilateral<8>::key:
      INTREPID2_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): Cell topology not supported. ");
      break;

      // Base and Extended Line, Beam and Shell topologies  
    case shards::Line<3>::key:
    case shards::Beam<2>::key:
    case shards::Beam<3>::key:
    case shards::ShellLine<2>::key:
    case shards::ShellLine<3>::key:
    case shards::ShellTriangle<3>::key:
    case shards::ShellTriangle<6>::key:
    case shards::ShellQuadrilateral<4>::key:
    case shards::ShellQuadrilateral<8>::key:
    case shards::ShellQuadrilateral<9>::key:
      INTREPID2_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): Cell topology not supported. ");
      break;
    default:
      INTREPID2_TEST_FOR_EXCEPTION( (true), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::mapToPhysicalFrame): Cell topology not supported.");        
    }// switch  
    // Temp (F,P) array for the values of nodal basis functions at the reference points
    int basisCardinality = HGRAD_Basis -> getCardinality();
    FieldContainer<Scalar> basisVals(basisCardinality, numPoints);

    // Initialize physPoints
    if(getrank(physPoints)==3){
      for(index_type i = 0; i < static_cast<index_type>(physPoints.dimension(0)); i++) {
        for(index_type j = 0; j < static_cast<index_type>(physPoints.dimension(1)); j++){
          for(index_type k = 0; k < static_cast<index_type>(physPoints.dimension(2)); k++){ 
            physPointsWrap(i,j,k) = 0.0;
          }
        }
      }
    }else if(getrank(physPoints)==2){
      for(index_type i = 0; i < static_cast<index_type>(physPoints.dimension(0)); i++){
	for(index_type j = 0; j < static_cast<index_type>(physPoints.dimension(1)); j++){ 
          physPointsWrap(i,j) = 0.0;
        }
      }
    }
    //#else
    //   Kokkos::deep_copy(physPoints.get_kokkos_view(), Scalar(0.0));  
    //#endif
    // handle separately rank-2 (P,D) and rank-3 (C,P,D) cases of refPoints
    switch(getrank(refPoints)) {
    
      // refPoints is (P,D): single set of ref. points is mapped to one or multiple physical cells
    case 2:
      {
        // getValues requires rank-2 (P,D) input array, but refPoints cannot be passed directly as argument because they are a user type
        FieldContainer<Scalar> tempPoints( static_cast<index_type>(refPoints.dimension(0)), static_cast<index_type>(refPoints.dimension(1)) );
        // Copy point set corresponding to this cell oridinal to the temp (P,D) array
        for(index_type pt = 0; pt < static_cast<index_type>(refPoints.dimension(0)); pt++){
          for(index_type dm = 0; dm < static_cast<index_type>(refPoints.dimension(1)) ; dm++){
            tempPoints(pt, dm) = refPointsWrap(pt, dm);
          }//dm
        }//pt
        HGRAD_Basis -> getValues(basisVals, tempPoints, OPERATOR_VALUE);

        // If whichCell = -1, ref pt. set is mapped to all cells, otherwise, the set is mapped to one cell only
        index_type cellLoop = (whichCell == -1) ? numCells : 1 ;

        // Compute the map F(refPoints) = sum node_coordinate*basis(refPoints)
        for(index_type cellOrd = 0; cellOrd < cellLoop; cellOrd++) {
          for(index_type pointOrd = 0; pointOrd < numPoints; pointOrd++) {
            for(index_type dim = 0; dim < spaceDim; dim++){
              for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
                
                if(whichCell == -1){
                  physPointsWrap(cellOrd, pointOrd, dim) += cellWorksetWrap(cellOrd, bfOrd, dim)*basisVals(bfOrd, pointOrd);
                }
                else{
                  physPointsWrap(pointOrd, dim) += cellWorksetWrap(whichCell, bfOrd, dim)*basisVals(bfOrd, pointOrd);
                }
              } // bfOrd
            }// dim
          }// pointOrd
        }//cellOrd
      }// case 2
  
      break;
      
      // refPoints is (C,P,D): multiple sets of ref. points are mapped to matching number of physical cells.  
    case 3:
      {
        // getValues requires rank-2 (P,D) input array, refPoints cannot be used as argument: need temp (P,D) array
        FieldContainer<Scalar> tempPoints( static_cast<index_type>(refPoints.dimension(1)), static_cast<index_type>(refPoints.dimension(2)) );
        
        // Compute the map F(refPoints) = sum node_coordinate*basis(refPoints)
        for(index_type cellOrd = 0; cellOrd < numCells; cellOrd++) {
          
          // Copy point set corresponding to this cell oridinal to the temp (P,D) array
          for(index_type pt = 0; pt < static_cast<index_type>(refPoints.dimension(1)); pt++){
            for(index_type dm = 0; dm < static_cast<index_type>(refPoints.dimension(2)) ; dm++){
              tempPoints(pt, dm) = refPointsWrap(cellOrd, pt, dm);
            }//dm
          }//pt
          
          // Compute basis values for this set of ref. points
          HGRAD_Basis -> getValues(basisVals, tempPoints, OPERATOR_VALUE);
          
          for(index_type pointOrd = 0; pointOrd < numPoints; pointOrd++) {
            for(index_type dim = 0; dim < spaceDim; dim++){
              for(int bfOrd = 0; bfOrd < basisCardinality; bfOrd++){
                
                physPointsWrap(cellOrd, pointOrd, dim) += cellWorksetWrap(cellOrd, bfOrd, dim)*basisVals(bfOrd, pointOrd);
                
              } // bfOrd
            }// dim
          }// pointOrd
        }//cellOrd        
      }// case 3
      break;
      
   
    }
  }	

  template<class Scalar>
  template<class ArraySubcellPoint, class ArrayParamPoint>
  void CellTools<Scalar>::mapToReferenceSubcell(ArraySubcellPoint     &       refSubcellPoints,
                                                const ArrayParamPoint &       paramPoints,
                                                const int                     subcellDim,
                                                const int                     subcellOrd,
                                                const shards::CellTopology &  parentCell){
  
    int cellDim = parentCell.getDimension();
    index_type numPts  = static_cast<index_type>(paramPoints.dimension(0));

#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( !(hasReferenceCell(parentCell) ), std::invalid_argument, 
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): the specified cell topology does not have a reference cell.");
  
    INTREPID2_TEST_FOR_EXCEPTION( !( (1 <= subcellDim) && (subcellDim <= 2 ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): method defined only for 1 and 2-dimensional subcells.");  
  
    INTREPID2_TEST_FOR_EXCEPTION( !( (0 <= subcellOrd) && (subcellOrd < (int)parentCell.getSubcellCount(subcellDim) ) ), std::invalid_argument,
                                  ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): subcell ordinal out of range.");
  
    // refSubcellPoints is rank-2 (P,D1), D1 = cell dimension
    std::string errmsg = ">>> ERROR (Intrepid2::mapToReferenceSubcell):";
    INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, refSubcellPoints, 2,2), std::invalid_argument, errmsg);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, refSubcellPoints, 1, cellDim, cellDim), std::invalid_argument, errmsg);
                    
    // paramPoints is rank-2 (P,D2) with D2 = subcell dimension
    INTREPID2_TEST_FOR_EXCEPTION( !requireRankRange(errmsg, paramPoints, 2,2), std::invalid_argument, errmsg);
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionRange(errmsg, paramPoints, 1, subcellDim, subcellDim), std::invalid_argument, errmsg);    
  
    // cross check: refSubcellPoints and paramPoints: dimension 0 must match
    INTREPID2_TEST_FOR_EXCEPTION( !requireDimensionMatch(errmsg, refSubcellPoints, 0,  paramPoints, 0), std::invalid_argument, errmsg);      
#endif
  
  
    // Get the subcell map, i.e., the coefficients of the parametrization function for the subcell
    const FieldContainer<double>& subcellMap = getSubcellParametrization(subcellDim, parentCell);

    // Apply the parametrization map to every point in parameter domain
    if(subcellDim == 2) {
      for(index_type pt = 0; pt < numPts; pt++){
        double u = paramPoints(pt,0);
        double v = paramPoints(pt,1);
      
        // map_dim(u,v) = c_0(dim) + c_1(dim)*u + c_2(dim)*v because both Quad and Tri ref faces are affine!
        for(int  dim = 0; dim < cellDim; dim++){
          refSubcellPoints(pt, dim) = subcellMap(subcellOrd, dim, 0) + \
            subcellMap(subcellOrd, dim, 1)*u + \
            subcellMap(subcellOrd, dim, 2)*v;
        }
      }
    }
    else if(subcellDim == 1) {    
      for(index_type pt = 0; pt < numPts; pt++){
        for(int dim = 0; dim < cellDim; dim++) {
          refSubcellPoints(pt, dim) = subcellMap(subcellOrd, dim, 0) + subcellMap(subcellOrd, dim, 1)*paramPoints(pt,0);
        }
      }
    }
    else{
      INTREPID2_TEST_FOR_EXCEPTION( !( (subcellDim == 1) || (subcellDim == 2) ), std::invalid_argument, 
                                    ">>> ERROR (Intrepid2::CellTools::mapToReferenceSubcell): method defined only for 1 and 2-subcells");
    }
  }
  
}

#endif
