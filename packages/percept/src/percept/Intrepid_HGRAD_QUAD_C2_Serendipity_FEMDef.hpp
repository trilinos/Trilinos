// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef INTREPID_HGRAD_QUAD_C2_SERENDIPITY_FEMDEF_HPP
#define INTREPID_HGRAD_QUAD_C2_SERENDIPITY_FEMDEF_HPP

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

/** \file   Intrepid_HGRAD_QUAD_C2_Serendipity_FEMDef.hpp
    \brief  Definition file for bi-linear FEM basis functions for H(grad) functions on QUAD cells.
    \author Created by P. Bochev and D. Ridzal. (serendipity version by S. R. Kennon, srkenno@sandia.gov)
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_QUAD_C2_Serendipity_FEM<Scalar, ArrayScalar>::Basis_HGRAD_QUAD_C2_Serendipity_FEM()
  {
    this -> basisCardinality_  = 8;
    this -> basisDegree_       = 2;    
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }
    
  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_QUAD_C2_Serendipity_FEM<Scalar, ArrayScalar>::initializeTags() {
  
    // Basis-dependent intializations
    int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

    // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
    int tags[]  = { 0, 0, 0, 1,
                    0, 1, 0, 1,
                    0, 2, 0, 1,
                    0, 3, 0, 1,
                    // edge midpoints
                    1, 0, 0, 1,
                    1, 1, 0, 1,
                    1, 2, 0, 1,
                    1, 3, 0, 1
                  
                    //// quad center
                    // ,2, 0, 0, 1
    };
  
    // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
    Intrepid::setOrdinalTagData(this -> tagToOrdinal_,
                                this -> ordinalToTag_,
                                tags,
                                this -> basisCardinality_,
                                tagSize,
                                posScDim,
                                posScOrd,
                                posDfOrd);
  }



  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_QUAD_C2_Serendipity_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
                                                                           const ArrayScalar &  inputPoints,
                                                                           const EOperator      operatorType) const {
  
    // Verify arguments
#ifdef HAVE_INTREPID_DEBUG
    Intrepid::getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
                                                        inputPoints,
                                                        operatorType,
                                                        this -> getBaseCellTopology(),
                                                        this -> getCardinality() );
#endif
  
    static double xc[] = {-1, 1, 1, -1};
    static double yc[] = {-1, -1, 1, 1};
    static double xm[] = {0, 1, 0, -1};
    static double ym[] = {-1, 0, 1, 0};

    // Number of evaluation points = dim 0 of inputPoints
    int dim0 = inputPoints.dimension(0);  
  
    // Temporaries: (x,y) coordinates of the evaluation point
    Scalar x = 0.0;                                    
    Scalar y = 0.0;                                    
  
    switch (operatorType) {
    
    case OPERATOR_VALUE:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0, 0);
        y = inputPoints(i0, 1);
        
        // corner nodes
        // 1/4 (1 + x x_i) (1 + y y_i) ( x x_i + y y_i - 1.0)
        // mid nodes on x
        // 1/2 (1 - x^2)(1 + y y_i)
        // mid nodes on y
        // 1/2 (1 - y^2)(1 + x x_i)

        // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
        
        for (int ic=0; ic < 4; ic++)
          {
            outputValues(ic, i0) =  (1.0 + x * xc[ic]) * (1.0 + y * yc[ic]) * ( x * xc[ic] + y * yc[ic] - 1.0) / 4.0;
          }

        for (int im=0; im < 3; im += 2)
          {
            outputValues(im+4, i0)   =  (1.0 - x*x)*(1.0 + y * ym[im]) / 2.0;     // 4,6 = 0, 2
            outputValues(im+4+1, i0) =  (1.0 - y*y)*(1.0 + x * xm[im+1]) / 2.0;   // 5,7 = 1, 3
          }
      }
      break;
      
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        for (int ic=0; ic < 4; ic++)
          {
            //outputValues(ic, i0) =  (1.0 + x * xc[ic]) * (1.0 + y * yc[ic]) * ( x * xc[ic] + y * yc[ic] - 1.0) / 4.0;
            outputValues(ic, i0, 0) =  
              (          xc[ic]) * (1.0 + y * yc[ic]) * ( x * xc[ic] + y * yc[ic] - 1.0) / 4.0 +
              (1.0 + x * xc[ic]) * (1.0 + y * yc[ic]) * (     xc[ic]                   ) / 4.0;
            outputValues(ic, i0, 1) =  
              (1.0 + x * xc[ic]) * (          yc[ic]) * ( x * xc[ic] + y * yc[ic] - 1.0) / 4.0 +
              (1.0 + x * xc[ic]) * (1.0 + y * yc[ic]) * (                  yc[ic]      ) / 4.0;
          }

        for (int im=0; im < 3; im += 2)
          {
            //outputValues(im+4, i0)   =  (1.0 -   x*x)*(1.0 + y * ym[im]) / 2.0;     // 4,6 = 0, 2
            outputValues(im+4, i0, 0)  =  (    - 2.0*x)*(1.0 + y * ym[im]) / 2.0;
            outputValues(im+4, i0, 1)  =  (1.0 -   x*x)*(          ym[im]) / 2.0;    

            //outputValues(im+4+1, i0)  =  (1.0 -   y*y)*(1.0 + x * xm[im+1]) / 2.0;   // 5,7 = 1, 3
            outputValues(im+4+1, i0, 0) =  (1.0 -   y*y)*(          xm[im+1]) / 2.0;   
            outputValues(im+4+1, i0, 1) =  (    - 2.0*y)*(1.0 + x * xm[im+1]) / 2.0;   
          }


      }
      break;
      
    case OPERATOR_CURL:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        // CURL(u) = (u_y, -u_x), is rotated GRAD
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        for (int ic=0; ic < 4; ic++)
          {
            //outputValues(ic, i0) =  (1.0 + x * xc[ic]) * (1.0 + y * yc[ic]) * ( x * xc[ic] + y * yc[ic] - 1.0) / 4.0;
            outputValues(ic, i0, 1) = -( 
                                        (          xc[ic]) * (1.0 + y * yc[ic]) * ( x * xc[ic] + y * yc[ic] - 1.0) / 4.0 +
                                        (1.0 + x * xc[ic]) * (1.0 + y * yc[ic]) * (     xc[ic]                   ) / 4.0);
            outputValues(ic, i0, 0) =  
              (1.0 + x * xc[ic]) * (          yc[ic]) * ( x * xc[ic] + y * yc[ic] - 1.0) / 4.0 +
              (1.0 + x * xc[ic]) * (1.0 + y * yc[ic]) * (                  yc[ic]      ) / 4.0;
          }

        for (int im=0; im < 3; im += 2)
          {
            //outputValues(im+4, i0)   =   (1.0 -   x*x)*(1.0 + y * ym[im]) / 2.0;     // 4,6 = 0, 2
            outputValues(im+4, i0, 1)  = - (    - 2.0*x)*(1.0 + y * ym[im]) / 2.0;
            outputValues(im+4, i0, 0)  =   (1.0 -   x*x)*(          ym[im]) / 2.0;    

            //outputValues(im+4+1, i0)  =   (1.0 -   y*y)*(1.0 + x * xm[im+1]) / 2.0;   // 5,7 = 1, 3
            outputValues(im+4+1, i0, 1) = - (1.0 -   y*y)*(          xm[im+1]) / 2.0;   
            outputValues(im+4+1, i0, 0) =   (    - 2.0*y)*(1.0 + x * xm[im+1]) / 2.0;   
          }

      }
      break;
      
    case OPERATOR_DIV:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_QUAD_C2_Serendipity_FEM): DIV is invalid operator for rank-0 (scalar) functions in 2D");
      break;
      
    case OPERATOR_D2:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
       
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, D2Cardinality=3) 
        for (int ic=0; ic < 4; ic++)
          {
            //outputValues(ic, i0) =  (1.0 + x * xc[ic]) * (1.0 + y * yc[ic]) * ( x * xc[ic] + y * yc[ic] - 1.0) / 4.0;
            outputValues(ic, i0, 0) =  
              (          xc[ic]) * (1.0 + y * yc[ic]) * (     xc[ic]                   ) / 4.0 +
              (          xc[ic]) * (1.0 + y * yc[ic]) * (     xc[ic]                   ) / 4.0;
            outputValues(ic, i0, 1) =  
              (          xc[ic]) * (          yc[ic]) * ( x * xc[ic] + y * yc[ic] - 1.0) / 4.0 +
              (          xc[ic]) * (1.0 + y * yc[ic]) * (                  yc[ic]      ) / 4.0 +
              (1.0 + x * xc[ic]) * (          yc[ic]) * (     xc[ic]                   ) / 4.0;

            outputValues(ic, i0, 2) =  
              (1.0 + x * xc[ic]) * (          yc[ic]) * (                  yc[ic]      ) / 4.0 +
              (1.0 + x * xc[ic]) * (          yc[ic]) * (                  yc[ic]      ) / 4.0;
          }

        for (int im=0; im < 3; im += 2)
          {
            //outputValues(im+4, i0)   =  (1.0 -   x*x)*(1.0 + y * ym[im]) / 2.0;     // 4,6 = 0, 2
            outputValues(im+4, i0, 0)  =  (    - 2.0  )*(1.0 + y * ym[im]) / 2.0;
            outputValues(im+4, i0, 1)  =  (    - 2.0*x)*(          ym[im]) / 2.0;

            outputValues(im+4, i0, 2)  =  0.0;

            //outputValues(im+4+1, i0)  =  (1.0 -   y*y)*(1.0 + x * xm[im+1]) / 2.0;   // 5,7 = 1, 3
            outputValues(im+4+1, i0, 0) =  0.0;
            outputValues(im+4+1, i0, 1) =  (    - 2.0*y)*(          xm[im+1]) / 2.0;   

            outputValues(im+4+1, i0, 2) =  (    - 2.0  )*(1.0 + x * xm[im+1]) / 2.0;   
          }
      }
      break;
      
    case OPERATOR_D3:
      // x
      // y

      // xx, xy
      // yy

      // xxx, xxy
      // xyy
      // yyy

      // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, D3Cardinality=4) 
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);

        for (int ic=0; ic < 4; ic++)
          {
            //outputValues(ic, i0) =  (1.0 + x * xc[ic]) * (1.0 + y * yc[ic]) * ( x * xc[ic] + y * yc[ic] - 1.0) / 4.0;
            outputValues(ic, i0, 0) =  0.0;

            outputValues(ic, i0, 1) =  
              (          xc[ic]) * (          yc[ic]) * (     xc[ic]                   ) / 4.0 +
              (          xc[ic]) * (          yc[ic]) * (     xc[ic]                   ) / 4.0;

            outputValues(ic, i0, 2) =  
              (          xc[ic]) * (          yc[ic]) * (                  yc[ic]      ) / 4.0 +
              (          xc[ic]) * (          yc[ic]) * (                  yc[ic]      ) / 4.0;


            outputValues(ic, i0, 3) =  0.0;
          }

        for (int im=0; im < 3; im += 2)
          {
            //outputValues(im+4, i0)   =  (1.0 -   x*x)*(1.0 + y * ym[im]) / 2.0;     // 4,6 = 0, 2
            outputValues(im+4, i0, 0)  =  0.0;
            outputValues(im+4, i0, 1)  =  (    - 2.0  )*(          ym[im]) / 2.0;

            outputValues(im+4, i0, 2)  =  0.0;

            outputValues(im+4, i0, 3)  =  0.0;

            //outputValues(im+4+1, i0)  =  (1.0 -   y*y)*(1.0 + x * xm[im+1]) / 2.0;   // 5,7 = 1, 3
            outputValues(im+4+1, i0, 0) =  0.0;
            outputValues(im+4+1, i0, 1) =  0.0;

            outputValues(im+4+1, i0, 2) =  (    - 2.0  )*(          xm[im+1]) / 2.0;   

            outputValues(im+4+1, i0, 3) =  0.0;
          }
      }
      break;
    
    case OPERATOR_D4:
      // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, D4Cardinality=5) 
      for (int i0 = 0; i0 < dim0; i0++) {
        
        for (int ib = 0; ib < 8; ib++)
          {
            for (int dc = 0; dc < 5; dc++)
              {
                outputValues(ib, i0, dc) = 0.0;
              }
          }
      }
      break;
      
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      {
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, DkCardinality)
        int DkCardinality = Intrepid::getDkCardinality(operatorType, 
                                                       this -> basisCellTopology_.getDimension() );
        for(int dofOrd = 0; dofOrd < this -> basisCardinality_; dofOrd++) {
          for (int i0 = 0; i0 < dim0; i0++) {
            for(int dkOrd = 0; dkOrd < DkCardinality; dkOrd++){
              outputValues(dofOrd, i0, dkOrd) = 0.0;
            }
          }
        }
      }
      break;
      
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( !( Intrepid::isValidOperator(operatorType) ), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_QUAD_C2_Serendipity_FEM): Invalid operator type");
    }
  }


  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_QUAD_C2_Serendipity_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                                           const ArrayScalar &    inputPoints,
                                                                           const ArrayScalar &    cellVertices,
                                                                           const EOperator        operatorType) const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                        ">>> ERROR (Basis_HGRAD_QUAD_C2_Serendipity_FEM): FEM Basis calling an FVD member function");
  }



  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_QUAD_C2_Serendipity_FEM<Scalar, ArrayScalar>::getDofCoords(ArrayScalar & DofCoords) const {
#ifdef HAVE_INTREPID_DEBUG
    // Verify rank of output array.
    TEUCHOS_TEST_FOR_EXCEPTION( !(DofCoords.rank() == 2), std::invalid_argument,
                        ">>> ERROR: (Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM::getDofCoords) rank = 2 required for DofCoords array");
    // Verify 0th dimension of output array.
    TEUCHOS_TEST_FOR_EXCEPTION( !( DofCoords.dimension(0) == this -> basisCardinality_ ), std::invalid_argument,
                        ">>> ERROR: (Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM::getDofCoords) mismatch in number of DoF and 0th dimension of DofCoords array");
    // Verify 1st dimension of output array.
    TEUCHOS_TEST_FOR_EXCEPTION( !( DofCoords.dimension(1) == (int)(this -> basisCellTopology_.getDimension()) ), std::invalid_argument,
                        ">>> ERROR: (Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM::getDofCoords) incorrect reference cell (1st) dimension in DofCoords array");
#endif

    DofCoords(0,0) = -1.0;   DofCoords(0,1) = -1.0;
    DofCoords(1,0) =  1.0;   DofCoords(1,1) = -1.0;
    DofCoords(2,0) =  1.0;   DofCoords(2,1) =  1.0;
    DofCoords(3,0) = -1.0;   DofCoords(3,1) =  1.0;

    DofCoords(4,0) =  0.0;   DofCoords(4,1) = -1.0;
    DofCoords(5,0) =  1.0;   DofCoords(5,1) =  0.0;
    DofCoords(6,0) =  0.0;   DofCoords(6,1) =  1.0;
    DofCoords(7,0) = -1.0;   DofCoords(7,1) =  0.0;

    //DofCoords(8,0) =  0.0;   DofCoords(8,1) =  0.0;

  }

}// namespace Intrepid
#endif
