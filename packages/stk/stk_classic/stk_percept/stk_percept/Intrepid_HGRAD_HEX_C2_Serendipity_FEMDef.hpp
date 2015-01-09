#ifndef INTREPID_HGRAD_HEX_C2_SERENDIPITY_FEMDEF_HPP
#define INTREPID_HGRAD_HEX_C2_SERENDIPITY_FEMDEF_HPP
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

/** \file   Intrepid_HGRAD_HEX_C2_Serendipity_FEMDef.hpp
    \brief  Definition file for bi-linear FEM basis functions for H(grad) functions on Hexahedron cells.
    \author Created by P. Bochev and D. Ridzal. (serendipity version by S. R. Kennon, srkenno@sandia.gov)
*/

namespace Intrepid {

  
  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_HEX_C2_Serendipity_FEM<Scalar,ArrayScalar>::Basis_HGRAD_HEX_C2_Serendipity_FEM()
  {
    this -> basisCardinality_  = 20;
    this -> basisDegree_       = 2;    
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }
  
  
  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_HEX_C2_Serendipity_FEM<Scalar, ArrayScalar>::initializeTags() {
  
    // Basis-dependent intializations
    int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

    // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
    int tags[]  = { 0, 0, 0, 1,     // Nodes 0 to 7 follow vertex order of the topology
                    0, 1, 0, 1,
                    0, 2, 0, 1,
                    0, 3, 0, 1,
                    0, 4, 0, 1,
                    0, 5, 0, 1,
                    0, 6, 0, 1,
                    0, 7, 0, 1,
                    1, 0, 0, 1,      // Node 8  -> edge 0
                    1, 1, 0, 1,      // Node 9  -> edge 1
                    1, 2, 0, 1,      // Node 10 -> edge 2
                    1, 3, 0, 1,      // Node 11 -> edge 3
                    1, 8, 0, 1,      // Node 12 -> edge 8
                    1, 9, 0, 1,      // Node 13 -> edge 9
                    1,10, 0, 1,      // Node 14 -> edge 10
                    1,11, 0, 1,      // Node 15 -> edge 11
                    1, 4, 0, 1,      // Node 16 -> edge 4
                    1, 5, 0, 1,      // Node 17 -> edge 5
                    1, 6, 0, 1,      // Node 18 -> edge 6
                    1, 7, 0, 1      // Node 19 -> edge 7

                    //                   ,3, 0, 0, 1,      // Node 20 -> Hexahedron
                    //                   2, 4, 0, 1,      // Node 21 -> face 4
                    //                   2, 5, 0, 1,      // Node 22 -> face 5
                    //                   2, 3, 0, 1,      // Node 23 -> face 3
                    //                   2, 1, 0, 1,      // Node 24 -> face 1
                    //                   2, 0, 0, 1,      // Node 25 -> face 0
                    //                   2, 2, 0, 1,      // Node 26 -> face 2

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
  void Basis_HGRAD_HEX_C2_Serendipity_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
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
  
    // Number of evaluation points = dim 0 of inputPoints
    int dim0 = inputPoints.dimension(0);  
  
    // Temporaries: (x,y,z) coordinates of the evaluation point
    Scalar x = 0.0;                                    
    Scalar y = 0.0;                                    
    Scalar z = 0.0;                                    
    Scalar one = 1.0;
    Scalar two = 2.0;

    Scalar one8th = 1./8.;
    Scalar one4th = 1./4.;
    Scalar half = 1./2.;

    switch (operatorType) {
    
    case OPERATOR_VALUE:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0, 0);
        y = inputPoints(i0, 1);
        z = inputPoints(i0, 2);
        Scalar s = x;
        Scalar t = y;
        Scalar u = z;

        //Scalar stu = s * t * u;
        //Scalar st  = s * t;
        //Scalar su  = s * u;
        //Scalar tu  = t * u;

        Scalar one_m_s = one - s;
        Scalar one_p_s = one + s;
        Scalar one_m_t = one - t;
        Scalar one_p_t = one + t;
        Scalar one_m_u = one - u;
        Scalar one_p_u = one + u;

        Scalar um2 = u - two;
        Scalar up2 = u + two;

        Scalar one_m_ss = one - s * s;
        Scalar one_m_tt = one - t * t;
        Scalar one_m_uu = one - u * u;


        // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
        outputValues( 0, i0)= one8th * one_m_s  * one_m_t  * one_m_u * (-s-t-up2);
        outputValues( 1, i0)= one8th * one_p_s  * one_m_t  * one_m_u * ( s-t-up2);
        outputValues( 2, i0)= one8th * one_p_s  * one_p_t  * one_m_u * ( s+t-up2);
        outputValues( 3, i0)= one8th * one_m_s  * one_p_t  * one_m_u * (-s+t-up2);
        outputValues( 4, i0)= one8th * one_m_s  * one_m_t  * one_p_u * (-s-t+um2);
        outputValues( 5, i0)= one8th * one_p_s  * one_m_t  * one_p_u * ( s-t+um2);
        outputValues( 6, i0)= one8th * one_p_s  * one_p_t  * one_p_u * ( s+t+um2);
        outputValues( 7, i0)= one8th * one_m_s  * one_p_t  * one_p_u * (-s+t+um2);
        outputValues( 8, i0)= one4th * one_m_ss * one_m_t  * one_m_u ;
        outputValues( 9, i0)= one4th * one_p_s  * one_m_tt * one_m_u ;
        outputValues(10, i0)= one4th * one_m_ss * one_p_t  * one_m_u ;
        outputValues(11, i0)= one4th * one_m_s  * one_m_tt * one_m_u ;
        outputValues(12, i0)= one4th * one_m_s  * one_m_t  * one_m_uu ;
        outputValues(13, i0)= one4th * one_p_s  * one_m_t  * one_m_uu ;
        outputValues(14, i0)= one4th * one_p_s  * one_p_t  * one_m_uu ;
        outputValues(15, i0)= one4th * one_m_s  * one_p_t  * one_m_uu ;
        outputValues(16, i0)= one4th * one_m_ss * one_m_t  * one_p_u ;
        outputValues(17, i0)= one4th * one_p_s  * one_m_tt * one_p_u ;
        outputValues(18, i0)= one4th * one_m_ss * one_p_t  * one_p_u ;
        outputValues(19, i0)= one4th * one_m_s  * one_m_tt * one_p_u ;
        
      }
      break;
      
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        z = inputPoints(i0,2);

        Scalar s = x;
        Scalar t = y;
        Scalar u = z;

        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)


        double one_m_s = one - s ;
        double one_p_s = one + s ;
        double one_m_t = one - t ;
        double one_p_t = one + t ;
        double one_m_u = one - u ;
        double one_p_u = one + u;

        double s2_m1 = two * s - 1;
        double t2_m1 = two * t - 1;
        double u2_m1 = two * u - 1;
 
        double s2_p1 = two * s + 1;
        double t2_p1 = two * t + 1;
        double u2_p1 = two * u + 1 ;

        double one_m_ss = one - s * s ;
        double one_m_tt = one - t * t ;
        double one_m_uu = one - u * u ;
          
        //double one_m_2s = one - two * s ;
        //double one_m_2t = one - two * t ;
        //double one_m_2u = one - two * u ;

        // shape function derivative in the s direction

        outputValues( 0, i0, 0)= -one8th * one_m_t  * one_m_u  * (-s2_p1-t-u);
        outputValues( 1, i0, 0)=  one8th * one_m_t  * one_m_u  * ( s2_m1-t-u);
        outputValues( 2, i0, 0)=  one8th * one_p_t  * one_m_u  * ( s2_m1+t-u);
        outputValues( 3, i0, 0)= -one8th * one_p_t  * one_m_u  * (-s2_p1+t-u);
        outputValues( 4, i0, 0)= -one8th * one_m_t  * one_p_u  * (-s2_p1-t+u);
        outputValues( 5, i0, 0)=  one8th * one_m_t  * one_p_u  * ( s2_m1-t+u);
        outputValues( 6, i0, 0)=  one8th * one_p_t  * one_p_u  * ( s2_m1+t+u);
        outputValues( 7, i0, 0)= -one8th * one_p_t  * one_p_u  * (-s2_p1+t+u);
        outputValues( 8, i0, 0)=   -half * one_m_t  * one_m_u  * s;
        outputValues( 9, i0, 0)=  one4th * one_m_tt * one_m_u;
        outputValues(10, i0, 0)=   -half * one_p_t  * one_m_u  * s;
        outputValues(11, i0, 0)= -one4th * one_m_tt * one_m_u;
        outputValues(12, i0, 0)= -one4th * one_m_t  * one_m_uu;
        outputValues(13, i0, 0)=  one4th * one_m_t  * one_m_uu;
        outputValues(14, i0, 0)=  one4th * one_p_t  * one_m_uu;
        outputValues(15, i0, 0)= -one4th * one_p_t  * one_m_uu;
        outputValues(16, i0, 0)=   -half * one_m_t  * one_p_u  * s;
        outputValues(17, i0, 0)=  one4th * one_m_tt * one_p_u;
        outputValues(18, i0, 0)=   -half * one_p_t  * one_p_u  * s;
        outputValues(19, i0, 0)= -one4th * one_m_tt * one_p_u;

        // shape function derivative in the t direction
         
        outputValues( 0, i0, 1)= -one8th * one_m_s  * one_m_u  * (-s-t2_p1-u);
        outputValues( 1, i0, 1)= -one8th * one_p_s  * one_m_u  * ( s-t2_p1-u);
        outputValues( 2, i0, 1)=  one8th * one_p_s  * one_m_u  * ( s+t2_m1-u);
        outputValues( 3, i0, 1)=  one8th * one_m_s  * one_m_u  * (-s+t2_m1-u);
        outputValues( 4, i0, 1)= -one8th * one_m_s  * one_p_u  * (-s-t2_p1+u);
        outputValues( 5, i0, 1)= -one8th * one_p_s  * one_p_u  * ( s-t2_p1+u);
        outputValues( 6, i0, 1)=  one8th * one_p_s  * one_p_u  * ( s+t2_m1+u);
        outputValues( 7, i0, 1)=  one8th * one_m_s  * one_p_u  * (-s+t2_m1+u);
        outputValues( 8, i0, 1)= -one4th * one_m_ss * one_m_u;
        outputValues( 9, i0, 1)=   -half * one_p_s  * one_m_u  * t;
        outputValues(10, i0, 1)=  one4th * one_m_ss * one_m_u;
        outputValues(11, i0, 1)=   -half * one_m_s  * one_m_u  * t;
        outputValues(12, i0, 1)= -one4th * one_m_s  * one_m_uu;
        outputValues(13, i0, 1)= -one4th * one_p_s  * one_m_uu;
        outputValues(14, i0, 1)=  one4th * one_p_s  * one_m_uu;
        outputValues(15, i0, 1)=  one4th * one_m_s  * one_m_uu;
        outputValues(16, i0, 1)= -one4th * one_m_ss * one_p_u;
        outputValues(17, i0, 1)=   -half * one_p_s  * one_p_u  * t;
        outputValues(18, i0, 1)=  one4th * one_m_ss * one_p_u;
        outputValues(19, i0, 1)=   -half * one_m_s  * one_p_u  * t;

        // shape function derivative in the u direction
         
        outputValues( 0, i0, 2)= -one8th * one_m_s  * one_m_t  * (-s-t-u2_p1);
        outputValues( 1, i0, 2)= -one8th * one_p_s  * one_m_t  * ( s-t-u2_p1);
        outputValues( 2, i0, 2)= -one8th * one_p_s  * one_p_t  * ( s+t-u2_p1);
        outputValues( 3, i0, 2)= -one8th * one_m_s  * one_p_t  * (-s+t-u2_p1);
        outputValues( 4, i0, 2)=  one8th * one_m_s  * one_m_t  * (-s-t+u2_m1);
        outputValues( 5, i0, 2)=  one8th * one_p_s  * one_m_t  * ( s-t+u2_m1);
        outputValues( 6, i0, 2)=  one8th * one_p_s  * one_p_t  * ( s+t+u2_m1);
        outputValues( 7, i0, 2)=  one8th * one_m_s  * one_p_t  * (-s+t+u2_m1);
        outputValues( 8, i0, 2)= -one4th * one_m_ss * one_m_t;
        outputValues( 9, i0, 2)= -one4th * one_p_s  * one_m_tt;
        outputValues(10, i0, 2)= -one4th * one_m_ss * one_p_t;
        outputValues(11, i0, 2)= -one4th * one_m_s  * one_m_tt;
        outputValues(12, i0, 2)=   -half * one_m_s  * one_m_t  * u;
        outputValues(13, i0, 2)=   -half * one_p_s  * one_m_t  * u;
        outputValues(14, i0, 2)=   -half * one_p_s  * one_p_t  * u;
        outputValues(15, i0, 2)=   -half * one_m_s  * one_p_t  * u;
        outputValues(16, i0, 2)=  one4th * one_m_ss * one_m_t;
        outputValues(17, i0, 2)=  one4th * one_p_s  * one_m_tt;
        outputValues(18, i0, 2)=  one4th * one_m_ss * one_p_t;
        outputValues(19, i0, 2)=  one4th * one_m_s  * one_m_tt;


      }
      break;
      
    case OPERATOR_CURL:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_CURL), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_HEX_C2_Serendipity_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
      break;
      
    case OPERATOR_DIV:
      TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_HEX_C2_Serendipity_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
      break;
      
    case OPERATOR_D2:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        z = inputPoints(i0,2);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, D2Cardinality = 6) 
        TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_D2), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_HEX_C2_Serendipity_FEM): OPERATOR_D2 not yet coded");
        
      }
      break;
      
    case OPERATOR_D3:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        z = inputPoints(i0,2);
                
        TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_D3), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_HEX_C2_Serendipity_FEM): OPERATOR_D3 not yet coded");

      }
      break;
      
    case OPERATOR_D4:
      {
        // Non-zero entries have Dk (derivative cardinality) indices {3,4,5,7,8,12}, all other entries are 0.
        // Intitialize array by zero and then fill only non-zero entries.

        TEUCHOS_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_D4), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_HEX_C2_Serendipity_FEM): OPERATOR_D4 not yet coded");

      }
      break;
      
    case OPERATOR_D5:
    case OPERATOR_D6:
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_HEX_C2_Serendipity_FEM): operator not supported");
      break;
      
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
                          ">>> ERROR (Basis_HGRAD_HEX_C2_Serendipity_FEM): Invalid operator type");
    }
  }



  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_HEX_C2_Serendipity_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                                          const ArrayScalar &    inputPoints,
                                                                          const ArrayScalar &    cellVertices,
                                                                          const EOperator        operatorType) const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                        ">>> ERROR (Basis_HGRAD_HEX_C2_Serendipity_FEM): FEM Basis calling an FVD member function");
  }

  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_HEX_C2_Serendipity_FEM<Scalar, ArrayScalar>::getDofCoords(ArrayScalar & DofCoords) const {
#ifdef HAVE_INTREPID_DEBUG
    // Verify rank of output array.
    TEUCHOS_TEST_FOR_EXCEPTION( !(DofCoords.rank() == 2), std::invalid_argument,
                        ">>> ERROR: (Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM::getDofCoords) rank = 2 required for DofCoords array");
    // Verify 0th dimension of output array.
    TEUCHOS_TEST_FOR_EXCEPTION( !( DofCoords.dimension(0) == this -> basisCardinality_ ), std::invalid_argument,
                        ">>> ERROR: (Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM::getDofCoords) mismatch in number of DoF and 0th dimension of DofCoords array");
    // Verify 1st dimension of output array.
    TEUCHOS_TEST_FOR_EXCEPTION( !( DofCoords.dimension(1) == (int)(this -> basisCellTopology_.getDimension()) ), std::invalid_argument,
                        ">>> ERROR: (Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM::getDofCoords) incorrect reference cell (1st) dimension in DofCoords array");
#endif

    DofCoords(0,0) = -1.0;   DofCoords(0,1) = -1.0; DofCoords(0,2) = -1.0;  
    DofCoords(1,0) =  1.0;   DofCoords(1,1) = -1.0; DofCoords(1,2) = -1.0;  
    DofCoords(2,0) =  1.0;   DofCoords(2,1) =  1.0; DofCoords(2,2) = -1.0;  
    DofCoords(3,0) = -1.0;   DofCoords(3,1) =  1.0; DofCoords(3,2) = -1.0;  
    DofCoords(4,0) = -1.0;   DofCoords(4,1) = -1.0; DofCoords(4,2) =  1.0;  
    DofCoords(5,0) =  1.0;   DofCoords(5,1) = -1.0; DofCoords(5,2) =  1.0;  
    DofCoords(6,0) =  1.0;   DofCoords(6,1) =  1.0; DofCoords(6,2) =  1.0;  
    DofCoords(7,0) = -1.0;   DofCoords(7,1) =  1.0; DofCoords(7,2) =  1.0;  

    DofCoords(8,0) =   0.0;   DofCoords(8,1) =  -1.0; DofCoords(8,2) =  -1.0;  
    DofCoords(9,0) =   1.0;   DofCoords(9,1) =   0.0; DofCoords(9,2) =  -1.0;  
    DofCoords(10,0) =  0.0;   DofCoords(10,1) =  1.0; DofCoords(10,2) = -1.0;  
    DofCoords(11,0) = -1.0;   DofCoords(11,1) =  0.0; DofCoords(11,2) = -1.0;  
    DofCoords(12,0) = -1.0;   DofCoords(12,1) = -1.0; DofCoords(12,2) =  0.0;  
    DofCoords(13,0) =  1.0;   DofCoords(13,1) = -1.0; DofCoords(13,2) =  0.0;  
    DofCoords(14,0) =  1.0;   DofCoords(14,1) =  1.0; DofCoords(14,2) =  0.0;  
    DofCoords(15,0) = -1.0;   DofCoords(15,1) =  1.0; DofCoords(15,2) =  0.0;  
    DofCoords(16,0) =  0.0;   DofCoords(16,1) = -1.0; DofCoords(16,2) =  1.0;  
    DofCoords(17,0) =  1.0;   DofCoords(17,1) =  0.0; DofCoords(17,2) =  1.0;  
    DofCoords(18,0) =  0.0;   DofCoords(18,1) =  1.0; DofCoords(18,2) =  1.0;  
    DofCoords(19,0) = -1.0;   DofCoords(19,1) =  0.0; DofCoords(19,2) =  1.0;  

    //     DofCoords(20,0) =  0.0;   DofCoords(20,1) =  0.0; DofCoords(20,2) =  0.0;  

    //     DofCoords(21,0) =  0.0;   DofCoords(21,1) =  0.0; DofCoords(21,2) = -1.0;  
    //     DofCoords(22,0) =  0.0;   DofCoords(22,1) =  0.0; DofCoords(22,2) =  1.0;  
    //     DofCoords(23,0) = -1.0;   DofCoords(23,1) =  0.0; DofCoords(23,2) =  0.0;  
    //     DofCoords(24,0) =  1.0;   DofCoords(24,1) =  0.0; DofCoords(24,2) =  0.0;  
    //     DofCoords(25,0) =  0.0;   DofCoords(25,1) = -1.0; DofCoords(25,2) =  0.0;  
    //     DofCoords(26,0) =  0.0;   DofCoords(26,1) =  1.0; DofCoords(26,2) =  0.0;  
  }

}// namespace Intrepid
#endif
