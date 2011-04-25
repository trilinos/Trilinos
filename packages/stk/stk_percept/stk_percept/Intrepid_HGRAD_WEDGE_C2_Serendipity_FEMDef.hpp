#ifndef INTREPID_HGRAD_WEDGE_C2_SERENDIPITY_FEMDEF_HPP
#define INTREPID_HGRAD_WEDGE_C2_SERENDIPITY_FEMDEF_HPP
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

/** \file   Intrepid_HGRAD_WEDGE_C2_Serendipity_FEMDef.hpp
    \brief  Definition file for incomplete bi-quadratic FEM basis functions for H(grad) functions on WEDGE cells.
    \author Created by P. Bochev and D. Ridzal. (serendipity version by S. R. Kennon, srkenno@sandia.gov)
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_WEDGE_C2_Serendipity_FEM<Scalar, ArrayScalar>::Basis_HGRAD_WEDGE_C2_Serendipity_FEM()
  {
    this -> basisCardinality_  = 15;
    this -> basisDegree_       = 2;    
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Wedge<6> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }
  
  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_WEDGE_C2_Serendipity_FEM<Scalar, ArrayScalar>::initializeTags() {
  
    // Basis-dependent intializations
    int tagSize  = 4;        // size of DoF tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

    // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
    int tags[]  = { 0, 0, 0, 1,
                    0, 1, 0, 1,
                    0, 2, 0, 1,
                    0, 3, 0, 1,
                    0, 4, 0, 1,
                    0, 5, 0, 1,
                    1, 0, 0, 1,
                    1, 1, 0, 1,
                    1, 2, 0, 1,
                    1, 6, 0, 1,
                    1, 7, 0, 1,
                    1, 8, 0, 1,
                    1, 3, 0, 1,
                    1, 4, 0, 1,
                    1, 5, 0, 1
                    // ,2, 0, 0, 1,
                    //  2, 1, 0, 1,
                    //  2, 2, 0, 1
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
  void Basis_HGRAD_WEDGE_C2_Serendipity_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &    outputValues,
                                                                            const ArrayScalar &  inputPoints,
                                                                            const EOperator      operatorType) const {
  
    // Verify arguments
    //#ifdef HAVE_INTREPID_DEBUG
    Intrepid::getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
                                                        inputPoints,
                                                        operatorType,
                                                        this -> getBaseCellTopology(),
                                                        this -> getCardinality() );
    //#endif
  
    // Number of evaluation points = dim 0 of inputPoints
    int dim0 = inputPoints.dimension(0);  
  
    // Temporaries: (x,y,z) coordinates of the evaluation point
    Scalar x = 0.0;                                    
    Scalar y = 0.0;   
    Scalar z = 0.0;
  
    switch (operatorType) {
    
    case OPERATOR_VALUE:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0, 0);
        y = inputPoints(i0, 1);
        z = inputPoints(i0, 2);
        
        // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
        //
        // triangle C2
        //         outputValues(0, i0) = (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        //         outputValues(1, i0) = x*(2.0*x - 1.0);
        //         outputValues(2, i0) = y*(2.0*y - 1.0);

        //         outputValues(3, i0) = -4.0*x*(x + y - 1.0);
        //         outputValues(4, i0) =  4.0*x*y;
        //         outputValues(5, i0) = -4.0*y*(x + y - 1.0);

        double z0 = (1.0 - z)/2.0;
        double z1 = (1.0 + z)/2.0;
        double zq0 = (-z)*z0;
        double zq1 = ( z)*z1;
        double zqh = 4.0*z0*z1;

        //--- 0,1,2
        outputValues(0, i0) = zq0* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        outputValues(1, i0) = zq0* x*(2.0*x - 1.0);
        outputValues(2, i0) = zq0* y*(2.0*y - 1.0);

        //--- 3,4,5
        // same as 0,1,2 with zq1 for zq0
        outputValues(3, i0) = zq1* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        outputValues(4, i0) = zq1* x*(2.0*x - 1.0);
        outputValues(5, i0) = zq1* y*(2.0*y - 1.0);

        //--- 6,7,8
        outputValues(6, i0) = z0* (-4.0*x*(x + y - 1.0));
        outputValues(7, i0) = z0* ( 4.0*x*y);
        outputValues(8, i0) = z0* (-4.0*y*(x + y - 1.0));


        //--- 9,10,11
        // same as 0,1,2 with zqh for zq0
        outputValues(9, i0)  = zqh* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        outputValues(10, i0) = zqh* x*(2.0*x - 1.0);
        outputValues(11, i0) = zqh* y*(2.0*y - 1.0);

        //--- 12,13,14
        // same as 6,7,8 with z1 for z0
        outputValues(12, i0) = z1* (-4.0*x*(x + y - 1.0));
        outputValues(13, i0) = z1* ( 4.0*x*y);
        outputValues(14, i0) = z1* (-4.0*y*(x + y - 1.0));

      }
      break;
      
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        z = inputPoints(i0,2);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        double z0 = (1.0 - z)/2.0;
        double z1 = (1.0 + z)/2.0;
        double zq0 = (-z)*z0;
        double zq1 = ( z)*z1;
        double zqh = 4.0*z0*z1;

        double dz0 = -0.5;
        double dz1 = 0.5; 
        double dzq0 = (-1.0)*z0 + (-z)*dz0;
        double dzq1 = ( 1.0)*z1 + ( z)*dz1;
        double dzqh = 4.0*(dz0*z1 + z0*dz1);
        
        //outputValues(0, i0) = zq0* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        outputValues(0, i0, 0) = zq0* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        outputValues(0, i0, 1) = zq0* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        outputValues(0, i0, 2) = dzq0* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);

        //outputValues(1, i0) = zq0* x*(2.0*x - 1.0);
        outputValues(1, i0, 0) = zq0* ( 4.0*x - 1.0 );
        outputValues(1, i0, 1) = 0.0;
        outputValues(1, i0, 2) = dzq0* x*(2.0*x - 1.0);

        //outputValues(2, i0) = zq0* y*(2.0*y - 1.0);
        outputValues(2, i0, 0) = 0.0;
        outputValues(2, i0, 1) = zq0* (4.0*y - 1.0);
        outputValues(2, i0, 2) = dzq0* y*(2.0*y - 1.0);

        //outputValues(3, i0) = zq1* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        outputValues(3, i0, 0) = zq1* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        outputValues(3, i0, 1) = zq1* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        outputValues(3, i0, 2) = dzq1* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);

        //outputValues(4, i0) = zq1* x*(2.0*x - 1.0);
        outputValues(4, i0, 0) = zq1* ( 4.0*x - 1.0);
        outputValues(4, i0, 1) = 0.0;
        outputValues(4, i0, 2) = dzq1* x*(2.0*x - 1.0);

        //outputValues(5, i0) = zq1* y*(2.0*y - 1.0);
        outputValues(5, i0, 0) = 0.0;
        outputValues(5, i0, 1) = zq1* (4.0*y - 1.0);
        outputValues(5, i0, 2) = dzq1* y*(2.0*y - 1.0);

        //outputValues(6, i0) = z0* (-4.0*x*(x + y - 1.0));
        outputValues(6, i0, 0) = z0* (-4.0*(2.0*x + y - 1.0) );
        outputValues(6, i0, 1) = z0* (-4.0*x*(1.0));
        outputValues(6, i0, 2) = dz0* (-4.0*x*(x + y - 1.0));

        //outputValues(7, i0) = z0* ( 4.0*x*y);
        outputValues(7, i0, 0) = z0* ( 4.0*y);
        outputValues(7, i0, 1) = z0* ( 4.0*x);
        outputValues(7, i0, 2) = dz0* ( 4.0*x*y);

        //outputValues(8, i0) = z0* (-4.0*y*(x + y - 1.0));
        outputValues(8, i0, 0) = z0* (-4.0*y*(1.0));
        outputValues(8, i0, 1) = z0* (-4.0*(x + 2.0*y - 1.0));
        outputValues(8, i0, 2) = dz0* (-4.0*y*(x + y - 1.0));

        //outputValues(9, i0) = zqh* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        outputValues(9, i0, 0) = zqh* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        outputValues(9, i0, 1) = zqh* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        outputValues(9, i0, 2) = dzqh* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);

        //outputValues(10, i0) = zqh* x*(2.0*x - 1.0);
        outputValues(10, i0, 0) = zqh* ( 1.0*(2.0*x - 1.0) + x*(2.0) );
        outputValues(10, i0, 1) = 0.0;
        outputValues(10, i0, 2) = dzqh* x*(2.0*x - 1.0);

        //outputValues(11, i0) = zqh* y*(2.0*y - 1.0);
        outputValues(11, i0, 0) = 0.0;
        outputValues(11, i0, 1) = zqh* (4.0*y - 1.0);
        outputValues(11, i0, 2) = dzqh* y*(2.0*y - 1.0);

        //outputValues(12, i0) = z1* (-4.0*x*(x + y - 1.0));
        outputValues(12, i0, 0) = z1* (-4.0*(2.0*x + y - 1.0));
        outputValues(12, i0, 1) = z1* (-4.0*x*(1.0));
        outputValues(12, i0, 2) = dz1* (-4.0*x*(x + y - 1.0));

        //outputValues(13, i0) = z1* ( 4.0*x*y);
        outputValues(13, i0, 0) = z1* ( 4.0*y);
        outputValues(13, i0, 1) = z1* ( 4.0*x);
        outputValues(13, i0, 2) = dz1* ( 4.0*x*y);

        //outputValues(14, i0) = z1* (-4.0*y*(x + y - 1.0));
        outputValues(14, i0, 0) = z1* (-4.0*y*(1.0));
        outputValues(14, i0, 1) = z1* (-4.0*(x + 2.0*y - 1.0));
        outputValues(14, i0, 2) = dz1* (-4.0*y*(x + y - 1.0));
        
      }
      break;
      
    case OPERATOR_CURL:
      TEST_FOR_EXCEPTION( (operatorType == OPERATOR_CURL), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_WEDGE_C2_Serendipity_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
      break;
      
    case OPERATOR_DIV:
      TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_WEDGE_C2_Serendipity_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
      break;
      
    case OPERATOR_D2:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        z = inputPoints(i0,2);

        double z0 = (1.0 - z)/2.0;
        double z1 = (1.0 + z)/2.0;
        double zq0 = (-z)*z0;
        double zq1 = ( z)*z1;
        double zqh = 4.0*z0*z1;

        double dz0 = -0.5;
        double dz1 = 0.5; 
        double dzq0 = (-1.0)*z0 + (-z)*dz0;
        double dzq1 = ( 1.0)*z1 + ( z)*dz1;
        double dzqh = 4.0*(dz0*z1 + z0*dz1);

        double d2z0 = 0.0;
        double d2z1 = 0.0;
        double d2zq0 = (-1.0)*dz0 + (-1.0)*dz0;
        double d2zq1 = ( 1.0)*dz1 + ( 1.0)*dz1;
        double d2zqh = 8.0*( dz0*dz1);
        
        //---- 0,1,2
        //outputValues(0, i0) = zq0* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        //         outputValues(0, i0, 0) = zq0* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        //         outputValues(0, i0, 1) = zq0* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        //         outputValues(0, i0, 2) = dzq0* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        outputValues(0, i0, 0) = zq0* ( 4.0 );
        outputValues(0, i0, 1) = zq0* ( 4.0 );
        outputValues(0, i0, 2) = dzq0* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        outputValues(0, i0, 3) = zq0* ( (1.0)*( 2.0) + (1.0)*(2.0) );
        outputValues(0, i0, 4) = dzq0* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        outputValues(0, i0, 5) = d2zq0* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);

        //outputValues(1, i0) = zq0* x*(2.0*x - 1.0);
        //         outputValues(1, i0, 0) = zq0* ( 4.0*x - 1.0 );
        //         outputValues(1, i0, 1) = 0.0;
        //         outputValues(1, i0, 2) = dzq0* x*(2.0*x - 1.0);
        outputValues(1, i0, 0) = zq0* ( 4.0 );
        outputValues(1, i0, 1) = 0.0;
        outputValues(1, i0, 2) = dzq0* ( 4.0*x - 1.0 );

        outputValues(1, i0, 3) = 0.0;
        outputValues(1, i0, 4) = 0.0;

        outputValues(1, i0, 5) = d2zq0* x*(2.0*x - 1.0);

        //outputValues(2, i0) = zq0* y*(2.0*y - 1.0);
        //         outputValues(2, i0, 0) = 0.0;
        //         outputValues(2, i0, 1) = zq0* (4.0*y - 1.0);
        //         outputValues(2, i0, 2) = dzq0* y*(2.0*y - 1.0);

        outputValues(2, i0, 0) = 0.0;
        outputValues(2, i0, 1) = 0.0;
        outputValues(2, i0, 2) = 0.0;

        outputValues(2, i0, 3) = zq0* (4.0);
        outputValues(2, i0, 4) = dzq0* (4.0*y - 1.0);

        outputValues(2, i0, 5) = d2zq0* y*(2.0*y - 1.0);

        //---- 3,4,5
        // same as 0,1,2 with zq1 for zq0
        //outputValues(3, i0) = zq1* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        //         outputValues(3, i0, 0) = zq1* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        //         outputValues(3, i0, 1) = zq1* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        //         outputValues(3, i0, 2) = dzq1* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        outputValues(3, i0, 0) = zq1* ( 4.0 );
        outputValues(3, i0, 1) = zq1* ( 4.0 );
        outputValues(3, i0, 2) = dzq1* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        outputValues(3, i0, 3) = zq1* ( (1.0)*( 2.0) + (1.0)*(2.0) );
        outputValues(3, i0, 4) = dzq1* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        outputValues(3, i0, 5) = d2zq1* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);

        //outputValues(4, i0) = zq1* x*(2.0*x - 1.0);
        //         outputValues(4, i0, 0) = zq1* ( 4.0*x - 1.0 );
        //         outputValues(4, i0, 1) = 0.0;
        //         outputValues(4, i0, 2) = dzq1* x*(2.0*x - 1.0);
        outputValues(4, i0, 0) = zq1* ( 4.0 );
        outputValues(4, i0, 1) = 0.0;
        outputValues(4, i0, 2) = dzq1* ( 4.0*x - 1.0 );

        outputValues(4, i0, 3) = 0.0;
        outputValues(4, i0, 4) = 0.0;

        outputValues(4, i0, 5) = d2zq1* x*(2.0*x - 1.0);

        //outputValues(5, i0) = zq1* y*(2.0*y - 1.0);
        //         outputValues(5, i0, 0) = 0.0;
        //         outputValues(5, i0, 1) = zq1* (4.0*y - 1.0);
        //         outputValues(5, i0, 2) = dzq1* y*(2.0*y - 1.0);

        outputValues(5, i0, 0) = 0.0;
        outputValues(5, i0, 1) = 0.0;
        outputValues(5, i0, 2) = 0.0;

        outputValues(5, i0, 3) = zq1* (4.0);
        outputValues(5, i0, 4) = dzq1* (4.0*y - 1.0);

        outputValues(5, i0, 5) = d2zq1* y*(2.0*y - 1.0);


        //--- 6,7,8
        //outputValues(6, i0) = z0* (-4.0*x*(x + y - 1.0));
        //         outputValues(6, i0, 0) = z0* (-4.0*(2.0*x + y - 1.0) );
        //         outputValues(6, i0, 1) = z0* (-4.0*x*(1.0));
        //         outputValues(6, i0, 2) = dz0* (-4.0*x*(x + y - 1.0));
        
        outputValues(6, i0, 0) = z0* (-8.0);
        outputValues(6, i0, 1) = z0* (-4.0);
        outputValues(6, i0, 2) = dz0* (-4.0*(2.0*x + y - 1.0) );

        outputValues(6, i0, 3) = 0.0;
        outputValues(6, i0, 4) = dz0* (-4.0*x*(1.0));

        outputValues(6, i0, 5) = d2z0* (-4.0*x*(x + y - 1.0));

        //outputValues(7, i0) = z0* ( 4.0*x*y);
        //         outputValues(7, i0, 0) = z0* ( 4.0*y);
        //         outputValues(7, i0, 1) = z0* ( 4.0*x);
        //         outputValues(7, i0, 2) = dz0* ( 4.0*x*y);

        outputValues(7, i0, 0) = 0.0;
        outputValues(7, i0, 1) = z0* ( 4.0 );
        outputValues(7, i0, 2) = dz0* ( 4.0*y);

        outputValues(7, i0, 3) = 0.0;
        outputValues(7, i0, 4) = dz0* ( 4.0*x);

        outputValues(7, i0, 5) = d2z0* ( 4.0*x*y);

        //outputValues(8, i0) = z0* (-4.0*y*(x + y - 1.0));
        //         outputValues(8, i0, 0) = z0* (-4.0*y*(1.0));
        //         outputValues(8, i0, 1) = z0* (-4.0*(x + 2.0*y - 1.0));
        //         outputValues(8, i0, 2) = dz0* (-4.0*y*(x + y - 1.0));

        outputValues(8, i0, 0) = 0.0;
        outputValues(8, i0, 1) = z0* (-4.0);
        outputValues(8, i0, 2) = dz0* (-4.0*y*(1.0));

        outputValues(8, i0, 3) = z0* (-8.0);
        outputValues(8, i0, 4) = dz0* (-4.0*(x + 2.0*y - 1.0));

        outputValues(8, i0, 5) = d2z0* (-4.0*y*(x + y - 1.0));

        //---- 9,10,11
        // same as 0,1,2 with zqh for zq0

        //outputValues(9, i0) = zqh* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        //         outputValues(9, i0, 0) = zqh* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        //         outputValues(9, i0, 1) = zqh* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        //         outputValues(9, i0, 2) = dzqh* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);

        outputValues(9, i0, 0) = zqh* ( 4.0 );
        outputValues(9, i0, 1) = zqh* ( 4.0 );
        outputValues(9, i0, 2) = dzqh* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        outputValues(9, i0, 3) = zqh* ( (1.0)*( 2.0) + (1.0)*(2.0) );
        outputValues(9, i0, 4) = dzqh* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        outputValues(9, i0, 5) = d2zqh* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);

        //outputValues(10, i0) = zqh* x*(2.0*x - 1.0);
        //         outputValues(10, i0, 0) = zqh* ( 4.0*x - 1.0 );
        //         outputValues(10, i0, 1) = 0.0;
        //         outputValues(10, i0, 2) = dzqh* x*(2.0*x - 1.0);
        outputValues(10, i0, 0) = zqh* ( 4.0 );
        outputValues(10, i0, 1) = 0.0;
        outputValues(10, i0, 2) = dzqh* ( 4.0*x - 1.0 );

        outputValues(10, i0, 3) = 0.0;
        outputValues(10, i0, 4) = 0.0;

        outputValues(10, i0, 5) = d2zqh* x*(2.0*x - 1.0);

        //outputValues(11, i0) = zqh* y*(2.0*y - 1.0);
        //         outputValues(11, i0, 0) = 0.0;
        //         outputValues(11, i0, 1) = zqh* (4.0*y - 1.0);
        //         outputValues(11, i0, 2) = dzqh* y*(2.0*y - 1.0);

        outputValues(11, i0, 0) = 0.0;
        outputValues(11, i0, 1) = 0.0;
        outputValues(11, i0, 2) = 0.0;

        outputValues(11, i0, 3) = zqh* (4.0);
        outputValues(11, i0, 4) = dzqh* (4.0*y - 1.0);

        outputValues(11, i0, 5) = d2zqh* y*(2.0*y - 1.0);

        //--- 12,13,14
        // same as 6,7,8 with z1 for z0
        //outputValues(12, i0) = z1* (-4.0*x*(x + y - 1.0));
        //         outputValues(12, i0, 0) = z1* (-4.0*(2.0*x + y - 1.0) );
        //         outputValues(12, i0, 1) = z1* (-4.0*x*(1.0));
        //         outputValues(12, i0, 2) = dz1* (-4.0*x*(x + y - 1.0));
        
        outputValues(12, i0, 0) = z1* (-8.0);
        outputValues(12, i0, 1) = z1* (-4.0);
        outputValues(12, i0, 2) = dz1* (-4.0*(2.0*x + y - 1.0) );

        outputValues(12, i0, 3) = 0.0;
        outputValues(12, i0, 4) = dz1* (-4.0*x*(1.0));

        outputValues(12, i0, 5) = d2z1* (-4.0*x*(x + y - 1.0));

        //outputValues(13, i0) = z1* ( 4.0*x*y);
        //         outputValues(13, i0, 0) = z1* ( 4.0*y);
        //         outputValues(13, i0, 1) = z1* ( 4.0*x);
        //         outputValues(13, i0, 2) = dz1* ( 4.0*x*y);

        outputValues(13, i0, 0) = 0.0;
        outputValues(13, i0, 1) = z1* ( 4.0 );
        outputValues(13, i0, 2) = dz1* ( 4.0*y);

        outputValues(13, i0, 3) = 0.0;
        outputValues(13, i0, 4) = dz1* ( 4.0*x);

        outputValues(13, i0, 5) = d2z1* ( 4.0*x*y);

        //outputValues(14, i0) = z1* (-4.0*y*(x + y - 1.0));
        //         outputValues(14, i0, 0) = z1* (-4.0*y*(1.0));
        //         outputValues(14, i0, 1) = z1* (-4.0*(x + 2.0*y - 1.0));
        //         outputValues(14, i0, 2) = dz1* (-4.0*y*(x + y - 1.0));

        outputValues(14, i0, 0) = 0.0;
        outputValues(14, i0, 1) = z1* (-4.0);
        outputValues(14, i0, 2) = dz1* (-4.0*y*(1.0));

        outputValues(14, i0, 3) = z1* (-8.0);
        outputValues(14, i0, 4) = dz1* (-4.0*(x + 2.0*y - 1.0));

        outputValues(14, i0, 5) = d2z1* (-4.0*y*(x + y - 1.0));


      }
      break;
      
    case OPERATOR_D3:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        z = inputPoints(i0,2);

        // 1st
        // x
        // y
        // z

        // 2nd
        // xx, xy, xz
        //     yy, yz
        //         zz

        // 3rd
        // xxx, xxy, [xxz],z, (2)
        // xyy, [xyz],z, (4)
        // xzz

        //  yyy, [yyz],z (7)
        //  yzz

        //  zzz


        // 4th
        // xxxx, xxxy, xxxz
        // xxyy, xxyz
        // xxzz, (5)

        // xyyy, xyyz,
        // xyzz, (8)

        // xzzz

        //  yyyy, yyyz,
        //  yyzz, (12)

        //  yzzz

        //  zzzz

        double z0 = (1.0 - z)/2.0;
        double z1 = (1.0 + z)/2.0;
        //double zq0 = (-z)*z0;
        //double zq1 = ( z)*z1;
        //double zqh = 4.0*z0*z1;

        double dz0 = -0.5;
        double dz1 = 0.5; 
        double dzq0 = (-1.0)*z0 + (-z)*dz0;
        double dzq1 = ( 1.0)*z1 + ( z)*dz1;
        double dzqh = 4.0*(dz0*z1 + z0*dz1);

        double d2z0 = 0.0;
        double d2z1 = 0.0;
        double d2zq0 = (-1.0)*dz0 + (-1.0)*dz0;
        double d2zq1 = ( 1.0)*dz1 + ( 1.0)*dz1;
        double d2zqh = 8.0*( dz0*dz1);

        double d3z0 = 0.0;
        double d3z1 = 0.0;
        double d3zq0 = 0.0;
        double d3zq1 = 0.0;
        double d3zqh = 0.0;
        
        //---- 0,1,2
        //outputValues(0, i0) = zq0* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        //         outputValues(0, i0, 0) = zq0* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        //         outputValues(0, i0, 1) = zq0* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        //         outputValues(0, i0, 2) = dzq0* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        //         outputValues(0, i0, 0) = zq0* ( 4.0 );
        //         outputValues(0, i0, 1) = zq0* ( 4.0 );
        //         outputValues(0, i0, 2) = dzq0* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        //         outputValues(0, i0, 3) = zq0* ( (1.0)*( 2.0) + (1.0)*(2.0) );
        //         outputValues(0, i0, 4) = dzq0* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        //         outputValues(0, i0, 5) = d2zq0* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        outputValues(0, i0, 0) = 0.0;
        outputValues(0, i0, 1) = 0.0;
        outputValues(0, i0, 2) = dzq0* ( 4.0 );

        outputValues(0, i0, 3) = 0.0;
        outputValues(0, i0, 4) = dzq0* ( 4.0 );

        outputValues(0, i0, 5) = d2zq0* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        outputValues(0, i0, 6) = 0.0;
        outputValues(0, i0, 7) = dzq0* ( (1.0)*( 2.0) + (1.0)*(2.0) );

        outputValues(0, i0, 8) = d2zq0* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        outputValues(0, i0, 9) = d3zq0* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);


        //outputValues(1, i0) = zq0* x*(2.0*x - 1.0);
        //         outputValues(1, i0, 0) = zq0* ( 4.0*x - 1.0 );
        //         outputValues(1, i0, 1) = 0.0;
        //         outputValues(1, i0, 2) = dzq0* x*(2.0*x - 1.0);
        //         outputValues(1, i0, 0) = zq0* ( 4.0 );
        //         outputValues(1, i0, 1) = 0.0;
        //         outputValues(1, i0, 2) = dzq0* ( 4.0*x - 1.0 );

        //         outputValues(1, i0, 3) = 0.0;
        //         outputValues(1, i0, 4) = 0.0;

        //         outputValues(1, i0, 5) = d2zq0* x*(2.0*x - 1.0);

        outputValues(1, i0, 0) = 0.0;
        outputValues(1, i0, 1) = 0.0;
        outputValues(1, i0, 2) = dzq0* ( 4.0 );

        outputValues(1, i0, 3) = 0.0;
        outputValues(1, i0, 4) = 0.0;

        outputValues(1, i0, 5) = d2zq0* ( 4.0*x - 1.0 );


        outputValues(1, i0, 6) = 0.0;
        outputValues(1, i0, 7) = 0.0;

        outputValues(1, i0, 8) = 0.0;

        outputValues(1, i0, 9) = d3zq0* x*(2.0*x - 1.0);

        //outputValues(2, i0) = zq0* y*(2.0*y - 1.0);
        //         outputValues(2, i0, 0) = 0.0;
        //         outputValues(2, i0, 1) = zq0* (4.0*y - 1.0);
        //         outputValues(2, i0, 2) = dzq0* y*(2.0*y - 1.0);

        //         outputValues(2, i0, 0) = 0.0;
        //         outputValues(2, i0, 1) = 0.0;
        //         outputValues(2, i0, 2) = 0.0;

        //         outputValues(2, i0, 3) = zq0* (4.0);
        //         outputValues(2, i0, 4) = dzq0* (4.0*y - 1.0);

        //         outputValues(2, i0, 5) = d2zq0* y*(2.0*y - 1.0);

        outputValues(2, i0, 0) = 0.0;
        outputValues(2, i0, 1) = 0.0;
        outputValues(2, i0, 2) = 0.0;

        outputValues(2, i0, 3) = 0.0;
        outputValues(2, i0, 4) = 0.0;

        outputValues(2, i0, 5) = 0.0;


        outputValues(2, i0, 6) = 0.0;
        outputValues(2, i0, 7) = dzq0* (4.0);

        outputValues(2, i0, 8) = d2zq0* (4.0*y - 1.0);

        outputValues(2, i0, 9) = d3zq0* y*(2.0*y - 1.0);

        //---- 3,4,5
        // same as 0,1,2 with zq1 for zq0
        //outputValues(3, i0) = zq1* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        //         outputValues(3, i0, 0) = zq1* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        //         outputValues(3, i0, 1) = zq1* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        //         outputValues(3, i0, 2) = dzq1* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        //         outputValues(3, i0, 0) = zq1* ( 4.0 );
        //         outputValues(3, i0, 1) = zq1* ( 4.0 );
        //         outputValues(3, i0, 2) = dzq1* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        //         outputValues(3, i0, 3) = zq1* ( (1.0)*( 2.0) + (1.0)*(2.0) );
        //         outputValues(3, i0, 4) = dzq1* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        //         outputValues(3, i0, 5) = d2zq1* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        outputValues(3, i0, 0) = 0.0;
        outputValues(3, i0, 1) = 0.0;
        outputValues(3, i0, 2) = dzq1* ( 4.0 );

        outputValues(3, i0, 3) = 0.0;
        outputValues(3, i0, 4) = dzq1* ( 4.0 );

        outputValues(3, i0, 5) = d2zq1* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        outputValues(3, i0, 6) = 0.0;
        outputValues(3, i0, 7) = dzq1* ( (1.0)*( 2.0) + (1.0)*(2.0) );

        outputValues(3, i0, 8) = d2zq1* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        outputValues(3, i0, 9) = d3zq1* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);


        //outputValues(4, i0) = zq1* x*(2.0*x - 1.0);
        //         outputValues(4, i0, 0) = zq1* ( 4.0*x - 1.0 );
        //         outputValues(4, i0, 1) = 0.0;
        //         outputValues(4, i0, 2) = dzq1* x*(2.0*x - 1.0);
        //         outputValues(4, i0, 0) = zq1* ( 4.0 );
        //         outputValues(4, i0, 1) = 0.0;
        //         outputValues(4, i0, 2) = dzq1* ( 4.0*x - 1.0 );

        //         outputValues(4, i0, 3) = 0.0;
        //         outputValues(4, i0, 4) = 0.0;

        //         outputValues(4, i0, 5) = d2zq1* x*(2.0*x - 1.0);

        outputValues(4, i0, 0) = 0.0;
        outputValues(4, i0, 1) = 0.0;
        outputValues(4, i0, 2) = dzq1* ( 4.0 );

        outputValues(4, i0, 3) = 0.0;
        outputValues(4, i0, 4) = 0.0;

        outputValues(4, i0, 5) = d2zq1* ( 4.0*x - 1.0 );


        outputValues(4, i0, 6) = 0.0;
        outputValues(4, i0, 7) = 0.0;

        outputValues(4, i0, 8) = 0.0;

        outputValues(4, i0, 9) = d3zq1* x*(2.0*x - 1.0);

        //outputValues(5, i0) = zq1* y*(2.0*y - 1.0);
        //         outputValues(5, i0, 0) = 0.0;
        //         outputValues(5, i0, 1) = zq1* (4.0*y - 1.0);
        //         outputValues(5, i0, 2) = dzq1* y*(2.0*y - 1.0);

        //         outputValues(5, i0, 0) = 0.0;
        //         outputValues(5, i0, 1) = 0.0;
        //         outputValues(5, i0, 2) = 0.0;

        //         outputValues(5, i0, 3) = zq1* (4.0);
        //         outputValues(5, i0, 4) = dzq1* (4.0*y - 1.0);

        //         outputValues(5, i0, 5) = d2zq1* y*(2.0*y - 1.0);

        outputValues(5, i0, 0) = 0.0;
        outputValues(5, i0, 1) = 0.0;
        outputValues(5, i0, 2) = 0.0;

        outputValues(5, i0, 3) = 0.0;
        outputValues(5, i0, 4) = 0.0;

        outputValues(5, i0, 5) = 0.0;


        outputValues(5, i0, 6) = 0.0;
        outputValues(5, i0, 7) = dzq1* (4.0);

        outputValues(5, i0, 8) = d2zq1* (4.0*y - 1.0);

        outputValues(5, i0, 9) = d3zq1* y*(2.0*y - 1.0);



        //--- 6,7,8
        //outputValues(6, i0) = z0* (-4.0*x*(x + y - 1.0));
        //         outputValues(6, i0, 0) = z0* (-4.0*(2.0*x + y - 1.0) );
        //         outputValues(6, i0, 1) = z0* (-4.0*x*(1.0));
        //         outputValues(6, i0, 2) = dz0* (-4.0*x*(x + y - 1.0));
        
        //         outputValues(6, i0, 0) = z0* (-8.0);
        //         outputValues(6, i0, 1) = z0* (-4.0);
        //         outputValues(6, i0, 2) = dz0* (-4.0*(2.0*x + y - 1.0) );

        //         outputValues(6, i0, 3) = 0.0;
        //         outputValues(6, i0, 4) = dz0* (-4.0*x*(1.0));

        //         outputValues(6, i0, 5) = d2z0* (-4.0*x*(x + y - 1.0));

        outputValues(6, i0, 0) = 0.0;
        outputValues(6, i0, 1) = 0.0;
        outputValues(6, i0, 2) = dz0* (-8.0);

        outputValues(6, i0, 3) = 0.0;
        outputValues(6, i0, 4) = dz0* (-4.0);

        outputValues(6, i0, 5) = d2z0* (-4.0*(2.0*x + y - 1.0) );


        outputValues(6, i0, 6) = 0.0;
        outputValues(6, i0, 7) = 0.0;

        outputValues(6, i0, 8) = d2z0* (-4.0*x*(1.0));

        outputValues(6, i0, 9) = d3z0* (-4.0*x*(x + y - 1.0));

        //outputValues(7, i0) = z0* ( 4.0*x*y);
        //         outputValues(7, i0, 0) = z0* ( 4.0*y);
        //         outputValues(7, i0, 1) = z0* ( 4.0*x);
        //         outputValues(7, i0, 2) = dz0* ( 4.0*x*y);

        //         outputValues(7, i0, 0) = 0.0;
        //         outputValues(7, i0, 1) = z0* ( 4.0 );
        //         outputValues(7, i0, 2) = dz0* ( 4.0*y);

        //         outputValues(7, i0, 3) = 0.0;
        //         outputValues(7, i0, 4) = dz0* ( 4.0*x);

        //         outputValues(7, i0, 5) = d2z0* ( 4.0*x*y);

        outputValues(7, i0, 0) = 0.0;
        outputValues(7, i0, 1) = 0.0;
        outputValues(7, i0, 2) = 0.0;

        outputValues(7, i0, 3) = 0.0;
        outputValues(7, i0, 4) = dz0* ( 4.0 );

        outputValues(7, i0, 5) = d2z0* ( 4.0*y);


        outputValues(7, i0, 6) = 0.0;
        outputValues(7, i0, 7) = 0.0;

        outputValues(7, i0, 8) = d2z0* ( 4.0*x);

        outputValues(7, i0, 9) = d3z0* ( 4.0*x*y);

        //outputValues(8, i0) = z0* (-4.0*y*(x + y - 1.0));
        //         outputValues(8, i0, 0) = z0* (-4.0*y*(1.0));
        //         outputValues(8, i0, 1) = z0* (-4.0*(x + 2.0*y - 1.0));
        //         outputValues(8, i0, 2) = dz0* (-4.0*y*(x + y - 1.0));

        //         outputValues(8, i0, 0) = 0.0;
        //         outputValues(8, i0, 1) = z0* (-4.0);
        //         outputValues(8, i0, 2) = dz0* (-4.0*y*(1.0));

        //         outputValues(8, i0, 3) = z0* (-8.0);
        //         outputValues(8, i0, 4) = dz0* (-4.0*(x + 2.0*y - 1.0));

        //         outputValues(8, i0, 5) = d2z0* (-4.0*y*(x + y - 1.0));

        outputValues(8, i0, 0) = 0.0;
        outputValues(8, i0, 1) = 0.0;
        outputValues(8, i0, 2) = 0.0;

        outputValues(8, i0, 3) = 0.0;
        outputValues(8, i0, 4) = dz0* (-4.0);

        outputValues(8, i0, 5) = d2z0* (-4.0*y*(1.0));


        outputValues(8, i0, 6) = 0.0;
        outputValues(8, i0, 7) = dz0* (-8.0);

        outputValues(8, i0, 8) = d2z0* (-4.0*(x + 2.0*y - 1.0));

        outputValues(8, i0, 9) = d3z0* (-4.0*y*(x + y - 1.0));

        //---- 9,10,11
        // same as 0,1,2 with zqh for zq0

        //outputValues(9, i0) = zqh* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        //         outputValues(9, i0, 0) = zqh* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        //         outputValues(9, i0, 1) = zqh* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );
        //         outputValues(9, i0, 2) = dzqh* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        //         outputValues(9, i0, 0) = zqh* ( 4.0 );
        //         outputValues(9, i0, 1) = zqh* ( 4.0 );
        //         outputValues(9, i0, 2) = dzqh* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        //         outputValues(9, i0, 3) = zqh* ( (1.0)*( 2.0) + (1.0)*(2.0) );
        //         outputValues(9, i0, 4) = dzqh* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        //         outputValues(9, i0, 5) = d2zqh* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
        outputValues(9, i0, 0) = 0.0;
        outputValues(9, i0, 1) = 0.0;
        outputValues(9, i0, 2) = dzqh* ( 4.0 );

        outputValues(9, i0, 3) = 0.0;
        outputValues(9, i0, 4) = dzqh* ( 4.0 );

        outputValues(9, i0, 5) = d2zqh* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        outputValues(9, i0, 6) = 0.0;
        outputValues(9, i0, 7) = dzqh* ( (1.0)*( 2.0) + (1.0)*(2.0) );

        outputValues(9, i0, 8) = d2zqh* ( (1.0)*(2.0*x + 2.0*y - 1.0) + (x + y - 1.0)*(2.0) );

        outputValues(9, i0, 9) = d3zqh* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);


        //outputValues(10, i0) = zqh* x*(2.0*x - 1.0);
        //         outputValues(10, i0, 0) = zqh* ( 4.0*x - 1.0 );
        //         outputValues(10, i0, 1) = 0.0;
        //         outputValues(10, i0, 2) = dzqh* x*(2.0*x - 1.0);
        //         outputValues(10, i0, 0) = zqh* ( 4.0 );
        //         outputValues(10, i0, 1) = 0.0;
        //         outputValues(10, i0, 2) = dzqh* ( 4.0*x - 1.0 );

        //         outputValues(10, i0, 3) = 0.0;
        //         outputValues(10, i0, 4) = 0.0;

        //         outputValues(10, i0, 5) = d2zqh* x*(2.0*x - 1.0);

        outputValues(10, i0, 0) = 0.0;
        outputValues(10, i0, 1) = 0.0;
        outputValues(10, i0, 2) = dzqh* ( 4.0 );

        outputValues(10, i0, 3) = 0.0;
        outputValues(10, i0, 4) = 0.0;

        outputValues(10, i0, 5) = d2zqh* ( 4.0*x - 1.0 );


        outputValues(10, i0, 6) = 0.0;
        outputValues(10, i0, 7) = 0.0;

        outputValues(10, i0, 8) = 0.0;

        outputValues(10, i0, 9) = d3zqh* x*(2.0*x - 1.0);

        //outputValues(11, i0) = zqh* y*(2.0*y - 1.0);
        //         outputValues(11, i0, 0) = 0.0;
        //         outputValues(11, i0, 1) = zqh* (4.0*y - 1.0);
        //         outputValues(11, i0, 2) = dzqh* y*(2.0*y - 1.0);

        //         outputValues(11, i0, 0) = 0.0;
        //         outputValues(11, i0, 1) = 0.0;
        //         outputValues(11, i0, 2) = 0.0;

        //         outputValues(11, i0, 3) = zqh* (4.0);
        //         outputValues(11, i0, 4) = dzqh* (4.0*y - 1.0);

        //         outputValues(11, i0, 5) = d2zqh* y*(2.0*y - 1.0);

        outputValues(11, i0, 0) = 0.0;
        outputValues(11, i0, 1) = 0.0;
        outputValues(11, i0, 2) = 0.0;

        outputValues(11, i0, 3) = 0.0;
        outputValues(11, i0, 4) = 0.0;

        outputValues(11, i0, 5) = 0.0;


        outputValues(11, i0, 6) = 0.0;
        outputValues(11, i0, 7) = dzqh* (4.0);

        outputValues(11, i0, 8) = d2zqh* (4.0*y - 1.0);

        outputValues(11, i0, 9) = d3zqh* y*(2.0*y - 1.0);


        //--- 12,13,14
        // same as 6,7,8 with z1 for z0

        //outputValues(12, i0) = z1* (-4.0*x*(x + y - 1.0));
        //         outputValues(12, i0, 0) = z1* (-4.0*(2.0*x + y - 1.0) );
        //         outputValues(12, i0, 1) = z1* (-4.0*x*(1.0));
        //         outputValues(12, i0, 2) = dz1* (-4.0*x*(x + y - 1.0));
        
        //         outputValues(12, i0, 0) = z1* (-8.0);
        //         outputValues(12, i0, 1) = z1* (-4.0);
        //         outputValues(12, i0, 2) = dz1* (-4.0*(2.0*x + y - 1.0) );

        //         outputValues(12, i0, 3) = 0.0;
        //         outputValues(12, i0, 4) = dz1* (-4.0*x*(1.0));

        //         outputValues(12, i0, 5) = d2z1* (-4.0*x*(x + y - 1.0));

        outputValues(12, i0, 0) = 0.0;
        outputValues(12, i0, 1) = 0.0;
        outputValues(12, i0, 2) = dz1* (-8.0);

        outputValues(12, i0, 3) = 0.0;
        outputValues(12, i0, 4) = dz1* (-4.0);

        outputValues(12, i0, 5) = d2z1* (-4.0*(2.0*x + y - 1.0) );


        outputValues(12, i0, 6) = 0.0;
        outputValues(12, i0, 7) = 0.0;

        outputValues(12, i0, 8) = d2z1* (-4.0*x*(1.0));

        outputValues(12, i0, 9) = d3z1* (-4.0*x*(x + y - 1.0));

        //outputValues(13, i0) = z1* ( 4.0*x*y);
        //         outputValues(13, i0, 0) = z1* ( 4.0*y);
        //         outputValues(13, i0, 1) = z1* ( 4.0*x);
        //         outputValues(13, i0, 2) = dz1* ( 4.0*x*y);

        //         outputValues(13, i0, 0) = 0.0;
        //         outputValues(13, i0, 1) = z1* ( 4.0 );
        //         outputValues(13, i0, 2) = dz1* ( 4.0*y);

        //         outputValues(13, i0, 3) = 0.0;
        //         outputValues(13, i0, 4) = dz1* ( 4.0*x);

        //         outputValues(13, i0, 5) = d2z1* ( 4.0*x*y);

        outputValues(13, i0, 0) = 0.0;
        outputValues(13, i0, 1) = 0.0;
        outputValues(13, i0, 2) = 0.0;

        outputValues(13, i0, 3) = 0.0;
        outputValues(13, i0, 4) = dz1* ( 4.0 );

        outputValues(13, i0, 5) = d2z1* ( 4.0*y);


        outputValues(13, i0, 6) = 0.0;
        outputValues(13, i0, 7) = 0.0;

        outputValues(13, i0, 8) = d2z1* ( 4.0*x);

        outputValues(13, i0, 9) = d3z1* ( 4.0*x*y);

        //outputValues(14, i0) = z1* (-4.0*y*(x + y - 1.0));
        //         outputValues(14, i0, 0) = z1* (-4.0*y*(1.0));
        //         outputValues(14, i0, 1) = z1* (-4.0*(x + 2.0*y - 1.0));
        //         outputValues(14, i0, 2) = dz1* (-4.0*y*(x + y - 1.0));

        //         outputValues(14, i0, 0) = 0.0;
        //         outputValues(14, i0, 1) = z1* (-4.0);
        //         outputValues(14, i0, 2) = dz1* (-4.0*y*(1.0));

        //         outputValues(14, i0, 3) = z1* (-8.0);
        //         outputValues(14, i0, 4) = dz1* (-4.0*(x + 2.0*y - 1.0));

        //         outputValues(14, i0, 5) = d2z1* (-4.0*y*(x + y - 1.0));

        outputValues(14, i0, 0) = 0.0;
        outputValues(14, i0, 1) = 0.0;
        outputValues(14, i0, 2) = 0.0;

        outputValues(14, i0, 3) = 0.0;
        outputValues(14, i0, 4) = dz1* (-4.0);

        outputValues(14, i0, 5) = d2z1* (-4.0*y*(1.0));


        outputValues(14, i0, 6) = 0.0;
        outputValues(14, i0, 7) = dz1* (-8.0);

        outputValues(14, i0, 8) = d2z1* (-4.0*(x + 2.0*y - 1.0));

        outputValues(14, i0, 9) = d3z1* (-4.0*y*(x + y - 1.0));

        
      }
      break;
      
    case OPERATOR_D4:
      {
        // There are only few constant non-zero entries. Initialize by zero and then assign non-zero entries.
        int DkCardinality = Intrepid::getDkCardinality(operatorType, this -> basisCellTopology_.getDimension() );
        for(int dofOrd = 0; dofOrd < this -> basisCardinality_; dofOrd++) {
          for (int i0 = 0; i0 < dim0; i0++) {
            for(int dkOrd = 0; dkOrd < DkCardinality; dkOrd++){
              outputValues(dofOrd, i0, dkOrd) = 0.0;
            }
          }
        }    
        
        for (int i0 = 0; i0 < dim0; i0++) {

          //double z0 = (1.0 - z)/2.0;
          //double z1 = (1.0 + z)/2.0;
          //double zq0 = (-z)*z0;
          //double zq1 = ( z)*z1;
          //double zqh = 4.0*z0*z1;

          double dz0 = -0.5;
          double dz1 = 0.5; 
          //double dzq0 = (-1.0)*z0 + (-z)*dz0;
          //double dzq1 = ( 1.0)*z1 + ( z)*dz1;
          //double dzqh = 4.0*(dz0*z1 + z0*dz1);

          double d2z0 = 0.0;
          double d2z1 = 0.0;
          double d2zq0 = (-1.0)*dz0 + (-1.0)*dz0;
          double d2zq1 = ( 1.0)*dz1 + ( 1.0)*dz1;
          double d2zqh = 8.0*( dz0*dz1);

          //           double d3z0 = 0.0;
          //           double d3z1 = 0.0;
          //           double d3zq0 = 0.0;
          //           double d3zq1 = 0.0;
          //           double d3zqh = 0.0;
        
          //        (2,4,7)-->(5,8,12)

          //---- 0,1,2
          //outputValues(0, i0) = zq0* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
          outputValues(0, i0, 5) = d2zq0* ( 4.0 );
          outputValues(0, i0, 8) = d2zq0* ( 4.0 );
          outputValues(0, i0, 12) = d2zq0* ( (1.0)*( 2.0) + (1.0)*(2.0) );

          //outputValues(1, i0) = zq0* x*(2.0*x - 1.0);
          outputValues(1, i0, 5) = d2zq0* ( 4.0 );
          outputValues(1, i0, 8) = 0.0;
          outputValues(1, i0, 12) = 0.0;

          //outputValues(2, i0) = zq0* y*(2.0*y - 1.0);
          outputValues(2, i0, 5) = 0.0;
          outputValues(2, i0, 8) = 0.0;
          outputValues(2, i0, 12) = d2zq0* (4.0);

          //---- 3,4,5
          // same as 0,1,2 with zq1 for zq0
          //outputValues(3, i0) = zq1* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
          outputValues(3, i0, 5) = d2zq1* ( 4.0 );
          outputValues(3, i0, 8) = d2zq1* ( 4.0 );
          outputValues(3, i0, 12) = d2zq1* ( (1.0)*( 2.0) + (1.0)*(2.0) );

          //outputValues(4, i0) = zq1* x*(2.0*x - 1.0);
          outputValues(4, i0, 5) = d2zq1* ( 4.0 );
          outputValues(4, i0, 8) = 0.0;
          outputValues(4, i0, 12) = 0.0;

          //outputValues(5, i0) = zq1* y*(2.0*y - 1.0);
          outputValues(5, i0, 5) = 0.0;
          outputValues(5, i0, 8) = 0.0;
          outputValues(5, i0, 12) = d2zq1* (4.0);

          //--- 6,7,8
          //outputValues(6, i0) = z0* (-4.0*x*(x + y - 1.0));
          outputValues(6, i0, 5) = d2z0* (-8.0);
          outputValues(6, i0, 8) = d2z0* (-4.0);
          outputValues(6, i0, 12) = 0.0;

          //outputValues(7, i0) = z0* ( 4.0*x*y);
          outputValues(7, i0, 5) = 0.0;
          outputValues(7, i0, 8) = d2z0* ( 4.0 );
          outputValues(7, i0, 12) = 0.0;

          //outputValues(8, i0) = z0* (-4.0*y*(x + y - 1.0));
          outputValues(8, i0, 5) = 0.0;
          outputValues(8, i0, 8) = d2z0* (-4.0);
          outputValues(8, i0, 12) = d2z0* (-8.0);

          //---- 9,10,11
          // same as 0,1,2 with zqh for zq0
          //outputValues(9, i0) = zqh* (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
          outputValues(9, i0, 5) = d2zqh* ( 4.0 );
          outputValues(9, i0, 8) = d2zqh* ( 4.0 );
          outputValues(9, i0, 12) = d2zqh* ( (1.0)*( 2.0) + (1.0)*(2.0) );

          //outputValues(10, i0) = zqh* x*(2.0*x - 1.0);
          outputValues(10, i0, 5) = d2zqh* ( 4.0 );
          outputValues(10, i0, 8) = 0.0;
          outputValues(10, i0, 12) = 0.0;

          //outputValues(11, i0) = zqh* y*(2.0*y - 1.0);
          outputValues(11, i0, 5) = 0.0;
          outputValues(11, i0, 8) = 0.0;
          outputValues(11, i0, 12) = d2zqh* (4.0);

          //--- 12,13,14
          // same as 6,7,8 with z1 for z0

          //outputValues(12, i0) = z1* (-4.0*x*(x + y - 1.0));
          outputValues(12, i0, 5) = d2z1* (-8.0);
          outputValues(12, i0, 8) = d2z1* (-4.0);
          outputValues(12, i0, 12) = 0.0;

          //outputValues(13, i0) = z1* ( 4.0*x*y);
          outputValues(13, i0, 5) = 0.0;
          outputValues(13, i0, 8) = d2z1* ( 4.0 );
          outputValues(13, i0, 12) = 0.0;

          //outputValues(14, i0) = z1* (-4.0*y*(x + y - 1.0));
          outputValues(14, i0, 5) = 0.0;
          outputValues(14, i0, 8) = d2z1* (-4.0);
          outputValues(14, i0, 12) = d2z1* (-8.0);

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
      TEST_FOR_EXCEPTION( !( Intrepid::isValidOperator(operatorType) ), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_WEDGE_C2_Serendipity_FEM): Invalid operator type");
    }
  }


  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_WEDGE_C2_Serendipity_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&       outputValues,
                                                                            const ArrayScalar &    inputPoints,
                                                                            const ArrayScalar &    cellVertices,
                                                                            const EOperator        operatorType) const {
    TEST_FOR_EXCEPTION( (true), std::logic_error,
                        ">>> ERROR (Basis_HGRAD_WEDGE_C2_Serendipity_FEM): FEM Basis calling an FVD member function");
  }
}// namespace Intrepid
#endif
