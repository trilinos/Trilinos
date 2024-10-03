// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  JacobianFlopEstimate.hpp
//  Trilinos
//
//  Created by Roberts, Nathan V on 7/7/21.
//

#ifndef Intrepid2_JacobianFlopEstimate_hpp
#define Intrepid2_JacobianFlopEstimate_hpp

inline double flopsPerJacobian(const int &spaceDim, const int &numPoints, const int &numGeometryNodes)
{
  // implementation looks like:
//  for (ordinal_type i=0;i<dim;++i)
//  for (ordinal_type j=0;j<dim;++j) {
//    _jacobian(cell, point, i, j) = 0;
//    for (ordinal_type bf=0;bf<cardinality;++bf)
//      _jacobian(cell, point, i, j) += _worksetCells(cell+_startCell, bf, i) * _basisGrads(bf, point, j); // 2 flops: one multiply, one add
//  }
  return 2.0 * spaceDim * spaceDim * numPoints * numGeometryNodes;
}

inline double flopsPerJacobianDet(const int &spaceDim, const int &numPoints)
{
  //implementation in RealSpaceTools:
  /*value_type r_val = 0.0;
  switch (dim) {
  case 3:
    r_val = ( inMat(0,0) * inMat(1,1) * inMat(2,2) + // 3 flops: 2 mults, 1 add
              inMat(1,0) * inMat(2,1) * inMat(0,2) + // 3 flops: 2 mults, 1 add
              inMat(2,0) * inMat(0,1) * inMat(1,2) - // 3 flops: 2 mults, 1 subtract
              inMat(2,0) * inMat(1,1) * inMat(0,2) - // 3 flops: 2 mults, 1 subtract
              inMat(0,0) * inMat(2,1) * inMat(1,2) - // 3 flops: 2 mults, 1 subtract
              inMat(1,0) * inMat(0,1) * inMat(2,2) ); // 2 flops: 2 mults
    break;
  case 2:
    r_val = ( inMat(0,0) * inMat(1,1) -
              inMat(0,1) * inMat(1,0) );
    break;
  case 1:
    r_val = ( inMat(0,0) );
    break;
  }
  return r_val;*/
  int r_val;
  switch (spaceDim) {
    case 3: r_val = 17.0 * numPoints; break;
    case 2: r_val = 3.0 * numPoints; break;
    case 1: r_val = 0.0; break;
    default: INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unhandled spaceDim");
  }
  return r_val;
}

inline double flopsPerJacobianInverse(const int &spaceDim, const int &numPoints)
{
  // implementation looks like:
  // const value_type val = RealSpaceTools<>::Serial::det(mat);
  double totalFlops = flopsPerJacobianDet(spaceDim, numPoints);
  
  /*
   case 1: {
     inv(0,0) = 1.0/mat(0,0);
     break;
   }
   case 2: {
     inv(0,0) =   mat(1,1)/val;
     inv(1,1) =   mat(0,0)/val;

     inv(1,0) = - mat(1,0)/val;
     inv(0,1) = - mat(0,1)/val;
     break;
   }
   */
  
  if (spaceDim == 1)
  {
    // inv(0,0) = 1.0/mat(0,0); // 1 flop: 1 division
    totalFlops += 1 * numPoints;
  }
  else if (spaceDim == 2)
  {
//    inv(0,0) =   mat(1,1)/val; // 1 flop: 1 division
//    inv(1,1) =   mat(0,0)/val; // 1 flop: 1 division
//
//    inv(1,0) = - mat(1,0)/val; // 2 flops: 1 negation, 1 division
//    inv(0,1) = - mat(0,1)/val; // 2 flops: 1 negation, 1 division
    totalFlops += 6 * numPoints;
  }
  else if (spaceDim == 3)
  {
  //  val0 =   mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2); // 3 flops: 2 mults, 1 subtraction
  //  val1 = - mat(1,0)*mat(2,2) + mat(2,0)*mat(1,2); // 4 flops: 2 mults, 1 negation, 1 add
  //  val2 =   mat(1,0)*mat(2,1) - mat(2,0)*mat(1,1); // 3 flops: 2 mults, 1 subtraction
  //
  //  inv(0,0) = val0/val; // 1 flop
  //  inv(1,0) = val1/val; // 1 flop
  //  inv(2,0) = val2/val; // 1 flop
  //
  //  val0 =   mat(2,1)*mat(0,2) - mat(0,1)*mat(2,2); // 3
  //  val1 =   mat(0,0)*mat(2,2) - mat(2,0)*mat(0,2); // 3
  //  val2 = - mat(0,0)*mat(2,1) + mat(2,0)*mat(0,1); // 4
  //
  //  inv(0,1) = val0/val; // 1
  //  inv(1,1) = val1/val; // 1
  //  inv(2,1) = val2/val; // 1
  //
  //  val0 =   mat(0,1)*mat(1,2) - mat(1,1)*mat(0,2); // 3
  //  val1 = - mat(0,0)*mat(1,2) + mat(1,0)*mat(0,2); // 4
  //  val2 =   mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1); // 3
  //
  //  inv(0,2) = val0/val; // 1
  //  inv(1,2) = val1/val; // 1
  //  inv(2,2) = val2/val; // 1
    totalFlops += 36.0 * numPoints;
  }
  else
  {
    INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unhandled spaceDim");
  }
  
  return totalFlops;
}

#endif /* Intrepid2_JacobianFlopEstimate_hpp */
