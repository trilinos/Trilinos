// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  Intrepid2_DataDef.hpp
//  Created by Roberts, Nathan V on 6/1/23.
//

#ifndef Intrepid2_DataDef_h
#define Intrepid2_DataDef_h

/** \file  Intrepid2_DataDef.hpp
   \brief  Defines implementations for the Data class that are not present in the declaration.
   \author Created by N.V. Roberts.
*/

#include "Intrepid2_DataCombiners.hpp"

namespace Intrepid2 {
  // forward declaration of a class that assists in in-place sums/products, etc. for Data containers.  Defined in Intrepid2_DataCombiners.hpp.
  template <class DataScalar,typename DeviceType, class BinaryOperator>
  class DataCombiner;

  template<class DataScalar,typename DeviceType>
  template<class BinaryOperator>
  void Data<DataScalar,DeviceType>::storeInPlaceCombination(const Data<DataScalar,DeviceType> &A, const Data<DataScalar,DeviceType> &B, BinaryOperator binaryOperator)
  {
    using DC = DataCombiner<DataScalar,DeviceType,BinaryOperator>;
    DC::storeInPlaceCombination(*this,A,B,binaryOperator);
  }
}
#endif /* Intrepid2_DataDef_h */
