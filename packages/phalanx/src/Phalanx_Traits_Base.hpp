// @HEADER
// ************************************************************************
// 
//            Phalanx: A Partial Differential Equation Assembly 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PHX_TRAITS_BASE_HPP
#define PHX_TRAITS_BASE_HPP

#include "Phalanx_ConfigDefs.hpp"

namespace PHX {
  
  struct TraitsBase {
    
    template<typename ScalarT> 
    struct TypeString
    { static const std::string value; };

    /*  Alternative way to handle connector types
    template<typename DataT> 
    struct DataTypeInfo
    { 
      void Error_You_Must_Specialize_This_Struct_For_All_Data_Types();
      // Please implement:
      // (1) a typedef for "scalar_type" associated with data type
      // (2) a typedef for "entity_type" associated with data type
    };
    */
    

  };
  
  template<typename ScalarT>
  const std::string TraitsBase::TypeString<ScalarT>::value = 
    "WARNING: Undefined Type! Please specialize the PHX::TraitsBase::TypeString struct in your traits class for all possible Scalar types (ScalarT) and Data types (DataT) that are defined in your traits class.";
  
}
   
#endif
