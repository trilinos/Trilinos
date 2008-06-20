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
