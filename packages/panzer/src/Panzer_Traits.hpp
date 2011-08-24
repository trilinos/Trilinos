#ifndef PANZER_TRAITS_HPP
#define PANZER_TRAITS_HPP

#include "Panzer_config.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// add embedded UQ
#ifdef HAVE_STOKHOS
   #include "Stokhos_Sacado.hpp"
#endif

// mpl (Meta Programming Library) templates
#include "Sacado_mpl_vector.hpp"
#include "Sacado_mpl_find.hpp"
#include "boost/mpl/map.hpp"
#include "boost/mpl/find.hpp"

// Scalar types
#include "Sacado.hpp"

// traits Base Class
#include "Phalanx_Traits_Base.hpp"

// Include User Data Types
#include "Phalanx_Allocator_Contiguous.hpp"
#include "Panzer_Workset.hpp"

// Debugging information
#include "Phalanx_TypeStrings.hpp"

// add embedded UQ
#ifdef HAVE_STOKHOS
   #include "Stokhos_StandardStorage.hpp"
   #include "Sacado_PCE_OrthogPoly.hpp"
   #include "Stokhos_OrthogPolyExpansion.hpp"
#endif

namespace panzer {

  struct Traits : public PHX::TraitsBase {

    // ******************************************************************
    // *** Scalar Types
    // ******************************************************************
    
    // Scalar types we plan to use
    typedef double RealType;
    typedef Sacado::Fad::DFad<double> FadType;

    #ifdef HAVE_STOKHOS
       typedef Stokhos::StandardStorage<int,RealType> SGStorageType;
       typedef Sacado::PCE::OrthogPoly<RealType,SGStorageType> SGType;
       typedef Sacado::Fad::DFad<SGType> SGFadType;
    #endif
    
    // ******************************************************************
    // *** Evaluation Types
    // ******************************************************************
    struct Residual { typedef RealType ScalarT; };
    struct Jacobian { typedef FadType ScalarT;  };
    #ifdef HAVE_STOKHOS
       struct SGResidual { typedef SGType ScalarT; };
       struct SGJacobian { typedef SGFadType ScalarT;  };
    #endif
    typedef Sacado::mpl::vector<Residual, Jacobian
                                #ifdef HAVE_STOKHOS
                                   , SGResidual, SGJacobian
                                #endif
                               > EvalTypes;

    // ******************************************************************
    // *** Response Types
    // ******************************************************************

    struct Value {      typedef RealType ScalarT; typedef Residual EvalType; };
    struct Derivative { typedef FadType ScalarT;  typedef Jacobian EvalType; };

    typedef Sacado::mpl::vector<Value, Derivative 
                               > RespTypes;

    // ******************************************************************
    // *** Data Types
    // ******************************************************************
    
    // Create the data types for each evaluation type
    
    // Residual (default scalar type is RealType)
    typedef Sacado::mpl::vector< RealType > ResidualDataTypes;
  
    // Jacobian (default scalar type is Fad<double, double>)
    typedef Sacado::mpl::vector< FadType > JacobianDataTypes;

    #ifdef HAVE_STOKHOS
       typedef Sacado::mpl::vector< SGType > SGResidualDataTypes;
       typedef Sacado::mpl::vector< SGFadType > SGJacobianDataTypes;
    #endif

    // Maps the key EvalType a vector of DataTypes
    typedef boost::mpl::map<
      boost::mpl::pair<Residual, ResidualDataTypes>,
      boost::mpl::pair<Jacobian, JacobianDataTypes>
      #ifdef HAVE_STOKHOS
         , boost::mpl::pair<SGResidual, SGResidualDataTypes>
         , boost::mpl::pair<SGJacobian, SGJacobianDataTypes>
      #endif
    >::type EvalToDataMap;

    // ******************************************************************
    // *** Allocator Type
    // ******************************************************************
    typedef PHX::ContiguousAllocator<double> Allocator;

    // ******************************************************************
    // *** User Defined Object Passed in for Evaluation Method
    // ******************************************************************
    struct SD { 
      Teuchos::RCP< std::vector<panzer::Workset> > worksets_;
    };
    typedef SD SetupData;
    typedef panzer::Workset& EvalData;
    typedef void* PreEvalData;
    typedef void* PostEvalData;

  };
 
  // ******************************************************************
  // ******************************************************************
  // Debug strings.  Specialize the Evaluation and Data types for the
  // TypeString object in the phalanx/src/Phalanx_TypeStrings.hpp file.
  // ******************************************************************
  // ******************************************************************

}

namespace PHX {
  
  // Evaluation Types
  template<> struct TypeString<panzer::Traits::Residual> 
  { static const std::string value; };

  template<> struct TypeString<panzer::Traits::Jacobian> 
  { static const std::string value; };

  #ifdef HAVE_STOKHOS
     template<> struct TypeString<panzer::Traits::SGResidual> 
     { static const std::string value; };
   
     template<> struct TypeString<panzer::Traits::SGJacobian> 
     { static const std::string value; };
  #endif 

  // Data Types
  template<> struct TypeString<double> 
  { static const std::string value; };

  template<> struct TypeString< Sacado::Fad::DFad<double> > 
  { static const std::string value; };

  #ifdef HAVE_STOKHOS
     template<> struct TypeString<panzer::Traits::SGType> 
     { static const std::string value; };
   
     template<> struct TypeString< panzer::Traits::SGFadType> 
     { static const std::string value; };
  #endif 

  // Response Types
  template<> struct TypeString<panzer::Traits::Value> 
  { static const std::string value; };

  template<> struct TypeString<panzer::Traits::Derivative> 
  { static const std::string value; };
}

#endif
