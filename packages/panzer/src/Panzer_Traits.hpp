// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

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
#include "Sacado_CacheFad_DFad.hpp"
#include "Sacado_ELRFad_DFad.hpp"
#include "Sacado_ELRCacheFad_DFad.hpp"

// traits Base Class
#include "Phalanx_Traits_Base.hpp"

// Include User Data Types
#include "Phalanx_Allocator_Contiguous.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

// Debugging information
#include "Phalanx_TypeStrings.hpp"

// add embedded UQ
#ifdef HAVE_STOKHOS
   #include "Stokhos_StandardStorage.hpp"
   #include "Sacado_PCE_OrthogPoly.hpp"
   #include "Stokhos_OrthogPolyExpansion.hpp"
#endif

namespace panzer {

  class LinearObjContainer;

  struct Traits : public PHX::TraitsBase {

    // ******************************************************************
    // *** Scalar Types
    // ******************************************************************
    
    // Scalar types we plan to use
    typedef double RealType;
    typedef Sacado::Fad::DFad<double> FadType;
    //typedef Sacado::CacheFad::DFad<double> FadType;
    //typedef Sacado::ELRFad::DFad<double> FadType;
    //typedef Sacado::ELRCacheFad::DFad<double> FadType;

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

/*
    struct PED {
       struct DirichletData {
          //! Stores information about which global indices were set as dirichlet conditions
          Teuchos::RCP<LinearObjContainer> ghostedCounter;
       } dirichletData;
    };
    typedef PED& PreEvalData;
*/
    typedef GlobalEvaluationDataContainer& PreEvalData;

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

  template<> struct TypeString< Sacado::CacheFad::DFad<double> > 
  { static const std::string value; };

  template<> struct TypeString< Sacado::ELRFad::DFad<double> > 
  { static const std::string value; };

  template<> struct TypeString< Sacado::ELRCacheFad::DFad<double> > 
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
