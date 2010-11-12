// @HEADER
// @HEADER

// ******************************************************************
// ******************************************************************
// Debug strings.  Specialize the Evaluation and Data types for the
// TypeString object in the phalanx/src/Phalanx_TypeStrings.hpp file.
// ******************************************************************
// ******************************************************************

#include <string>
#include "Panzer_Traits.hpp"

const std::string PHX::TypeString<panzer::Traits::Residual>::value = "Residual";

const std::string PHX::TypeString<panzer::Traits::Jacobian>::value = "Jacobian";

const std::string PHX::TypeString<double>::value = "double";

const std::string PHX::TypeString< VectorT<double> >::value = 
  "VectorT<double>";

const std::string PHX::TypeString< TensorT<double> >::value = 
  "TensorT<double>";

const std::string PHX::TypeString< Sacado::Fad::DFad<double> >::value = 
  "Sacado::Fad::DFad<double>";

const std::string PHX::TypeString< VectorT<Sacado::Fad::DFad<double> > >::value = 
  "Sacado::Fad::DFad< VectorT<double> >";

const std::string PHX::TypeString< TensorT<Sacado::Fad::DFad<double> > >::value = 
  "Sacado::Fad::DFad< TensorT<double> >";

// Extra declarations not used in field manager

const std::string PHX::TypeString<Vector>::value = 
  "Vector";

const std::string PHX::TypeString<Tensor>::value = 
  "Tensor";
