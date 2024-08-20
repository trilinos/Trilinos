// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EXAMPLE_FACTORY_TRAITS_HPP
#define EXAMPLE_FACTORY_TRAITS_HPP

// mpl (Meta Programming Library) templates
#include "Sacado_mpl_vector.hpp"

// User Defined Evaluator Types
#include "Evaluator_Constant.hpp"
#include "Evaluator_Density.hpp"
#include "Evaluator_EnergyFlux_Fourier.hpp"
#include "Evaluator_FEInterpolation.hpp"
#include "Evaluator_NonlinearSource.hpp"


#include "Sacado_mpl_placeholders.hpp"
using namespace Sacado::mpl::placeholders;

/*! \brief Struct to define Evaluator objects for the EvaluatorFactory.
    
    Preconditions:
    - You must provide a Sacado::mpl::vector named EvaluatorTypes that contain all Evaluator objects that you wish the factory to build.  Do not confuse evaluator types (concrete instances of evaluator objects) with evaluation types (types of evaluations to perform, i.e., Residual, Jacobian). 

*/
template<typename Traits>
struct MyFactoryTraits {
  
  static const int id_constant = 0;
  static const int id_density = 1;
  static const int id_fourier = 2;
  static const int id_feinterpolation = 3;
  static const int id_nonlinearsource = 4;

  typedef Sacado::mpl::vector< Constant<_,Traits>,             // 0
			       Density<_,Traits>,              // 1
			       Fourier<_,Traits>,              // 2
			       FEInterpolation<_,Traits>,      // 3
			       NonlinearSource<_,Traits>       // 4
  > EvaluatorTypes;
  
};

#endif
