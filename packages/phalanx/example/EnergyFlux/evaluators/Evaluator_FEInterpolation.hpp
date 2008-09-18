// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
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

#ifndef PHX_EXAMPLE_VP_FE_INTERPOLATION_HPP
#define PHX_EXAMPLE_VP_FE_INTERPOLATION_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_Field.hpp"

/** \brief Finite Element Interpolation Evaluator

    This object evaluates a scalar field and it's gradient at the
    quadrature points for a specific variable.

*/
template<typename EvalT, typename Traits>
class FEInterpolation : public PHX::EvaluatorWithBaseImpl<Traits>,
			public PHX::EvaluatorDerived<EvalT, Traits>  {
  
public:
  
  FEInterpolation(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);
  
private:

  typedef typename EvalT::ScalarT ScalarT;

  //! Values at nodes
  PHX::Field< ScalarT > val_node;

  //! Values at quadrature points
  PHX::Field< ScalarT > val_qp;

  //! Gradient values at quadrature points
  PHX::Field< MyVector<ScalarT> > val_grad_qp;
  
};

#include "Evaluator_FEInterpolation_Def.hpp"

#endif
