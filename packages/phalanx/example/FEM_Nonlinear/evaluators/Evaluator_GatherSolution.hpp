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

#ifndef PHX_EXAMPLE_GATHER_SOLUTION_HPP
#define PHX_EXAMPLE_GATHER_SOLUTION_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

#include "Epetra_Vector.h"

/** \brief Gathers solution values from the Newton solution vector into 
    the nodal fields of the field manager

    Currently makes an assumption that the stride is constant for dofs
    and that the nmber of dofs is equal to the size of the solution
    names vector.

*/
template<typename EvalT, typename Traits> class GatherSolution;

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Residual 
// **************************************************************
template<typename Traits>
class GatherSolution<PHX::MyTraits::Residual,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<PHX::MyTraits::Residual, Traits>  {
  
public:
  
  GatherSolution(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);
  
private:

  typedef typename PHX::MyTraits::Residual::ScalarT ScalarT;

  std::vector< PHX::MDField<ScalarT,PHX::NaturalOrder,Cell,Node> > val;
 
  Teuchos::RCP<Epetra_Vector> x;

  std::size_t num_nodes;
};

// **************************************************************
// Jacobian
// **************************************************************
template<typename Traits>
class GatherSolution<PHX::MyTraits::Jacobian,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<PHX::MyTraits::Jacobian, Traits>  {
  
public:
  
  GatherSolution(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);
  
private:

  typedef typename PHX::MyTraits::Jacobian::ScalarT ScalarT;

  std::vector< PHX::MDField<ScalarT,PHX::NaturalOrder,Cell,Node> > val;
 
  Teuchos::RCP<Epetra_Vector> x;

  std::size_t num_nodes;
};

// **************************************************************

#include "Evaluator_GatherSolution_Def.hpp"

// **************************************************************
#endif
