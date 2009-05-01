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

#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename Traits>
GatherSolution<PHX::MyTraits::Residual, Traits>::
GatherSolution(const Teuchos::ParameterList& p)
{ 
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Solution Names"));

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");

  x = p.get< Teuchos::RCP<Epetra_Vector> >("Solution Vector");

  val.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    PHX::MDField<ScalarT,Cell,Node> f(names[eq],dl);
    val[eq] = f;
    this->addEvaluatedField(val[eq]);
  }

  this->setName("Gather Solution");
}

// **********************************************************************
template<typename Traits> 
void GatherSolution<PHX::MyTraits::Residual, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& fm)
{
  for (std::size_t eq = 0; eq < val.size(); ++eq)
    this->utils.setFieldData(val[eq],fm);

  num_nodes = val[0].dimension(1);
}

// **********************************************************************
template<typename Traits>
void GatherSolution<PHX::MyTraits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  const std::size_t num_eq = val.size();
  std::vector<Element_Linear2D>::iterator element = workset.begin;

  std::size_t cell = 0;
  for (; element != workset.end; ++element,++cell) {
    
    for (std::size_t node = 0; node < num_nodes; node++) {
      int node_GID = element->globalNodeId(node);
      //int firstDOF = x->Map().LID(node_GID * num_eq);
      int firstDOF = x->Map().LID(node_GID) * num_eq;
      for (std::size_t eq = 0; eq < val.size(); eq++)
	(val[eq])(cell,node) = (*x)[firstDOF + eq];
    }

  }

}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename Traits>
GatherSolution<PHX::MyTraits::Jacobian, Traits>::
GatherSolution(const Teuchos::ParameterList& p)
{ 
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Solution Names"));

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");

  x = p.get< Teuchos::RCP<Epetra_Vector> >("Solution Vector");

  val.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    PHX::MDField<ScalarT,Cell,Node> f(names[eq],dl);
    val[eq] = f;
    this->addEvaluatedField(val[eq]);
  }

  this->setName("Gather Solution");
}

// **********************************************************************
template<typename Traits> 
void GatherSolution<PHX::MyTraits::Jacobian, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& fm)
{
  for (std::size_t eq = 0; eq < val.size(); ++eq)
    this->utils.setFieldData(val[eq],fm);

  num_nodes = val[0].dimension(1);
}

// **********************************************************************
template<typename Traits>
void GatherSolution<PHX::MyTraits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  const std::size_t num_eq = val.size();
  const std::size_t num_dof = num_eq * num_nodes;
  std::vector<Element_Linear2D>::iterator element = workset.begin;

  std::size_t cell = 0;
  for (; element != workset.end; ++element,++cell) {
    
    for (std::size_t node = 0; node < num_nodes; node++) {
      int node_GID = element->globalNodeId(node);
      int firstDOF = x->Map().LID(node_GID) * num_eq;

      for (std::size_t eq = 0; eq < val.size(); eq++) {
	(val[eq])(cell,node) = 
	  Sacado::Fad::DFad<double>(num_dof, (*x)[firstDOF + eq]);
	(val[eq])(cell,node).fastAccessDx(num_eq * node + eq) = 1.0;

	// std::cout << "val[" << eq << "](" << cell << "," << node << ") = "
// 		  << (val[eq])(cell,node) << std::endl;
      }
    }

  }

}

// **********************************************************************
