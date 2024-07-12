// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_SCALAR_CONTAINER_BASE_DEF_HPP
#define PHX_SCALAR_CONTAINER_BASE_DEF_HPP

#include "Teuchos_Assert.hpp"

// **************************************************************************
template<typename Traits>
PHX::EvaluationContainerBase<Traits>::EvaluationContainerBase()
{

}

// **************************************************************************
template<typename Traits>
PHX::EvaluationContainerBase<Traits>::~EvaluationContainerBase()
{

}

// **************************************************************************
template<typename Traits>
void PHX::EvaluationContainerBase<Traits>::
requireField(const PHX::FieldTag& f) 
{ 
  dag_manager_.requireField(f);
}

// **************************************************************************
template<typename Traits>
void PHX::EvaluationContainerBase<Traits>::
registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& e) 
{ 
  dag_manager_.registerEvaluator(e);
}

// **************************************************************************
template<typename Traits>
void PHX::EvaluationContainerBase<Traits>::
writeGraphvizFile(const std::string filename,
		  bool writeEvaluatedFields,
		  bool writeDependentFields,
		  bool debugRegisteredEvaluators) const 
{ 
  dag_manager_.writeGraphvizFile(filename, 
				 writeEvaluatedFields, 
				 writeDependentFields, 
				 debugRegisteredEvaluators);

  {
    std::ofstream ofs;
    std::string ext_filename = filename+".txt";
    ofs.open(ext_filename.c_str());
    this->print(ofs);
  }
}
    
// **************************************************************************
template<typename Traits>
std::ostream&
PHX::operator<<(std::ostream& os, const PHX::EvaluationContainerBase<Traits>& sc)
{ 
  sc.print(os);
  return os;
}

// **************************************************************************

#endif
