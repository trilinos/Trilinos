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

#ifndef PHX_FIELD_MANAGER_HPP
#define PHX_FIELD_MANAGER_HPP

#include <cstddef>
#include <iostream>
#include <vector>
#include <algorithm>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_Field.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_EvaluationContainer_TemplateManager.hpp"

namespace PHX {

  template<typename Traits>
  class FieldManager {
    
  public:

    typedef typename PHX::EvaluationContainer_TemplateManager<Traits>::iterator iterator;

    FieldManager();

    ~FieldManager();
    
    void requireFieldForAllEvaluationTypes(const PHX::FieldTag& t);
    
    template<typename EvalT>
    void requireField(const PHX::FieldTag& t);

    void registerEvaluatorForAllEvaluationTypes(const Teuchos::RCP< PHX::Evaluator<Traits> >& e);
    
    template<typename EvalT>
    void registerEvaluator(const Teuchos::RCP< PHX::Evaluator<Traits> >& e);

    void registerEvaluator(typename PHX::FieldManager<Traits>::iterator it,
			   const Teuchos::RCP< PHX::Evaluator<Traits> >& e);
    
    template<typename DataT, typename EvalT> 
    void getFieldData(PHX::Field<DataT>& f);
    
    template<typename DataT, typename EvalT, PHX::ArrayOrder Order, 
	     typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	     typename Tag4, typename Tag5, typename Tag6, typename Tag7> 
    void getFieldData(PHX::MDField<DataT,Order,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,
		      Tag6,Tag7>& f);
    
    template<typename DataT, typename EvalT> 
    void getFieldData(const PHX::FieldTag& t, Teuchos::ArrayRCP<DataT>& d);
    
    //! Allows for different size worksets for each evaluation type
    template<typename EvalT>
    void postRegistrationSetupForType(std::size_t max_num_cells);

    //! Forces the same size workset for all evaluation types
    void postRegistrationSetup(std::size_t max_num_cells);

    template<typename EvalT>
    void evaluateFields(typename Traits::EvalData d);

    template<typename EvalT>
    void preEvaluate(typename Traits::PreEvalData d);

    template<typename EvalT>
    void postEvaluate(typename Traits::PostEvalData d);

    template<typename EvalT>
    std::size_t getMaxNumCells() const;

    //! Return iterator to first EvaluationContainer
    typename FieldManager::iterator begin();

    //! Return iterator to last EvaluationContainer
    typename FieldManager::iterator end();

    void print(std::ostream& os) const;

  private:

    typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;

    std::size_t m_num_evaluation_types;

    PHX::EvaluationContainer_TemplateManager<Traits> m_eval_containers;

    std::vector<std::size_t> m_max_num_cells;

  };

  template<typename Traits>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::FieldManager<Traits>& vm);

} 

#include "Phalanx_FieldManager_Def.hpp"

#endif 
