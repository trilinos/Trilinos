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

#ifndef PHX_EVALUATOR_UTILITIES_H
#define PHX_EVALUATOR_UTILITIES_H

#include <vector>

#include "Phalanx_Field.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_FieldManager.hpp"

namespace PHX {

  /*! @brief Utilities to hide templating in concrete Evaluators.
   
  */
  template<typename EvalT, typename Traits> 
  struct EvaluatorUtilities {
    
    template <typename DataT>
    void setFieldData(PHX::Field<DataT>& f, PHX::FieldManager<Traits>& fm) 
    {
      fm.template getFieldData<DataT,EvalT>(f);
    }

    template <typename DataT,
	      typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	      typename Tag4, typename Tag5, typename Tag6, typename Tag7>
    void setFieldData(PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,
		      Tag6,Tag7>& f, PHX::FieldManager<Traits>& fm) 
    {
      fm.template getFieldData<DataT,EvalT>(f);
    }

    template <typename DataT>
    std::size_t getWorksetSize(PHX::Field<DataT>& f, 
			       PHX::FieldManager<Traits>& fm) 
    {
      return fm.template 
	getWorksetSize<EvalT>(f.fieldTag().dataLayout().worksetType());
    }

    template <typename DataT,
	      typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	      typename Tag4, typename Tag5, typename Tag6, typename Tag7>
    std::size_t getWorksetSize(PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,
			       Tag4,Tag5,Tag6,Tag7>& f, 
			       PHX::FieldManager<Traits>& fm) 
    {
      return fm.template
	getWorksetSize<EvalT>(f.fieldTag().dataLayout().worksetType());
    }

  };
}

#endif
