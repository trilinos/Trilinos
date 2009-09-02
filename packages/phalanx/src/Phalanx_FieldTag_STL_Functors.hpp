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

#ifndef PHX_FIELDTAG_STL_FUNCTORS_HPP
#define PHX_FIELDTAG_STL_FUNCTORS_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Phalanx_FieldTag.hpp"

namespace PHX {

  //! Functor for RCP<FieldTag> comparisons in stl std::map.  This must follow strict weak ordering rules.
  struct FTComp
  {
    bool operator()(const Teuchos::RCP<const PHX::FieldTag> a, 
		    const Teuchos::RCP<const PHX::FieldTag> b) const
    {
      return strcmp(a->identifier().c_str(), b->identifier().c_str()) < 0;
    }
  };

  //! Predicate for searches of stl std::vector< RCP<FieldTag> > using find_if algorithm.
  struct FTPred {
    const Teuchos::RCP<FieldTag> m_a;
    FTPred(const Teuchos::RCP<FieldTag> a) : m_a(a) {}
    bool operator()(const Teuchos::RCP<FieldTag> b) 
    {
      if (*m_a == *b)
	return true;
      else 
	return false;
    }
  };

  //! Predicate for searches of stl std::vector< RCP<FieldTag> > using find_if algorithm.
  struct FTPredRef {
    const FieldTag& m_a;
    FTPredRef(const FieldTag& a) : m_a(a) {}
    bool operator()(const Teuchos::RCP<FieldTag> b) 
    {
      if (m_a == *b)
	return true;
      else 
	return false;
    }
  };

} 

#endif 
