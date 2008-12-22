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

#ifndef PHX_FIELDTAG_HPP
#define PHX_FIELDTAG_HPP

#include <string>
#include <typeinfo>
#include <iostream>
#include "Teuchos_RCP.hpp"
#include "Phalanx_ConfigDefs.hpp"

namespace PHX {

  class DataLayout;

  class FieldTag {

  public:

    FieldTag() {}
    
    virtual ~FieldTag() {}

    virtual Teuchos::RCP<FieldTag> clone() const = 0;

    virtual bool operator==(const FieldTag& t) const = 0;
    
    virtual bool operator!=(const FieldTag& t) const
    { return !(*this == t); };
    
    virtual const std::string& name() const = 0;

    virtual const PHX::DataLayout& dataLayout() const = 0;

    virtual const std::type_info& dataTypeInfo() const = 0;
    
    //! Unique name identifier that can be used for strict weak ordering in stl std::map keys.
    virtual const std::string identifier() const = 0;

    virtual void print(std::ostream& os) const = 0;

  };

  std::ostream& operator<<(std::ostream& os, const PHX::FieldTag& t);
  
} 

#endif 
