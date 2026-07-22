// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_FIELDTAG_HPP
#define PHX_FIELDTAG_HPP

#include <string>
#include <typeinfo>
#include <iostream>
#include "Teuchos_RCP.hpp"
#include "Phalanx_config.hpp"

namespace PHX {

  class DataLayout;

  class FieldTag {

  public:

    FieldTag() = default;
    
    virtual ~FieldTag() = default;

    virtual Teuchos::RCP<FieldTag> clone() const = 0;

    virtual bool operator==(const FieldTag& t) const = 0;
    
    virtual bool operator!=(const FieldTag& t) const
    { return !(*this == t); };
    
    virtual const std::string& name() const = 0;

    virtual const PHX::DataLayout& dataLayout() const = 0;

    virtual PHX::DataLayout& nonConstDataLayout() = 0;

    virtual const std::type_info& dataTypeInfo() const = 0;
    
    //! Unique name identifier that can be used for strict weak ordering in stl std::map keys.
    virtual const std::string identifier() const = 0;

    virtual void print(std::ostream& os) const = 0;

  };

  std::ostream& operator<<(std::ostream& os, const PHX::FieldTag& t);
  
} 

#endif 
