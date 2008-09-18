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

#ifndef PHX_FIELD_H
#define PHX_FIELD_H

#include <iostream>
#include <string>
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_FieldTag_Tag.hpp"

namespace PHX {

  template<typename DataT>
  class Field {
    
  public:

    typedef DataT value_type;
    
    Field(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& t);
    
    Field(const PHX::Tag<DataT>& v);
    
    Field();
    
    ~Field();
    
    const PHX::FieldTag& fieldTag() const;

    DataT& operator[](int index);

    typename Teuchos::ArrayRCP<DataT>::Ordinal size() const;

    void setFieldTag(const PHX::Tag<DataT>& t);
    
    void setFieldData(const Teuchos::ArrayRCP<DataT>& d);
    
    void print(std::ostream& os) const;

  private:
    
    PHX::Tag<DataT> m_tag;
    
    Teuchos::ArrayRCP<DataT> m_field_data;

#ifdef PHX_DEBUG
    bool m_tag_set;
    bool m_data_set;
    static const std::string m_field_tag_error_msg;
    static const std::string m_field_data_error_msg;
#endif

  };
  
  template<typename DataT>
  std::ostream& operator<<(std::ostream& os, const PHX::Field<DataT>& h);
  
} 

#include "Phalanx_Field_Def.hpp"

#endif 
