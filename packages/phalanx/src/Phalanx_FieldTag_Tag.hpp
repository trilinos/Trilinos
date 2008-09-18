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

#ifndef PHX_FIELDTAG_TAG_HPP
#define PHX_FIELDTAG_TAG_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_DataLayout.hpp"

namespace PHX {

  /*! \brief Typed Field Tag

      This class is a concrete implementation of the FieldTag base
      class that is templated on the data type to determine type
      information.

  */
  template<typename DataT>
  class Tag : public PHX::FieldTag {

  public:

    typedef DataT value_type;

    Tag(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& dl);
    
    ~Tag();

    Teuchos::RCP<FieldTag> clone() const;

    void operator=(const PHX::Tag<DataT>& t);
    
    bool operator==(const FieldTag& t) const;
    
    const std::string& name() const;

    const PHX::DataLayout& dataLayout() const;

    const std::type_info& dataTypeInfo() const;

    const std::string identifier() const;

    void print(std::ostream& os) const;

  protected:

    std::string m_name;
    
    Teuchos::RCP<PHX::DataLayout> m_data_layout;

  };

  template<typename DataT>
  std::ostream& operator<<(std::ostream& os, const PHX::Tag<DataT>& t);
  
} 

#include "Phalanx_FieldTag_Tag_Def.hpp"

#endif 
