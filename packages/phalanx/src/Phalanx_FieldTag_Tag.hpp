// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#ifndef PHX_FIELDTAG_TAG_HPP
#define PHX_FIELDTAG_TAG_HPP

#include <string>
#include <typeinfo>
#include <iostream>
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
