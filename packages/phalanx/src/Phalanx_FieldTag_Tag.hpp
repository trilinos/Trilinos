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
#include <type_traits>
#include "Phalanx_config.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_DataLayout.hpp"

namespace PHX {

  /*! \brief Typed Field Tag

      This class is a concrete implementation of the FieldTag base
      class that is templated on the data type to determine type
      information. NOTE: The constness on the DataT is ignored for
      object comparisons and identifier string creation.

  */
  template<typename DataT>
  class Tag : public PHX::FieldTag {

  public:

    typedef DataT value_type;

    // Default constructor
    Tag();

    Tag(const std::string& name, const Teuchos::RCP<PHX::DataLayout>& dl);

    // Use SFINAE to remove this ctor if the data types are not the
    // same (excluding const differences). This works but the error
    // message is cryptic. Instead use the static_assert version below
    // template<typename InDataT, typename T = DataT>
    // Tag(const Tag<InDataT>& t,
    //     typename std::enable_if<std::is_same<typename std::remove_const<InDataT>::type,typename std::remove_const<T>::type>::value>::type * = nullptr)
    //   : Tag(t.m_name,t.m_data_layout) {}

    template<typename InDataT, typename T = DataT>
    Tag(const Tag<InDataT>& t) : Tag(t.m_name,t.m_data_layout)
    {
      using type1 = typename std::remove_const<InDataT>::type;
      using type2 = typename std::remove_const<T>::type;
      static_assert(std::is_same<type1,type2>::value,
                    "** ERROR ** PHX::Tag ctor: tag data types are not the same (excluding const)!");
    }

    ~Tag() noexcept;

    Teuchos::RCP<FieldTag> clone() const override;

    void operator=(const PHX::Tag<const DataT>& t);
    
    bool operator==(const FieldTag& t) const override;
    
    const std::string& name() const override;

    const PHX::DataLayout& dataLayout() const override;

    PHX::DataLayout& nonConstDataLayout() override;

    const std::type_info& dataTypeInfo() const override;

    const std::string identifier() const override;

    void print(std::ostream& os) const override;

    template<typename> friend class Tag;

  protected:

    std::string m_name;
    
    Teuchos::RCP<PHX::DataLayout> m_data_layout;

  };

  template<typename DataT>
  std::ostream& operator<<(std::ostream& os, const PHX::Tag<DataT>& t);
  
} 

#include "Phalanx_FieldTag_Tag_Def.hpp"

#endif 
