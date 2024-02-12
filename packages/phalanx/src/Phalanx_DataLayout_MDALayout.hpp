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


#ifndef PHX_DATA_LAYOUT_MDALAYOUT
#define PHX_DATA_LAYOUT_MDALAYOUT

#include "Phalanx_DataLayout.hpp"
#include <string>
#include <iostream>
#include <vector>
#include <type_traits>

namespace PHX {

  // ***********************
  // Pack counting template
  // ***********************
  template<typename...ParamPack> struct PackSize;

  template<>
  struct PackSize<> : std::integral_constant<int,0> {};

  template<typename NextMember, typename...ParamPack>
  struct PackSize<NextMember,ParamPack...> : std::integral_constant<int,1+PackSize<ParamPack...>::value>{};

  // ******************************************************************
  /*! \brief A concrete implementation of the DataLayout class for compile time checked multidimensional arrays.
  */
  template<typename... Tags>
  class MDALayout : public DataLayout {

  public:
    static constexpr int Rank = PHX::PackSize<Tags...>::value;
    static_assert(Rank > 0);

    template<typename... Extents>
    MDALayout(Extents... extents);

    template<typename... Extents>
    MDALayout(const std::string& prefix,Extents... extents);

    template<typename... Extents>
    MDALayout(const char* prefix,Extents... extents);

    virtual ~MDALayout() noexcept {}

    virtual bool operator==(const DataLayout& right) const override;

    virtual PHX::Device::size_type rank() const override;

    virtual PHX::Device::size_type dimension(size_type ordinal) const override;

    virtual PHX::Device::size_type extent(size_type ordinal) const override;

    virtual int extent_int(size_type ordinal) const override;

    virtual void dimensions(std::vector<PHX::Device::size_type>& dim) const override;

    virtual std::string name(size_type ordinal) const override;

    virtual void names(std::vector<std::string>& names) const override;

    virtual PHX::Device::size_type size() const override;

    virtual PHX::DataLayout::KokkosLayoutType kokkosLayout() const override;

    virtual std::string identifier() const override;

    virtual void print(std::ostream& os, int offset) const override;

  private:

    std::string createIdentifier(const std::string& prefix = "");

    template<typename IndexType>
    typename std::enable_if<std::is_signed<IndexType>::value>::type
    checkForValidRank(const IndexType& ordinal) const;

    /** \brief Specialization to remove compiler warnings about
	pointless comparison of unsigned integer with zero.
     */
    template<typename IndexType>
    typename std::enable_if<std::is_unsigned<IndexType>::value>::type
    checkForValidRank(const IndexType& ordinal) const;

  protected:

    virtual void
    setExtentsOnDerivedClass(const std::vector<PHX::Device::size_type>& extents) override;

  private:

    std::vector<std::string> m_dim_name;

    std::vector<PHX::Device::size_type> m_dim_size;

    PHX::Device::size_type m_size;

    std::string m_identifier;

  };

}

#include "Phalanx_DataLayout_MDALayout_Def.hpp"

#endif
