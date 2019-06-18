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


#ifndef PHX_DATA_LAYOUT_LAYOUT
#define PHX_DATA_LAYOUT_LAYOUT

#include "Phalanx_DataLayout.hpp"
#include <string>
#include <iostream>
#include <vector>

namespace PHX {

  //! Default DataLayout implementation that allows for runtime sizing.
  class Layout : public DataLayout {

  public:
    using KokkosLayout = PHX::DataLayout::KokkosLayoutType;

    Layout(const std::string& id = "");

    template<typename... extent_pack>
    //Layout(extent_pack... extents) : m_extents{extents...}
    Layout(const std::string& id, extent_pack... extents) :
      m_identifier(id),
      m_extents(sizeof...(extents)),
      m_kokkos_layout_type(KokkosLayout::Default)
    {
      static_assert(sizeof...(extents) > 0,
                    "Error - PHX::Layout - rank must be greater than zero!");
      static_assert(sizeof...(extents) < 9,
                    "Error - PHX::Layout - rank must be less than 9!");

      PHX::SetExtentsImpl<extent_pack...>::setExtents(0,m_extents,extents...);
    }

    template<typename... extent_pack>
    void setExtents(extent_pack... extents)
    {
      static_assert(sizeof...(extents) > 0,
                    "Error - PHX::Layout - rank must be greater than zero!");
      static_assert(sizeof...(extents) < 9,
                    "Error - PHX::Layout - rank must be less than 9!");

      //m_extents = {extents...};
      m_extents.resize(sizeof...(extents));
      PHX::SetExtentsImpl<extent_pack...>::setExtents(0,m_extents,extents...);
    }

    virtual void setKokkosLayout(const PHX::DataLayout::KokkosLayoutType& klt);

    virtual ~Layout() noexcept {}

    virtual bool operator==(const DataLayout& src) const override;

    virtual PHX::Device::size_type rank() const override;

    virtual PHX::Device::size_type dimension(size_type ordinal) const override;

    virtual PHX::Device::size_type extent(size_type ordinal) const override;

    virtual int extent_int(size_type ordinal) const override;

    virtual void dimensions(std::vector<PHX::Device::size_type>& dim) const override;

    virtual PHX::Device::size_type size() const override;

    virtual std::string name(size_type ordinal) const override;

    virtual void names(std::vector<std::string>& names) const override;

    virtual PHX::DataLayout::KokkosLayoutType kokkosLayout() const override;

    virtual std::string identifier() const override;

    virtual void print(std::ostream& os, int offset) const override;

  protected:

    virtual void
    setExtentsOnDerivedClass(const std::vector<PHX::Device::size_type>& extents) override;

  private:

    std::string m_identifier;
    std::vector<PHX::Device::size_type> m_extents;
    PHX::DataLayout::KokkosLayoutType m_kokkos_layout_type;
  };

  std::ostream& operator<<(std::ostream& os, const PHX::Layout& t);

}

#endif
