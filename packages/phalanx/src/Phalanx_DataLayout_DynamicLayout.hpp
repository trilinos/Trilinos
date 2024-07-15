// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
