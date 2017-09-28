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


#ifndef PHX_DATA_LAYOUT
#define PHX_DATA_LAYOUT

#include <iostream>
#include <vector>
#include <string>
#include "Phalanx_KokkosDeviceTypes.hpp"

namespace PHX{

  template<typename... extent_pack>
  struct SetExtentsImpl;

  //! Used to set the extents from a parameter pack. Can't use a simple initializer list due to narrowing from in size_type.
  template<typename T, typename... extent_pack>
  struct SetExtentsImpl<T,extent_pack...> {
    static void setExtents(PHX::Device::size_type index,
                           std::vector<PHX::Device::size_type>& extents,
                           T val,
                           extent_pack... pack)
    {
      extents[index] = val;
      PHX::SetExtentsImpl<extent_pack...>::setExtents(index+1,extents,pack...);
    }
  };

  //! Used to set the extents from a parameter pack. this implementation ends the recursion.
  template<>
  struct SetExtentsImpl<> {
    static void setExtents(PHX::Device::size_type ,
                           std::vector<PHX::Device::size_type>& )
    {}
  };

  // ******************************************************************
  /*! \brief A pure virtual class to provide size and rank information
      and a unique identifier for a fields.

      The DataLayout class is used to (1) specify the rank and extents
      of fields, and (2) to provide a unique identifier that can be
      used to differentiate fields.  For example suppose we want to
      store density at both the basis points and the quadrature points
      in a cell.  If we use the same string name for the field in the
      FieldTag, the DataLayout could be used to differentiate the objects.

      NOTE: We could probably just use an enumerated type here, but
      the DataLayout class allows users to derive and pass in
      auxiliary data via the tag.
  */
  class DataLayout {

  public:

    // typedef long unsigned int size_type;
    typedef PHX::Device::size_type size_type;

    DataLayout() = default;

    virtual ~DataLayout() = default;

    virtual PHX::Device::size_type rank() const = 0; 

    virtual PHX::Device::size_type dimension(size_type ordinal) const = 0; 

    virtual PHX::Device::size_type extent(size_type ordinal) const = 0; 

    virtual int extent_int(size_type ordinal) const = 0; 

    virtual void 
    dimensions(std::vector<PHX::Device::size_type>& dim) const = 0; 

    virtual PHX::Device::size_type size() const = 0;

    virtual std::string name(size_type ordinal) const = 0;

    virtual void names(std::vector<std::string>& names) const = 0; 

    virtual bool operator==(const DataLayout& left) const = 0;

    virtual bool operator!=(const DataLayout& left) const
    { return !(*this == left); }

    //! Unique name identifier that can be used for strict weak ordering in stl std::map keys.
    virtual std::string identifier() const = 0;

    virtual void print(std::ostream& os, int indent = 0) const = 0;

    template<typename... extent_pack>
    void setExtents(extent_pack... extents)
    {
      static_assert(sizeof...(extents) > 0,
                    "Error - PHX::Layout - rank must be greater than zero!");
      static_assert(sizeof...(extents) < 9,
                    "Error - PHX::Layout - rank must be less than 9!");

      std::vector<PHX::Device::size_type> extents_vec(sizeof...(extents));
      PHX::SetExtentsImpl<extent_pack...>::setExtents(0,extents_vec,extents...);
      this->setExtentsOnDerivedClass(extents_vec);
    }

  protected:

    virtual void setExtentsOnDerivedClass(const std::vector<PHX::Device::size_type>& extents) = 0;

  };

  std::ostream& operator<<(std::ostream& os, const PHX::DataLayout& t);
  
}

#endif
