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


#ifndef PHX_DATA_LAYOUT_DYNAMIC_LAYOUT_DEF
#define PHX_DATA_LAYOUT_DYNAMIC_LAYOUT_DEF

#include "Phalanx_DataLayout_DynamicLayout.hpp"
#include "Teuchos_Assert.hpp"

//**********************************************************************
PHX::Layout::Layout(const std::string& name) : m_identifier(name) {}

//**********************************************************************
bool PHX::Layout::operator==(const PHX::DataLayout& src) const
{
  const auto* src_cast = dynamic_cast<const PHX::Layout*>(&src);
  if (src_cast == nullptr)
    return false;

  // Check identifier
  if (m_identifier != src_cast->m_identifier)
    return false;

  // Check rank
  if (m_extents.size() != src_cast->m_extents.size())
    return false;

  // Check dimensions
  for (size_t i=0; i < m_extents.size(); ++i)
    if (m_extents[i] != src_cast->m_extents[i])
      return false;

  return true;
}

//**********************************************************************
PHX::Device::size_type PHX::Layout::rank() const
{ return static_cast<PHX::Device::size_type>(m_extents.size()); }

//**********************************************************************
void PHX::Layout::dimensions(std::vector<PHX::Device::size_type>& dim) const
{
  dim.resize(m_extents.size());
  for(std::size_t i=0; i < m_extents.size(); ++i)
    dim[i] = m_extents[i];
}

//**********************************************************************
PHX::Device::size_type PHX::Layout::size() const
{
  PHX::Device::size_type my_size = 0;
  if (m_extents.size() > 0) {
    my_size = 1;
    for (const auto& i : m_extents) {
      my_size *= i;
    }
  }
  return my_size;
}

//**********************************************************************
std::string PHX::Layout::identifier() const
{
  return m_identifier;
}

//**********************************************************************
PHX::Device::size_type PHX::Layout::dimension(size_type ordinal) const
{
  return m_extents[ordinal];
}

//**********************************************************************
PHX::Device::size_type PHX::Layout::extent(size_type ordinal) const
{
  return m_extents[ordinal];
}

//**********************************************************************
int PHX::Layout::extent_int(size_type ordinal) const
{
  return static_cast<int>(m_extents[ordinal]);
}

//**********************************************************************
std::string PHX::Layout::name(size_type ordinal) const
{
  return "EXT"+std::to_string(ordinal);
}

//**********************************************************************
void PHX::Layout::names(std::vector<std::string>& names) const
{
  names.resize(m_extents.size());
  for (std::size_t i=0; i < m_extents.size(); ++i)
    names[i] = "EXT"+std::to_string(i);
}

//**********************************************************************
void PHX::Layout::print(std::ostream& os, int /* offset */) const
{
  os << m_identifier << "(";
  for (size_t i=0; i < m_extents.size(); ++i) {
    if (i != 0) os << ",";
    os << m_extents[i];
  }
  os << ")";
}

//**********************************************************************
void
PHX::Layout::
setExtentsOnDerivedClass(const std::vector<PHX::Device::size_type>& extents)
{
  m_extents = std::move(extents);
}

//**********************************************************************
std::ostream& PHX::operator<<(std::ostream& os, const PHX::Layout& v)
{
  v.print(os,0);
  return os;
}

//**********************************************************************

#endif
