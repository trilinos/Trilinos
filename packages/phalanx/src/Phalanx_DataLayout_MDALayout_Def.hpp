// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_DATA_LAYOUT_MDALAYOUT_DEF
#define PHX_DATA_LAYOUT_MDALAYOUT_DEF

#include <iostream>
#include <sstream>
#include <typeinfo>
#include "Teuchos_Assert.hpp"
#include "Phalanx_ExtentTraits.hpp"

namespace PHX {

//**********************************************************************
template<int N, typename CurrentExtent, typename... Extents>
void setExtentsVariadic(std::vector<PHX::Device::size_type>& e, CurrentExtent ce, Extents... extents)
{
  e[N] = ce;

  if constexpr (sizeof...(extents) != 0) {
    setExtentsVariadic<N+1>(e,extents...);
  }
}

//**********************************************************************
template<int N, typename CurrentTag, typename... Tags>
void setTagNames(std::vector<std::string>& names)
{
  names.push_back(PHX::print<CurrentTag>());

  if constexpr (sizeof...(Tags) != 0) {
    setTagNames<N+1,Tags...>(names);
  }
}

} // namespace PHX

//**********************************************************************
template<typename... Tags>
template<typename... Extents>
PHX::MDALayout<Tags...>::MDALayout(Extents... extents)
{
  static_assert(Rank == PackSize<Extents...>::value);

  m_dim_size.resize(Rank);
  setExtentsVariadic<0>(m_dim_size,extents...);

  m_dim_name.clear();
  setTagNames<0,Tags...>(m_dim_name);

  m_size = 1;
  for (int i=0; i < Rank; ++i)
    m_size *= m_dim_size[i];

  m_identifier = this->createIdentifier();
}

//**********************************************************************
template<typename... Tags>
template<typename... Extents>
PHX::MDALayout<Tags...>::
MDALayout(const std::string& prefix, Extents... extents)
  : MDALayout<Tags...>::MDALayout(extents...)
{
  m_identifier = this->createIdentifier(prefix);
}

//**********************************************************************
template<typename... Tags>
template<typename... Extents>
PHX::MDALayout<Tags...>::
MDALayout(const char* prefix, Extents... extents)
  : MDALayout<Tags...>::MDALayout(extents...)
{
  m_identifier = this->createIdentifier(std::string(prefix));
}

//**********************************************************************
template<typename... Tags>
bool PHX::MDALayout<Tags...>::
operator==(const PHX::DataLayout& right) const
{
  const PHX::MDALayout<Tags...>* tmp = 0;
  tmp = dynamic_cast< const PHX::MDALayout<Tags...>* >(&right);

  if (tmp == 0)
    return false;

  for (size_type i=0; i < Rank; ++i)
    if (m_dim_size[i] != tmp->m_dim_size[i])
      return false;

  return (this->size() == tmp->size());
}

//**********************************************************************
template<typename... Tags>
PHX::Device::size_type
PHX::MDALayout<Tags...>::rank() const
{ return Rank; }

//**********************************************************************
template<typename... Tags>
void PHX::MDALayout<Tags...>::
dimensions(std::vector<PHX::Device::size_type>& dim) const
{
  dim.resize(Rank);
  for(std::size_t i=0; i < Rank; ++i)
    dim[i] = m_dim_size[i];
}

//**********************************************************************
template<typename... Tags>
void PHX::MDALayout<Tags...>::
names(std::vector<std::string>& names) const
{
  names.resize(Rank);
  for(std::size_t i=0; i < names.size(); ++i)
    names[i] = m_dim_name[i];
}

//**********************************************************************
template<typename... Tags>
PHX::Device::size_type PHX::MDALayout<Tags...>::
size() const
{ return m_size; }

//**********************************************************************
template<typename... Tags>
PHX::DataLayout::KokkosLayoutType
PHX::MDALayout<Tags...>::
kokkosLayout() const
{
  return PHX::DataLayout::KokkosLayoutType::Default;
}

//**********************************************************************
template<typename... Tags>
std::string PHX::MDALayout<Tags...>::
identifier() const
{
  return m_identifier;
}

//**********************************************************************
template<typename... Tags>
PHX::Device::size_type
PHX::MDALayout<Tags...>::
dimension(size_type ordinal) const
{
#ifdef PHX_DEBUG
  this->checkForValidRank(ordinal);
#endif
  return m_dim_size[ordinal];
}

//**********************************************************************
template<typename... Tags>
PHX::Device::size_type
PHX::MDALayout<Tags...>::
extent(size_type ordinal) const
{
#ifdef PHX_DEBUG
  this->checkForValidRank(ordinal);
#endif
  return m_dim_size[ordinal];
}

//**********************************************************************
template<typename... Tags>
int
PHX::MDALayout<Tags...>::
extent_int(size_type ordinal) const
{
#ifdef PHX_DEBUG
  this->checkForValidRank(ordinal);
#endif
  return static_cast<int>(m_dim_size[ordinal]);
}

//**********************************************************************
template<typename... Tags>
std::string
PHX::MDALayout<Tags...>::
name(size_type ordinal) const
{
#ifdef PHX_DEBUG
  this->checkForValidRank(ordinal);
#endif
  return m_dim_name[ordinal];
}

//**********************************************************************
template<typename... Tags>
void PHX::MDALayout<Tags...>::
print(std::ostream& os, int /* offset */) const
{
  os << m_identifier;
}

//**********************************************************************
template<typename... Tags>
std::string PHX::MDALayout<Tags...>::
createIdentifier(const std::string& prefix)
{
  std::ostringstream os;
  os << prefix << "<";
  for (std::size_t i=0; i < m_dim_name.size(); ++i) {
    if (i > 0)
      os << ",";
    os << std::string(m_dim_name[i]);
  }
  os << ">(";
  for (size_type i=0; i < Rank; ++i) {
    if (i > 0)
      os << ",";
    os << m_dim_size[i];
  }
  os << ")";

  return os.str();
}

//**********************************************************************
template<typename... Tags>
void
PHX::MDALayout<Tags...>::
setExtentsOnDerivedClass(const std::vector<PHX::Device::size_type>& extents)
{
#ifdef PHX_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(extents.size() != Rank, std::runtime_error,
                             "ERROR - MDALayout::setExtentsOnDerivedClass()"
                             << " - Rank mismatch while setting extents!");
#endif
  for (std::size_t i=0; i < extents.size(); ++i)
    m_dim_size[i] = extents[i];
}

//**********************************************************************
template<typename... Tags>
template<typename IndexType>
typename std::enable_if<std::is_signed<IndexType>::value>::type
PHX::MDALayout<Tags...>::
checkForValidRank(const IndexType& ordinal) const
{
  if ( (ordinal >= Rank) || (ordinal < 0) ) {
    std::ostringstream os;
    os << "Requested Ordinal " << ordinal
       << " is outside the valid range of 0 - " << Rank - 1
       << " in DataLayout object:\n"
       << m_identifier << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(ordinal >= Rank || ordinal < 0,
			       std::runtime_error, os.str());
  }
}

//**********************************************************************
template<typename... Tags>
template<typename IndexType>
typename std::enable_if<std::is_unsigned<IndexType>::value>::type
PHX::MDALayout<Tags...>::
checkForValidRank(const IndexType& ordinal) const
{
  if (ordinal >= Rank) {
    std::ostringstream os;
    os << "Requested Ordinal " << ordinal
       << " is outside the valid range of 0 - " << Rank - 1
       << " in DataLayout object:\n"
       << m_identifier << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(ordinal >= Rank,
			       std::runtime_error, os.str());
  }
}

//**********************************************************************

#endif
