// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_MaterialModelEntry.hpp"

//=======================================================================
//=======================================================================
panzer::MaterialModelEntry::MaterialModelEntry()
{ }

//=======================================================================
//=======================================================================
panzer::MaterialModelEntry::
MaterialModelEntry(const std::string factory_name) :
  m_factory_name(factory_name)
{ }

//=======================================================================
//=======================================================================
panzer::MaterialModelEntry::
MaterialModelEntry(const std::string factory_name,
		   const Teuchos::ParameterList& p) :
  m_factory_name(factory_name),
  m_params(p)
{ }

//=======================================================================
//=======================================================================
std::string panzer::MaterialModelEntry::factoryName() const
{
  return m_factory_name;
}

//=======================================================================
//=======================================================================
const Teuchos::ParameterList& panzer::MaterialModelEntry::params() const
{
  return m_params;
}

//=======================================================================
//=======================================================================
void panzer::MaterialModelEntry::operator=(const MaterialModelEntry& e)
{
  m_factory_name = e.m_factory_name;
  m_params = e.m_params;
}

//=======================================================================
//=======================================================================
bool panzer::MaterialModelEntry::operator==(const MaterialModelEntry& e) const
{
  return ( (m_factory_name == e.m_factory_name) &&
	   (this->m_params == e.m_params) );
}

//=======================================================================
//=======================================================================
bool panzer::MaterialModelEntry::operator!=(const MaterialModelEntry& e) const
{
  return (!(this->operator==(e)));
}

//=======================================================================
//=======================================================================
void panzer::MaterialModelEntry::print(std::ostream& os) const
{
  os << "Material Model Entry: " << m_factory_name;
}

//=======================================================================
//=======================================================================
std::ostream& 
panzer::operator<<(std::ostream& os, panzer::MaterialModelEntry& m)
{
  m.print(os);
  return os;
}
