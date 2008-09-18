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

#include <sstream>
#include <typeinfo>
#include "Phalanx_DataLayout_FlatLayout.hpp"

//**********************************************************************
PHX::FlatLayout::FlatLayout(const std::string& unique_identifier, 
			    std::size_t size) :
  m_name(unique_identifier),
  m_size(size)
{ }

//**********************************************************************
PHX::FlatLayout::~FlatLayout()
{ }

//**********************************************************************
bool PHX::FlatLayout::operator==(const PHX::DataLayout& right) const
{
  const PHX::FlatLayout* tmp = 0;
  tmp = dynamic_cast< const PHX::FlatLayout* >(&right);
  
  if (tmp == 0)
    return false;
  
  return (  (this->name() == tmp->name()) &&
	    (this->size() == tmp->size()) );
}

//**********************************************************************
const std::string& PHX::FlatLayout::name() const
{ return m_name; }

//**********************************************************************
PHX::DataLayout::size_type PHX::FlatLayout::rank() const
{ return 1; }

//**********************************************************************
void PHX::FlatLayout::
dimensions(std::vector<PHX::DataLayout::size_type>& dim) const
{ 
  dim.resize(1);
  dim[0] = m_size;
}

//**********************************************************************
PHX::DataLayout::size_type PHX::FlatLayout::size() const
{ return m_size; }

//**********************************************************************
std::string PHX::FlatLayout::identifier() const
{ 
  std::ostringstream ost;
  ost << this->name() << this->size();
  return ost.str(); 
}

//**********************************************************************
void PHX::FlatLayout::print(std::ostream& os, int indent) const
{
  std::ostringstream s;
  for (int i = 0; i < indent; i++)
    s << " ";

  os << s.str() << m_name << ", size = " << m_size;
}

//**********************************************************************
std::ostream& PHX::operator<<(std::ostream& os, const PHX::FlatLayout& v)
{
  v.print(os);
  return os;
}

//**********************************************************************
