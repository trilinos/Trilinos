// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_XMLParameterListWriter.hpp"	

using namespace Teuchos;

XMLParameterListWriter::XMLParameterListWriter()
{;}


XMLObject XMLParameterListWriter::toXML(const ParameterList& p) const
{
  XMLObject rtn("ParameterList");
  
  for (ParameterList::ConstIterator i=p.begin(); i!=p.end(); ++i)
    {
      const ParameterEntry& val = p.entry(i);
      const string& name = p.name(i);
      XMLObject child = toXML(val);
      child.addAttribute("name", name);
      rtn.addChild(child);
    }

  return rtn;
}

XMLObject XMLParameterListWriter::toXML(const ParameterEntry& entry) const
{
  if (entry.isList())
    {
      return toXML(getValue<ParameterList>(entry));
    }

  XMLObject rtn("Parameter");
  string type;
  string value;

  if (entry.isType<int>())
    {
      type = "int";
      value = toString(getValue<int>(entry));
    }
  else if (entry.isType<double>())
    {
      type = "double";
      value = toString(getValue<double>(entry));
    }
  else if (entry.isType<float>())
    {
      type = "float";
      value = toString(getValue<float>(entry));
    }
  else if (entry.isType<string>())
    {
      type = "string";
      value = getValue<string>(entry);
    }
  else if (entry.isType<char>())
    {
      type = "char";
      value = toString(getValue<char>(entry));
    }
  else if (entry.isType<bool>())
    {
      type = "bool";
      value = toString(getValue<bool>(entry));
    }
  else
    {
      type = "any";
      TeuchosOStringStream ss;
      ss << entry;
      value = TEUCHOS_OSTRINGSTREAM_GET_C_STR(ss);
    }
  

  rtn.addAttribute("type", type);
  rtn.addAttribute("value", value);
  
  if (entry.isDefault())
    {
      rtn.addAttribute("isDefault", "true");
    }

  return rtn;
}
