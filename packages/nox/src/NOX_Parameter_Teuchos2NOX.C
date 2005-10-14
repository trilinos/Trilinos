//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "NOX_Parameter_Teuchos2NOX.H"

// Included for XML stuff
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_XMLParameterListReader.hpp"

// Teuchos typedefs
typedef Teuchos::map<string, Teuchos::ParameterEntry> tMap;
typedef tMap::iterator tIterator;
typedef tMap::const_iterator tConstIterator;

// NOX typedefs
typedef map<string, NOX::Parameter::Entry> nMap;
typedef nMap::iterator nIterator;
typedef nMap::const_iterator nConstIterator;

using namespace NOX::Parameter;
using namespace Teuchos;

RefCountPtr<NOX::Parameter::List> Teuchos2NOX::toNOX(const ParameterList& p) const
{
  RefCountPtr<NOX::Parameter::List> rtn = rcp(new NOX::Parameter::List());

  for (ParameterList::ConstIterator i=p.begin(); i!=p.end(); ++i)
    {
      const ParameterEntry& val = p.entry(i);
      const string& name = p.name(i);
      
      if (val.isList())
        {
          rtn->sublist(name) = *(toNOX(getValue<ParameterList>(val)));
        }
      else if (val.isType<int>())
        {
          rtn->setParameter(name, getValue<int>(val));
        }
      else if (val.isType<double>())
        {
          rtn->setParameter(name, getValue<double>(val));
        }
      else if (val.isType<bool>())
        {
          rtn->setParameter(name, getValue<bool>(val));
        }
      else if (val.isType<string>())
        {
          rtn->setParameter(name, getValue<string>(val));
        }
      else
        {
          any data = val.getAny();
          rtn->setParameter(name, AnyPtr(data));
        }
    }

  return rtn;
}

ParameterList Teuchos2NOX::toTeuchos(const NOX::Parameter::List& npl) const
{
  ParameterList rtn;

  for (nConstIterator i=npl.begin(); i!=npl.end(); ++i)
    {
      const Entry& val = npl.entry(i);
      const string& name = npl.name(i);
      
      if (val.isList())
        {
          rtn.sublist(name) = toTeuchos(val.getListValue());
        }
      else if (val.isInt())
        {
	  rtn.set(name, val.getIntValue());
        }
      else if (val.isDouble())
        {
	  rtn.set(name, val.getDoubleValue());
        }
      else if (val.isBool())
        {
	  rtn.set(name, val.getBoolValue());
        }
      else if (val.isString())
        {
	  rtn.set(name, val.getStringValue());
        }
      else if (val.isArbitrary())
        {
	  // crj 9/29/05
	  // WARNING - Converting arbitrary values to Teuchos::ParameterList
	  //	has not been tested.  I don't know if it makes sense or not.
	  cerr << "Teuchos2NOX::toTeuchos - Warning: detected arbitrary value"
		<< " in NOX::Parameter::List.  This is untested and might not." 
		<< "  Make sense in what needs to be done in Teuchos."  << endl;
	  rtn.set(name, val.getArbitraryValue().clone());
        }
      else
        {
	  cerr << "NOX::Parameter::Teuchos2NOX::toTeuchos - unrecognized type "
		<< "name " << name << " value " << val << endl;
	  throw "NOX Error";
        }
    }

  return rtn;
}

#ifdef HAVE_TEUCHOS_EXPAT

void Teuchos2NOX::SaveToXMLFile(const string filename, const Parameter::List& npl) const
{
  /////////////////////////////////////////////////
  //  Convert the Parameter::List to an XML file
  ///////////////////////////////////////////////// 

  // create a parameter list converter
  NOX::Parameter::Teuchos2NOX pl_converter;

  // Convert the nox parameter list to a teuchos parameter list
  Teuchos::ParameterList tpl = pl_converter.toTeuchos(npl);

  // create a xml converter
  Teuchos::XMLParameterListWriter xml_converter;

  // Convert the teuchos parameter list to an XMLObject
  Teuchos::XMLObject xml_pl = xml_converter.toXML(tpl);

  // Write the xml to a file
  ofstream of(filename.c_str()); 
  of << xml_pl << endl;
  of.close();
}

RefCountPtr<NOX::Parameter::List> Teuchos2NOX::ReadFromXMLFile(const string filename) const
{
  // read in a file using teuchos
  Teuchos::FileInputSource fileSrc(filename);

  // Convert the file data into an xml object
  Teuchos::XMLObject xml_obj = fileSrc.getObject();

  // Create a xml to teuchos::parameterlist converter
  Teuchos::XMLParameterListReader xml_to_pl_converter;

  // convert the teuchos xml object to a parameter list
  Teuchos::ParameterList pl = xml_to_pl_converter.toParameterList(xml_obj);

  // create a parameter list to list converter
  NOX::Parameter::Teuchos2NOX pl_converter;

  // convert the parameter list to a nox list
  RefCountPtr<NOX::Parameter::List> l = pl_converter.toNOX(pl);

  return l;
}

#endif //#ifdef HAVE_TEUCHOS_EXPAT

