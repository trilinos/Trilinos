//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_PARAMETER_TEUCHOS2NOX_H
#define NOX_PARAMETER_TEUCHOS2NOX_H

#include "NOX_Common.H" // for NOX_Config.h
#include "NOX_Parameter_List.H"		
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace NOX {

  namespace Parameter {

    /**
     * Teuchos2NOX has methods to convert a Teuchos parameter list 
     * to a NOX parameter list, as well as the reverse. There are
     * also methods to read and write XML files using the converters
     * and the Teuchos parameter list utilities.
     *
     * @author Kevin Long
     */
    class Teuchos2NOX
    {
    public:
      /** Construct a Teuchos2NOX converter */
      Teuchos2NOX(){;}
      
      /** Utility to Convert a Teuchos::ParameterList to a NOX::Parameter::List */
      Teuchos::RefCountPtr<NOX::Parameter::List> toNOX(const Teuchos::ParameterList& tpl) const ;

      /** Utility to Convert a NOX::Parameter::List to a Teuchos::ParameterList */
      Teuchos::ParameterList toTeuchos(const NOX::Parameter::List& npl) const ;

#ifdef HAVE_TEUCHOS_EXPAT
      /** Utility to write a NOX::Parameter::List to an XML file (via Teuchos) */
      void SaveToXMLFile(const string filename, const NOX::Parameter::List& npl) const;

      /** Utility to read a NOX::Parameter::List from an XML file (via Teuchos) */
      Teuchos::RefCountPtr<NOX::Parameter::List> ReadFromXMLFile(const string filename) const;
#endif //#ifdef HAVE_TEUCHOS_EXPAT

    protected:
   
    }; // class Teuchos2NOX
  } // namespace Parameter
} // namespace NOX

#endif // NOX_PARAMETER_TEUCHOS2NOX_H
