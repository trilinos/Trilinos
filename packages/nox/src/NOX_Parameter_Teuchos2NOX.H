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

#ifdef HAVE_NOX_ANY

namespace NOX {

  namespace Parameter {

    /**
     * Teuchos2NOX converts a Teuchos parameter list to a NOX parameter list
     *
     * @author Kevin Long
     */
    class Teuchos2NOX
    {
    public:
      /** Construct a Teuchos2NOX converter */
      Teuchos2NOX(){;}
      
      /** */
      Parameter::List toNOX(const Teuchos::ParameterList& tpl) const ;
    protected:
      
   
    }; // class Teuchos2NOX
  } // namespace Parameter
} // namespace NOX


#endif // HAVE_NOX_ANY


#endif // NOX_PARAMETER_TEUCHOS2NOX_H
