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

#include "Teuchos_Dependency.hpp"


namespace Teuchos{


Dependency::Dependency(DependeeList& dependees, DependeeList& dependents):
  dependees_(dependees), dependents_(dependents)
{}

Dependency::Dependency(DependeeList& dependees, RCP<ParameterEntry> dependent):
  dependees_(dependees), dependents(DependentList(1, dependent))
{}

Dependency::Dependency(RCP<ParameterEntry> dependee, DependentList dependents):
  dependees_(DependeeList(1, dependee)), dependents_(dependents)
{}
  
Dependency::Dependency(RCP<ParameterEntry> dependee, RCP<ParameterEntry> dependent)
  dependees_(DependeeList(1, dependee)), dependents_(DependentList(1, dependent))
{}


} //namespace Teuchos

