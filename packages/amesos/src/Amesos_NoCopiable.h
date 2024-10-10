/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef AMESOS_NOCOPIABLE_H
#define AMESOS_NOCOPIABLE_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

/*!
 \class Amesos_NoCopiable
 
 \brief Amesos_NoCopiable: Simple class to prevent the usage of
     copy constructor and operator =

 \author Marzio Sala, SNL 9214

 \date Last updated on 24-May-05 (Champions' League Final day)
*/
class Amesos_NoCopiable
{
public:
  //! Default constructor
  Amesos_NoCopiable() {}

  //! Default destructor
  ~Amesos_NoCopiable() {}

private:
  //! Copy constructor
  Amesos_NoCopiable(const Amesos_NoCopiable& rhs);

  //! Copy destructor
  Amesos_NoCopiable& operator=(const Amesos_NoCopiable& rhs);
};
#endif
