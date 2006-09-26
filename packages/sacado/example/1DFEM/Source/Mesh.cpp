// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Mesh.hpp"

Mesh::Mesh() :
  elements()
{
}

Mesh::~Mesh()
{
}

void
Mesh::addElement(const Teuchos::RefCountPtr<AbstractElement>& element)
{
  elements.push_back(element);
}

unsigned int
Mesh::numElements() const
{
  return elements.size();
}

Mesh::iterator
Mesh::begin()
{
  return elements.begin();
}

Mesh::const_iterator
Mesh::begin() const
{
  return elements.begin();
}

Mesh::iterator
Mesh::end()
{
  return elements.end();
}

Mesh::const_iterator
Mesh::end() const
{
  return elements.end();
}
