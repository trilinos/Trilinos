// @HEADER
// ************************************************************************
//
//                  Galeri Matrix Generation Package
//                 Copyright (2006) Sandia Corporation
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
// ************************************************************************
// @HEADER

#ifndef PHX_GRID_ELEMENT_H
#define PHX_GRID_ELEMENT_H

#include "Galeri_ConfigDefs.h"
#include "phx_core_Object.h"
#include <vector>

namespace phx {
namespace grid {

class Element : public phx::core::Object
{
public: 
  // @{ \name constructor and destructor
  //! Default constructor.
  Element() :
    numVertices_(0),
    numComponents_(0)
  {}

  //! Copy constructor.
  Element(const Element& rhs)
  {
    setNumVertices(rhs.getNumVertices());
    setNumComponents(rhs.getNumComponents());
    for (int i = 0; i < rhs.getNumComponents(); ++i)
      setComponent(i, rhs.getComponent(i));
    setLabel(rhs.getLabel());
  }

  Element& operator=(const Element& rhs)
  {
    setNumVertices(rhs.getNumVertices());
    setNumComponents(rhs.getNumComponents());
    for (int i = 0; i < rhs.getNumComponents(); ++i)
      setComponent(i, rhs.getComponent(i));
    setLabel(rhs.getLabel());

    return(*this);
  }

  //! destructor.
  ~Element() {}

  // @}
  // @{ \name Get methods

  int getNumVertices() const
  {
    return(numVertices_);
  }

  int getNumComponents() const 
  {
    return(numComponents_);
  }

  const Element& getComponent(const int which) const
  {
    return(components_[which]);
  }

  // @}
  // @{ \name Set methods

  void setNumVertices(const int numVertices)
  {
    numVertices_ = numVertices;
  }

  void setNumComponents(const int numComponents)
  {
    numComponents_ = numComponents;
    components_.resize(numComponents);
  }

  void setComponent(const int which, const Element& what)
  {
    components_[which] = what;
  }

  virtual void print(ostream & os) const
  {
    cout << endl;
    cout << "Number of vertices = " << getNumVertices() << endl;
    cout << "Number of components = " << getNumComponents() << endl;
    for (int i = 0; i < getNumComponents(); ++i)
      cout << "Label of component " << i << " = " <<
        getComponent(i).getLabel() << endl;
    cout << endl;
  }

  // @}
private:
  // @{ \name Private data and methods

  int numVertices_;
  int numComponents_;
  vector<Element> components_;

  // @}
}; // class Element

} // namespace grid
} // namespace phx
#endif // PHX_BASE_ELEMENT_H
