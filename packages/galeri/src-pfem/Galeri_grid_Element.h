// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

/*! 
 * \file Galeri_grid_Element.h
 *
 * \brief Base class for all grid elements.
 *
 * \author Marzio Sala, ETHZ
 *
 * \date Last modified on Aug-06
 */

#ifndef GALERI_GRID_ELEMENT_H
#define GALERI_GRID_ELEMENT_H

#include "Galeri_ConfigDefs.h"
#include "Galeri_core_Object.h"
#include "Teuchos_TestForException.hpp"
#include <vector>

namespace Galeri {
namespace grid {

/*! \class Galeri::grid::Element  
 *
 * \brief Basic class for grid elements.
 *
 * Class Galeri::grid::Element specifies the interface for a generic 1D, 2D or
 * 3D grid element. This class is a semi-virtual class.
 *
 * In Galeri/pfem, an element is defined as a geometric object, composed by
 * <i>vertices</i> and <i>components</i>. The former are easy to define; the
 * latter, instead, define the sub-entities of the element, and they are
 * either segments (for 2D elements) or faces (for 3D elements). 
 *
 * The idea is that you can define an object by deriving this class, and set
 * the correct number of vertices and components. You also have to decide 
 * the local ID of each of the components. Mixed components are perfectly
 * legal. Since a component is nothing but an Element-derived class,
 * you can recursively assemble components and create the object you need for
 * your discretization.
 *
 * This class is derived by Galeri::grid::Point, Galeri::grid::Segment,
 * Galeri::grid::Triangle, Galeri::grid::Quad, 
 * Galeri::grid::Tet and Galeri::grid::Hex. New elements can be easily if
 * required.
 */
class Element : public Galeri::core::Object
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

  //! operator =
  virtual Element& operator=(const Element& rhs)
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
  // @{ \name Get/Set methods

  //! Sets the number of vertices in \c this object.
  void setNumVertices(const int numVertices)
  {
    numVertices_ = numVertices;
  }

  //! Gets the number of vertices associated with \c this object.
  int getNumVertices() const
  {
    return(numVertices_);
  }

  //! Sets the number of components in \c this object.
  void setNumComponents(const int numComponents)
  {
    numComponents_ = numComponents;
    components_.resize(numComponents);
  }

  //! Gets the number of components associated with \c this object.
  int getNumComponents() const 
  {
    return(numComponents_);
  }

  //! Sets the element type for component \c which.
  void setComponent(const int which, const Element& what)
  {
    TEST_FOR_EXCEPTION(which >= components_.size(), std::logic_error,
                       "accessing element " << which << ", but the "
                       << " element only has " <<  components_.size()
                       << " components.");

    components_[which] = what;
  }

  //! Gets the element type of the specified component.
  const Element& getComponent(const int which) const
  {
    return(components_[which]);
  }

  // @}
  // @{ Other methods
  
  //! Prints the output on \c os.
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
  // @{ \name Private data and methods
private:

  //! Number of vertices.
  int numVertices_;
  //! Number of components.
  int numComponents_;
  //! Vector containing the element of each component.
  vector<Element> components_;

  // @}
}; // class Element

} // namespace grid
} // namespace Galeri

#endif // GALERI_BASE_ELEMENT_H
