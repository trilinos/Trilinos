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

#ifndef PHX_DATA_LAYOUT
#define PHX_DATA_LAYOUT

#include <iostream>
#include <vector>
#include <string>

namespace PHX{


  /*! \brief A pure virtual class to distinguish a unique data layout in a cell.

      The DataLayout class is used to (1) specify the array size of a
      an algebraic type in a single cell, and (2) to differentiate
      FieldTags that have the same name, but have different
      DataLayouts in the FieldManager.  For example suppose we want to
      store density at both the nodes and the quadrature points in a
      cell.  If we use the same string name for the FieldTag, the
      DataLayout will differentiate the objects.  We could probably
      just use an enumerated type here, but the DataLayout class
      allows users to derive and pass in auxiliary data via the tag.

  */
  class DataLayout {

  public:

    typedef int size_type;

    DataLayout() {}

    virtual ~DataLayout() {}

    virtual size_type rank() const = 0; 

    virtual size_type dimension(size_type ordinal) const = 0; 

    virtual void dimensions(std::vector<size_type>& dim) const = 0; 

    //! Returns the name of the input ordinal
    virtual std::string name(size_type ordinal) const = 0;

    //! Returns the names of all ordinals in a vector
    virtual void names(std::vector<std::string>& names) const = 0; 

    virtual size_type size() const = 0;

    virtual bool operator==(const DataLayout& left) const = 0;

    virtual bool operator!=(const DataLayout& left) const
    { return !(*this == left); }

    //! Unique name identifier that can be used for strict weak ordering in stl std::map keys.
    virtual std::string identifier() const = 0;

    virtual void print(std::ostream& os, int indent = 0) const = 0;

  };

  std::ostream& operator<<(std::ostream& os, const PHX::DataLayout& t);
  
}

#endif
