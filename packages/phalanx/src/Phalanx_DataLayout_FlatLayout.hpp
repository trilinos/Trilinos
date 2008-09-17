// @HEADER
// ************************************************************************
// 
//            Phalanx: A Partial Differential Equation Assembly 
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

#ifndef PHX_DATA_LAYOUT_FLAT_LAYOUT
#define PHX_DATA_LAYOUT_FLAT_LAYOUT

#include "Phalanx_DataLayout.hpp"
#include <string>

namespace PHX{

  class FlatLayout : public DataLayout {

  public:

    FlatLayout(const std::string& unique_identifier, std::size_t i);

    virtual ~FlatLayout();

    virtual bool operator==(const DataLayout& right) const;

    virtual const std::string& name() const;

    virtual size_type rank() const; 

    virtual void dimensions(std::vector<size_type>& dim) const; 

    virtual size_type size() const;

    virtual std::string identifier() const;

    virtual void print(std::ostream& os, int indent = 0) const;

  private:

    const std::string m_name;

    const size_type m_size;

  };

  std::ostream& operator<<(std::ostream& os, const PHX::FlatLayout& t);

}

#endif
