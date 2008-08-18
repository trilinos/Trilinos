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

#ifndef PHX_DIMTAG_HPP
#define PHX_DIMTAG_HPP

#include "Phalanx_ConfigDefs.hpp"

namespace PHX {

  struct DimTag {
    virtual const char * name() const = 0 ;
    
    virtual std::string to_string( size_t , int ) const ;
    virtual int          to_index( size_t , const std::string & ) const ;
    
    // A derived type must provide the following static method:
    //
    // static const Derived_DimTag_Type & descriptor();
    
  protected:
    virtual ~DimTag();
    DimTag() {}
  private:
    DimTag( const DimTag & );
    DimTag & operator = ( const DimTag & );
  };

} 

#endif 
