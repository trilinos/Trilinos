// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_UNIT_TEST_HELPERS_HPP
#define STOKHOS_UNIT_TEST_HELPERS_HPP

namespace Stokhos {

  template<class OrdinalType, class ValueType>
  bool comparePCEs(const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& a1, 
		   const std::string& a1_name,
		   const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>&a2, 
		   const std::string& a2_name,
		   const ValueType& rel_tol, const ValueType& abs_tol,
		   Teuchos::FancyOStream& out)
  {
    bool success = true;
    
    out << "Comparing " << a1_name << " == " << a2_name << " ... ";
    
    const OrdinalType n = a1.size();
    
    // Compare sizes
    if (a2.size() != n) {
      out << "\nError, "<<a1_name<<".size() = "<<a1.size()<<" == " 
	  << a2_name<<".size() = "<<a2.size()<<" : failed!\n";
      return false;
    }
    
    // Compare elements
    for( OrdinalType i = 0; i < n; ++i ) {
      ValueType err = std::abs(a1[i] - a2[i]);
      ValueType tol = 
	abs_tol + rel_tol*std::max(std::abs(a1[i]),std::abs(a2[i]));
      if (err  > tol) {
	out
	  <<"\nError, relErr("<<a1_name<<"["<<i<<"],"
	  <<a2_name<<"["<<i<<"]) = relErr("<<a1[i]<<","<<a2[i]<<") = "
	  <<err<<" <= tol = "<<tol<<": failed!\n";
	success = false;
      }
    }
    if (success) {
      out << "passed\n";
    }
    else {
      out << std::endl 
	  << a1_name << " = " << a1 << std::endl
	  << a2_name << " = " << a2 << std::endl;
    }
    
    return success;
  }
}

#endif // STOKHOS_UNIT_TEST_HELPERS_HPP
