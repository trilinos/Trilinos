// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef STOKHOS_IFPACK2_HPP
#define STOKHOS_IFPACK2_HPP

// PCE includes
#include "Stokhos_Sacado.hpp"

// Belos-PCE-Tpetra adapter
#include "BelosPCETpetraAdapter.hpp"

// Ifpack2 includes
#include "Ifpack2_RILUK.hpp"
#include "Ifpack2_ILUT.hpp"
#include "Ifpack2_Chebyshev.hpp"
#include "Ifpack2_Diagonal.hpp"
#include "Ifpack2_Relaxation.hpp"
#include "Ifpack2_Krylov.hpp"

namespace Ifpack2 {
  
  //! Specialization of BelosScalarType to PCE types
  template <typename T, typename S>
  struct BelosScalarType< Sacado::PCE::OrthogPoly<T,S> > {
    typedef T type;
  };

}

#endif
