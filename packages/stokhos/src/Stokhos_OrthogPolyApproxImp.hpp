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

template <typename T> 
Stokhos::OrthogPolyApprox<T>::
OrthogPolyApprox() :
  coeff_(1)
{
}

template <typename T> 
Stokhos::OrthogPolyApprox<T>::
OrthogPolyApprox(const T& x) :
  coeff_(1)
{
  coeff_[0] = x;
}

template <typename T> 
Stokhos::OrthogPolyApprox<T>::
OrthogPolyApprox(unsigned int sz, const T& x) :
  coeff_(sz)
{
  coeff_[0] = x;
}

template <typename T> 
Stokhos::OrthogPolyApprox<T>::
OrthogPolyApprox(unsigned int sz) :
  coeff_(sz)
{
}

template <typename T> 
Stokhos::OrthogPolyApprox<T>::
OrthogPolyApprox(unsigned int sz, unsigned int l) :
  coeff_(l)
{
  coeff_.resize(sz);
}

template <typename T> 
Stokhos::OrthogPolyApprox<T>::
OrthogPolyApprox(const Stokhos::OrthogPolyApprox<T>& x) :
  coeff_(x.coeff_) 
{
}

template <typename T> 
Stokhos::OrthogPolyApprox<T>::
~OrthogPolyApprox()
{
}

template <typename T> 
Stokhos::OrthogPolyApprox<T>&
Stokhos::OrthogPolyApprox<T>::
operator=(const Stokhos::OrthogPolyApprox<T>& x) 
{
  if (this != &x) {
    coeff_ = x.coeff_;
  }

  return *this;
}

template <typename T> 
void
Stokhos::OrthogPolyApprox<T>::
resize(unsigned int sz)
{
  coeff_.resize(sz);
}

template <typename T> 
void
Stokhos::OrthogPolyApprox<T>::
reserve(unsigned int sz)
{
  coeff_.reserve(sz);
}

template <typename T> 
template <typename BasisT>
Stokhos::Polynomial<T>
Stokhos::OrthogPolyApprox<T>::
toStandardBasis(const BasisT& basis) const
{
  return basis.toStandardBasis(&coeff_[0], coeff_.size());
}

template <typename T>
std::ostream& 
Stokhos::operator << (std::ostream& os, 
		      const Stokhos::OrthogPolyApprox<T>& a)
{
  os << "[ ";
      
  for (unsigned int i=0; i<a.size(); i++) {
    os << a[i] << " ";
  }

  os << "]\n";
  return os;
}
