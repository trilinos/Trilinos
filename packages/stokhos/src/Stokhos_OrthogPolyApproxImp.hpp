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
T&
Stokhos::OrthogPolyApprox<T>::
term(const BasisT& basis, 
     int i0, int i1, int i2, int i3, int i4,
     int i5, int i6, int i7, int i8, int i9)
{
  std::vector<unsigned int> trm;
  unsigned int d = basis.dimension();
  if (i0 >= 0 && d >= 1)
    trm.push_back(i0);
  if (i1 >= 0 && d >= 2)
    trm.push_back(i1);
  if (i2 >= 0 && d >= 3)
    trm.push_back(i2);
  if (i3 >= 0 && d >= 4)
    trm.push_back(i3);
  if (i4 >= 0 && d >= 5)
    trm.push_back(i4);
  if (i5 >= 0 && d >= 6)
    trm.push_back(i5);
  if (i6 >= 0 && d >= 7)
    trm.push_back(i6);
  if (i7 >= 0 && d >= 8)
    trm.push_back(i7);
  if (i8 >= 0 && d >= 9)
    trm.push_back(i8);
  if (i9 >= 0 && d >= 10)
    trm.push_back(i9);

  unsigned int index = basis.getIndex(trm);
  return coeff_[index];
}

template <typename T> 
template <typename BasisT>
const T&
Stokhos::OrthogPolyApprox<T>::
term(const BasisT& basis, 
     int i0, int i1, int i2, int i3, int i4,
     int i5, int i6, int i7, int i8, int i9) const
{
  std::vector<unsigned int> trm;
  unsigned int d = basis.dimension();
  if (i0 >= 0 && d >= 1)
    trm.push_back(i0);
  if (i1 >= 0 && d >= 2)
    trm.push_back(i1);
  if (i2 >= 0 && d >= 3)
    trm.push_back(i2);
  if (i3 >= 0 && d >= 4)
    trm.push_back(i3);
  if (i4 >= 0 && d >= 5)
    trm.push_back(i4);
  if (i5 >= 0 && d >= 6)
    trm.push_back(i5);
  if (i6 >= 0 && d >= 7)
    trm.push_back(i6);
  if (i7 >= 0 && d >= 8)
    trm.push_back(i7);
  if (i8 >= 0 && d >= 9)
    trm.push_back(i8);
  if (i9 >= 0 && d >= 10)
    trm.push_back(i9);

  unsigned int index = basis.getIndex(trm);
  return coeff_[index];
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
template <typename BasisT>
T
Stokhos::OrthogPolyApprox<T>::
evaluate(const BasisT& basis, const std::vector<T>& point) const
{
  std::vector<T> basis_pts(coeff_.size());
  basis.evaluateBases(point, basis_pts);
  T val = T(0.0);
  for (unsigned int i=0; i<coeff_.size(); i++)
    val += coeff_[i]*basis_pts[i];

  return val;
}

template <typename T> 
template <typename BasisT>
std::ostream&
Stokhos::OrthogPolyApprox<T>::
print(const BasisT& basis, std::ostream& os) const
{
  std::vector<unsigned int> trm;
  os << "Stokhos::OrthogPolyApprox of size " << coeff_.size() << " in basis "
     << "\n\t" << basis.getName() << ":" << std::endl;
  for (unsigned int i=0; i<coeff_.size(); i++) {
    trm = basis.getTerm(i);
    os << "\t\t(";
    for (unsigned int j=0; j<trm.size()-1; j++)
      os << trm[j] << ", ";
    os << trm[trm.size()-1] << ") = " << coeff_[i] << std::endl;
  }

  return os;
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
