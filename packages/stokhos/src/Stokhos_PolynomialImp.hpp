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
Stokhos::Polynomial<T>::
Polynomial(unsigned int deg) :
  coeffs(deg+1, T(0.))
{
}

template <typename T>
Stokhos::Polynomial<T>::
Polynomial(const Teuchos::Array<T>& coefficients) :
  coeffs(coefficients)
{
}

template <typename T>
Stokhos::Polynomial<T>::
Polynomial(const Stokhos::Polynomial<T>& p) :
  coeffs(p.coeffs)
{
}

template <typename T>
Stokhos::Polynomial<T>::
~Polynomial()
{
}

template <typename T>
Stokhos::Polynomial<T>&
Stokhos::Polynomial<T>::
operator=(const Stokhos::Polynomial<T>& p)
{
  if (this != &p)
    coeffs = p.coeffs;
  return *this;
}

template <typename T>
unsigned int
Stokhos::Polynomial<T>::
degree() const
{
  return coeffs.size()-1;
}

template <typename T>
const T&
Stokhos::Polynomial<T>::
coeff(unsigned int i) const
{
  return coeffs[i];
}

template <typename T>
T&
Stokhos::Polynomial<T>::
coeff(unsigned int i)
{
  return coeffs[i];
}

template <typename T>
const T&
Stokhos::Polynomial<T>::
operator[](unsigned int i) const
{
  return coeffs[i];
}

template <typename T>
T&
Stokhos::Polynomial<T>::
operator[](unsigned int i)
{
  return coeffs[i];
}

template <typename T>
void
Stokhos::Polynomial<T>::
multiply(const T& alpha, 
	 const Stokhos::Polynomial<T>& a,
	 const Stokhos::Polynomial<T>& b, 
	 const T& beta)
{
  const unsigned int d = this->degree();
  const unsigned int da = a.degree();
  const unsigned int db = b.degree();

  if (alpha == T(0.0))
    for (unsigned int k=0; k<=d; k++)
      coeffs[k] = beta*coeffs[k];
  else if (da >= d && db >=d)
    for (unsigned int k=0; k<=d; k++) {
      T t = 0.0;
      for (unsigned int i=0; i<=k; i++)
	t += a.coeffs[i]*b.coeffs[k-i];
      coeffs[k] = beta*coeffs[k] + alpha*t;
    }
  else if (da >= d) {
    for (unsigned int k=0; k<=db; k++) {
      T t = 0.0;
      for (unsigned int i=0; i<=k; i++)
	t += a.coeffs[k-i]*b.coeffs[i];
      coeffs[k] = beta*coeffs[k] + alpha*t;
    }
    for (unsigned int k=db+1; k<=d; k++) {
      T t = 0.0;
      for (unsigned int i=k-db; i<=k; i++)
	t += a.coeffs[i]*b.coeffs[k-i];
      coeffs[k] = beta*coeffs[k] + alpha*t;
    }
  }
  else if (db >= d) {
    for (unsigned int k=0; k<=da; k++) {
      T t = 0.0;
      for (unsigned int i=0; i<=k; i++)
	t += a.coeffs[i]*b.coeffs[k-i];
      coeffs[k] = beta*coeffs[k] + alpha*t;
    }
    for (unsigned int k=da+1; k<=d; k++) {
      T t = 0.0;
      for (unsigned int i=k-da; i<=k; i++)
	t += a.coeffs[k-i]*b.coeffs[i];
      coeffs[k] = beta*coeffs[k] + alpha*t;
    }
  }
  else if (da >= db) {
    for (unsigned int k=0; k<=db; k++) {
      T t = 0.0;
      for (unsigned int i=0; i<=k; i++)
	t += a.coeffs[i]*b.coeffs[k-i];
      coeffs[k] = beta*coeffs[k] + alpha*t;
    }
    for (unsigned int k=db+1; k<=da; k++) {
      T t = 0.0;
      for (unsigned int i=k-db; i<=k; i++)
	t += a.coeffs[i]*b.coeffs[k-i];
      coeffs[k] = beta*coeffs[k] + alpha*t;
    }
    for (unsigned int k=da+1; k<=d; k++) {
      T t = 0.0;
      for (unsigned int i=k-db; i<=da; i++)
	t += a.coeffs[i]*b.coeffs[k-i];
      coeffs[k] = beta*coeffs[k] + alpha*t;
    }
  }
  else {
    for (unsigned int k=0; k<=da; k++) {
      T t = 0.0;
      for (unsigned int i=0; i<=k; i++)
	t += a.coeffs[k-i]*b.coeffs[i];
      coeffs[k] = beta*coeffs[k] + alpha*t;
    }
    for (unsigned int k=da+1; k<=db; k++) {
      T t = 0.0;
      for (unsigned int i=k-da; i<=k; i++)
	t += a.coeffs[k-i]*b.coeffs[i];
      coeffs[k] = beta*coeffs[k] + alpha*t;
    }
    for (unsigned int k=db+1; k<=d; k++) {
      T t = 0.0;
      for (unsigned int i=k-da; i<=db; i++)
	t += a.coeffs[k-i]*b.coeffs[i];
      coeffs[k] = beta*coeffs[k] + alpha*t;
    }
  }
}

template <typename T>
void
Stokhos::Polynomial<T>::
add(const T& alpha, 
    const Stokhos::Polynomial<T>& a,
    const T& gamma)
{
  const unsigned int d = this->degree();
  const unsigned int da = a.degree();

  if (da >= d)
    for (unsigned int i=0; i<=d; i++)
      coeffs[i] = alpha*a.coeffs[i] + gamma*coeffs[i];
  else {
    for (unsigned int i=0; i<=da; i++)
      coeffs[i] = alpha*a.coeffs[i] + gamma*coeffs[i];
    for (unsigned int i=da+1; i<=d; i++)
      coeffs[i] = gamma*coeffs[i];
  }
}

template <typename T>
void
Stokhos::Polynomial<T>::
add(const T& alpha, 
    const Stokhos::Polynomial<T>& a,
    const T& beta,
    const Stokhos::Polynomial<T>& b, 
    const T& gamma)
{
  const unsigned int d = this->degree();
  const unsigned int da = a.degree();
  const unsigned int db = b.degree();

  if (da >= d && db >= d)
    for (unsigned int i=0; i<=d; i++)
      coeffs[i] = alpha*a.coeffs[i] + beta*b.coeffs[i] + gamma*coeffs[i];
  else if (da < d) {
    for (unsigned int i=0; i<=da; i++)
      coeffs[i] = alpha*a.coeffs[i] + beta*b.coeffs[i] + gamma*coeffs[i];
    if (db <= d) {
      for (unsigned int i=da+1; i<=db; i++)
	coeffs[i] = beta*b.coeffs[i] + gamma*coeffs[i];
      for (unsigned int i=db+1; i<=d; i++)
	coeffs[i] = gamma*coeffs[i];
    }
    else
      for (unsigned int i=da+1; i<=d; i++)
	coeffs[i] = beta*b.coeffs[i] + gamma*coeffs[i];
  }
  else {
    for (unsigned int i=0; i<=db; i++)
      coeffs[i] = alpha*a.coeffs[i] + beta*b.coeffs[i] + gamma*coeffs[i];
    if (da <= d) {
      for (unsigned int i=db+1; i<=da; i++)
	coeffs[i] = alpha*a.coeffs[i] + gamma*coeffs[i];
      for (unsigned int i=da+1; i<=d; i++)
	coeffs[i] = gamma*coeffs[i];
    }
    else
      for (unsigned int i=db+1; i<=d; i++)
	coeffs[i] = alpha*a.coeffs[i] + gamma*coeffs[i];
  }
}

template <typename T>
void
Stokhos::Polynomial<T>::
print(std::ostream& os) const
{
  os << "[";
  for (unsigned int i=0; i<coeffs.size(); i++)
    os << " " << coeffs[i];
  os << " ]\n";
}
