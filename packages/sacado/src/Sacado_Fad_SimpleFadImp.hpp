// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Sacado_ConfigDefs.h"

template <typename ValueT, typename ScalarT> 
inline Sacado::Fad::SimpleFad<ValueT,ScalarT>& 
Sacado::Fad::SimpleFad<ValueT,ScalarT>::
operator += (const Sacado::Fad::SimpleFad<ValueT,ScalarT>& x)
{
  int xsz = x.size(), sz = this->s_.size();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for (int i=0; i<sz; ++i)
	  this->s_.dx_[i] += x.fastAccessDx(i);
      else
	for (int i=0; i<sz; ++i)
	  this->s_.dx_[i] += x.dx(i);
    }
    else {
      this->s_.resize(xsz);
      if (x.hasFastAccess())
	for (int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = x.dx(i);
    }
  }

  this->s_.val_ += x.val();

  return *this;
}

template <typename ValueT, typename ScalarT> 
inline Sacado::Fad::SimpleFad<ValueT,ScalarT>& 
Sacado::Fad::SimpleFad<ValueT,ScalarT>::
operator -= (const Sacado::Fad::SimpleFad<ValueT,ScalarT>& x)
{
  int xsz = x.size(), sz = this->s_.size();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  this->s_.dx_[i] -= x.fastAccessDx(i);
      else
	for (int i=0; i<sz; ++i)
	  this->s_.dx_[i] -= x.dx(i);
    }
    else {
      this->s_.resize(xsz);
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = -x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = -x.dx(i);
    }
  }

  this->s_.val_ -= x.val();


  return *this;
}

template <typename ValueT, typename ScalarT> 
inline Sacado::Fad::SimpleFad<ValueT,ScalarT>& 
Sacado::Fad::SimpleFad<ValueT,ScalarT>::
operator *= (const Sacado::Fad::SimpleFad<ValueT,ScalarT>& x)
{
  int xsz = x.size(), sz = this->s_.size();
  ValueT xval = x.val();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  this->s_.dx_[i] = this->s_.val_ * x.fastAccessDx(i) + this->s_.dx_[i] * xval;
      else
	for (int i=0; i<sz; ++i)
	  this->s_.dx_[i] = this->s_.val_ * x.dx(i) + this->s_.dx_[i] * xval;
    }
    else {
      this->s_.resize(xsz);
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = this->s_.val_ * x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = this->s_.val_ * x.dx(i);
    }
  }
  else {
    if (sz) {
      for (int i=0; i<sz; ++i)
	this->s_.dx_[i] *= xval;
    }
  }

  this->s_.val_ *= xval;

  return *this;
}

template <typename ValueT, typename ScalarT>
inline Sacado::Fad::SimpleFad<ValueT,ScalarT>& 
Sacado::Fad::SimpleFad<ValueT,ScalarT>::
operator /= (const Sacado::Fad::SimpleFad<ValueT,ScalarT>& x)
{
  int xsz = x.size(), sz = this->s_.size();
  ValueT xval = x.val();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  this->s_.dx_[i] = ( this->s_.dx_[i]*xval - this->s_.val_*x.fastAccessDx(i) )/ (xval*xval);
      else
	for (int i=0; i<sz; ++i)
	  this->s_.dx_[i] = ( this->s_.dx_[i]*xval - this->s_.val_*x.dx(i) )/ (xval*xval);
    }
    else {
      this->s_.resize(xsz);
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = - this->s_.val_*x.fastAccessDx(i) / (xval*xval);
      else
	for (int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = -this->s_.val_ * x.dx(i) / (xval*xval);
    }
  }
  else {
    if (sz) {
      for (int i=0; i<sz; ++i)
	this->s_.dx_[i] /= xval;
    }
  }

  this->s_.val_ /= xval;

  return *this;
}

template <typename ValueT> 
inline Sacado::Fad::SimpleFad<ValueT,ValueT>& 
Sacado::Fad::SimpleFad<ValueT,ValueT>::
operator += (const Sacado::Fad::SimpleFad<ValueT,ValueT>& x)
{
  int xsz = x.size(), sz = this->s_.size();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for (int i=0; i<sz; ++i)
	  this->s_.dx_[i] += x.fastAccessDx(i);
      else
	for (int i=0; i<sz; ++i)
	  this->s_.dx_[i] += x.dx(i);
    }
    else {
      this->s_.resize(xsz);
      if (x.hasFastAccess())
	for (int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = x.dx(i);
    }
  }

  this->s_.val_ += x.val();

  return *this;
}

template <typename ValueT> 
inline Sacado::Fad::SimpleFad<ValueT,ValueT>& 
Sacado::Fad::SimpleFad<ValueT,ValueT>::
operator -= (const Sacado::Fad::SimpleFad<ValueT,ValueT>& x)
{
  int xsz = x.size(), sz = this->s_.size();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  this->s_.dx_[i] -= x.fastAccessDx(i);
      else
	for (int i=0; i<sz; ++i)
	  this->s_.dx_[i] -= x.dx(i);
    }
    else {
      this->s_.resize(xsz);
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = -x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = -x.dx(i);
    }
  }

  this->s_.val_ -= x.val();


  return *this;
}

template <typename ValueT> 
inline Sacado::Fad::SimpleFad<ValueT,ValueT>& 
Sacado::Fad::SimpleFad<ValueT,ValueT>::
operator *= (const Sacado::Fad::SimpleFad<ValueT,ValueT>& x)
{
  int xsz = x.size(), sz = this->s_.size();
  ValueT xval = x.val();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  this->s_.dx_[i] = this->s_.val_ * x.fastAccessDx(i) + this->s_.dx_[i] * xval;
      else
	for (int i=0; i<sz; ++i)
	  this->s_.dx_[i] = this->s_.val_ * x.dx(i) + this->s_.dx_[i] * xval;
    }
    else {
      this->s_.resize(xsz);
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = this->s_.val_ * x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = this->s_.val_ * x.dx(i);
    }
  }
  else {
    if (sz) {
      for (int i=0; i<sz; ++i)
	this->s_.dx_[i] *= xval;
    }
  }

  this->s_.val_ *= xval;

  return *this;
}

template <typename ValueT>
inline Sacado::Fad::SimpleFad<ValueT,ValueT>& 
Sacado::Fad::SimpleFad<ValueT,ValueT>::
operator /= (const Sacado::Fad::SimpleFad<ValueT,ValueT>& x)
{
  int xsz = x.size(), sz = this->s_.size();
  ValueT xval = x.val();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  this->s_.dx_[i] = ( this->s_.dx_[i]*xval - this->s_.val_*x.fastAccessDx(i) )/ (xval*xval);
      else
	for (int i=0; i<sz; ++i)
	  this->s_.dx_[i] = ( this->s_.dx_[i]*xval - this->s_.val_*x.dx(i) )/ (xval*xval);
    }
    else {
      this->s_.resize(xsz);
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = - this->s_.val_*x.fastAccessDx(i) / (xval*xval);
      else
	for (int i=0; i<xsz; ++i)
	  this->s_.dx_[i] = -this->s_.val_ * x.dx(i) / (xval*xval);
    }
  }
  else {
    if (sz) {
      for (int i=0; i<sz; ++i)
	this->s_.dx_[i] /= xval;
    }
  }

  this->s_.val_ /= xval;

  return *this;
}

