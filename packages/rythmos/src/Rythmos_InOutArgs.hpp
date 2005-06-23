//
// @HEADER
// ***********************************************************************
// 
//                           Rythmos Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef Rythmos_INOUT_ARGS_H
#define Rythmos_INOUT_ARGS_H

#include "Thyra_VectorBase.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Rythmos {

template<class Scalar> 
class InArgs
{
  public:
    
    // Destructor
    ~InArgs();

    // Constructor
    InArgs(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &x, Scalar t);
    InArgs();

    // Specify solution vector
    void set_x(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &x);
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > get_x() const;

    // Specify time
    void set_t(Scalar t);
    Scalar get_t() const;

  protected:
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x_;
    Scalar t_;

};

template<class Scalar> 
InArgs<Scalar>::~InArgs()
{
}

template<class Scalar> 
InArgs<Scalar>::InArgs(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &x, Scalar t)
{ 
  x_ = x; 
  t_ = t; 
}

template<class Scalar> 
InArgs<Scalar>::InArgs()
{
}

template<class Scalar> 
void InArgs<Scalar>::set_x(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &x)
{ 
  x_ = x; 
}

template<class Scalar> 
Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > InArgs<Scalar>::get_x() const
{ 
  return(x_); 
}

template<class Scalar> 
void InArgs<Scalar>::set_t(Scalar t)
{ 
  t_ = t; 
}

template<class Scalar> 
Scalar InArgs<Scalar>::get_t() const
{ 
  return(t_); 
}

template<class Scalar> 
class OutArgs
{
  public:
    
    // Destructor
    ~OutArgs();

    // Constructor
    OutArgs();

    // Request residual:
    void request_F(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &F);
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > get_F() const;

  protected:
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > F_;
    
};

template<class Scalar> 
OutArgs<Scalar>::~OutArgs()
{
}

template<class Scalar> 
OutArgs<Scalar>::OutArgs()
{
}

template<class Scalar> 
void OutArgs<Scalar>::request_F(const Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &F) 
{ 
  F_ = F; 
}

template<class Scalar> 
Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > OutArgs<Scalar>::get_F() const
{ 
  return(F_); 
}

}// namespace Rythmos 




#endif // Rythmos_INOUT_ARGS_H
