// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(Intrepid_MiniTensor_TensorBase_i_h)
#define Intrepid_MiniTensor_TensorBase_i_h

namespace Intrepid
{

//
// Default constructor.
//
template<typename T, typename ST>
inline
TensorBase<T, ST>::TensorBase() :
dimension_(0)
{
  if (ST::IS_DYNAMIC == true) {
    set_number_components(0);
  } else {
    fill(NANS);
  }
  return;
}

//
// Construction that initializes to NaNs
//
template<typename T, typename ST>
inline
TensorBase<T, ST>::TensorBase(Index const dimension, Index const order) :
dimension_(0)
{
  set_dimension(dimension, order);

  fill(NANS);

  return;
}

//
// Create with specified value
//
template<typename T, typename ST>
inline
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    ComponentValue const value) :
    dimension_(0)
{
  set_dimension(dimension, order);

  fill(value);

  return;
}

//
// Construction from a scalar
//
template<typename T, typename ST>
inline
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    T const & s) :
    dimension_(0)
{
  set_dimension(dimension, order);

  fill(s);

  return;
}

//
// Construction from array
//Kokkos data Types:
template<typename T, typename ST>
template<class ArrayT, typename iType>
inline
TensorBase<T, ST>::TensorBase(
   Index const dimension,
   Index const order,
   ArrayT &data,
   iType indx1):
   dimension_(0)
{
 set_dimension(dimension, order);
 fill(data,indx1);
  return;
}



template<typename T, typename ST>
template<class ArrayT, typename iType>
inline
TensorBase<T, ST>::TensorBase(
   Index const dimension,
   Index const order,
   ArrayT &data,
   iType indx1,
   iType indx2):
   dimension_(0)
{
 set_dimension(dimension, order);
 fill(data,indx1,indx2);
  return;
}

template<typename T, typename ST>
template<class ArrayT, typename iType>
inline
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    ArrayT &data,
    iType indx1, 
    iType indx2,
    iType indx3) :
    dimension_(0)
{
  set_dimension(dimension, order);

  fill(data,indx1,indx2,indx3);

  return;
}

template<typename T, typename ST>
template<class ArrayT, typename iType>
inline
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    ArrayT &data,
    iType indx1,
    iType indx2,
    iType indx3,
    iType indx4) :
    dimension_(0)
{
  set_dimension(dimension, order);

  fill(data,indx1,indx2,indx3, indx4);

  return;
}

template<typename T, typename ST>
template<class ArrayT, typename iType>
inline
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    ArrayT &data,
    iType indx1,
    iType indx2,
    iType indx3,
    iType indx4,
    iType indx5) :
    dimension_(0)
{
  set_dimension(dimension, order);

  fill(data,indx1,indx2,indx3, indx4, indx5);

  return;
}

template<typename T, typename ST>
template<class ArrayT, typename iType>
inline
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    ArrayT &data,
    iType indx1,
    iType indx2,
    iType indx3,
    iType indx4,
    iType indx5, 
    iType indx6) :
    dimension_(0)
{
  set_dimension(dimension, order);

  fill(data,indx1,indx2,indx3, indx4, indx5, indx6);

  return;
}



template<typename T, typename ST>
inline
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    T const * data_ptr) :
    dimension_(0)
{
  set_dimension(dimension, order);

  fill(data_ptr);

  return;
}
//
// Copy constructor
//
template<typename T, typename ST>
inline
TensorBase<T, ST>::TensorBase(TensorBase<T, ST> const & X) :
dimension_(X.dimension_)
{
  Index const
  number_components = X.get_number_components();

  set_number_components(number_components);

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] = X[i];
  }

  return;
}

//
// Copy assignment
//
template<typename T, typename ST>
inline
TensorBase<T, ST> &
TensorBase<T, ST>::operator=(TensorBase<T, ST> const & X)
{
  if (this == &X) return *this;

  dimension_ = X.dimension_;

  Index const
  number_components = X.get_number_components();

  set_number_components(number_components);

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] = X[i];
  }

  return *this;
}

//
// Simple destructor
//
template<typename T, typename ST>
inline
TensorBase<T, ST>::~TensorBase()
{
  return;
}

//
// Get dimension
//
template<typename T, typename ST>
inline
Index
TensorBase<T, ST>::get_dimension() const
{
  assert(ST::IS_DYNAMIC == true);

  return dimension_;
}

//
// Set dimension
//
template<typename T, typename ST>
inline
void
TensorBase<T, ST>::set_dimension(Index const dimension, Index const order)
{
  if (ST::IS_STATIC == true) return;

  dimension_ = dimension;

  Index const
  number_components = integer_power(dimension, order);

  set_number_components(number_components);

  return;
}

//
// Linear access to components
//
template<typename T, typename ST>
inline
T const &
TensorBase<T, ST>::operator[](Index const i) const
{
  return components_[i];
}

//
// Linear access to components
//
template<typename T, typename ST>
inline
T &
TensorBase<T, ST>::operator[](Index const i)
{
  return components_[i];
}

//
// Get total number of components
//
template<typename T, typename ST>
inline
Index
TensorBase<T, ST>::get_number_components() const
{
  return components_.size();
}

//
// Allocate space for components
//
template<typename T, typename ST>
inline
void
TensorBase<T, ST>::set_number_components(Index const number_components)
{
  components_.resize(number_components);

  return;
}

//
// Fill components with value.
//
template<typename T, typename ST>
inline
void
TensorBase<T, ST>::fill(ComponentValue const value)
{
  Index const
  number_components = get_number_components();

  switch (value) {

    case ZEROS:
      for (Index i = 0; i < number_components; ++i) {
        (*this)[i] = 0;
      }
      break;

    case ONES:
      for (Index i = 0; i < number_components; ++i) {
        (*this)[i] = 1;
      }
      break;

    case SEQUENCE:
      for (Index i = 0; i < number_components; ++i) {
        (*this)[i] = static_cast<T>(i);
      }
      break;

    case RANDOM_UNIFORM:
      for (Index i = 0; i < number_components; ++i) {
        (*this)[i] = random_uniform<T>();
      }
      break;

    case RANDOM_NORMAL:
      for (Index i = 0; i < number_components; ++i) {
        (*this)[i] = random_normal<T>();
      }
      break;

    case NANS:
      for (Index i = 0; i < number_components; ++i) {
        (*this)[i] = not_a_number<T>();
      }
      break;

    default:
      std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
      std::cerr << std::endl;
      std::cerr << "Unknown specification of value for filling components.";
      std::cerr << std::endl;
      exit(1);
      break;
  }

  return;
}

//
// Fill components from argument
//
template<typename T, typename ST>
inline
void
TensorBase<T, ST>::fill(T const & s)
{
  Index const
  number_components = get_number_components();

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] = s;
  }

  return;
}

//
// Fill components from array defined by pointer.
//#ifdef HAVE_INTREPID_KOKKOSCORE
template<typename T, typename ST>
template<class ArrayT, typename iType>
inline
void
TensorBase<T, ST>::fill(ArrayT  &data, iType indx1)
{

  Index const number_components = get_number_components();
  Index const rank=number_components/data.dimension(0);

 if (indx1==0)
 {
  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] = data(i);
  }
 }
 else
    TEUCHOS_TEST_FOR_EXCEPTION( !( (indx1 == 0)  ), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): las input parameter should always be 0");
  return;
}

template<typename T, typename ST>
template<class ArrayT, typename iType>
inline
void
TensorBase<T, ST>::fill(ArrayT  &data, iType indx1,iType indx2)
{

  Index const  number_components = get_number_components();

  Index rank=0;
  Index temp=number_components;
  Index const dim=data.dimension(1);

  while (temp!=1){
   temp =temp/dim;
   rank=rank +1;
   if (temp<1) TEUCHOS_TEST_FOR_EXCEPTION( ( (temp<1)  ), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): rank calculation is not correct");
  }

   TEUCHOS_TEST_FOR_EXCEPTION( !( (indx2 == 0)  ), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): las input parameter should always be 0"); 

 
if (rank==1){
  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] = data(indx1,i);
  }
}
else if (rank==2){
  for (Index i = 0; i <dim; ++i)
     for (Index j = 0; j <dim; ++j)
      (*this)[dim*i+j] =data(i,j);
}
else  TEUCHOS_TEST_FOR_EXCEPTION( !( rank==1 || rank==2), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): method doesn't comput correct rank");

  return;
}

template<typename T, typename ST>
template<class ArrayT, typename iType>
inline
void
TensorBase<T, ST>::fill(ArrayT  &data, iType indx1, iType indx2, iType indx3)
{

  Index const number_components = get_number_components();
  Index const dim = data.dimension(2);
  
  Index rank=0;
  Index temp=number_components;

  while (temp!=1){
   temp =temp/dim;
   rank=rank +1;
   if (temp<1) TEUCHOS_TEST_FOR_EXCEPTION( ( (temp<1)  ), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): rank calculation is not correct");
  } 

 TEUCHOS_TEST_FOR_EXCEPTION( !( (indx3 == 0)  ), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): las input parameter should always be 0");
  if (rank == 1){
   for (Index i = 0; i < dim; ++i) {
          (*this)[i] = data(indx1,indx2,i);
   }
  }
  else if (rank==2)
  {
    for (Index i2 = 0; i2 < dim; ++i2) { 
     for (Index i3 = 0; i3 < dim; ++i3)      
      (*this)[i2*dim+i3] = data(indx1,i2,i3);
     }
  }
  else if (rank==3)
 {
  for (Index i1 = 0; i1 < dim; ++i1) {
   for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
      (*this)[i1*dim*dim+i2*dim+i3] = data(i1,i2,i3);
     }
   }
 } 
  else  TEUCHOS_TEST_FOR_EXCEPTION( !( rank==1)&&!( rank==2)&&!( rank==3), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): method doesn't comput correct rank");

  return;
}

template<typename T, typename ST>
template<class ArrayT, typename iType>
inline
void
TensorBase<T, ST>::fill(ArrayT  &data, iType  indx1,  iType  indx2, iType indx3, iType indx4)
{

  Index const
  number_components = get_number_components();
  Index const dim = data.dimension(3);

  Index rank=0;
  Index temp=number_components;

  while (temp!=1){
   temp =temp/dim;
   rank=rank +1;
   if (temp<1) TEUCHOS_TEST_FOR_EXCEPTION( ( (temp<1)  ), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): rank calculation is not correct");
  }


  TEUCHOS_TEST_FOR_EXCEPTION( !( (indx4 == 0)  ), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): las input parameter should always be 0");
   if (rank == 1){
   for (Index i = 0; i < dim; ++i) {
          (*this)[i] = data(indx1,indx2,indx3,i);
   }
  }
  else if (rank==2)
  {
    for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
      (*this)[i2*dim+i3] = data(indx1,indx2,i2,i3);
     }
  }
  else if (rank==3)
 {
  for (Index i1 = 0; i1 < dim; ++i1) {
   for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
      (*this)[i1*dim*dim+i2*dim+i3] = data(indx1,i1,i2,i3);
     }
   }
 }
 else if (rank==4)
 {
  for (Index i1 = 0; i1 < dim; ++i1) {
   for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
       for (Index i4 = 0; i4 < dim; ++i4)
         (*this)[i1*dim*dim*dim+i2*dim*dim+i3*dim+i4] = data(i1,i2,i3,i4);
     }
   }
 }
  else  TEUCHOS_TEST_FOR_EXCEPTION( !( rank==1)&&!( rank==2)&&!( rank==3)&&!( rank==4), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): method doesn't comput correct rank");

  return;
}

template<typename T, typename ST>
template<class ArrayT, typename iType>
inline
void
TensorBase<T, ST>::fill(ArrayT  &data, iType indx1, iType indx2, iType indx3, iType indx4, iType indx5)
{

  Index const
  number_components = get_number_components();

  TEUCHOS_TEST_FOR_EXCEPTION( !( (indx5 == 0)  ), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): las input parameter should always be 0");

  Index rank=0;
  Index temp=number_components;
  Index const dim = data.dimension(4);

  while (temp!=1){
   temp =temp/dim;
   rank=rank +1;
   if (temp<1) TEUCHOS_TEST_FOR_EXCEPTION( ( (temp<1)  ), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): rank calculation is not correct");
  }

   if (rank == 1){
   for (Index i = 0; i < dim; ++i) {
          (*this)[i] = data(indx1,indx2,indx3,indx4,i);
   }
  }
  else if (rank==2)
  {
    for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
      (*this)[i2*dim+i3] = data(indx1,indx2,indx3,i2,i3);
     }
  }
  else if (rank==3)
 {
  for (Index i1 = 0; i1 < dim; ++i1) {
   for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
      (*this)[i1*dim*dim+i2*dim+i3] = data(indx1,indx2,i1,i2,i3);
     }
   }
 }
 else if (rank==4)
 {
  for (Index i1 = 0; i1 < dim; ++i1) {
   for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
       for (Index i4 = 0; i4 < dim; ++i4)
         (*this)[i1*dim*dim*dim+i2*dim*dim+i3*dim+i4] = data(indx1,i1,i2,i3,i4);
     }
   }
 }
 else if (rank==5)
 {
  for (Index i1 = 0; i1 < dim; ++i1) {
   for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
       for (Index i4 = 0; i4 < dim; ++i4)
         for (Index i5 = 0; i5 < dim; ++i5)
           (*this)[i1*dim*dim*dim*dim+i2*dim*dim*dim+i3*dim*dim+i4*dim+i5] = data(i1,i2,i3,i4,i5);
     }
   }
 }
  else  TEUCHOS_TEST_FOR_EXCEPTION( !( rank==1)&&!( rank==2)&&!( rank==3)&&!( rank==4)&&!( rank==5), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): method doesn't comput correct rank");



  return;
}

template<typename T, typename ST>
template<class ArrayT, typename iType>
inline
void
TensorBase<T, ST>::fill(ArrayT  &data, iType indx1, iType indx2, iType indx3, iType indx4, iType indx5, iType indx6)
{

  Index const
  number_components = get_number_components();
  Index const dim = data.dimension(5);

  TEUCHOS_TEST_FOR_EXCEPTION( !( (indx6 == 0)  ), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): las input parameter should always be 0");

  Index rank=0;
  Index temp=number_components;

  while (temp!=1){
   temp =temp/dim;
   rank=rank +1;
   if (temp<1) TEUCHOS_TEST_FOR_EXCEPTION( ( (temp<1)  ), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): rank calculation is not correct");
  }

  if (rank == 1){
   for (Index i = 0; i < dim; ++i) {
          (*this)[i] = data(indx1,indx2,indx3,indx4,indx5,i);
   }
  }
  else if (rank==2)
  {
    for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
      (*this)[i2*dim+i3] = data(indx1,indx2,indx3,indx4,i2,i3);
     }
  }
  else if (rank==3)
 {
  for (Index i1 = 0; i1 < dim; ++i1) {
   for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
      (*this)[i1*dim*dim+i2*dim+i3] = data(indx1,indx2,indx3,i1,i2,i3);
     }
   }
 }
 else if (rank==4)
 {
  for (Index i1 = 0; i1 < dim; ++i1) {
   for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
       for (Index i4 = 0; i4 < dim; ++i4)
         (*this)[i1*dim*dim*dim+i2*dim*dim+i3*dim+i4] = data(indx1,indx2,i1,i2,i3,i4);
     }
   }
 }
 else if (rank==5)
 {
  for (Index i1 = 0; i1 < dim; ++i1) {
   for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
       for (Index i4 = 0; i4 < dim; ++i4)
         for (Index i5 = 0; i5 < dim; ++i5)
           (*this)[i1*dim*dim*dim*dim+i2*dim*dim*dim+i3*dim*dim+i4*dim+i5] = data(indx1,i1,i2,i3,i4,i5);
     }
   }
 }
 else if (rank==6)
 {
  for (Index i1 = 0; i1 < dim; ++i1) {
   for (Index i2 = 0; i2 < dim; ++i2) {
     for (Index i3 = 0; i3 < dim; ++i3)
       for (Index i4 = 0; i4 < dim; ++i4)
         for (Index i5 = 0; i5 < dim; ++i5)
           for (Index i6 = 0; i6 < dim; ++i6)
           (*this)[i1*dim*dim*dim*dim*dim+i2*dim*dim*dim*dim+i3*dim*dim*dim+i4*dim*dim+i5*dim+i6] = data(i1,i2,i3,i4,i5,i6);
     }
   }
 }
  else  TEUCHOS_TEST_FOR_EXCEPTION( !( rank==1)&&!( rank==2)&&!( rank==3)&&!( rank==4)&&!( rank==5)&&!( rank==6), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): method doesn't comput correct rank");

 return;
}
template<typename T, typename ST>
inline
void
TensorBase<T, ST>::fill(T const * data_ptr)
{
  assert(data_ptr != NULL);

  Index const
  number_components = get_number_components();

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] = data_ptr[i];
  }

  return;
}


//
// Component increment
//
template<typename T, typename ST>
template<typename S, typename SS>
inline
TensorBase<T, ST> &
TensorBase<T, ST>::operator+=(TensorBase<S, SS> const & X)
{
  Index const
  number_components = get_number_components();

  assert(number_components == X.get_number_components());

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] += X[i];
  }

  return *this;
}

//
// Component decrement
//
template<typename T, typename ST>
template<typename S, typename SS>
inline
TensorBase<T, ST> &
TensorBase<T, ST>::operator-=(TensorBase<S, SS> const & X)
{
  Index const
  number_components = get_number_components();

  assert(number_components == X.get_number_components());

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] -= X[i];
  }

  return *this;
}

//
// Component scale
//
template<typename T, typename ST>
template<typename S>
inline
TensorBase<T, ST> &
TensorBase<T, ST>::operator*=(S const & X)
{
  Index const
  number_components = get_number_components();

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] *= X;
  }
  return *this;
}

//
// Component divide
//
template<typename T, typename ST>
template<typename S>
inline
TensorBase<T, ST> &
TensorBase<T, ST>::operator/=(S const & X)
{
  Index const
  number_components = get_number_components();

  for (Index i = 0; i < number_components; ++i) {
    (*this)[i] /= X;
  }
  return *this;
}

//
// Fill with zeros
//
template<typename T, typename ST>
inline
void
TensorBase<T, ST>::clear()
{
  fill(ZEROS);
  return;
}

//
// Square of Frobenius norm
//
template<typename T, typename ST>
T
norm_f_square(TensorBase<T, ST> const & X)
{
  T
  s = 0.0;

  for (Index i = 0; i < X.get_number_components(); ++i) {
    s += X[i] * X[i];
  }

  return s;
}

//
// Frobenius norm
//
template<typename T, typename ST>
T
norm_f(TensorBase<T, ST> const & X)
{
  return std::sqrt(norm_f_square(X));
}

//
// Base addition
//
template<typename R, typename S, typename T, typename SR, typename SS,
typename ST>
void
add(
    TensorBase<R, SR> const & A,
    TensorBase<S, SS> const & B,
    TensorBase<T, ST> & C
)
{
  Index const
  number_components = A.get_number_components();

  assert(B.get_number_components() == number_components);
  assert(C.get_number_components() == number_components);

  for (Index i = 0; i < number_components; ++i) {
    C[i] = A[i] + B[i];
  }

  return;
}

//
// Base subtraction
//
template<typename R, typename S, typename T, typename SR, typename SS,
typename ST>
void
subtract(
    TensorBase<R, SR> const & A,
    TensorBase<S, SS> const & B,
    TensorBase<T, ST> & C)
{
  Index const
  number_components = A.get_number_components();

  assert(B.get_number_components() == number_components);
  assert(C.get_number_components() == number_components);

  for (Index i = 0; i < number_components; ++i) {
    C[i] = A[i] - B[i];
  }

  return;
}

//
// Base minus
//
template<typename T, typename ST>
void
minus(TensorBase<T, ST> const & A, TensorBase<T, ST> & B)
{
  Index const
  number_components = A.get_number_components();

  assert(B.get_number_components() == number_components);

  for (Index i = 0; i < number_components; ++i) {
    B[i] = - A[i];
  }

  return;
}

//
// Base equality
//
template<typename T, typename ST>
bool
equal(TensorBase<T, ST> const & A, TensorBase<T, ST> const & B)
{
  Index const
  number_components = A.get_number_components();

  assert(B.get_number_components() == number_components);

  for (Index i = 0; i < number_components; ++i) {
    if (A[i] != B[i]) return false;
  }

  return true;
}

//
// Base not equality
//
template<typename T, typename ST>
bool
not_equal(TensorBase<T, ST> const & A, TensorBase<T, ST> const & B)
{
  return !(equal(A, B));
}

//
// Base scaling
//
template<typename R, typename S, typename T, typename SR, typename ST>
void
scale(TensorBase<R, SR> const & A, S const & s, TensorBase<T, ST> & B)
{
  Index const
  number_components = A.get_number_components();

  assert(B.get_number_components() == number_components);

  for (Index i = 0; i < number_components; ++i) {
    B[i] = s * A[i];
  }

  return;
}

//
// Base division
//
template<typename R, typename S, typename T, typename SR, typename ST>
void
divide(TensorBase<R, SR> const & A, S const & s, TensorBase<T, ST> & B)
{
  Index const
  number_components = A.get_number_components();

  assert(B.get_number_components() == number_components);

  for (Index i = 0; i < number_components; ++i) {
    B[i] = A[i] / s;
  }

  return;
}

//
// Base split (scalar divided by tensor)
//
template<typename R, typename S, typename T, typename SR, typename ST>
void
split(TensorBase<R, SR> const & A, S const & s, TensorBase<T, ST> & B)
{
  Index const
  number_components = A.get_number_components();

  assert(B.get_number_components() == number_components);

  for (Index i = 0; i < number_components; ++i) {
    B[i] = s / A[i];
  }

  return;
}

} // namespace Intrepid

#endif // Intrepid_MiniTensor_TensorBase_i_h
