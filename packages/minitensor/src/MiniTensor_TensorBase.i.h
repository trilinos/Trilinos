// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_TensorBase_i_h)
#define MiniTensor_TensorBase_i_h

namespace minitensor
{

//
// Default constructor.
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
TensorBase<T, ST>::TensorBase()
{
  Index const
  static_size = ST::static_size();

  set_number_components(static_size);
  fill(Filler::NANS);

  return;
}

//
// Construction that initializes to NaNs
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
TensorBase<T, ST>::TensorBase(Index const dimension, Index const order)
{
  set_dimension(dimension, order);
  fill(Filler::NANS);
  return;
}

//
// Create with specified value
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    Filler const value)
{
  set_dimension(dimension, order);
  fill(value);
  return;
}

//
// Construction from a scalar
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    T const & s)
{
  set_dimension(dimension, order);
  fill(s);
  return;
}

//
// Construction from array
//Kokkos data Types:
template<typename T, typename ST>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    ArrayT & data,
    Index index1)
{
  set_dimension(dimension, order);
  fill(data, index1);
  return;
}

template<typename T, typename ST>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    ArrayT & data,
    Index index1,
    Index index2)
{
  set_dimension(dimension, order);
  fill(data, index1, index2);
  return;
}

template<typename T, typename ST>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3)
{
  set_dimension(dimension, order);
  fill(data, index1, index2, index3);
  return;
}

template<typename T, typename ST>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3,
    Index index4)
{
  set_dimension(dimension, order);
  fill(data, index1, index2, index3, index4);
  return;
}

template<typename T, typename ST>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3,
    Index index4,
    Index index5)
{
  set_dimension(dimension, order);
  fill(data, index1, index2, index3, index4, index5);
  return;
}

template<typename T, typename ST>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3,
    Index index4,
    Index index5,
    Index index6)
{
  set_dimension(dimension, order);
  fill(data, index1, index2, index3, index4, index5, index6);
  return;
}

template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
TensorBase<T, ST>::TensorBase(
    Index const dimension,
    Index const order,
    T const * data_ptr)
{
  set_dimension(dimension, order);
  fill(data_ptr);
  return;
}

//
// Copy constructor
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
// Get dimension
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
Index
TensorBase<T, ST>::get_dimension(Index const order) const
{
  return dimension_;
}

//
// Set dimension
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
void
TensorBase<T, ST>::set_dimension(Index const dimension, Index const order)
{
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
KOKKOS_INLINE_FUNCTION
T const &
TensorBase<T, ST>::operator[](Index const i) const
{
  return components_[i];
}

//
// Linear access to components
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
T &
TensorBase<T, ST>::operator[](Index const i)
{
  return components_[i];
}

//
// Get total number of components
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
Index
TensorBase<T, ST>::get_number_components() const
{
  return components_.size();
}

//
// Allocate space for components
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
void
TensorBase<T, ST>::set_number_components(Index const number_components)
{
  using S = typename Sacado::ScalarType<T>::type;

  Index const
  old_size = get_number_components();

  Index const
  new_size = number_components;

  if (new_size < old_size) {
    for (auto i = new_size; i < old_size; ++i) {
      auto & entry = (*this)[i];
      fill_AD<T>(entry, not_a_number<S>());
      entry = not_a_number<T>();
    }
  }

  components_.resize(number_components);

  if (new_size > old_size) {
    for (auto i = old_size; i < new_size; ++i) {
      auto & entry = (*this)[i];
      fill_AD<T>(entry, not_a_number<S>());
      entry = not_a_number<T>();
    }
  }

  return;
}

//
// Fill components with value.
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
void
TensorBase<T, ST>::fill(Filler const value)
{
  using S = typename Sacado::ScalarType<T>::type;

  Index const
  number_components = get_number_components();

  switch (value) {

  case Filler::ZEROS:
    for (Index i = 0; i < number_components; ++i) {
      auto & entry = (*this)[i];
      fill_AD<T>(entry, Kokkos::ArithTraits<S>::zero());
      entry = Kokkos::ArithTraits<S>::zero();
    }
    break;

  case Filler::ONES:
    for (Index i = 0; i < number_components; ++i) {
      auto & entry = (*this)[i];
      fill_AD<T>(entry, Kokkos::ArithTraits<S>::zero());
      entry = Kokkos::ArithTraits<S>::one();
    }
    break;

  case Filler::SEQUENCE:
    for (Index i = 0; i < number_components; ++i) {
      auto & entry = (*this)[i];
      fill_AD<T>(entry, Kokkos::ArithTraits<S>::zero());
      entry = static_cast<S>(i);
    }
    break;

  case Filler::NANS:
    for (Index i = 0; i < number_components; ++i) {
      auto & entry = (*this)[i];
      fill_AD<T>(entry, not_a_number<S>());
      entry = not_a_number<S>();
    }
    break;

  case Filler::RANDOM:
    KOKKOS_IF_ON_HOST((
    for (Index i = 0; i < number_components; ++i) {
      auto & entry = (*this)[i];
      fill_AD<T>(entry, Kokkos::ArithTraits<S>::zero());
      entry = random<S>();
    }
    break;
    ))
    KOKKOS_IF_ON_DEVICE((
    [[fallthrough]];
    ))

  case Filler::RANDOM_UNIFORM:
    KOKKOS_IF_ON_HOST((
    for (Index i = 0; i < number_components; ++i) {
      auto & entry = (*this)[i];
      fill_AD<T>(entry, Kokkos::ArithTraits<S>::zero());
      entry = random_uniform<S>();
    }
    break;
    ))
    KOKKOS_IF_ON_DEVICE((
    [[fallthrough]];
    ))

  case Filler::RANDOM_NORMAL:
    KOKKOS_IF_ON_HOST((
    for (Index i = 0; i < number_components; ++i) {
      auto & entry = (*this)[i];
      fill_AD<T>(entry, Kokkos::ArithTraits<S>::zero());
      entry = random_normal<S>();
    }
    break;
    ))
    KOKKOS_IF_ON_DEVICE((
    [[fallthrough]];
    ))

  default:
    MT_ERROR_EXIT("Unknown or undefined (in execution space) specification of "
                  "value for filling components.");
    break;
  }

  return;
}

//
// Fill components from argument
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
void
TensorBase<T, ST>::fill(T const & s)
{
  using S = typename Sacado::ScalarType<T>::type;

  Index const
  number_components = get_number_components();

  for (Index i = 0; i < number_components; ++i) {
    auto & entry = (*this)[i];
    fill_AD<T>(entry, Kokkos::ArithTraits<S>::zero());
    entry = s;
  }

  return;
}

//
// Fill components from array defined by pointer.
//
template<typename T, typename ST>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
void
TensorBase<T, ST>::fill(
    ArrayT & data,
    Index index1)
{
  assert(index1 == 0);

  Index const
  number_components = get_number_components();

  Index const
  rank = number_components / data.extent(0);

  switch (rank) {

  default:
    MT_ERROR_EXIT("Invalid rank.");
    break;

  case 1:
    for (Index i = 0; i < number_components; ++i) {
      (*this)[i] = data(i);
    }
    break;
  }

  return;
}

template<typename T, typename ST>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
void
TensorBase<T, ST>::fill(
    ArrayT & data,
    Index index1,
    Index index2)
{
  assert(index2 == 0);

  Index const
  number_components = get_number_components();

  Index
  rank = 0;

  Index
  sub_dimension = number_components;

  Index const
  dim = data.extent(1);

  while (sub_dimension != 1) {

    sub_dimension /= dim;
    ++rank;

    assert(sub_dimension >= 1);

  }

  switch (rank) {

  default:
    MT_ERROR_EXIT("Invalid rank.");
    break;

  case 1:
    for (Index j = 0; j < number_components; ++j) {
      (*this)[j] = data(index1, j);
    }
    break;

  case 2:
    for (Index i = 0; i < dim; ++i) {
      for (Index j = 0; j < dim; ++j) {
        (*this)[dim * i + j] = data(i, j);
      }
    }
    break;
  }

  return;
}

template<typename T, typename ST>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
void
TensorBase<T, ST>::fill(
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3)
{
  assert(index3 == 0);

  Index const
  number_components = get_number_components();

  Index const
  dim = data.extent(2);

  Index
  rank = 0;

  Index
  sub_dimension = number_components;

  while (sub_dimension != 1) {

    sub_dimension /= dim;
    ++rank;

    assert(sub_dimension >= 1);

  }

  switch (rank) {

  default:
    MT_ERROR_EXIT("Invalid rank.");
    break;

  case 1:
    for (Index k = 0; k < number_components; ++k) {
      (*this)[k] = data(index1, index2, k);
    }
    break;

  case 2:
    for (Index j = 0; j < dim; ++j) {
      for (Index k = 0; k < dim; ++k) {
        (*this)[dim * j + k] = data(index1, j, k);
      }
    }
    break;

  case 3:
    for (Index i = 0; i < dim; ++i) {
      for (Index j = 0; j < dim; ++j) {
        for (Index k = 0; k < dim; ++k) {
          (*this)[dim * (dim * i + j) + k] = data(i, j, k);
        }
      }
    }
    break;
  }

  return;
}

template<typename T, typename ST>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
void
TensorBase<T, ST>::fill(
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3,
    Index index4)
{
  assert(index4 == 0);

  Index const
  number_components = get_number_components();

  Index const
  dim = data.extent(2);

  Index
  rank = 0;

  Index
  sub_dimension = number_components;

  while (sub_dimension != 1) {

    sub_dimension /= dim;
    ++rank;

    assert(sub_dimension >= 1);

  }

  switch (rank) {

  default:
    MT_ERROR_EXIT("Invalid rank.");
    break;

  case 1:
    for (Index l = 0; l < number_components; ++l) {
      (*this)[l] = data(index1, index2, index3, l);
    }
    break;

  case 2:
    for (Index k = 0; k < dim; ++k) {
      for (Index l = 0; l < dim; ++l) {
        (*this)[dim * k + l] = data(index1, index2, k, l);
      }
    }
    break;

  case 3:
    for (Index j = 0; j < dim; ++j) {
      for (Index k = 0; k < dim; ++k) {
        for (Index l = 0; l < dim; ++l) {
          (*this)[dim * (dim * j + k) + l] = data(index1, j, k, l);
        }
      }
    }
    break;

  case 4:
    for (Index i = 0; i < dim; ++i) {
      for (Index j = 0; j < dim; ++j) {
        for (Index k = 0; k < dim; ++k) {
          for (Index l = 0; l < dim; ++l) {
            (*this)[dim * (dim * (dim * i + j) + k) + l] =
                data(i, j, k, l);
          }
        }
      }
    }
    break;
  }

  return;
}

template<typename T, typename ST>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
void
TensorBase<T, ST>::fill(
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3,
    Index index4,
    Index index5)
{
  assert(index5 == 0);

  Index const
  number_components = get_number_components();

  Index const
  dim = data.extent(2);

  Index
  rank = 0;

  Index
  sub_dimension = number_components;

  while (sub_dimension != 1) {

    sub_dimension /= dim;
    ++rank;

    assert(sub_dimension >= 1);

  }

  switch (rank) {

  default:
    MT_ERROR_EXIT("Invalid rank.");
    break;

  case 1:
    for (Index m = 0; m < number_components; ++m) {
      (*this)[m] = data(index1, index2, index3, index4, m);
    }
    break;

  case 2:
    for (Index l = 0; l < dim; ++l) {
      for (Index m = 0; m < dim; ++m) {
        (*this)[dim * l + m] = data(index1, index2, index3, l, m);
      }
    }
    break;

  case 3:
    for (Index k = 0; k < dim; ++k) {
      for (Index l = 0; l < dim; ++l) {
        for (Index m = 0; m < dim; ++m) {
          (*this)[dim * (dim * k + l) + m] = data(index1, index2, k, l, m);
        }
      }
    }
    break;

  case 4:
    for (Index j = 0; j < dim; ++j) {
      for (Index k = 0; k < dim; ++k) {
        for (Index l = 0; l < dim; ++l) {
          for (Index m = 0; m < dim; ++m) {
            (*this)[dim * (dim * (dim * j + k) + l) + m] =
                data(index1, j, k, l, m);
          }
        }
      }
    }
    break;

  case 5:
    for (Index i = 0; i < dim; ++i) {
      for (Index j = 0; j < dim; ++j) {
        for (Index k = 0; k < dim; ++k) {
          for (Index l = 0; l < dim; ++l) {
            for (Index m = 0; m < dim; ++m) {
              (*this)[dim * (dim * (dim * (dim * i + j) + k) + l) + m] =
                  data(i, j, k, l, m);
            }
          }
        }
      }
    }
    break;
  }

  return;
}

template<typename T, typename ST>
template<class ArrayT>
KOKKOS_INLINE_FUNCTION
void
TensorBase<T, ST>::fill(
    ArrayT & data,
    Index index1,
    Index index2,
    Index index3,
    Index index4,
    Index index5,
    Index index6)
{
  assert(index6 == 0);

  Index const
  number_components = get_number_components();

  Index const
  dim = data.extent(2);

  Index
  rank = 0;

  Index
  sub_dimension = number_components;

  while (sub_dimension != 1) {

    sub_dimension /= dim;
    ++rank;

    assert(sub_dimension >= 1);

  }

  switch (rank) {

  default:
    MT_ERROR_EXIT("Invalid rank.");
    break;

  case 1:
    for (Index n = 0; n < number_components; ++n) {
      (*this)[n] = data(index1, index2, index3, index4, index5, n);
    }
    break;

  case 2:
    for (Index m = 0; m < dim; ++m) {
      for (Index n = 0; n < dim; ++n) {
        (*this)[dim * m + n] = data(index1, index2, index3, index4, m, n);
      }
    }
    break;

  case 3:
    for (Index l = 0; l < dim; ++l) {
      for (Index m = 0; m < dim; ++m) {
        for (Index n = 0; n < dim; ++n) {
          (*this)[dim * (dim * l + m) + n] =
              data(index1, index2, index3, l, m, n);
        }
      }
    }
    break;

  case 4:
    for (Index k = 0; k < dim; ++k) {
      for (Index l = 0; l < dim; ++l) {
        for (Index m = 0; m < dim; ++m) {
          for (Index n = 0; n < dim; ++n) {
            (*this)[dim * (dim * (dim * k + l) + m) + n] =
                data(index1, index2, k, l, m, n);
          }
        }
      }
    }
    break;

  case 5:
    for (Index j = 0; j < dim; ++j) {
      for (Index k = 0; k < dim; ++k) {
        for (Index l = 0; l < dim; ++l) {
          for (Index m = 0; m < dim; ++m) {
            for (Index n = 0; n < dim; ++n) {
              (*this)[dim * (dim * (dim * (dim * j + k) + l) + m) + n] =
                  data(index1, j, k, l, m, n);
            }
          }
        }
      }
    }
    break;

  case 6:
    for (Index i = 0; i < dim; ++i) {
      for (Index j = 0; j < dim; ++j) {
        for (Index k = 0; k < dim; ++k) {
          for (Index l = 0; l < dim; ++l) {
            for (Index m = 0; m < dim; ++m) {
              for (Index n = 0; n < dim; ++n) {
                (*this)[dim * (dim * (dim * (dim * (dim *
                    i + j) + k) + l) + m) + n] = data(i, j, k, l, m, n);
              }
            }
          }
        }
      }
    }
    break;
  }

  return;
}
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
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
// Fill components from array defined by pointer.
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
void
TensorBase<T, ST>::fill(
    T const * data_ptr,
    ComponentOrder const component_order)
{
  assert(data_ptr != NULL);

  TensorBase<T, ST> &
  self = (*this);

  Index const
  number_components = self.get_number_components();

  switch (number_components) {

  default:
    self.fill(data_ptr);
    break;

  case 9:

    switch (component_order) {

    case ComponentOrder::CANONICAL:
      self.fill(data_ptr);
      break;

    case ComponentOrder::SIERRA_FULL:
      //  0  1  2  3  4  5  6  7  8
      // XX YY ZZ XY YZ ZX YX ZY XZ
      //  0  4  8  1  5  6  3  7  2
      self[0] = data_ptr[0];
      self[4] = data_ptr[1];
      self[8] = data_ptr[2];

      self[1] = data_ptr[3];
      self[5] = data_ptr[4];
      self[6] = data_ptr[5];

      self[3] = data_ptr[6];
      self[7] = data_ptr[7];
      self[2] = data_ptr[8];
      break;

    case ComponentOrder::SIERRA_SYMMETRIC:
      self[0] = data_ptr[0];
      self[4] = data_ptr[1];
      self[8] = data_ptr[2];

      self[1] = data_ptr[3];
      self[5] = data_ptr[4];
      self[6] = data_ptr[5];

      self[3] = data_ptr[3];
      self[7] = data_ptr[4];
      self[2] = data_ptr[5];
      break;

    default:
      MT_ERROR_EXIT("Unknown component order.");
      break;

    }

    break;
  }

  return;
}

//
// Component increment
//
template<typename T, typename ST>
template<typename S, typename SS>
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
void
TensorBase<T, ST>::clear()
{
  fill(Filler::ZEROS);
  return;
}

//
// Square of Frobenius norm
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
T
norm_f(TensorBase<T, ST> const & X)
{
  T const
  s = norm_f_square(X);

  if (s > 0.0) return std::sqrt(s);

  return 0.0;
}

//
// Base addition
//
template<typename R, typename S, typename T, typename SR, typename SS,
    typename ST>
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
void
minus(TensorBase<T, ST> const & A, TensorBase<T, ST> & B)
{
  Index const
  number_components = A.get_number_components();

  assert(B.get_number_components() == number_components);

  for (Index i = 0; i < number_components; ++i) {
    B[i] = -A[i];
  }

  return;
}

//
// Base equality
//
template<typename T, typename ST>
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
bool
not_equal(TensorBase<T, ST> const & A, TensorBase<T, ST> const & B)
{
  return !(equal(A, B));
}

//
// Base scaling
//
template<typename R, typename S, typename T, typename SR, typename ST>
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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

} // namespace minitensor

#endif // MiniTensor_TensorBase_i_h
