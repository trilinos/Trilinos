// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Geometry_i_h)
#define MiniTensor_Geometry_i_h


namespace minitensor {

//
// Helper functions for determining the type of element
//
namespace {

KOKKOS_INLINE_FUNCTION
ELEMENT::Type
find_type_1D(Index const nodes)
{
  switch (nodes) {
    case 2:   return ELEMENT::SEGMENTAL;
    default:  return ELEMENT::UNKNOWN;
  }
}

KOKKOS_INLINE_FUNCTION
ELEMENT::Type
find_type_2D(Index const nodes)
{
  switch (nodes) {
    case 3:   return ELEMENT::TRIANGULAR;
    case 4:   return ELEMENT::QUADRILATERAL;
    default:  return ELEMENT::UNKNOWN;
  }
}

KOKKOS_INLINE_FUNCTION
ELEMENT::Type
find_type_3D(Index const nodes)
{
  switch (nodes) {
    case 4:   return ELEMENT::TETRAHEDRAL;
    case 8:   return ELEMENT::HEXAHEDRAL;
    default:  return ELEMENT::UNKNOWN;
  }
}

} // anonymous namespace


//
//
//
KOKKOS_INLINE_FUNCTION
ELEMENT::Type
find_type(Index const dimension, Index const number_nodes)
{

  ELEMENT::Type
  type = ELEMENT::UNKNOWN;

  switch (dimension) {

    case 1:
      type = find_type_1D(number_nodes);
      break;

    case 2:
      type = find_type_2D(number_nodes);
      break;

    case 3:
      type = find_type_3D(number_nodes);
      break;

    default:
      break;

  }

  if (type == ELEMENT::UNKNOWN) {
    MT_ERROR_EXIT("Unknown element type.");
  }

  return type;
}

//
// Constructor for SphericalParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
SphericalParametrization<T, N>::SphericalParametrization(
    Tensor4<T, N> const & A) : tangent_(A)
{
#if defined(KOKKOS_ENABLE_CUDA)
  minimum_ = DBL_MAX;
  maximum_ = DBL_MIN;
#else
  minimum_ = std::numeric_limits<T>::max();
  maximum_ = std::numeric_limits<T>::min();
#endif
  return;
}

//
// Normal vector for SphericalParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
SphericalParametrization<T, N>::get_normal(
    Vector<T, dimension_const<N, 2>::value> const & parameters
) const
{
  T const &
  phi = parameters(0);

  T const &
  theta = parameters(1);

  Vector<T, N>
  normal(sin(phi) * sin(theta), cos(phi), sin(phi) * cos(theta));

  return normal;
}

//
// Evaluation for SphericalParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
SphericalParametrization<T, N>::operator()(
    Vector<T, dimension_const<N, 2>::value> const & parameters
)
{
  assert(parameters.get_dimension() == 2);

  Vector<T, N> const
  normal = get_normal(parameters);

  // Localization tensor
  Tensor<T, N> const
  Q = dot2(normal, dot(tangent_, normal));

  T const
  determinant = det(Q);

  if (determinant < minimum_) {
    minimum_ = determinant;
    arg_minimum_ = parameters;
  }

  if (determinant > maximum_) {
    maximum_ = determinant;
    arg_maximum_ = parameters;
  }

  return;
}

//
// Constructor for StereographicParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
StereographicParametrization<T, N>::StereographicParametrization(
    Tensor4<T, N> const & A) : tangent_(A)
{
#if defined(KOKKOS_ENABLE_CUDA)
  minimum_ = DBL_MAX;
  maximum_ = DBL_MIN;
#else
  minimum_ = std::numeric_limits<T>::max();
  maximum_ = std::numeric_limits<T>::min();
#endif
  return;
}

//
// Normal vector for StereographicParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
StereographicParametrization<T, N>::get_normal(
    Vector<T, dimension_const<N, 2>::value> const & parameters
) const
{
  T const &
  x = parameters(0);

  T const &
  y = parameters(1);

  T const
  r2 = x * x + y * y;

  Vector<T, N>
  normal(2.0 * x, 2.0 * y, r2 - 1.0);

  normal /= (r2 + 1.0);

  return normal;
}

//
// Evaluation for StereographicParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
StereographicParametrization<T, N>::operator()(
    Vector<T, dimension_const<N, 2>::value> const & parameters
)
{
  assert(parameters.get_dimension() == 2);

  Vector<T, N> const
  normal = get_normal(parameters);

  // Localization tensor
  Tensor<T, N> const
  Q = dot2(normal, dot(tangent_, normal));

  T const
  determinant = det(Q);

  if (determinant < minimum_) {
    minimum_ = determinant;
    arg_minimum_ = parameters;
  }

  if (determinant > maximum_) {
    maximum_ = determinant;
    arg_maximum_ = parameters;
  }

  return;
}

//
// Constructor for ProjectiveParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
ProjectiveParametrization<T, N>::ProjectiveParametrization(
    Tensor4<T, N> const & A) : tangent_(A)
{
#if defined(KOKKOS_ENABLE_CUDA)
  minimum_ = DBL_MAX;
  maximum_ = DBL_MIN;
#else
  minimum_ = std::numeric_limits<T>::max();
  maximum_ = std::numeric_limits<T>::min();
#endif
  return;
}

//
// Normal vector for ProjectiveParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
ProjectiveParametrization<T, N>::get_normal(
    Vector<T, dimension_const<N, 3>::value> const & parameters
) const
{
  T const &
  x = parameters(0);

  T const &
  y = parameters(1);

  T const &
  z = parameters(2);

  Vector<T, N>
  normal(x, y, z);

  T const
  n = norm(normal);

  if (n > 0.0) {
    normal /= n;
  } else {
    normal = Vector<T, N>(1.0, 1.0, 1.0);
  }

  return normal;
}

//
// Evaluation for ProjectiveParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
ProjectiveParametrization<T, N>::operator()(
    Vector<T, dimension_const<N, 3>::value> const & parameters
)
{
  assert(parameters.get_dimension() == 3);

  Vector<T, N> const
  normal = get_normal(parameters);

  // Localization tensor
  Tensor<T, N> const
  Q = dot2(normal, dot(tangent_, normal));

  T const
  determinant = det(Q);

  if (determinant < minimum_) {
    minimum_ = determinant;
    arg_minimum_ = parameters;
  }

  if (determinant > maximum_) {
    maximum_ = determinant;
    arg_maximum_ = parameters;
  }

  return;
}

//
// Constructor for TangentParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
TangentParametrization<T, N>::TangentParametrization(
    Tensor4<T, N> const & A) : tangent_(A)
{
#if defined(KOKKOS_ENABLE_CUDA)
  minimum_ = DBL_MAX;
  maximum_ = DBL_MIN;
#else
  minimum_ = std::numeric_limits<T>::max();
  maximum_ = std::numeric_limits<T>::min();
#endif
  return;
}

//
// Normal vector for TangentParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
TangentParametrization<T, N>::get_normal(
    Vector<T, dimension_const<N, 2>::value> const & parameters
) const
{
  T const &
  x = parameters(0);

  T const &
  y = parameters(1);

  T const
  r = std::sqrt(x * x + y * y);

  Vector<T, N>
  normal(3, Filler::ZEROS);

  if (r > 0.0) {
    normal(0) = x * sin(r) / r;
    normal(1) = y * sin(r) / r;
    normal(2) = cos(r);
  } else {
    normal(2) = 1.0;
  }

  return normal;
}

//
// Evaluation for TangentParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
TangentParametrization<T, N>::operator()(
    Vector<T, dimension_const<N, 2>::value> const & parameters
)
{
  assert(parameters.get_dimension() == 2);

  Vector<T, N> const
  normal = get_normal(parameters);

  // Localization tensor
  Tensor<T, N> const
  Q = dot2(normal, dot(tangent_, normal));

  T const
  determinant = det(Q);

  if (determinant < minimum_) {
    minimum_ = determinant;
    arg_minimum_ = parameters;
  }

  if (determinant > maximum_) {
    maximum_ = determinant;
    arg_maximum_ = parameters;
  }

  return;
}

//
// Constructor for CartesianParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
CartesianParametrization<T, N>::CartesianParametrization(
    Tensor4<T, N> const & A) : tangent_(A)
{
#if defined(KOKKOS_ENABLE_CUDA)
  minimum_ = DBL_MAX;
  maximum_ = DBL_MIN;
#else
  minimum_ = std::numeric_limits<T>::max();
  maximum_ = std::numeric_limits<T>::min();
  return;
#endif
}

//
// Normal vector for CartesianParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
CartesianParametrization<T, N>::get_normal(
    Vector<T, dimension_const<N, 3>::value> const & parameters
) const
{
  T const &
  x = parameters(0);

  T const &
  y = parameters(1);

  T const
  z = parameters(2);

  Vector<T, N> const
  normal(x, y, z);

  return normal;
}

//
// Evaluation for CartesianParametrization
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
void
CartesianParametrization<T, N>::operator()(
    Vector<T, dimension_const<N, 3>::value> const & parameters
)
{
  Vector<T, N>
  normal = get_normal(parameters);

  // Localization tensor
  Tensor<T, N> const
  Q = dot2(normal, dot(tangent_, normal));

  T const
  determinant = det(Q);

  if (determinant < minimum_) {
    minimum_ = determinant;
    arg_minimum_ = parameters;
  }

  if (determinant > maximum_) {
    maximum_ = determinant;
    arg_maximum_ = parameters;
  }

  return;
}

//
// Constructor for ParametricGrid
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
ParametricGrid<T, N>::ParametricGrid(
    Vector<T, N> const & lower,
    Vector<T, N> const & upper,
    Vector<Index, N> const & points_per_dimension)
{
  assert(lower.get_dimension() == upper.get_dimension());
  assert(lower.get_dimension() == points_per_dimension.get_dimension());

  lower_ = lower;
  upper_ = upper;
  points_per_dimension_ = points_per_dimension;

  return;
}

//
// Traverse the grid and apply the visitor to each point.
//
template<typename T, Index N>
template<typename Visitor>
KOKKOS_INLINE_FUNCTION
void
ParametricGrid<T, N>::traverse(Visitor & visitor) const
{
  // Loop over the grid
  Index const
  number_parameters = lower_.get_dimension();

  LongIndex
  total_number_points = 1;

  for (Index dimension = 0; dimension < number_parameters; ++dimension) {
    total_number_points *= points_per_dimension_(dimension);
  }

  Vector<LongIndex, N> steps(number_parameters, Filler::ONES);

  for (Index dimension = 1; dimension < number_parameters; ++dimension) {
    steps(dimension) =
        steps(dimension - 1) * points_per_dimension_(dimension - 1);
  }

  Vector<Index, N>
  indices(number_parameters, Filler::ZEROS);

  Vector<T, N>
  position_in_grid(number_parameters, Filler::ZEROS);

  Vector<T, N> const
  span = upper_ - lower_;

  for (LongIndex point = 1;  point <= total_number_points; ++point) {

    //std::cout << "Indices : ";

    for (Index dimension = 0; dimension < number_parameters; ++dimension) {

/*           if ( points_per_dimension_(dimension) == 1 ) {     

             position_in_grid(dimension) = lower_(dimension);
        
          } else {
      
            position_in_grid(dimension) = indices(dimension) * span(dimension) /
                    (points_per_dimension_(dimension) - 1) + lower_(dimension);
          }
*/
        position_in_grid(dimension) = indices(dimension) * span(dimension) /
              (points_per_dimension_(dimension) - 1) + lower_(dimension);
        visitor(position_in_grid);

      //std::cout << indices(dimension) << " ";

      // Check if index needs to be increased or rolled back
      if (point % steps(dimension) == 0) {
        ++indices(dimension);
      }
      if (indices(dimension) == points_per_dimension_(dimension)) {
        indices(dimension) = 0;
      }

    }
     visitor(position_in_grid);
    //std::cout << std::endl;
    //std::cout << "Position : " << position_in_grid << std::endl;

  }

  return;
}

} // namespace minitensor

#endif // MiniTensor_Geometry_i_h
