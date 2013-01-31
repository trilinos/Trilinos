#ifndef TEST_TILE_HPP
#define TEST_TILE_HPP

#include <KokkosArray_ParallelReduce.hpp>
#include <KokkosArray_View.hpp>

template < typename Device , typename TileLayout>
struct ReduceTileErrors
{
  typedef Device device_type;
  typedef typename Device::size_type size_type;

  typedef KokkosArray::View< ptrdiff_t**, TileLayout, Device> array_type;

  typedef ptrdiff_t value_type;

  ReduceTileErrors( array_type a )
    : m_array(a)
  {}


  KOKKOSARRAY_INLINE_FUNCTION
  static void init( value_type & errors )
  {
    errors = 0;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  static void join( volatile value_type & errors ,
                    const volatile value_type & src_errors )
  {
    errors += src_errors;
  }


  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( size_type global_i , value_type & errors ) const
  {
    typedef typename array_type::tile_type tile_type;

    size_t t0 = m_array.global_to_tile_index_0(global_i);
    size_t i = m_array.global_to_local_tile_index_0(global_i);

    for (size_t global_j = 0, dim_1 = m_array.dimension_1(); global_j < dim_1; ++global_j) {

      ptrdiff_t offset = &m_array(global_i,global_j) - &m_array(0,0);

      size_t t1 = m_array.global_to_tile_index_1(global_j);
      size_t j = m_array.global_to_local_tile_index_1(global_j);

      tile_type tile = m_array.tile(t0,t1);
      tile(i,j) = (tile.dimension_0() * tile.dimension_1()) * (t0 + m_array.tiles_in_dimension_0() * t1) + i + (j * tile.dimension_0());

      errors += (offset != m_array(global_i,global_j));
    }
  }


  array_type m_array;
};


#endif //TEST_TILE_HPP
