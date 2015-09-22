/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards
 */

#ifndef centroid_algorithm_hpp
#define centroid_algorithm_hpp

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <use_cases/UseCase_Common.hpp>

namespace stk_classic {
namespace mesh {
namespace use_cases {

enum {SpatialDim = 3};

//----------------------------------------------------------------------
// An example of a trivial element-block algorithm
// that just works on element field data
// written with a "generic programming" flavor.

template< class ElementTraits >
void centroid( unsigned number_elements ,
               double *  elem_centroid ,
               double ** elem_node_coordinates )
{
  enum { vertices_per_element = ElementTraits::vertex_count };
  enum { nodes_per_element    = ElementTraits::node_count };

  for ( unsigned j = 0 ; j < number_elements ; ++j ) {

    // Array is 'double * elem_node[ nodes_per_element * number_elements ]
    // so get the node coordinate pointer for this element.

    double ** node_coordinates = elem_node_coordinates + j * nodes_per_element ;

    double tmp[SpatialDim] = { 0 , 0 , 0 };

    for ( unsigned i = 0 ; i < vertices_per_element ; ++i ) {
      // Pointer to this node's coordinates,
      // accessing this field data in-place
      // as opposed to copying (a.k.a. gathering)
      // it into a local temporary array.

      double * coord = node_coordinates[i] ;

      // Sum for the mean.
      tmp[0] += coord[0] ;
      tmp[1] += coord[1] ;
      tmp[2] += coord[2] ;
    }

    tmp[0] /= vertices_per_element ;
    tmp[1] /= vertices_per_element ;
    tmp[2] /= vertices_per_element ;

    // The centroid field data for this element.

    double * centroid = elem_centroid + j * SpatialDim ;

    centroid[0] = tmp[0] ;
    centroid[1] = tmp[1] ;
    centroid[2] = tmp[2] ;
  }
}

//----------------------------------------------------------------------
// An example of a trivial element-loop algorithm
// written with a "generic programming" flavor.
//
// Loop over all elements in the given 'elem_part'.
// Access the element-node coordinates.

template< class ElementTraits >
void centroid_algorithm(
  const BulkData & bulkData ,
  const VectorFieldType             & elem_centroid ,
  const ElementNodePointerFieldType & elem_node_coord ,
  Part & elem_part,
  EntityRank element_rank )
{
  // Use the "homogeneous subset" concept (see the Domain Model document)
  // for field data storage.  A "homogeneous subset" is called
  // a 'Bucket'.

  // Iterate the set of element buckets:

  const std::vector<Bucket*> & buckets = bulkData.buckets( element_rank );

  for ( std::vector<Bucket*>::const_iterator
        k = buckets.begin() ; k != buckets.end() ; ++k ) {
    Bucket & bucket = **k ;

    // If this bucket is a subset of the given elem_part
    // then want to compute on it.

    if ( has_superset( bucket, elem_part ) ) {

      // Number of elements in the bucket:

      const unsigned size = bucket.size();

      // Aggressive "gather" field data for the elements
      // in the bucket.
      //   double * node_ptr[ nodes_per_element * number_of_elements ]

      double ** node_ptr = field_data( elem_node_coord , bucket.begin() );

      // Element centroid field data
      //   double elem_ptr[ SpatialDim * number_of_elements ]

      double *  elem_ptr = field_data( elem_centroid , bucket.begin() );

      // Call an element function to calculate centroid for
      // contiguous arrays of element field data.

      centroid< ElementTraits >( size , elem_ptr , node_ptr );
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template< class ElementTraits >
bool centroid_algorithm_unit_test_dimensions(
  const BulkData & bulkData ,
  const VectorFieldType             & elem_centroid ,
  const ElementNodePointerFieldType & elem_node_coord ,
  Part & elem_part,
  EntityRank element_rank )
{
  bool result = true;
  // Use the "homogeneous subset" concept (see the Domain Model document)
  // for field data storage.  A "homogeneous subset" is called
  // a 'Bucket'.

  // Iterate the set of element buckets:

  const std::vector<Bucket*> & buckets = bulkData.buckets( element_rank );

  for ( std::vector<Bucket*>::const_iterator
        k = buckets.begin() ; k != buckets.end() ; ++k ) {
     Bucket & bucket = **k ;

    // If this bucket is a subset of the given elem_part
    // then want to compute on it.

    if ( has_superset( bucket, elem_part ) ) {

      // Number of elements in the bucket:

      const unsigned size = bucket.size();

      // Unit testing the dimension feature
      {
        BucketArray< ElementNodePointerFieldType >
          array( elem_node_coord, bucket.begin(), bucket.end() );
        const unsigned n1 = array.template dimension<0>();
        const unsigned n2 = array.template dimension<1>();
        if ( n1 != ElementTraits::node_count ) {
          std::cerr << "Error!  n1 == " << n1 << " != " << ElementTraits::node_count << " == ElementTraits::node_count" << std::endl;
          result = false;
        }
        if ( n2 != size ) {
          std::cerr << "Error!  n2 == " << n2 << " != " << size << " == size" << std::endl;
          result = false;
        }
        if ( (unsigned) array.size() != n1 * n2 ) {
          std::cerr << "Error!  array.size() == " << array.size() << " != " << n1*n2 << " == n1*n2" << std::endl;
          result = false;
        }
      }

      {
        BucketArray< VectorFieldType > array( elem_centroid , bucket.begin(), bucket.end()  );
        const unsigned n1 = array.template dimension<0>();
        const unsigned n2 = array.template dimension<1>();
        if ( n1 != (unsigned) SpatialDim ) {
          std::cerr << "Error!  n1 == " << n1 << " != " << SpatialDim << " == SpatialDim" << std::endl;
          result = false;
        }
        if ( n2 != size ) {
          std::cerr << "Error!  n2 == " << n2 << " != " << size << " == size" << std::endl;
          result = false;
        }
        if ( (unsigned) array.size() != n1 * n2 ) {
          std::cerr << "Error!  array.size() == " << array.size() << " != " << n1*n2 << " == n1*n2" << std::endl;
          result = false;
        }
      }
    }
  }
  return result;
}

} // namespace use_cases
} // namespace mesh
} // namespace stk_classic

#endif // centroid_algorithm_hpp
