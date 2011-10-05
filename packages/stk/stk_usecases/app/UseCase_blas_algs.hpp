/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef USECASE_BLAS_ALGS_HPP
#define USECASE_BLAS_ALGS_HPP

#include <vector>
#include <algorithm>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_algsup/AlgorithmRunner.hpp>

//----------------------------------------------------------------------
/** AxpbyAlg calculates: y = alpha*x + beta*y
*/
template< class field_type >
class AxpbyAlg {
public:
  AxpbyAlg(double alpha, field_type& x,
           double beta,  field_type& y)
   : alpha_(alpha), beta_(beta), x_(x), y_(y)
  {
    // put a test here to insist on compatibility of x and y fields...
  }

  ~AxpbyAlg(){}

  enum { maximum_entity_count = 0 }; /* Don't slice buckets */

  void apply( stk::mesh::Bucket::iterator bucket_i ,
              stk::mesh::Bucket::iterator bucket_j ) const
  {
    typedef typename stk::mesh::FieldTraits<field_type>::data_type scalar ;

    stk::mesh::BucketArray<field_type> xa( x_ , bucket_i , bucket_j );
    stk::mesh::BucketArray<field_type> ya( y_ , bucket_i , bucket_j );

    scalar * xi = xa.contiguous_data();
    scalar * yi = ya.contiguous_data();
    scalar * ye = yi + ya.size();

    //is it worthwhile to optimize by adding special loops for
    //cases where alpha or beta equal 1.0 or 0.0, or where x==y?

    while(yi < ye) {
      *yi = alpha_* *xi++ + beta_* *yi;
      ++yi;
    }
  }

private:
  AxpbyAlg& operator=(const AxpbyAlg&);
  AxpbyAlg(const AxpbyAlg&);

  double alpha_, beta_;
  field_type& x_;
  field_type& y_;
};

//----------------------------------------------------------------------
/** FillAlg fills the specified field with the specified scalar...
*/
template< class field_type >
class FillAlg {
public:
  FillAlg(field_type& x, double scalar_value)
   : x_(x), scalar_(scalar_value) {}

  ~FillAlg(){}

  enum { maximum_entity_count = 0 }; /* Don't slice buckets */

  void apply( stk::mesh::Bucket::iterator bucket_i ,
              stk::mesh::Bucket::iterator bucket_j ) const
  {
    typedef typename stk::mesh::FieldTraits<field_type>::data_type scalar ;

    stk::mesh::BucketArray<field_type> xa( x_ , bucket_i , bucket_j );

    scalar * xi = xa.contiguous_data();
    scalar * xe = xi + xa.size();

    std::fill(xi, xe, scalar_);
  }

private:
  FillAlg& operator=(const FillAlg&);
  FillAlg(const FillAlg&);

  field_type& x_;
  double scalar_;
};


//----------------------------------------------------------------------
//Below are some helper functions for running the above algorithms.
//Are these necessary?
//(they are currently called in UseCase_blas.cpp)

template< class field_type >
inline
void axpby( const stk::AlgorithmRunnerInterface& alg_runner,
            stk::mesh::BulkData& bulk,
            stk::mesh::EntityRank entitytype,
            double alpha,
            field_type& x,
            double beta,
            field_type& y)
{
  AxpbyAlg<field_type> axpby_alg(alpha, x, beta, y);

  const std::vector<stk::mesh::Bucket*>
    & all_entity_buckets = bulk.buckets( entitytype );

  const stk::mesh::MetaData & meta = stk::mesh::MetaData::get(bulk);
  stk::mesh::Selector select_used =
    meta.locally_owned_part() |
    meta.globally_shared_part() ;

  stk::mesh::PartVector empty_union_vector;
  alg_runner.run( select_used , empty_union_vector, all_entity_buckets , axpby_alg );
}

template< class field_type >
inline
void fill( const stk::AlgorithmRunnerInterface& alg_runner,
           stk::mesh::BulkData& bulk,
           stk::mesh::EntityRank entitytype,
           field_type& x,
           double value)
{
  FillAlg<field_type> fillalg(x, value);

  const std::vector<stk::mesh::Bucket*>
    & all_entity_buckets = bulk.buckets( entitytype );

  const stk::mesh::MetaData & meta = stk::mesh::MetaData::get(bulk);
  stk::mesh::Selector select_used =
    meta.locally_owned_part() |
    meta.globally_shared_part() ;

  stk::mesh::PartVector empty_union_vector;
  alg_runner.run( select_used , empty_union_vector, all_entity_buckets , fillalg );
}

//----------------------------------------------------------------------

template< class field_type >
class DotAlg {
public:
  const field_type & field_x ;
  const field_type & field_y ;

  typedef typename stk::mesh::FieldTraits<field_type>::data_type scalar ;

  DotAlg( const field_type & x , const field_type & y )
    : field_x( x ), field_y( y ) {}

  /* Threaded-safe reduction */

  typedef double reduce_type ;
  enum { reduce_count = 1 };

  void init( reduce_type * out ) const { *out = 0 ; }
  void join( reduce_type * inout , const reduce_type * in ) const
    { *inout += *in ; }

  /* Primary thread-safe operation */

  enum { maximum_entity_count = 0 }; /* Don't slice buckets */

  void apply( stk::mesh::Bucket::iterator bucket_i ,
              stk::mesh::Bucket::iterator bucket_j ,
              reduce_type * inout ) const
  {
    stk::mesh::BucketArray<field_type> xa( field_x , bucket_i , bucket_j );
    stk::mesh::BucketArray<field_type> ya( field_y , bucket_i , bucket_j );

    scalar * xi = xa.contiguous_data();
    scalar * yi = ya.contiguous_data();
    scalar * ye = yi + ya.size();

    scalar tmp = 0 ;
    while ( yi != ye ) { tmp += *xi++ * *yi++ ; }
    *inout += tmp ;
  }

private:
  DotAlg& operator=(const DotAlg&);
  DotAlg(const DotAlg&);
};


template< class field_type >
typename stk::mesh::FieldTraits<field_type>::data_type
dot( const stk::AlgorithmRunnerInterface & alg_runner ,
     const stk::mesh::BulkData& bulk,
     const stk::mesh::EntityRank entitytype,
     const field_type & x ,
     const field_type & y )
{
  // Construct algorithm to hand to local parallel algorithm runner

  DotAlg<field_type> alg_dot(x,y);

  const std::vector<stk::mesh::Bucket*>
    & all_entity_buckets = bulk.buckets( entitytype );

  stk::mesh::Selector
    select_owned( stk::mesh::MetaData::get(bulk).locally_owned_part());

  double local_dot = 0 ;
  double global_dot = 0 ;

  // Local parallel execution of the algorithm

  stk::mesh::PartVector empty_union_vector;
  alg_runner.run( select_owned , empty_union_vector, all_entity_buckets , alg_dot , & local_dot );

  // Global sum

  MPI_Comm comm = bulk.parallel();

  MPI_Allreduce(& local_dot , & global_dot, 1, MPI_DOUBLE, MPI_SUM, comm);

  return global_dot ;
}

//----------------------------------------------------------------------

template< class field_type >
class Norm2Alg {
public:
  const field_type & field_x ;

  typedef typename stk::mesh::FieldTraits<field_type>::data_type scalar ;

  explicit Norm2Alg( const field_type & x ) : field_x( x ) {}

  /* Threaded-safe reduction */

  typedef double reduce_type ;
  enum { reduce_count = 1 };

  void init( reduce_type * out ) const { *out = 0 ; }
  void join( reduce_type * inout , const reduce_type * in ) const
    { *inout += *in ; }

  /* Primary thread-safe operation */

  enum { maximum_entity_count = 0 }; /* Don't slice buckets */

  void apply( stk::mesh::Bucket::iterator bucket_i ,
              stk::mesh::Bucket::iterator bucket_j ,
              reduce_type * inout ) const
  {
    stk::mesh::BucketArray<field_type> xa( field_x , bucket_i , bucket_j );

    scalar * xi = xa.contiguous_data();
    scalar * xe = xi + xa.size();

    scalar tmp = 0 ;
    while ( xi != xe ) { tmp += *xi * *xi ; ++xi; }

    *inout += tmp ;
  }

private:
  Norm2Alg& operator=(const Norm2Alg&);
  Norm2Alg(const Norm2Alg&);
};

template< class field_type >
typename stk::mesh::FieldTraits<field_type>::data_type
norm2( const stk::AlgorithmRunnerInterface & alg_runner ,
       const stk::mesh::BulkData           & bulk ,
       const stk::mesh::EntityRank           entitytype ,
       const field_type                    & x )
{
  // Construct algorithm to hand to local parallel algorithm runner

  Norm2Alg<field_type> alg_norm2(x);

  const std::vector<stk::mesh::Bucket*>
    & all_entity_buckets = bulk.buckets( entitytype );

  stk::mesh::Selector
    select_owned(stk::mesh::MetaData::get(bulk).locally_owned_part());

  double local_dot = 0 ;
  double global_dot = 0 ;

  // Local parallel execution of the algorithm

  stk::mesh::PartVector empty_union_vector;
  alg_runner.run( select_owned , empty_union_vector, all_entity_buckets , alg_norm2 , & local_dot );

  // Global sum

  MPI_Comm comm = bulk.parallel();

  MPI_Allreduce(& local_dot , & global_dot, 1, MPI_DOUBLE, MPI_SUM, comm);

  return std::sqrt( global_dot ) ;
}

#endif

