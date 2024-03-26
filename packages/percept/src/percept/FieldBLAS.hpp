// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_FieldBLAS_hpp
#define percept_FieldBLAS_hpp

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace percept
{

//: Scaled sum  Y = alpha * X + beta * Y
void field_axpby(
  const stk::mesh::BulkData &bulkdata ,
  const double     alpha ,
  stk::mesh::Field<double> *X,
  const double     beta ,
  stk::mesh::Field<double> *Y)
{
  const stk::mesh::MetaData & meta = bulkdata.mesh_meta_data();
  stk::mesh::Selector select_used =
    meta.locally_owned_part() |
    meta.globally_shared_part() ;

  stk::mesh::BucketVector const& entity_buckets =
    bulkdata.get_buckets( stk::topology::NODE_RANK, select_used );

  for ( stk::mesh::BucketVector::const_iterator ib = entity_buckets.begin();
        ib != entity_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib ;

    const double * x = stk::mesh::field_data(*X, b);
    double * y = stk::mesh::field_data(*Y, b);

    const stk::mesh::Bucket::size_type length = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      y[k] = alpha * x[k] + beta * y[k];
    }
  }
}

//: Y(ts) = X
void field_copy_component(
  const stk::mesh::BulkData &bulkdata ,
  stk::mesh::Field<double> *X,
  stk::mesh::Field<double> *Y,
  const int index)
{
  const stk::mesh::MetaData & meta = bulkdata.mesh_meta_data();
  const unsigned nDim = meta.spatial_dimension();

  const unsigned field_size_y = Y->max_size();

  stk::mesh::Selector select_used =
    meta.locally_owned_part() |
    meta.globally_shared_part() ;

  stk::mesh::BucketVector const& entity_buckets =
    bulkdata.get_buckets( stk::topology::NODE_RANK, select_used );

  for ( stk::mesh::BucketVector::const_iterator ib = entity_buckets.begin();
        ib != entity_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib ;

    const double * x = stk::mesh::field_data(*X, b);
    double * y = stk::mesh::field_data(*Y, b);

    const stk::mesh::Bucket::size_type length = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      const unsigned offset_y = k * field_size_y + index * nDim;
      const unsigned offset_x = k * nDim;

      for (unsigned idim=0; idim<nDim; idim++) {
	y[offset_y + idim] = x[offset_x + idim];
      }
    }
  }
}

// compute the MAC array for two basis sets
void field_compute_MAC(
  const stk::mesh::BulkData &bulkdata,
  const stk::mesh::MetaData &metadata,
  stk::mesh::Field<double> * X,
  stk::mesh::Field<double> * Y,
  Kokkos::View<double**, Kokkos::HostSpace> & mac_values)
{
  const unsigned field_size = Y->max_size();

  const unsigned nDim = metadata.spatial_dimension();
  const unsigned nVecs = field_size / nDim;

  Kokkos::View<double*, Kokkos::HostSpace> inner_product_x("inner_product_x",nVecs);
  Kokkos::View<double*, Kokkos::HostSpace> inner_product_y("inner_product_x",nVecs);

  stk::mesh::Selector select_used = metadata.locally_owned_part();

  stk::mesh::BucketVector const& entity_buckets =
    bulkdata.get_buckets( stk::topology::NODE_RANK, select_used );

  for ( stk::mesh::BucketVector::const_iterator ib = entity_buckets.begin();
        ib != entity_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib ;

    const double * x = stk::mesh::field_data(*X, b);
    const double * y = stk::mesh::field_data(*Y, b);

    const stk::mesh::Bucket::size_type length = b.size();

    unsigned offset_k, offset_mk, offset_nk;

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      offset_k = k * field_size;

      for (unsigned m=0; m<nVecs; m++) {

	offset_mk = offset_k + m * nDim;

	for (unsigned idim=0; idim<nDim; idim++) {
	  inner_product_x(m) += pow(x[offset_mk + idim], 2);
	  inner_product_y(m) += pow(y[offset_mk + idim], 2);
	}

	for (unsigned n=0; n<nVecs; n++) {

	  offset_nk = offset_k + n * nDim;

	  for (unsigned idim=0; idim<nDim; idim++) {
	    mac_values(m, n) += x[offset_mk + idim] * y[offset_nk + idim];
	  }
	}
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &*(mac_values.data()), mac_values.size(), MPI_DOUBLE, MPI_SUM, bulkdata.parallel());
  MPI_Allreduce(MPI_IN_PLACE, &*(inner_product_x.data()), inner_product_x.size(), MPI_DOUBLE, MPI_SUM, bulkdata.parallel());
  MPI_Allreduce(MPI_IN_PLACE, &*(inner_product_y.data()), inner_product_y.size(), MPI_DOUBLE, MPI_SUM, bulkdata.parallel());

  for (unsigned m=0; m<nVecs; m++) {
    for (unsigned n=0; n<nVecs; n++) {
      mac_values(m, n) = mac_values(m, n) / sqrt(inner_product_x(m) * inner_product_y(n));
    }
  }
}

}//namespace percept

#endif
