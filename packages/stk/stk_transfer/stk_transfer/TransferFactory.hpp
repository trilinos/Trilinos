/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef STK_Transfer_hpp
#define STK_Transfer_hpp

#include <Intrepid_FieldContainer.hpp>

#include <stk_util/diag/Timer.hpp>

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>


#include <iostream>

#include <boost/shared_ptr.hpp>

#include <stk_transfer/Transfer.hpp>

namespace stk {
namespace transfer {

struct InterpolationParams {
  InterpolationParams() : 
    Tolerance(0), ExpansionFactor(0), Timer(0), comm() {} 
  double                Tolerance;
  double          ExpansionFactor;
  stk::diag::Timer         *Timer;
  const stk::ParallelMachine comm;
};

boost::shared_ptr<TransferBase> NodeToNode(
              std::vector<stk::mesh::FieldBase*>         ToFields,
              stk::mesh::Selector                         ToNodes,
              const stk::mesh::VectorFieldType     &ToCoordinates,
              const std::vector<stk::mesh::FieldBase*> FromFields,
              const stk::mesh::Selector                 FromNodes,
              const stk::mesh::VectorFieldType   &FromCoordinates,
              const InterpolationParams               &InterpParam) {

  stk::mesh::EntityVec from_nodes;
  stk::mesh::BulkData &from_bulk_data = FromCoordinates.get_mesh();
  stk::mesh::MetaData &from_meta_data = stk::mesh::MetaData::get(from_bulk_data);
  stk::mesh::get_selected_entities(FromNodes,
                                   from_bulk_data.buckets(stk::mesh::MetaData::NODE_RANK),
                                   from_nodes);

  stk::mesh::EntityVec to_nodes;
  stk::mesh::BulkData &to_bulk_data = ToCoordinates.get_mesh();
  stk::mesh::MetaData &to_meta_data = stk::mesh::MetaData::get(to_bulk_data);
  stk::mesh::get_selected_entities(ToNodes,
                                   to_bulk_data.buckets(stk::mesh::MetaData::NODE_RANK),
                                   to_nodes);

  boost::shared_ptr<TransferBase> transfer_base;
  const unsigned spatial_dimension = FromCoordinates.max_size(stk::mesh::MetaData::NODE_RANK);
  switch (spatial_dimension) {
    case 1 :
      stk::STKMesh FromMesh<1>(from_nodes, FromCoordinates,  FromFields);
      stk::STKMesh   ToMesh<1>(  to_nodes,   ToCoordinates,    ToFields);
      transfer_base.reset(new Transfer<LinearInterpoate<STKMesh,STKMesh> >(FromMesh, ToMesh, InterpParam));
      break;
    case 2 :
      stk::STKMesh FromMesh<2>(from_nodes, FromCoordinates,  FromFields);
      stk::STKMesh   ToMesh<2>(  to_nodes,   ToCoordinates,    ToFields);
      transfer_base.reset(new Transfer<LinearInterpoate<STKMesh,STKMesh> >(FromMesh, ToMesh, InterpParam));
      break;
    case 3 :
      stk::STKMesh FromMesh<3>( from_nodes, FromCoordinates,  FromFields);
      stk::STKMesh   ToMesh<3>(   to_nodes,   ToCoordinates,    ToFields);
      transfer_base.reset(new Transfer<LinearInterpoate<STKMesh,STKMesh> >(FromMesh, ToMesh, InterpParam));
      break;
    }
  return transfer_base;
}

#endif
