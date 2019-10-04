// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef ToMeshDef_hpp
#define ToMeshDef_hpp

namespace percept
{

template<class T>
void
ToMesh<T>::bounding_boxes(std::vector<BoundingBox> &v_box) const
{
  stk::mesh::EntityRank toRank = toFields_[0]->entity_rank();
  // std::cout << "P[" << toBulkData_.parallel_rank() << "] toRank = " << toRank << " transferType_= " 
  //           << transferType_ << " THREED_TO_THREED= " << THREED_TO_THREED << std::endl;

  stk::mesh::Selector sel = toMetaData_.locally_owned_part();
  Adaptor adaptor;
  adaptor.modify_selector(*this, sel);

  stk::mesh::BucketVector const& entity_buckets =
    toBulkData_.get_buckets( toRank, sel);

  for ( stk::mesh::BucketVector::const_iterator ib = entity_buckets.begin();
        ib != entity_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      stk::mesh::Entity entity = b[k];

      std::vector<double> coords(3);
      if (toRank==stk::topology::ELEMENT_RANK) {
        compute_element_centroid(*tocoordinates_, entity, &coords[0]);
      }
      else if (toRank==stk::topology::NODE_RANK) {
        compute_nodal_coords(*tocoordinates_, entity, &coords[0]);
      }

      // setup ident
      EntityProc theIdent(toBulkData_.entity_key(entity), toBulkData_.parallel_rank());

      BoundingBox theBox;
      double r, z;

      switch(transferType_)
      {
      case THREED_TO_THREED:
        theBox = BoundingBox(stk::search::Sphere<double>( stk::search::Point<double>(coords[0],coords[1],coords[2]), radius_), theIdent );
        break;
      case TWOD_TO_TWOD:
        theBox = BoundingBox(stk::search::Sphere<double>( stk::search::Point<double>(coords[0],coords[1]), radius_), theIdent );
        break;
      case TWOD_AXI_TO_THREED:
        r = sqrt(coords[0]*coords[0]+coords[1]*coords[1]);
        z = coords[2];

        theBox = BoundingBox( stk::search::Sphere<double>( stk::search::Point<double>(r,z), radius_), theIdent );
        break;
      default:
        std::exit(EXIT_FAILURE);
      }

      v_box.push_back(theBox);
    }
  }
}

template<class T>
void
ToMesh<T>::update_values () {
  std::vector<const stk::mesh::FieldBase *> fields;
  for (unsigned ii=0; ii < toFields_.size(); ++ii)
    fields.push_back(toFields_[ii]);
  stk::mesh::copy_owned_to_shared(toBulkData_, fields);
}

} // namespace percept

#endif
