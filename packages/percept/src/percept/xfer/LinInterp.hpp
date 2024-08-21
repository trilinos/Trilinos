// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef LinInterp_h
#define LinInterp_h


#include <string>
#include <vector>
#include <utility>

#include <stk_mesh/base/Entity.hpp>

#include "Intrepid2_CellTools.hpp"
#include <percept/element/intrepid/BasisTable.hpp>

#include <percept/xfer/TransferHelper.hpp>

namespace percept{

template <class FROM, class TO>
class LinInterp {

public :

  typedef FROM                            MeshA;
  typedef TO                              MeshB;
  typedef typename MeshA::EntityKey       EntityKeyA;
  typedef typename MeshB::EntityKey       EntityKeyB;
  typedef typename MeshA::EntityProc      EntityProcA;
  typedef typename MeshB::EntityProc      EntityProcB;

  typedef std::pair<EntityProcB, EntityProcA>          EntityProcRelation;
  typedef std::vector<EntityProcRelation>              EntityProcRelationVec;

  typedef std::multimap<EntityKeyB, EntityKeyA>        EntityKeyMap;

  static void filter_to_nearest(EntityKeyMap    &RangeToDomain,
      const MeshA     &FromElem,
      MeshB     &ToPoints) ;

  static void apply (MeshB         &ToPoints,
      const MeshA         &FromElem,
      const EntityKeyMap &RangeToDomain) ;

 private:
  static void apply_from_elem_field (
    const MeshA     &FromElem,
    MeshB           &ToPoints,
    stk::mesh::Entity theElem,
    stk::mesh::Entity theNode);

  static void apply_from_nodal_field (
    const MeshA     &FromElem,
    MeshB           &ToPoints,
    stk::mesh::Entity theElem,
    stk::mesh::Entity theNode,
    const std::vector<double> &isoParCoords);
};

template <class FROM, class TO>
void
LinInterp<FROM,TO>::filter_to_nearest (
  EntityKeyMap    &RangeToDomain,
  const MeshA     &FromElem,
  MeshB           &ToPoints) {

  using DynRankView = Kokkos::DynRankView<double, Kokkos::HostSpace>;

  typename MeshB::Adaptor adaptor;
  if (adaptor.filter_to_nearest(RangeToDomain, FromElem, ToPoints))
    {
      return;
    }

  const double parametric_tolerance = 0.5;

  const stk::mesh::BulkData &fromBulkData = FromElem.fromBulkData_;
  stk::mesh::BulkData &toBulkData = ToPoints.toBulkData_;

  const stk::mesh::FieldBase *fromcoordinates = FromElem.fromcoordinates_;
  const stk::mesh::FieldBase *tocoordinates   = ToPoints.tocoordinates_;

  stk::mesh::EntityRank toRank = ToPoints.toFields_[0]->entity_rank();

  // FixMe: Check nDim against Dimension to make sure they are equal
  const unsigned nDim = FromElem.fromMetaData_.spatial_dimension();

  typedef typename EntityKeyMap::iterator iterator;
  typedef typename EntityKeyMap::const_iterator const_iterator;

  for (const_iterator current_key=RangeToDomain.begin(); current_key!=RangeToDomain.end(); ) {

    double best_dist = std::numeric_limits<double>::max();

    const stk::mesh::EntityKey thePt  = current_key->first;
    stk::mesh::Entity theNode = toBulkData.get_entity(thePt);

    std::vector<double> coord(3);
    if (toRank==stk::topology::ELEMENT_RANK) {
      compute_element_centroid(*tocoordinates, theNode, &coord[0]);
    }
    else if (toRank==stk::topology::NODE_RANK) {
      compute_nodal_coords(*tocoordinates, theNode, &coord[0]);
    }

    DynRankView inputPhysicalPoints("LI:inputPhysicalPoints",1,1,nDim);

    double rz_coord[2];
    switch(ToPoints.transferType_) {
    case THREED_TO_THREED:
    case TWOD_TO_TWOD:
      std::copy(&coord[0],&coord[0]+nDim, inputPhysicalPoints.data());
      break;
    case TWOD_AXI_TO_THREED:
      rz_coord[0] = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
      rz_coord[1] = coord[2];
      std::copy(rz_coord,rz_coord+nDim, inputPhysicalPoints.data());
      break;
    default:
      std::exit(EXIT_FAILURE);
    }

    std::pair<iterator, iterator> keys=RangeToDomain.equal_range(current_key->first);
    iterator nearest = keys.second;

    for (iterator ii=keys.first; ii != keys.second; ++ii) {

      const stk::mesh::EntityKey theBox = ii->second;
      stk::mesh::Entity theElem = fromBulkData.get_entity(theBox);

      const stk::mesh::Bucket &theBucket = fromBulkData.bucket(theElem);
      const CellTopologyData * const bucket_cell_topo_data = stk::mesh::get_cell_topology(theBucket.topology()).getCellTopologyData();

      stk::mesh::Entity const* elem_node_rels = fromBulkData.begin_nodes(theElem);
      const int num_nodes = fromBulkData.num_nodes(theElem);

      DynRankView cellWorkset("LI:cellWorkset", 1, num_nodes, nDim);

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];

        const double * fromcoords = static_cast<double*>(stk::mesh::field_data(*fromcoordinates, node ));
        for ( unsigned j = 0; j < nDim; ++j ) {
	  cellWorkset(0,ni,j) = fromcoords[j];
        }
      }

      std::vector<double> isoParCoords(nDim);

      DynRankView outputParametricPoints("LI:outputParametricPoints",1,1,nDim);
      shards::CellTopology topo(bucket_cell_topo_data);

      double dist = 0.0;

      if (topo.getKey()==shards::Particle::key) {
        dist = 0.0;
        for ( unsigned j = 0; j < nDim; ++j ) {
          dist += std::pow(cellWorkset(0,0,j) - inputPhysicalPoints(0,0,j), 2);
        }
        dist = std::sqrt(dist);
      }
      else if (topo.getKey()==shards::Beam<2>::key) {
        // do nothing
      }
      else {
        /*std::cout << outputParametricPoints.rank() << " "
                  << inputPhysicalPoints.rank() << ""
                  << std::endl;
        std::cout << outputParametricPoints.size() << " "
                  << inputPhysicalPoints.size() << " "
                  << std::endl;*/
        Intrepid2::CellTools<Kokkos::HostSpace>::mapToReferenceFrame(outputParametricPoints,
                                                         inputPhysicalPoints,
                                                         cellWorkset,
                                                         topo);
        
        dist = parametricDistanceToEntity(outputParametricPoints.data(), topo);
      }

      if ( dist < (1.0 + parametric_tolerance) && dist < best_dist ) {

        best_dist = dist;

	for ( unsigned j = 0; j < nDim; ++j ) {
	  isoParCoords[j] = outputParametricPoints(0,0,j);
	}

        ToPoints.TransferInfo_[thePt] = isoParCoords;
        nearest = ii;
      }
    }

    current_key = keys.second;
    if (nearest != keys.first ) RangeToDomain.erase(keys.first, nearest);
    if (nearest != keys.second) RangeToDomain.erase(++nearest, keys.second);
  }
}

template <class FROM, class TO>
void
LinInterp<FROM,TO>::apply(
  MeshB              &ToPoints,
  const MeshA        &FromElem,
  const EntityKeyMap &RangeToDomain) {

  typename MeshB::Adaptor adaptor;
  if (adaptor.apply(ToPoints, FromElem, RangeToDomain))
    {
      return;
    }

  const stk::mesh::BulkData &fromBulkData = FromElem.fromBulkData_;
  stk::mesh::BulkData         &toBulkData = ToPoints.toBulkData_;

  stk::mesh::EntityRank fromRank = FromElem.fromFields_[0]->entity_rank();

  typename EntityKeyMap::const_iterator ii;
  for(ii=RangeToDomain.begin(); ii!=RangeToDomain.end(); ++ii ) {

    const stk::mesh::EntityKey thePt  = ii->first;
    const stk::mesh::EntityKey theBox = ii->second;

    if (1 != ToPoints.TransferInfo_.count(thePt)) {
      if (0 == ToPoints.TransferInfo_.count(thePt))
        throw std::runtime_error("Key not found in database");
      else
        throw std::runtime_error("Too many Keys found in database");
    }

    const std::vector<double> &isoParCoords_ = ToPoints.TransferInfo_[thePt];
    stk::mesh::Entity theNode =   toBulkData.get_entity(thePt);
    stk::mesh::Entity theElem = fromBulkData.get_entity(theBox);

    if (fromRank==stk::topology::ELEMENT_RANK) {
      apply_from_elem_field(FromElem, ToPoints, theElem, theNode);
    }
    else if (fromRank==stk::topology::NODE_RANK) {
      apply_from_nodal_field(FromElem, ToPoints, theElem, theNode, isoParCoords_);
    }
  }
}

template <class FROM, class TO>
void
LinInterp<FROM,TO>::apply_from_elem_field (
  const MeshA     &FromElem,
  MeshB           &ToPoints,
  stk::mesh::Entity theElem,
  stk::mesh::Entity theNode)
{
  stk::mesh::EntityRank toRank = ToPoints.toFields_[0]->entity_rank();
  const unsigned from_field_size = FromElem.fromFields_[0]->max_size();
  const double * fromField_data = (double *) stk::mesh::field_data(*(FromElem.fromFields_[0]), theElem);
  double         * toField_data = (double *) stk::mesh::field_data(*(ToPoints.toFields_[0]),  theNode);

  std::vector<double> coord(3);
  if (toRank==stk::topology::ELEMENT_RANK) {
    compute_element_centroid(*(ToPoints.tocoordinates_), theNode, &coord[0]);
  }
  else if (toRank==stk::topology::NODE_RANK) {
    compute_nodal_coords(*(ToPoints.tocoordinates_), theNode, &coord[0]);
  }

  switch(ToPoints.transferType_) {
  case THREED_TO_THREED:
  case TWOD_TO_TWOD:
    {
      // also if field is not comp of vector field (thvec, rzvec)
      for ( unsigned fi = 0; fi<from_field_size; fi++) {
        toField_data[fi] = fromField_data[fi];
      }
    }
    break;
  case TWOD_AXI_TO_THREED:
    {
      switch(ToPoints.srcFieldType_) {
      case SRC_FIELD:
        {
          for ( unsigned fi = 0; fi<from_field_size; fi++) {
            toField_data[fi] = fromField_data[fi];
          }
        }
        break;
      case SRC_RZN_FIELD:
        {
          const double radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
          const double inv_radius = (radius < 2*std::numeric_limits<double>::min()) ? 1.0 : 1./radius;

          // fromField_data[0] is the angular component
          toField_data[0] = fromField_data[0] * (-coord[1]*inv_radius);
          toField_data[1] = fromField_data[0] * ( coord[0]*inv_radius);
          toField_data[2] = 0;
        }
        break;
      case SRC_RZP_FIELD:
        {
          const double radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
          const double inv_radius = (radius < 2*std::numeric_limits<double>::min()) ? 1.0 : 1./radius;

          // fromField_data[0] is the radial component
          toField_data[0] = fromField_data[0] * ( coord[0]*inv_radius);
          toField_data[1] = fromField_data[0] * ( coord[1]*inv_radius);
          toField_data[2] = fromField_data[1]; // z
        }
      }
    }
  }
}

template <class FROM, class TO>
void
LinInterp<FROM,TO>::apply_from_nodal_field (
  const MeshA     &FromElem,
  MeshB           &ToPoints,
  stk::mesh::Entity theElem,
  stk::mesh::Entity theNode,
  const std::vector<double> &isoParCoords)
{
  using DynRankView = Kokkos::DynRankView<double, Kokkos::HostSpace>;

  const unsigned nDim = FromElem.fromMetaData_.spatial_dimension();
  const unsigned from_field_size = FromElem.fromFields_[0]->max_size();

  const stk::mesh::BulkData &fromBulkData = FromElem.fromBulkData_;

  const CellTopologyData * const cell_topo_data =
    stk::mesh::get_cell_topology(fromBulkData.bucket(theElem).topology()).getCellTopologyData();
  shards::CellTopology topo(cell_topo_data);

  stk::mesh::Entity const* elem_node_rels = fromBulkData.begin_nodes(theElem);
  const int num_nodes = fromBulkData.num_nodes(theElem);

  DynRankView Coeff("LI:Coeff", num_nodes, from_field_size);

  // now load the elemental values for future interpolation; fill in connected nodes
  for ( int ni = 0; ni < num_nodes; ++ni ) {
    stk::mesh::Entity node = elem_node_rels[ni];
    const double * fromField_data = (double *) stk::mesh::field_data(*(FromElem.fromFields_[0]), node);

    for ( unsigned fi = 0; fi<from_field_size; fi++) {
      Coeff(ni,fi) = fromField_data[fi];
    }
  }

  DynRankView inputParametricPoints("LI:inputParametricPoints", 1, nDim);

  std::copy(&isoParCoords[0],&isoParCoords[0]+nDim, inputParametricPoints.data());

  DynRankView basisVals("LI:basisVals", num_nodes, 1);

  if (topo.getKey()==shards::Particle::key) {
    basisVals(0,0) = 1.0;
  }
  else if(topo.getKey()==shards::Beam<2>::key) {
    basisVals(0,0) = 0.5;
    basisVals(1,0) = 0.5;    
  }
  else {
    auto HGRAD_Basis = BasisTable::getInstance()->getBasis(topo);

    HGRAD_Basis->getValues(basisVals, inputParametricPoints, Intrepid2::OPERATOR_VALUE);
  }

  // TODO: save off these basis values for re-use

  double * toField_data = (double *) stk::mesh::field_data(*(ToPoints.toFields_[0]), theNode);
  const double * coord = (double *) stk::mesh::field_data(*(ToPoints.tocoordinates_), theNode);

  switch(ToPoints.transferType_) {
  case THREED_TO_THREED:
  case TWOD_TO_TWOD:
  {
    for ( unsigned fi = 0; fi<from_field_size; fi++) {
      toField_data[fi] = 0.0;
    }
    for (int ni = 0; ni < num_nodes; ni++) {
      for ( unsigned fi = 0; fi<from_field_size; fi++) {
	toField_data[fi] += Coeff(ni,fi) * basisVals(ni, 0);
      }
    }
  }
  break;
  case TWOD_AXI_TO_THREED:
  {
    DynRankView from_field_values("LI:from_field_values", from_field_size);
    for ( unsigned fi = 0; fi<from_field_size; fi++) {
      from_field_values(fi) = 0;
      for (int ni = 0; ni < num_nodes; ni++) {
	      from_field_values(fi) += Coeff(ni,fi) * basisVals(ni, 0);
      }
    }
    
    switch(ToPoints.srcFieldType_) {
    case SRC_FIELD:
      {
        for ( unsigned fi = 0; fi<from_field_size; fi++) {
          toField_data[fi] = 0.0;
        }
        for (int ni = 0; ni < num_nodes; ni++) {
          for ( unsigned fi = 0; fi<from_field_size; fi++) {
            toField_data[fi] += Coeff(ni,fi) * basisVals(ni, 0);
          }
        }
      }
      break;
    case SRC_RZN_FIELD:
      {
        const double radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
        const double inv_radius = (radius < 2*std::numeric_limits<double>::min()) ? 1.0 : 1./radius;

        // fromField_data[0] is the angular component
        toField_data[0] = from_field_values(0) * (-coord[1]*inv_radius);
        toField_data[1] = from_field_values(0) * ( coord[0]*inv_radius);
        toField_data[2] = 0;
      }
      break;
      case SRC_RZP_FIELD:
        {
          const double radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
          const double inv_radius = (radius < 2*std::numeric_limits<double>::min()) ? 1.0 : 1./radius;

          // fromField_data[0] is the radial component
          toField_data[0] = from_field_values(0) * ( coord[0]*inv_radius);
          toField_data[1] = from_field_values(0) * ( coord[1]*inv_radius);
          toField_data[2] = from_field_values(1); // z
        }
      }
    }
  }
}

} // namespace percept

#endif
