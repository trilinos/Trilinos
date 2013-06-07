/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef  STK_TRANSFER_HPP
#define  STK_TRANSFER_HPP

#include <boost/smart_ptr/shared_ptr.hpp>

#include <stk_util/util/StaticAssert.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_mesh/base/Comm.hpp>



#include <stk_search/BoundingBox.hpp>

namespace stk {

template <class INTERPOLATE> class GeometricTransfer : public class GeometricTransferBase {


typedef INTERPOLATE::MeshA       MeshA;
typedef INTERPOLATE::MeshB       MeshB;
typedef INTERPOLATE::EntityKeyMap EntityKeyMap;
typedef MeshA::EntityKey    EntityKeyA;
typedef MeshB::EntityKey    EntityKeyB;
typedef MeshA::EntityKeySet EntityKeySetA;
typedef MeshB::EntityKeySet EntityKeySetB;


typedef stk::search::ident::IdentProc<EntityKeyA,unsigned> IdentProcA;
typedef stk::search::ident::IdentProc<EntityKeyB,unsigned> IdentProcB;
typedef std::pair<IdentProcA, IdentProcB>                  IdentProcRelation;
typedef std::vector<IdentProcRelation>                     IdentProcRelationVec;

enum {spatial_dimension = MeshA::Dimension};

Transfer(MeshA &mesha, MeshB &meshb);

virtual void initialize();
virtual void apply();

private :

MeshA        &m_mesha;
MeshB        &m_meshb;
const std::string m_name;

EntityKeyMap TotalRangeToDomain;
template <class MESH>
void boundingbox_vector(std::vector<BoundingBox> &vector,
                        const MESH::EntityKeySet &keys_to_match,
                        const MESH               &mesh,
                        const BoundingBox::Data   radius,
                        const stk::ParallelMachine comm);

EntityKeyMap copy_domain_to_range_processors(MeshA                     &mesh,
                                     IdentProcRelationVec      &RangeToDomain,
                                     const std::string         &transfer_name);

void point_to_point_coarse_search(IdentProcRelationVec &RangeToDomain,
                                  const EntityKeySetA  &keys_to_matcha,
                                  const EntityKeySetB  &keys_to_matchb,
                                  const MeshB     &meshb,
                                  const MeshA     &mesha,
                                  const BoundingBox::Data radius,
                                  const stk::ParallelMachine comm);

void delete_range_points_found(EntityKeySetA            &to_keys,
                               const EntityKeyMap   &entity_keys);

enum { dim_eq = stk::StaticAssert<MeshB::Dimension==MeshA::Dimension::OK };

};




template <class INTERPOLATE> GeometricTransfer<INTERPOLATE>::Transfer(
       MeshA &mesha, MeshB &meshb, const std::string &name) :
  m_mesha(mesha),
  m_meshb(meshb),
  m_name (name) {}

template <class INTERPOLATE> void GeometricTransfer<INTERPOLATE>::initialize() {

  EntityKeySetA keys_to_matcha;
  EntityKeySetB keys_to_matchb;
  MeshA.keys(keys_to_matcha);;
  MeshB.keys(keys_to_matchb);;

  while (keys_to_matcha.size()) { // Keep going until all range points are processed.
    IdentProcRelationVec RangeToDomain;
    point_to_point_coarse_search(RangeToDomain, 
                                 keys_to_matcha, 
                                 keys_to_matchb,
                                 mesha, meshb, 
                                 radius);
    EntityKeyMap entity_key_map;
    {   
      stk::diag::TimeBlock __timer_ghosting(GhostingTimer);
      if (mesha.has_communication_capabilities()) {
        entity_key_map = copy_domain_to_range_processors(mesha, RangeToDomain, m_name);
      } else if (meshb.has_communication_capabilities()) {
        //copy_range_to_domain_processors(meshb, RangeToDomain, m_name);
        ThrowRequireMsg (mesha.has_communication_capabilities() || mesha.has_communication_capabilities(),
          __FILE__<<":"<<__LINE__<<" Expected "<<span<<" relations."<<" Found:"<<num_relations);
      } else {
        ThrowRequireMsg (mesha.has_communication_capabilities() || mesha.has_communication_capabilities(),
          __FILE__<<":"<<__LINE__<<" Expected "<<span<<" relations."<<" Found:"<<num_relations);
      }
    }   
    INTERPOLATE::filter_with_fine_search(entity_key_map, mesha, meshb); 
    TotalRangeToDomain.insert(entity_key_map.begin(), entity_key_map.end());

    delete_range_points_found(keys_to_matcha, entity_key_map); 
    radius *= ExpansionFactor; // If points were missed, increase search radius.
  }
}

template <class INTERPOLATE> void GeometricTransfer<INTERPOLATE>::delete_range_points_found(
                               EntityKeySetA            &to_keys,
                               const EntityKeyMap   &entity_keys) {
  for (EntityKeyMap::const_iterator i=entity_keys.begin(); i!=entity_keys.end(); ++i) to_keys.erase(i->first);
}



template <class INTERPOLATE> void GeometricTransfer<NDIM,INTERPOLATE>::apply(){
  INTERPOLATE::apply(MeshB, MeshA, TotalRangeToDomain, comm);
}

template <class INTERPOLATE> temlate <class MESH> void GeometricTransfer<NDIM,INTERPOLATE>::boundingbox_vector(
                         std::vector<BoundingBox> &vector,
                        const MESH::EntityKeyVec &keys_to_match,
                        const MESH            &mesh,
                        const BoundingBox::Data   radius,
                        const stk::ParallelMachine comm) {
  BoundingBox::Data center[DIM];

  for (EntityKeyVec::const_iterator k=keys_to_match.begin(); k!=keys_to_match.end(); ++k) {
    const stk::mesh::EntityKeyA Id = *k;
    const BoundingBox B = mesh.boundingbox(Id, radius);
    vector.push_back(B);
  }
}


template <class INTERPOLATE> 
EntityKeyMap GeometricTransfer<INTERPOLATE>::copy_domain_to_range_processors(MeshA &mesh,
                                     const IdentProcRelationVec         &RangeToDomain,
                                     const std::string         &transfer_name) {

  stk::ParallelMachine comm = mesh.comm();
  const unsigned my_rank = stk::parallel_machine_rank(comm);
  
  std::vector<MESH::EntityProc> entities_to_copy ;

  const IdentProcRelationVec::const_iterator end=RangeToDomain.end();
  for (IdentProcRelationVec::const_iterator i=RangeToDomain.begin(); i!=end; ++i) {
    const unsigned            domain_owning_rank = i->second.proc;
    const unsigned             range_owning_rank = i->first.proc;
    if (domain_owning_rank == my_rank && range_owning_rank != my_rank) {
      const EntityKeyA entity = i->second.ident;
      const unsigned owner_rank = mesh.owner_rank(entity);
      if (owner_rank == (int)my_rank) { // should always be true by construction.
        const MESH::EntityProc ep(entity, range_owning_rank);
        entities_to_copy.push_back(ep);
      }   
    }   
  }
 
  mesh.copy_entities(entities_to_copy, transfer_name);

  EntityKeyMap entity_key_map;
  for (IdentProcRelationVec::const_iterator i=RangeToDomain.begin(); i!=end; ++i) {
    const unsigned range_owning_rank = i->first.proc;
    if (range_owning_rank == my_rank) {
      const EntityKeyB range_entity  = i->first.ident;
      const EntityKeyA domain_entity = i->second.ident;
      std::pair<EntityKeyB,EntityKeyA> key_map(range_entity,domain_entity);
      entity_key_map.insert(key_map);
    }   
  }   
  return entity_key_map;
}

template <class INTERPOLATE>  void point_to_point_coarse_search<NDIM,INTERPOLATE>::(
                                  IdentProcRelationVec   &RangeToDomain,
                                  const EntityKeyVecA &keys_to_matcha,
                                  const EntityKeyVecB &keys_to_matchb,
                                  const MeshB     &meshb,
                                  const MeshA     &mesha,
                                  const BoundingBox::Data radius,
                                  const stk::ParallelMachine comm) {
  
  std::vector<BoundingBox> range_vector;
  std::vector<BoundingBox> domain_vector;

  boundingbox_vector<spatial_dimension>(domain_vector, keys_to_matcha, mesha, radius, comm);
  boundingbox_vector<spatial_dimension>(range_vector,  keys_to_matchb, meshb, radius, comm);

  stk::search::FactoryOrder order;
  order.m_communicator = comm;

  // Slightly confusing: coarse_search documentation has domain->range
  // relations sorted by domain key.  We want range->domain type relations
  // sorted on range key. It might appear we have the arguments revered
  // in coarse_search call, but really, this is what we want.
  stk::search::coarse_search(RangeToDomain, domain_vector, range_vector, order);
}




#endif

