/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef  STK_GEOMETRICTRANSFER_HPP
#define  STK_GEOMETRICTRANSFER_HPP

#include <boost/smart_ptr/shared_ptr.hpp>

#include <stk_util/util/StaticAssert.hpp>
#include <stk_util/environment/ReportHandler.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_transfer/TransferBase.hpp>

namespace stk {
namespace transfer {

template <class INTERPOLATE> class GeometricTransfer : public TransferBase {

public :
  
  typedef INTERPOLATE                                     InterpolateClass;
  typedef typename InterpolateClass::MeshA                MeshA;
  typedef typename InterpolateClass::MeshB                MeshB;
  typedef typename MeshA::EntityKey                       EntityKeyA;
  typedef typename MeshB::EntityKey                       EntityKeyB;
  typedef typename std::multimap<EntityKeyB, EntityKeyA>  EntityKeyMap;
  typedef typename std::set     <EntityKeyA>              EntityKeySetA;
  typedef typename std::set     <EntityKeyB>              EntityKeySetB;
  typedef typename MeshA::BoundingBox                     BoundingBoxA;
  typedef typename MeshA::BoundingBox                     BoundingBoxB;
  
  typedef typename MeshA::EntityProc                      EntityProcA;
  typedef typename MeshB::EntityProc                      EntityProcB;

  typedef std::pair<EntityProcB, EntityProcA>             EntityProcRelation;
  typedef std::vector<EntityProcRelation>                 EntityProcRelationVec;
  
  enum {Dimension = MeshA::Dimension};
  
  
  GeometricTransfer(MeshA &mesha, 
                    MeshB &meshb,
                    const double radius,
                    const double expansion_factor,
                    const std::string &name);
  virtual void initialize();
  virtual void apply();
  
  
private :
  
  MeshA        &m_mesha;
  MeshB        &m_meshb;
  const double      m_radius;
  const double      m_expansion_factor;
  const std::string m_name;
  
  EntityKeyMap TotalRangeToDomain;
  
  
  template <class MESH>
  void boundingboxes(std::vector<typename MESH::BoundingBox> &vector,
                     const typename MESH::EntityKeySet       &keys_to_match,
                     const MESH                              &mesh,
                     const double                             radius);
  
  EntityKeyMap copy_domain_to_range_processors(MeshA                            &mesh,
                                               const EntityProcRelationVec      &RangeToDomain,
                                               const std::string                &transfer_name);
  
  void coarse_search(EntityProcRelationVec   &RangeToDomain,
                     const EntityKeySetA     &keys_to_matcha,
                     const EntityKeySetB     &keys_to_matchb,
                     const MeshA             &mesha,
                     const MeshB             &meshb,
                     const double             radius);

  void delete_range_points_found(EntityKeySetB            &to_keys,
                                 const EntityKeyMap   &entity_keys);
  
  enum { dim_eq = StaticAssert<MeshB::Dimension==MeshA::Dimension>::OK };
  
};



template <class INTERPOLATE> GeometricTransfer<INTERPOLATE>::GeometricTransfer (MeshA            &mesha,
                                                                                MeshB            &meshb,
                                                                                const double      radius,
                                                                                const double      expansion_factor,
                                                                                const std::string &name) :
  m_mesha(mesha),
  m_meshb(meshb),
  m_radius(radius),
  m_expansion_factor(expansion_factor),
  m_name (name) {}

template <class INTERPOLATE> void GeometricTransfer<INTERPOLATE>::initialize() {

  EntityKeySetA keys_to_matcha;
  EntityKeySetB keys_to_matchb;
  m_mesha.keys(keys_to_matcha);;
  m_meshb.keys(keys_to_matchb);;

  double radius = m_radius;
  while (keys_to_matchb.size()) { // Keep going until all range points are processed.
    EntityProcRelationVec RangeToDomain;
    coarse_search(RangeToDomain, 
                  keys_to_matcha, 
                  keys_to_matchb,
                  m_mesha, 
                  m_meshb, 
                  radius);
    EntityKeyMap entity_key_map;
    {   
//      diag::TimeBlock __timer_ghosting(GhostingTimer);
      ParallelMachine comm = m_mesha.comm();
      const unsigned p_size = parallel_machine_size(comm);
      if (m_mesha.has_communication_capabilities()) {
        entity_key_map = copy_domain_to_range_processors(m_mesha, RangeToDomain, m_name);
      } else if (m_meshb.has_communication_capabilities()) {
        ThrowRequireMsg (m_mesha.has_communication_capabilities() || m_mesha.has_communication_capabilities(),
          __FILE__<<":"<<__LINE__<<" Still working on communicaiton capabilities");
      } else if (1==p_size) {
        const typename EntityProcRelationVec::const_iterator end=RangeToDomain.end();
        for (typename EntityProcRelationVec::const_iterator i=RangeToDomain.begin(); i!=end; ++i) {
          const EntityKeyB range_entity    = i->first.ident;
          const EntityKeyA domain_entity   = i->second.ident;
          std::pair<EntityKeyB,EntityKeyA> key_map(range_entity, domain_entity);
          entity_key_map.insert(key_map);
        }
      } else               {
        ThrowRequireMsg (m_mesha.has_communication_capabilities() || m_mesha.has_communication_capabilities(),
          __FILE__<<":"<<__LINE__<<" Still working on communicaiton capabilities");
      }
    }   
    INTERPOLATE::filter_to_nearest(entity_key_map, m_mesha, m_meshb); 
    TotalRangeToDomain.insert(entity_key_map.begin(), entity_key_map.end());

    delete_range_points_found(keys_to_matchb, entity_key_map); 
    radius *= m_expansion_factor; // If points were missed, increase search radius.
  }
}

template <class INTERPOLATE> void GeometricTransfer<INTERPOLATE>::delete_range_points_found(
                               EntityKeySetB            &to_keys,
                               const EntityKeyMap       &del) {
  for (typename EntityKeyMap::const_iterator i=del.begin(); i!=del.end(); ++i) {
    const EntityKeyB e = i->first;
    to_keys.erase(e);
  }
}



template <class INTERPOLATE> void GeometricTransfer<INTERPOLATE>::apply(){
  m_mesha.update_values();
  INTERPOLATE::apply(m_meshb, m_mesha, TotalRangeToDomain);
}

template <class INTERPOLATE> template <class MESH> void GeometricTransfer<INTERPOLATE>::boundingboxes(
  std::vector<typename MESH::BoundingBox>        &vector,
  const typename MESH::EntityKeySet              &keys_to_match,
  const MESH                                     &mesh,
  const double                                    radius) {

  typedef typename MESH::BoundingBox        BoundingBox;
  typedef typename MESH::BoundingBox::Data  Data;
  typedef typename MESH::EntityKey          EntityKey;
  typedef typename MESH::EntityKeySet       EntityKeySet;
  vector.clear();
  vector.reserve(keys_to_match.size());
  const Data r = radius;
  for (typename EntityKeySet::const_iterator k=keys_to_match.begin(); k!=keys_to_match.end(); ++k) {
    const EntityKey  Id = *k;
    const BoundingBox B = mesh.boundingbox(Id, r);
    vector.push_back(B);
  }
}


template <class INTERPOLATE>
typename GeometricTransfer<INTERPOLATE>::EntityKeyMap GeometricTransfer<INTERPOLATE>::copy_domain_to_range_processors(
  MeshA                               &mesh,
  const EntityProcRelationVec         &RangeToDomain,
  const std::string                   &transfer_name) {
  
  ParallelMachine comm = mesh.comm();
  const unsigned my_rank = parallel_machine_rank(comm);
  
  typename MeshA::EntityProcVec entities_to_copy ;

  const typename EntityProcRelationVec::const_iterator end=RangeToDomain.end();
  for (typename EntityProcRelationVec::const_iterator i=RangeToDomain.begin(); i!=end; ++i) {
    const unsigned            domain_owning_rank = i->second.proc;
    const unsigned             range_owning_rank = i->first.proc;
    if (domain_owning_rank == my_rank && range_owning_rank != my_rank) {
      const EntityKeyA entity = i->second.ident;
      const typename MeshA::EntityProc ep(entity, range_owning_rank);
      entities_to_copy.push_back(ep);
    }   
  }
  std::sort(entities_to_copy.begin(), entities_to_copy.end());
  typename MeshA::EntityProcVec::iterator del = std::unique(entities_to_copy.begin(), entities_to_copy.end());
  entities_to_copy.resize(std::distance(entities_to_copy.begin(), del));
 
  mesh.copy_entities(entities_to_copy, transfer_name);

  EntityKeyMap entity_key_map;
  for (typename EntityProcRelationVec::const_iterator i=RangeToDomain.begin(); i!=end; ++i) {
    const unsigned range_owning_rank = i->first.proc;
    if (range_owning_rank == my_rank) {
      const EntityKeyB range_entity  = i->first.ident;
      const EntityKeyA domain_entity = i->second.ident;
      std::pair<EntityKeyB,EntityKeyA> key_map(range_entity, domain_entity);
      entity_key_map.insert(key_map);
    }   
  }   
  return entity_key_map;
}

template <class INTERPOLATE>  void GeometricTransfer<INTERPOLATE>::coarse_search
(EntityProcRelationVec   &RangeToDomain,
 const EntityKeySetA     &keys_to_matcha,
 const EntityKeySetB     &keys_to_matchb,
 const MeshA             &mesha,
 const MeshB             &meshb,
 const double            radius) {
  
  std::vector<BoundingBoxB> range_vector;
  std::vector<BoundingBoxA> domain_vector;

  boundingboxes(domain_vector, keys_to_matcha, mesha, radius);
  boundingboxes(range_vector,  keys_to_matchb, meshb, radius);

  search::FactoryOrder order;
  order.m_communicator = mesha.comm();

  // Slightly confusing: coarse_search documentation has domain->range
  // relations sorted by domain key.  We want range->domain type relations
  // sorted on range key. It might appear we have the arguments revered
  // in coarse_search call, but really, this is what we want.
  search::coarse_search(RangeToDomain, domain_vector, range_vector, order);
}

}
}


#endif

