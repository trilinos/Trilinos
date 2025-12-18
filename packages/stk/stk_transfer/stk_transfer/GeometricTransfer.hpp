// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef  STK_GEOMETRICTRANSFER_HPP
#define  STK_GEOMETRICTRANSFER_HPP

#include <algorithm>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <stk_util/util/StaticAssert.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/SearchMethod.hpp>
#include <stk_transfer/GeometricTransferImpl.hpp>
#include <stk_transfer/TransferBase.hpp>

namespace stk {
namespace transfer {

template <class INTERPOLATE> class GeometricTransfer : public TransferBase {

public :

  typedef INTERPOLATE                                     InterpolateClass;
  typedef typename InterpolateClass::MeshA                MeshA;
  typedef typename InterpolateClass::MeshB                MeshB;
  typedef typename InterpolateClass::EntityKeyA           EntityKeyA;
  typedef typename InterpolateClass::EntityKeyB           EntityKeyB;
  typedef typename InterpolateClass::EntityKeyMap         EntityKeyMap;

  typedef typename InterpolateClass::EntityProcA          EntityProcA;
  typedef typename InterpolateClass::EntityProcB          EntityProcB;

  typedef typename InterpolateClass::EntityProcRelation       EntityProcRelation;
  typedef typename InterpolateClass::EntityProcRelationVec    EntityProcRelationVec;

  typedef typename std::set     <EntityKeyA>              EntityKeySetA;
  typedef typename std::set     <EntityKeyB>              EntityKeySetB;
  typedef typename MeshA::BoundingBox                     BoundingBoxA;
  typedef typename MeshB::BoundingBox                     BoundingBoxB;


  enum {Dimension = 3};

  GeometricTransfer(std::shared_ptr<MeshA> &mesha,
                    std::shared_ptr<MeshB> &meshb,
                    const std::string &name,
                    stk::ParallelMachine pm,
                    const double expansion_factor = 1.5,
                    const stk::search::SearchMethod search_method = stk::search::KDTREE)
  : m_mesha(mesha),
  m_meshb(meshb),
  m_name (name),
  m_has_parallel_machine(true),
  m_parallel_machine(pm),
  m_expansion_factor(expansion_factor),
  m_search_method(search_method)
  {
  }

  GeometricTransfer(std::shared_ptr<MeshA> &mesha,
                    std::shared_ptr<MeshB> &meshb,
                    const std::string &name,
                    const double expansion_factor = 1.5,
                    const stk::search::SearchMethod search_method = stk::search::KDTREE)
  : m_mesha(mesha),
  m_meshb(meshb),
  m_name (name),
  m_has_parallel_machine(false),
  m_expansion_factor(expansion_factor),
  m_search_method(search_method)
  {
  }

  virtual ~GeometricTransfer(){};
  virtual void coarse_search();
  virtual void communication();
  virtual void local_search();
  virtual void apply();

  void determine_entities_to_copy(typename MeshA::EntityProcVec &entities_to_copy) const;
  const std::shared_ptr<MeshA> meshA() const {return m_mesha;}
  const std::shared_ptr<MeshB> meshB() const {return m_meshb;}

protected :
  std::shared_ptr<MeshA>               m_mesha;
  std::shared_ptr<MeshB>               m_meshb;

  void copy_domain_to_range_processors();
  void localize_entity_key_map();
  void repeat_search_if_needed();

  const std::string     m_name;
  bool m_has_parallel_machine;
  stk::ParallelMachine m_parallel_machine;
  const double          m_expansion_factor;
  const stk::search::SearchMethod m_search_method;

  EntityProcRelationVec m_global_range_to_domain;
  EntityKeyMap          m_local_range_to_domain;
};

template <class INTERPOLATE>
void GeometricTransfer<INTERPOLATE>::coarse_search()
{
  if (!m_has_parallel_machine)  //in some cases we needed delayed construction since bulk data might not be set at transfer construction
  {
    m_parallel_machine = m_mesha->comm();
    m_has_parallel_machine = true;
  }

  impl::coarse_search_impl<INTERPOLATE>(m_global_range_to_domain,
                m_parallel_machine,
                m_mesha.get(),
                m_meshb.get(),
                m_search_method,
                m_expansion_factor);
}

template <class INTERPOLATE>
void GeometricTransfer<INTERPOLATE>::communication()
{
  const unsigned p_size = parallel_machine_size(m_parallel_machine);
  if (1 < p_size) {
    copy_domain_to_range_processors();
  }
}

template <class INTERPOLATE>
void GeometricTransfer<INTERPOLATE>::local_search()
{
  localize_entity_key_map();
  INTERPOLATE::filter_to_nearest(m_local_range_to_domain, *m_mesha, *m_meshb);
}

template <class INTERPOLATE>
void GeometricTransfer<INTERPOLATE>::apply()
{
  m_mesha->update_values();
  INTERPOLATE::apply(*m_meshb, *m_mesha, m_local_range_to_domain);
  m_meshb->update_values();
}

template <class INTERPOLATE> void GeometricTransfer<INTERPOLATE>::determine_entities_to_copy(
                       typename MeshA::EntityProcVec   &entities_to_copy) const {
  entities_to_copy.clear();

  const unsigned my_rank = parallel_machine_rank(m_parallel_machine);

  const typename EntityProcRelationVec::const_iterator end=m_global_range_to_domain.end();
  for (typename EntityProcRelationVec::const_iterator i=m_global_range_to_domain.begin(); i!=end; ++i) {
    const unsigned            domain_owning_rank = i->second.proc();
    const unsigned             range_owning_rank = i->first.proc();
    if (domain_owning_rank == my_rank && range_owning_rank != my_rank) {
      const EntityKeyA entity = i->second.id();
      const typename MeshA::EntityProc ep(entity, range_owning_rank);
      entities_to_copy.push_back(ep);
    }
  }
  std::sort(entities_to_copy.begin(), entities_to_copy.end());
  typename MeshA::EntityProcVec::iterator del = std::unique(entities_to_copy.begin(), entities_to_copy.end());
  entities_to_copy.resize(std::distance(entities_to_copy.begin(), del));
}


template <class INTERPOLATE> void
GeometricTransfer<INTERPOLATE>::localize_entity_key_map()
{
  impl::localize_entity_key_map<INTERPOLATE>(*m_meshb, m_global_range_to_domain, m_local_range_to_domain);
}

template <class INTERPOLATE> void
GeometricTransfer<INTERPOLATE>::copy_domain_to_range_processors()  {

  typename MeshA::EntityProcVec entities_to_copy ;

  determine_entities_to_copy(entities_to_copy);
  impl::call_copy_entities(*m_mesha, entities_to_copy, m_name);
}

}
}


#endif

