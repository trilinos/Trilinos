/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICTRANSFER_HPP_
#define STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICTRANSFER_HPP_

#include "stk_search/CoarseSearch.hpp"
#include "stk_search/DistanceComparison.hpp"
#include "stk_search_util/spmd/GeometricSearch.hpp"
#include <stk_transfer/GeometricTransfer.hpp>
#include <stk_transfer/GeometricTransferImpl.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
//#include <stk_transfer_util/spmd/GeometricTransferDispatch.hpp>
#include "stk_util/util/ReportHandler.hpp"
#include <type_traits>
#include <string>
#include <typeinfo>
#include <utility>
#include <iostream>
#include <type_traits>

namespace stk {
namespace transfer {
namespace spmd {

template <typename INTERPOLATE>
class GeometricTransfer : public stk::transfer::GeometricTransfer<INTERPOLATE>
{
 public:
  using SENDMESH              = typename INTERPOLATE::MeshA;
  using RECVMESH              = typename INTERPOLATE::MeshB;
  using RecvEntityKey         = typename INTERPOLATE::EntityKeyB;
  using RecvEntityKeyVector   = typename std::vector<RecvEntityKey>;
  using RecvBoundingBox       = typename INTERPOLATE::MeshB::BoundingBox;

  using MeshA                 =  typename INTERPOLATE::MeshA;
  using MeshB                 =  typename INTERPOLATE::MeshB;
  using EntityKeyA            =  typename INTERPOLATE::EntityKeyA;
  using EntityKeyB            =  typename INTERPOLATE::EntityKeyB;
  using EntityKeyMap          =  typename INTERPOLATE::EntityKeyMap;

  using EntityProcA           =  typename INTERPOLATE::EntityProcA;
  using EntityProcB           =  typename INTERPOLATE::EntityProcB;

  using EntityProcRelation    =  typename INTERPOLATE::EntityProcRelation;
  using EntityProcRelationVec =  typename INTERPOLATE::EntityProcRelationVec;

  using BoundingBoxA          =  typename MeshA::BoundingBox;
  using BoundingBoxB          =  typename MeshB::BoundingBox;

  GeometricTransfer(std::shared_ptr<SENDMESH> sendMesh, std::shared_ptr<RECVMESH> recvMesh, const std::string& name,
                    stk::ParallelMachine pm, const double expansionFactor, const double expansionSum,
                    const stk::search::SearchMethod searchMethod = stk::search::KDTREE)
    : stk::transfer::GeometricTransfer<INTERPOLATE>(sendMesh, recvMesh, name, pm, expansionFactor, searchMethod)
    , m_expansionSum(expansionSum)
    , m_initialized(false)
  {
    create_search_object(sendMesh, recvMesh, name, pm, expansionFactor, expansionSum, searchMethod);
  }

  GeometricTransfer(std::shared_ptr<SENDMESH> sendMesh, std::shared_ptr<RECVMESH> recvMesh, const std::string& name,
                       const double expansionFactor, const double expansionSum,
                       const stk::search::SearchMethod searchMethod = stk::search::KDTREE)
    : stk::transfer::GeometricTransfer<INTERPOLATE>(sendMesh, recvMesh, name, expansionFactor, searchMethod)
    , m_expansionSum(expansionSum)
    , m_initialized(false)
  {
    create_search_object(sendMesh, recvMesh, name, sendMesh->comm(), expansionFactor, expansionSum, searchMethod);
  }

  virtual ~GeometricTransfer() {}
  virtual void coarse_search() override;
  virtual void local_search() override;
  virtual void apply() override;

  virtual const std::shared_ptr<SENDMESH> send_mesh() const {return this->m_mesha;}
  virtual const std::shared_ptr<RECVMESH> recv_mesh() const {return this->m_meshb;}

  virtual void register_inspector(const std::string& fileName, const RecvEntityKeyVector& rangeEntities);
  virtual void inspect_user_defined_entities() { this->m_stkSearch->inspect_user_defined_entities(); }

  virtual void initialize_meshes() {
    this->m_mesha->initialize();
    this->m_meshb->initialize();
  }

  virtual void use_centroid_for_geometric_proximity(const bool option)
  {
    this->m_stkSearch->set_use_centroid_for_geometric_proximity(option);
  }
  virtual bool use_centroid_for_geometric_proximity() const
  {
    return this->m_stkSearch->get_use_centroid_for_geometric_proximity();
  }

  virtual void print_search_warnings(const bool option) { this->m_stkSearch->set_print_search_warnings(option); }
  virtual bool print_search_warnings() const { return this->m_stkSearch->get_print_search_warnings(); }

  virtual void closest_bounding_box_using_nearest_node(const bool option)
  {
    this->m_stkSearch->set_closest_bounding_box_using_nearest_node(option);
  }
  virtual bool closest_bounding_box_using_nearest_node() const
  {
    return this->m_stkSearch->get_closest_bounding_box_using_nearest_node();
  }

  virtual bool is_initialized() const { return m_initialized; }
  virtual void set_initialized(bool flag) { m_initialized = flag; }

  virtual EntityProcRelationVec& get_range_to_domain() { return this->m_stkSearch->get_range_to_domain(); }
  virtual const EntityProcRelationVec& get_range_to_domain() const { return this->m_stkSearch->get_range_to_domain(); }

  virtual void set_expansion_factor(const double expansionFactor)
  {
    this->m_stkSearch->set_expansion_factor(expansionFactor);
  }
  virtual double get_expansion_factor() const
  {
    return this->m_stkSearch->get_expansion_factor();
  }

  virtual void set_expansion_padding(const double expansionPadding)
  {
    m_expansionSum = expansionPadding;
    this->m_stkSearch->set_expansion_padding(expansionPadding);
  }
  virtual double get_expansion_padding() const
  {
    return this->m_stkSearch->get_expansion_padding();
  }

  virtual void set_search_method(const stk::search::SearchMethod searchMethod)
  {
    this->m_stkSearch->set_search_method(searchMethod);
  }
  virtual stk::search::SearchMethod get_search_method() const
  {
    return this->m_stkSearch->get_search_method();
  }

 protected:
  double m_expansionSum;
  bool m_initialized;

  std::vector<RecvBoundingBox> m_unpairedRecvEntities;

  std::shared_ptr<stk::search::spmd::GeometricSearch<SENDMESH, RECVMESH>> m_stkSearch;

  virtual void create_search_object(std::shared_ptr<SENDMESH>& sendMesh, std::shared_ptr<RECVMESH>& recvMesh, const std::string& name,
                                    stk::ParallelMachine pm, const double expansionFactor, const double expansionSum,
                                    const stk::search::SearchMethod searchMethod);
};

template <class INTERPOLATE>
void GeometricTransfer<INTERPOLATE>::create_search_object(std::shared_ptr<SENDMESH>& sendMesh, std::shared_ptr<RECVMESH>& recvMesh,
                                                          const std::string& name, stk::ParallelMachine pm,
                                                          const double expansionFactor, const double expansionSum,
                                                          const stk::search::SearchMethod searchMethod)
{
  m_stkSearch = std::make_shared<stk::search::spmd::GeometricSearch<SENDMESH,RECVMESH>>(sendMesh, recvMesh, name, pm,
                                                                                        expansionFactor, expansionSum, searchMethod);
}

template <class INTERPOLATE>
void GeometricTransfer<INTERPOLATE>::coarse_search()
{
  if(!this->m_has_parallel_machine) // in some cases we needed delayed construction since bulk data might not be set at
                                    // transfer construction
  {
    this->m_parallel_machine = this->m_mesha->comm();
    this->m_has_parallel_machine = true;
  }
  this->m_stkSearch->coarse_search();
}

template <class INTERPOLATE>
void GeometricTransfer<INTERPOLATE>::local_search()
{
  this->m_stkSearch->local_search();
}

template <class INTERPOLATE>
void GeometricTransfer<INTERPOLATE>::apply()
{
  this->m_mesha->update_values();
  INTERPOLATE::apply(*this->m_meshb, *this->m_mesha,
                     this->m_stkSearch->get_range_to_domain(),
                     this->m_stkSearch->get_search_filter_result());
  this->m_meshb->update_values();
}

template <class INTERPOLATE>
void GeometricTransfer<INTERPOLATE>::register_inspector(const std::string& fileName, const RecvEntityKeyVector& rangeEntities)
{
  this->m_stkSearch->register_inspector(fileName, rangeEntities);
}

} // namespace spmd
} // namespace transfer
} // namespace stk

#endif /* STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICTRANSFER_HPP_ */
