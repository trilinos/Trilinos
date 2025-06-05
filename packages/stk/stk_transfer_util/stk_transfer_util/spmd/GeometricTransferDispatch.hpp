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

#ifndef STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICTRANSFERDISPATCH_HPP_
#define STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICTRANSFERDISPATCH_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_transfer/TransferBase.hpp>
#include <stk_transfer_util/spmd/ElementRecvMesh.hpp>
#include <stk_transfer_util/spmd/ElementSendMesh.hpp>
#include <stk_transfer_util/spmd/NodeRecvMesh.hpp>
#include <stk_transfer_util/spmd/NodeSendMesh.hpp>
#include <stk_transfer_util/spmd/GeometricInterp.hpp>
#include <stk_transfer_util/spmd/GeometricTransfer.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <memory>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {
namespace spmd {

template <typename INTERPOLATE>
class GeometricTransfer;


class GeometricTransferDispatchBase {
 public:
  GeometricTransferDispatchBase() {};

  virtual ~GeometricTransferDispatchBase() { }

  virtual void inspect_user_defined_entities() = 0;

  virtual void destroy_ghosting() = 0;

  virtual void ghost_from_elements() = 0;

  virtual void initialize_meshes()  = 0;

  virtual void use_centroid_for_geometric_proximity(const bool option) = 0;
  virtual bool use_centroid_for_geometric_proximity() const = 0;

  virtual void print_search_warnings(const bool option) = 0;
  virtual bool print_search_warnings() const = 0;

  virtual void closest_bounding_box_using_nearest_node(const bool option) = 0;
  virtual bool closest_bounding_box_using_nearest_node() const = 0;

  virtual void set_initialized(bool flag) = 0;
  virtual bool is_initialized() const = 0;

  virtual void set_expansion_factor(const double expansionFactor) = 0;
  virtual double get_expansion_factor() const = 0;

  virtual void set_expansion_padding(const double expansionPadding) = 0;
  virtual double get_expansion_padding() const = 0;

  virtual void set_search_method(const stk::search::SearchMethod searchMethod) = 0;
  virtual stk::search::SearchMethod get_search_method() const = 0;

};

template <typename TRANSFER, typename INTERPOLATE, bool ENFORCE_STK_SPMD_MESH = false>
class GeometricTransferDispatch : public GeometricTransferDispatchBase {
 public:
  using SEND                  = typename INTERPOLATE::MeshA;
  using RECV                  = typename INTERPOLATE::MeshB;

  using SendEntityKey         = typename INTERPOLATE::EntityKeyA;
  using SendEntityKeyVector   = typename std::vector<SendEntityKey>;
  using SendBoundingBox       = typename INTERPOLATE::MeshA::BoundingBox;

  using RecvEntityKey         = typename INTERPOLATE::EntityKeyB;
  using RecvEntityKeyVector   = typename std::vector<RecvEntityKey>;
  using RecvBoundingBox       = typename INTERPOLATE::MeshB::BoundingBox;

  GeometricTransferDispatch(std::shared_ptr<stk::transfer::TransferBase> transfer)
  : m_transfer(std::dynamic_pointer_cast<TRANSFER>(transfer))
  {
    if(!transfer || !m_transfer) {
      STK_ThrowRequireMsg(false,"Input transfer object is not convertible to TRANSFER template");
    }

    if(ENFORCE_STK_SPMD_MESH) {
      static_assert(std::is_base_of<stk::transfer::spmd::GeometricTransfer<INTERPOLATE>, TRANSFER>::value ,
                    "TRANSFER must be a derived class of stk::transfer::spmd/GeometricTransfer.");

      static_assert(std::is_base_of<stk::transfer::spmd::GeometricInterp<SEND, RECV>, INTERPOLATE>::value ,
                    "INTERPOLATE must be a derived class of stk::transfer::spmd/GeometricInterp.");
    }
  }

  GeometricTransferDispatch() {};
  virtual ~GeometricTransferDispatch() { }

  virtual void inspect_user_defined_entities() override
  {
    m_transfer->inspect_user_defined_entities();
  }

  virtual void destroy_ghosting() override
  {
    auto mesha = m_transfer->send_mesh();
    mesha->destroy_ghosting();
  }

  virtual void ghost_from_elements() override
  {
    typename SEND::EntityProcVec entity_keys;
    m_transfer->determine_entities_to_copy(entity_keys);

    auto mesha = m_transfer->send_mesh();
    mesha->update_ghosting(entity_keys);
  }

  virtual bool is_initialized() const override
  {
    return m_transfer->is_initialized();
  }

  virtual void set_initialized(bool flag) override
  {
    m_transfer->set_initialized(flag);
  }

  virtual bool use_centroid_for_geometric_proximity() const override
  {
    return m_transfer->use_centroid_for_geometric_proximity();
  }

  virtual void use_centroid_for_geometric_proximity(bool flag) override
  {
    m_transfer->use_centroid_for_geometric_proximity(flag);
  }

  virtual bool closest_bounding_box_using_nearest_node() const override
  {
    return m_transfer->closest_bounding_box_using_nearest_node();
  }

  virtual void closest_bounding_box_using_nearest_node(bool flag) override
  {
    m_transfer->closest_bounding_box_using_nearest_node(flag);
  }

  virtual void initialize_meshes() override
  {
    m_transfer->initialize_meshes();
  }

  virtual void print_search_warnings(const bool option) override
  {
    m_transfer->print_search_warnings(option);
  }

  virtual bool print_search_warnings() const override
  {
    return m_transfer->print_search_warnings();
  }

  virtual void set_expansion_factor(const double expansionFactor) override
  {
    m_transfer->set_expansion_factor(expansionFactor);
  }

  virtual double get_expansion_factor() const override
  {
    return m_transfer->get_expansion_factor();
  }

  virtual void set_expansion_padding(const double expansionPadding) override
  {
    m_transfer->set_expansion_padding(expansionPadding);
  }

  virtual double get_expansion_padding() const override
  {
    return m_transfer->get_expansion_padding();
  }

  virtual void set_search_method(const stk::search::SearchMethod searchMethod) override
  {
    m_transfer->set_search_method(searchMethod);
  }

  virtual stk::search::SearchMethod get_search_method() const override
  {
    return m_transfer->get_search_method();
  }

  void register_inspector(const std::string& fileName,
                          const std::vector<RecvEntityKey>& rangeEntities) /*override*/
  {
    m_transfer->register_inspector(fileName, rangeEntities);
  }

 private:
  std::shared_ptr<TRANSFER> m_transfer;
};

} // namespace spmd
} // namespace transfer
} // namespace stk

#endif /* STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICTRANSFERDISPATCH_HPP_ */
