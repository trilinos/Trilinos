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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SPMD_GEOMETRICSEARCHDISPATCH_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SPMD_GEOMETRICSEARCHDISPATCH_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <ostream>   // for std::ostream
#include <limits>    // for std::numeric_limits
#include <algorithm> // for std::sort

#include "stk_search_util/spmd/GeometricSearch.hpp"
#include <stk_search_util/spmd/ElementRecvMesh.hpp>
#include <stk_search_util/spmd/ElementSendMesh.hpp>
#include <stk_search_util/spmd/NodeRecvMesh.hpp>
#include <stk_search_util/spmd/NodeSendMesh.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <stk_util/parallel/OutputStreams.hpp>
#include <stk_util/parallel/CouplingVersions.hpp>
#include "stk_util/util/ReportHandler.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {
namespace spmd {

class GeometricSearchDispatchBase
{
public:
  GeometricSearchDispatchBase() {}
  virtual ~GeometricSearchDispatchBase() {}

  virtual void set_closest_bounding_box_using_nearest_node(bool flag) = 0;
  virtual bool get_closest_bounding_box_using_nearest_node() const = 0;

  virtual void set_use_centroid_for_geometric_proximity(bool flag) = 0;
  virtual bool get_use_centroid_for_geometric_proximity() const = 0;

  virtual void set_do_initial_search_expansion(bool flag) = 0;
  virtual bool get_do_initial_search_expansion() const = 0;

  virtual void set_print_search_warnings(bool flag) = 0;
  virtual bool get_print_search_warnings() const = 0;

  virtual void set_output_stream(std::ostream& out) = 0;
  virtual std::ostream& get_output_stream() = 0;

  virtual void set_fractional_limit_for_objects_outside_domain(double f) = 0;
  virtual double get_fractional_limit_for_objects_outside_domain() const = 0;

  virtual void set_expansion_limit(double expansionLimit) = 0;
  virtual double get_expansion_limit() const = 0;

  virtual void set_expansion_factor(const double expansionFactor) = 0;
  virtual double get_expansion_factor() const = 0;

  virtual void set_expansion_padding(const double expansionPadding) = 0;
  virtual double get_expansion_padding() const = 0;

  virtual void set_search_method(const stk::search::SearchMethod searchMethod) = 0;
  virtual stk::search::SearchMethod get_search_method() const = 0;
};

template <typename SENDMESH, typename RECVMESH, bool ENFORCE_STK_SPMD_MESH = false>
class GeometricSearchDispatch : public GeometricSearchDispatchBase {
 public:
  using SendEntityKey         = typename SENDMESH::EntityKey;
  using SendEntityKeyVector   = typename std::vector<SendEntityKey>;
  using SendBoundingBox       = typename SENDMESH::BoundingBox;

  using RecvEntityKey         = typename RECVMESH::EntityKey;
  using RecvEntityKeyVector   = typename std::vector<RecvEntityKey>;
  using RecvBoundingBox       = typename RECVMESH::BoundingBox;

  using BaseSearchClass       = typename stk::search::spmd::GeometricSearch<SENDMESH, RECVMESH>;

  GeometricSearchDispatch(std::shared_ptr<stk::search::spmd::SearchBase> search)
  : m_search(std::dynamic_pointer_cast<BaseSearchClass>(search))
  {
    if(!search || !m_search) {
      STK_ThrowRequireMsg(false,
          "Input search object is not convertible to stk::search::spmd::GeometricSearch<SENDMESH, RECVMESH> template");
    }

    if(ENFORCE_STK_SPMD_MESH) {
      static_assert(std::is_base_of<stk::search::spmd::NodeSendMesh   , SENDMESH>::value ||
                    std::is_base_of<stk::search::spmd::ElementSendMesh, SENDMESH>::value,
                    "SENDMESH must be a derived class of stk::search::spmd::[Node/Element]SendMesh.");

      static_assert(std::is_base_of<stk::search::spmd::NodeRecvMesh   , RECVMESH>::value ||
                    std::is_base_of<stk::search::spmd::ElementRecvMesh, RECVMESH>::value,
                    "RECVMESH must be a derived class of stk::search::spmd::[Node/Element]RecvMesh.");
    }
  }

  GeometricSearchDispatch() {};
  virtual ~GeometricSearchDispatch() { }

  virtual void set_closest_bounding_box_using_nearest_node(bool flag) override {
    m_search->set_closest_bounding_box_using_nearest_node(flag);
  }
  virtual bool get_closest_bounding_box_using_nearest_node() const override {
    return m_search->get_closest_bounding_box_using_nearest_node();
  }

  virtual void set_use_centroid_for_geometric_proximity(bool flag) override {
    m_search->set_use_centroid_for_geometric_proximity(flag);
  }
  virtual bool get_use_centroid_for_geometric_proximity() const override {
    return m_search->get_use_centroid_for_geometric_proximity();
  }

  virtual void set_do_initial_search_expansion(bool flag) override {
    m_search->set_do_initial_search_expansion(flag);
  }
  virtual bool get_do_initial_search_expansion() const override {
    return m_search->get_do_initial_search_expansion();
  }

  virtual void set_print_search_warnings(bool flag) override {
    m_search->set_print_search_warnings(flag);
  }
  virtual bool get_print_search_warnings() const override {
    return m_search->get_print_search_warnings();
  }

  virtual void set_output_stream(std::ostream& out) override {
    m_search->set_output_stream(out);
  }
  virtual std::ostream& get_output_stream() override {
    return m_search->get_output_stream();
  }

  virtual void set_fractional_limit_for_objects_outside_domain(double f) override {
    m_search->set_fractional_limit_for_objects_outside_domain(f);
  }
  virtual double get_fractional_limit_for_objects_outside_domain() const override {
    return m_search->get_fractional_limit_for_objects_outside_domain();
  }

  virtual void set_expansion_limit(double expansionLimit) override {
    m_search->set_expansion_limit(expansionLimit);
  }
  virtual double get_expansion_limit() const override {
    return m_search->get_expansion_limit();
  }

  virtual void set_expansion_factor(const double expansionFactor) override {
    m_search->set_expansion_factor(expansionFactor);
  }
  virtual double get_expansion_factor() const override {
    return m_search->get_expansion_factor();
  }

  virtual void set_expansion_padding(const double expansionPadding) override {
    m_search->set_expansion_padding(expansionPadding);
  }
  virtual double get_expansion_padding() const override {
    return m_search->get_expansion_padding();
  }

  virtual void set_search_method(const stk::search::SearchMethod searchMethod) override {
    m_search->set_search_method(searchMethod);
  }
  virtual stk::search::SearchMethod get_search_method() const override {
    return m_search->get_search_method();
  }

 private:
  std::shared_ptr<BaseSearchClass> m_search;
};

} // namespace spmd
} // namespace search
} // namespace stk

#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SPMD_GEOMETRICSEARCHDISPATCH_HPP_ */
