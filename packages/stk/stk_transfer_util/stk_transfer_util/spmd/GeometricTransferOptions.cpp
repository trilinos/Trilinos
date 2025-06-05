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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_search_util/MeshUtility.hpp"
#include "stk_transfer_util/spmd/GeometricTransfer.hpp"
#include "stk_transfer_util/spmd/GeometricTransferDispatch.hpp"
#include "stk_transfer_util/spmd/GeometricTransferOptions.hpp"
#include "stk_transfer_util/spmd/GeometricTransferUtils.hpp"
#include <stk_transfer_util/spmd/ElementRecvMesh.hpp>
#include <stk_transfer_util/spmd/ElementSendMesh.hpp>
#include <stk_transfer_util/spmd/NodeRecvMesh.hpp>
#include <stk_transfer_util/spmd/NodeSendMesh.hpp>
#include "stk_search_util/ObjectCoordinates.hpp"
#include <stk_mesh/base/Entity.hpp>         // for Entity
#include <stk_mesh/base/FieldBase.hpp>      // for FieldBase
#include "stk_mesh/base/MetaData.hpp"       // for MetaData
#include "stk_mesh/base/Selector.hpp"       // for Selector
#include "stk_mesh/base/Bucket.hpp"         // for Bucket
#include <stk_mesh/base/BulkData.hpp>       // for BulkData
#include "stk_mesh/base/CompositeRank.hpp"                 // for CompositeRank
#include "stk_topology/topology.hpp"        // for topology, operator<<, top...
#include "stk_util/environment/RuntimeWarning.hpp"         // for RuntimeWar...
#include "stk_util/util/ReportHandler.hpp"  // for eval_test_condition, STK_...
#include "stk_util/util/SortAndUnique.hpp"

#include <algorithm>                        // for max, sort, fill
#include <cstddef>                          // for size_t
#include <memory>                           // for shared_ptr, __shared_ptr_...
#include <type_traits>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {
namespace spmd {

#define ThrowOrReturnOnError(throwOnError, expr, message, errorReturnValue)  \
  if (throwOnError) {                                                        \
    STK_ThrowRequireMsg(expr, message);                                      \
  } else if(!(expr)) {                                                       \
    stk::RuntimeWarning() << message;                                        \
    return errorReturnValue;                                                 \
  }

#define ThrowOrWarnOnError(throwOnError, expr, message)                      \
  if (throwOnError) {                                                        \
    STK_ThrowRequireMsg(expr, message);                                      \
  } else if(!(expr)) {                                                       \
    stk::RuntimeWarning() << message;                                        \
  }


namespace impl {

static constexpr double __default_parametric_tolerance = 0.00001;
static constexpr double __default_mesh_bounding_box_scale_factor = 1.0e-9;
static constexpr double __default_search_expansion_factor = 0.1;

void print_error(const std::string& function, const std::string& errorMessage, const std::string& transferName)
{
  stk::outputP0() << function << " -> " << errorMessage << " for transfer named: " << transferName << std::endl;
}

void add_field_spec(const FieldSpec& fieldSpec, std::vector<FieldSpec>& fieldSpecs)
{
  fieldSpecs.push_back(fieldSpec);
}

void add_field_spec(const std::string& fieldName, std::vector<FieldSpec>& fieldSpecs)
{
  FieldSpec fieldSpec(fieldName);
  fieldSpecs.push_back(fieldSpec);
}

void add_field_spec(const std::string& fieldName,
                    stk::mesh::FieldState fieldState,
                    std::vector<FieldSpec>& fieldSpecs)
{
  FieldSpec fieldSpec(fieldName, fieldState);
  fieldSpecs.push_back(fieldSpec);
}

void add_field_spec(const std::string& fieldName,
                    stk::mesh::FieldState fieldState,
                    unsigned int fieldIndex,
                    std::vector<FieldSpec>& fieldSpecs)
{
  FieldSpec fieldSpec(fieldName, fieldState, fieldIndex);
  fieldSpecs.push_back(fieldSpec);
}

void fill_parts(const stk::mesh::MetaData& meta,
                const std::vector<std::string>& partNames,
                const std::string& transferName,
                stk::mesh::PartVector& parts)
{
  parts.clear();

  for(const auto& partName : partNames) {
    stk::mesh::Part* meshPart = meta.get_part(partName);
    if(nullptr == meshPart) {
      print_error("fill_parts()", "Could not find part named: " + partName, transferName);
    } else {
      parts.push_back(meshPart);
    }
  }
}

bool make_send_mesh_type(TransferOptionHelper* helper)
{
  return helper->make_send_mesh_type();
}
bool make_recv_mesh_type(TransferOptionHelper* helper)
{
  return helper->make_recv_mesh_type();
}
bool make_coordinate_field_name(TransferOptionHelper* helper)
{
  return helper->make_coordinate_field_name();
}
bool make_coordinate_field_state(TransferOptionHelper* helper)
{
  return helper->make_coordinate_field_state();
}
bool make_coordinate_field(TransferOptionHelper* helper)
{
  return helper->make_coordinate_field();
}
bool make_interpolation_type(TransferOptionHelper* helper)
{
  return helper->make_interpolation_type();
}
bool make_extrapolate_option(TransferOptionHelper* helper)
{
  return helper->make_extrapolate_option();
}
bool make_patch_recovery_type(TransferOptionHelper* helper)
{
  return helper->make_patch_recovery_type();
}
bool make_point_evaluator(TransferOptionHelper* helper)
{
  return helper->make_point_evaluator();
}
bool make_external_point_handler(TransferOptionHelper* helper)
{
  return helper->make_external_point_handler();
}
bool make_master_element_provider(TransferOptionHelper* helper)
{
  return helper->make_master_element_provider();
}
bool make_parametric_coordinates_finder(TransferOptionHelper* helper)
{
  return helper->make_parametric_coordinates_finder();
}
bool make_field_interpolator(TransferOptionHelper* helper)
{
  return helper->make_field_interpolator();
}
bool make_parametric_tolerance(TransferOptionHelper* helper)
{
  return helper->make_parametric_tolerance();
}
bool make_geometric_tolerance(TransferOptionHelper* helper)
{
  return helper->make_geometric_tolerance();
}
bool make_search_expansion_factor(TransferOptionHelper* helper)
{
  return helper->make_search_expansion_factor();
}
bool make_search_expansion_padding(TransferOptionHelper* helper)
{
  return helper->make_search_expansion_padding();
}
bool make_active_selector(TransferOptionHelper* helper)
{
  return helper->make_active_selector();
}
bool make_default_part(TransferOptionHelper* helper)
{
  return helper->make_default_part();
}
bool make_node_send_mesh(TransferOptionHelper* helper)
{
  return helper->make_node_send_mesh();
}
bool make_element_send_mesh(TransferOptionHelper* helper)
{
  return helper->make_element_send_mesh();
}
bool make_node_recv_mesh(TransferOptionHelper* helper)
{
  return helper->make_node_recv_mesh();
}
bool make_element_recv_mesh(TransferOptionHelper* helper)
{
  return helper->make_element_recv_mesh();
}
bool make_parallel_machine(TransferOptionHelper* helper)
{
  return helper->make_parallel_machine();
}

class PartChecker {
public:
  PartChecker(const GeometricTransferOptions& transfer, MPI_Comm comm)
    : m_transferOptions(transfer),
      m_comm(comm),
      m_hasBadParts(false),
      m_skipChecks(skip_checks())
  {
    const stk::transfer::spmd::GeometricSendTransferOptions& sendOptions = m_transferOptions.get_send_options();
    const stk::transfer::spmd::GeometricRecvTransferOptions& recvOptions = m_transferOptions.get_recv_options();

    impl::fill_parts(sendOptions.get_bulk().mesh_meta_data(), sendOptions.get_part_names(), m_transferOptions.name(), m_sendParts);
    impl::fill_parts(recvOptions.get_bulk().mesh_meta_data(), recvOptions.get_part_names(), m_transferOptions.name(), m_recvParts);
  }

  void process_fields(const stk::mesh::FieldBase& sendField, const stk::mesh::FieldBase& recvField)
  {
    if (m_skipChecks) return;

    check_part_registration("send", m_sendParts, sendField);
    check_part_registration("recv", m_recvParts, recvField);
  }

  void verify_field_part_definitions() const
  {
    if (m_skipChecks) return;

    STK_ThrowRequireMsg(!m_hasBadParts, m_badPartStream.str());
  }

private:
  bool skip_checks()
  {
    return (m_transferOptions.check_interpolation_type() &&
            m_transferOptions.get_interpolation_type() != InterpolationType::PATCH_RECOVERY);
  }

  void check_part_registration(std::string partType, const stk::mesh::PartVector& parts, const stk::mesh::FieldBase& field)
  {
    const stk::mesh::PartVector badParts = parts_without_field(parts, field);

    if (stk::is_true_on_all_procs(m_comm, badParts.empty())) return;

    m_hasBadParts = true;
    m_badPartStream << "Transfer block " << m_transferOptions.name() << ": ";
    m_badPartStream << "the field " << field.name();
    m_badPartStream << " is not defined on " << partType << " parts";
    for (stk::mesh::Part* part : badParts) {
      m_badPartStream << " " << part->name();
    }
    m_badPartStream << '\n' << '\n';
  }

  stk::mesh::PartVector parts_without_field(const stk::mesh::PartVector& parts, const stk::mesh::FieldBase& field)
  {
    stk::mesh::PartVector badParts;

    for (stk::mesh::Part* part : parts) {
      if (!field_exists_on_all_part_buckets(field, part)) {
        badParts.push_back(part);
      }
    }

    return badParts;
  }

  bool field_exists_on_all_part_buckets(const stk::mesh::FieldBase& field, stk::mesh::Part* part)
  {
    const stk::mesh::Selector fieldSelector(field);
    const stk::mesh::Selector partSelector(*part);
    const stk::mesh::BucketVector& buckets = partSelector.get_buckets(field.entity_rank());

    for (stk::mesh::Bucket* bucket : buckets) {
      if (!fieldSelector(bucket)) return false;
    }

    return true;
  }

  const GeometricTransferOptions& m_transferOptions;

  MPI_Comm m_comm;
  stk::mesh::PartVector m_sendParts;
  stk::mesh::PartVector m_recvParts;

  std::ostringstream m_badPartStream;
  bool m_hasBadParts;
  bool m_skipChecks;
};

bool is_recv_node_mesh(const GeometricRecvTransferOptions& recvOptions)
{
  bool recvMeshTypeSet = recvOptions.has_option(TransferOptions::RECV_MESH_TYPE);
  return recvMeshTypeSet ? recvOptions.get_mesh_type() == RecvMeshType::NODE : false;
}

bool is_recv_gauss_point_mesh(const GeometricRecvTransferOptions& recvOptions)
{
  bool recvMeshTypeSet = recvOptions.has_option(TransferOptions::RECV_MESH_TYPE);
  return recvMeshTypeSet ? (recvOptions.get_mesh_type() == RecvMeshType::EDGE_GAUSS_POINT ||
                            recvOptions.get_mesh_type() == RecvMeshType::FACE_GAUSS_POINT ||
                            recvOptions.get_mesh_type() == RecvMeshType::ELEMENT_GAUSS_POINT)
                         : false;
}

bool is_recv_centroid_mesh(const GeometricRecvTransferOptions& recvOptions)
{
  bool recvMeshTypeSet = recvOptions.has_option(TransferOptions::RECV_MESH_TYPE);
  return recvMeshTypeSet ? (recvOptions.get_mesh_type() == RecvMeshType::EDGE_CENTROID ||
                            recvOptions.get_mesh_type() == RecvMeshType::FACE_CENTROID ||
                            recvOptions.get_mesh_type() == RecvMeshType::ELEMENT_CENTROID)
                         : false;
}

bool is_send_node_mesh(const GeometricSendTransferOptions& sendOptions)
{
  bool sendMeshTypeSet = sendOptions.has_option(TransferOptions::SEND_MESH_TYPE);
  return sendMeshTypeSet ? sendOptions.get_mesh_type() == SendMeshType::NODE : false;
}

bool is_send_element_mesh(const GeometricSendTransferOptions& sendOptions)
{
  bool sendMeshTypeSet = sendOptions.has_option(TransferOptions::SEND_MESH_TYPE);
  return sendMeshTypeSet ? (sendOptions.get_mesh_type() == SendMeshType::ELEMENT ||
                            sendOptions.get_mesh_type() == SendMeshType::FACE    ||
                            sendOptions.get_mesh_type() == SendMeshType::EDGE)
                         : false;
}

}

std::ostream & operator<<(std::ostream &out, SendMeshType type)
{
  switch (type)
  {
  case SendMeshType::NODE:    out << "NODE"; break;
  case SendMeshType::EDGE:    out << "EDGE"; break;
  case SendMeshType::FACE:    out << "FACE"; break;
  case SendMeshType::ELEMENT: out << "ELEMENT"; break;
  case SendMeshType::INVALID: out << "INVALID"; break;
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, RecvMeshType type)
{
  switch (type)
  {
  case RecvMeshType::NODE:                out << "NODE"; break;
  case RecvMeshType::EDGE_CENTROID:       out << "EDGE_CENTROID"; break;
  case RecvMeshType::EDGE_GAUSS_POINT:    out << "EDGE_GAUSS_POINT"; break;
  case RecvMeshType::FACE_CENTROID:       out << "FACE_CENTROID"; break;
  case RecvMeshType::FACE_GAUSS_POINT:    out << "FACE_GAUSS_POINT"; break;
  case RecvMeshType::ELEMENT_CENTROID:    out << "ELEMENT_CENTROID"; break;
  case RecvMeshType::ELEMENT_GAUSS_POINT: out << "ELEMENT_GAUSS_POINT"; break;
  case RecvMeshType::INVALID:             out << "INVALID"; break;
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, InterpolationType type)
{
  switch (type)
  {
  case InterpolationType::COPY:           out << "COPY"; break;
  case InterpolationType::SUM:            out << "SUM"; break;
  case InterpolationType::MASTER_ELEMENT: out << "MASTER_ELEMENT"; break;
  case InterpolationType::PATCH_RECOVERY: out << "PATCH_RECOVERY"; break;
  case InterpolationType::INVALID:        out << "INVALID"; break;
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, TransferType type)
{
  switch (type)
  {
  case TransferType::NODE_TO_NODE:       out << "NODE_TO_NODE"; break;
  case TransferType::NODE_TO_ELEMENT:    out << "NODE_TO_ELEMENT"; break;
  case TransferType::ELEMENT_TO_NODE:    out << "ELEMENT_TO_NODE"; break;
  case TransferType::ELEMENT_TO_ELEMENT: out << "ELEMENT_TO_ELEMENT"; break;
  case TransferType::INVALID:            out << "INVALID"; break;
  }
  return out;
}

std::ostream & operator<<(std::ostream &out, TransferOptions option)
{
  switch (option)
  {
  case TransferOptions::NONE:                         out << "NONE"; break;
  case TransferOptions::SEND_MESH_TYPE:               out << "SEND_MESH_TYPE"; break;
  case TransferOptions::RECV_MESH_TYPE:               out << "RECV_MESH_TYPE"; break;
  case TransferOptions::COORDINATE_FIELD_NAME:        out << "COORDINATE_FIELD_NAME"; break;
  case TransferOptions::COORDINATE_FIELD_STATE:       out << "COORDINATE_FIELD_STATE"; break;
  case TransferOptions::COORDINATE_FIELD:             out << "COORDINATE_FIELD"; break;
  case TransferOptions::INTERPOLATION_TYPE:           out << "INTERPOLATION_TYPE"; break;
  case TransferOptions::EXTRAPOLATION_TYPE:           out << "EXTRAPOLATION_TYPE"; break;
  case TransferOptions::PATCH_RECOVERY_TYPE:          out << "PATCH_RECOVERY_TYPE"; break;
  case TransferOptions::POINT_EVALUATOR:              out << "POINT_EVALUATOR"; break;
  case TransferOptions::EXTERNAL_POINT_HANDLER:       out << "EXTERNAL_POINT_HANDLER"; break;
  case TransferOptions::MASTER_ELEMENT_PROVIDER:      out << "MASTER_ELEMENT_PROVIDER"; break;
  case TransferOptions::PARAMETRIC_COORDINATE_FINDER: out << "PARAMETRIC_COORDINATE_FINDER"; break;
  case TransferOptions::FIELD_INTERPOLATOR:           out << "FIELD_INTERPOLATOR"; break;
  case TransferOptions::PARAMETRIC_TOLERANCE:         out << "PARAMETRIC_TOLERANCE"; break;
  case TransferOptions::GEOMETRIC_TOLERANCE:          out << "GEOMETRIC_TOLERANCE"; break;
  case TransferOptions::BOX_EXPANSION_FACTOR:         out << "BOX_EXPANSION_FACTOR"; break;
  case TransferOptions::BOX_EXPANSION_PADDING:        out << "BOX_EXPANSION_PADDING"; break;
  case TransferOptions::ACTIVE_SELECTOR:              out << "ACTIVE_SELECTOR"; break;
  case TransferOptions::DEFAULT_PART:                 out << "DEFAULT_PART"; break;
  case TransferOptions::NODE_SEND_MESH:               out << "NODE_SEND_MESH"; break;
  case TransferOptions::ELEMENT_SEND_MESH:            out << "ELEMENT_SEND_MESH"; break;
  case TransferOptions::NODE_RECV_MESH:               out << "NODE_RECV_MESH"; break;
  case TransferOptions::ELEMENT_RECV_MESH:            out << "ELEMENT_RECV_MESH"; break;
  case TransferOptions::PARALLEL_MACHINE:             out << "PARALLEL_MACHINE"; break;
  default:                                            out << "TRANSFER_OPTION_" << static_cast<unsigned>(option); break;
  }
  return out;
}

void TransferOptionHelper::register_prerequisite_map_entry(TransferOptions option, DispatchFunction dispatchFunction, const std::string& optionString)
{
  auto iter = m_prerequisiteMap.find(option);
  STK_ThrowRequireMsg(iter == m_prerequisiteMap.end(), "Option: " << option << " has already been registered");

  m_prerequisiteMap[option] = std::make_pair(dispatchFunction, optionString);
}

void TransferOptionHelper::initialize_prerequisite_map()
{
  if(m_mapInitialized) return;

  register_prerequisite_map_entry(TransferOptions::SEND_MESH_TYPE              , impl::make_send_mesh_type               , "send mesh type");
  register_prerequisite_map_entry(TransferOptions::RECV_MESH_TYPE              , impl::make_recv_mesh_type               , "recv mesh type");
  register_prerequisite_map_entry(TransferOptions::COORDINATE_FIELD_NAME       , impl::make_coordinate_field_name        , "coordinate field name");
  register_prerequisite_map_entry(TransferOptions::COORDINATE_FIELD_STATE      , impl::make_coordinate_field_state       , "coordinate field state");
  register_prerequisite_map_entry(TransferOptions::COORDINATE_FIELD            , impl::make_coordinate_field             , "coordinate field");
  register_prerequisite_map_entry(TransferOptions::INTERPOLATION_TYPE          , impl::make_interpolation_type           , "interpolation type");
  register_prerequisite_map_entry(TransferOptions::EXTRAPOLATION_TYPE          , impl::make_extrapolate_option           , "extrapolate option");
  register_prerequisite_map_entry(TransferOptions::PATCH_RECOVERY_TYPE         , impl::make_patch_recovery_type          , "patch recovery type");
  register_prerequisite_map_entry(TransferOptions::POINT_EVALUATOR             , impl::make_point_evaluator              , "point evaluator");
  register_prerequisite_map_entry(TransferOptions::EXTERNAL_POINT_HANDLER      , impl::make_external_point_handler       , "external point handler");
  register_prerequisite_map_entry(TransferOptions::MASTER_ELEMENT_PROVIDER     , impl::make_master_element_provider      , "master element provider");
  register_prerequisite_map_entry(TransferOptions::PARAMETRIC_COORDINATE_FINDER, impl::make_parametric_coordinates_finder, "parametric coordinates finder");
  register_prerequisite_map_entry(TransferOptions::FIELD_INTERPOLATOR          , impl::make_field_interpolator           , "field interpolator");
  register_prerequisite_map_entry(TransferOptions::PARAMETRIC_TOLERANCE        , impl::make_parametric_tolerance         , "parametric_tolerance");
  register_prerequisite_map_entry(TransferOptions::GEOMETRIC_TOLERANCE         , impl::make_geometric_tolerance          , "geometric tolerance");
  register_prerequisite_map_entry(TransferOptions::BOX_EXPANSION_FACTOR        , impl::make_search_expansion_factor      , "search expansion factor");
  register_prerequisite_map_entry(TransferOptions::BOX_EXPANSION_PADDING       , impl::make_search_expansion_padding     , "search expansion padding");
  register_prerequisite_map_entry(TransferOptions::ACTIVE_SELECTOR             , impl::make_active_selector              , "active selector");
  register_prerequisite_map_entry(TransferOptions::DEFAULT_PART                , impl::make_default_part                 , "default part");
  register_prerequisite_map_entry(TransferOptions::NODE_SEND_MESH              , impl::make_node_send_mesh               , "node send mesh");
  register_prerequisite_map_entry(TransferOptions::ELEMENT_SEND_MESH           , impl::make_element_send_mesh            , "element send mesh");
  register_prerequisite_map_entry(TransferOptions::NODE_RECV_MESH              , impl::make_node_recv_mesh               , "node recv mesh");
  register_prerequisite_map_entry(TransferOptions::ELEMENT_RECV_MESH           , impl::make_element_recv_mesh            , "element recv mesh");
  register_prerequisite_map_entry(TransferOptions::PARALLEL_MACHINE            , impl::make_parallel_machine             , "parallel machine");

  m_mapInitialized = true;
}

bool TransferOptionHelper::make_prerequisite(TransferOptions option, const std::string& preamble)
{
  STK_ThrowRequireMsg(m_mapInitialized, "Pre-requesites map has not been initialized");

  bool status{false};

  auto iter = m_prerequisiteMap.find(option);
  if(iter != m_prerequisiteMap.end()) {
    auto dispatchFunction = iter->second.first;
    const std::string& optionString = iter->second.second;

    status = (*dispatchFunction)(this);
    if(!status && m_verboseError) {
      impl::print_error(m_prefix + ":" + preamble, "Could not create pre-requisite " + optionString,  m_transferName);
    }
  } else {
    std::ostringstream oss;
    oss << "Could not find option: " << option << " in pre-requisite map";
    impl::print_error(m_prefix + ":make_prerequisite()", oss.str(), m_transferName);
  }

  return status;
}

bool TransferOptionHelper::check_prerequisite(TransferOptions option, const std::string& preamble)
{
  STK_ThrowRequireMsg(m_mapInitialized, "Pre-requesites map has not been initialized");

  bool status{false};

  auto iter = m_prerequisiteMap.find(option);
  if(iter != m_prerequisiteMap.end()) {
    const std::string& optionString = iter->second.second;

    status = has_option(option);
    if(!status && m_verboseError) {
      impl::print_error(m_prefix + ":" + preamble, optionString + " has not been set",  m_transferName);
    }
  } else {
    std::ostringstream oss;
    oss << "Could not find option: " << option << " in pre-requisite map";
    impl::print_error(m_prefix + ":" + "check_prerequisite()", oss.str(), m_transferName);
  }

  return status;
}

void TransferOptionHelper::set_send_mesh_type(SendMeshType meshType)
{
  m_sendMeshType = meshType;
  set_option(TransferOptions::SEND_MESH_TYPE);
}
SendMeshType TransferOptionHelper::get_send_mesh_type() const { return m_sendMeshType; }
bool TransferOptionHelper::check_send_mesh_type() const { return has_option(TransferOptions::SEND_MESH_TYPE); }
bool TransferOptionHelper::make_send_mesh_type()
{
  return check_prerequisite(TransferOptions::SEND_MESH_TYPE, "make_mesh_type()");
}

void TransferOptionHelper::set_recv_mesh_type(RecvMeshType meshType)
{
  m_recvMeshType = meshType;
  set_option(TransferOptions::RECV_MESH_TYPE);
}
RecvMeshType TransferOptionHelper::get_recv_mesh_type() const { return m_recvMeshType; }
bool TransferOptionHelper::check_recv_mesh_type() const { return has_option(TransferOptions::RECV_MESH_TYPE); }
bool TransferOptionHelper::make_recv_mesh_type()
{
  return check_prerequisite(TransferOptions::RECV_MESH_TYPE, "make_mesh_type()");
}

void TransferOptionHelper::set_coordinate_field_name(const std::string& coordinateFieldName)
{
  m_coordinateFieldName = coordinateFieldName;
  set_option(TransferOptions::COORDINATE_FIELD_NAME);
}
const std::string& TransferOptionHelper::get_coordinate_field_name() const { return m_coordinateFieldName; }
bool TransferOptionHelper::check_coordinate_field_name() const { return has_option(TransferOptions::COORDINATE_FIELD_NAME); }
bool TransferOptionHelper::make_coordinate_field_name()
{
  if(has_option(TransferOptions::COORDINATE_FIELD_NAME)) return true;
  set_coordinate_field_name(m_bulk.mesh_meta_data().coordinate_field_name());
  return has_option(TransferOptions::COORDINATE_FIELD_NAME);
}

void TransferOptionHelper::set_coordinate_field_state(stk::mesh::FieldState coordinateFieldState)
{
  m_coordinateFieldState = coordinateFieldState;
  set_option(TransferOptions::COORDINATE_FIELD_STATE);
}
stk::mesh::FieldState TransferOptionHelper::get_coordinate_field_state() const { return m_coordinateFieldState; }
bool TransferOptionHelper::check_coordinate_field_state() const { return has_option(TransferOptions::COORDINATE_FIELD_STATE); }
bool TransferOptionHelper::make_coordinate_field_state()
{
  if(has_option(TransferOptions::COORDINATE_FIELD_STATE)) return true;
  const stk::mesh::MetaData& meta = m_bulk.mesh_meta_data();
  const stk::mesh::FieldBase* coordinates = meta.get_field(stk::topology::NODE_RANK, meta.coordinate_field_name());
  if(nullptr != coordinates) {
    set_coordinate_field_state(coordinates->state());
  }
  return has_option(TransferOptions::COORDINATE_FIELD_STATE);
}

void TransferOptionHelper::set_coordinate_field(const stk::mesh::FieldBase* coordinateField)
{
  if(nullptr != coordinateField) {
    m_coordinateField = coordinateField;
    set_option(TransferOptions::COORDINATE_FIELD);
  }
}
const stk::mesh::FieldBase* TransferOptionHelper::get_coordinate_field() const { return m_coordinateField; }
bool TransferOptionHelper::check_coordinate_field() const { return has_option(TransferOptions::COORDINATE_FIELD); }
bool TransferOptionHelper::make_coordinate_field()
{
  static std::string preamble{"make_coordinate_field()"};
  if(has_option(TransferOptions::COORDINATE_FIELD)) return true;

  if(!make_prerequisite(TransferOptions::COORDINATE_FIELD_STATE, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::COORDINATE_FIELD_NAME, preamble)) { return false; }

  stk::mesh::FieldBase* coords = m_bulk.mesh_meta_data().get_field(stk::topology::NODE_RANK, m_coordinateFieldName);
  STK_ThrowRequireMsg(coords, "NULL nodal send coordinate field named: " << m_coordinateFieldName << " with state: " << m_coordinateFieldState);
  set_coordinate_field(coords->field_state(m_coordinateFieldState));

  return has_option(TransferOptions::COORDINATE_FIELD);
}

void TransferOptionHelper::add_field_spec(const FieldSpec& fieldSpec)
{
  impl::add_field_spec(fieldSpec, m_fieldSpecs);
}
void TransferOptionHelper::add_field_spec(const std::vector<FieldSpec>& fieldSpecs)
{
  for(const auto& fieldSpec : fieldSpecs) {
    impl::add_field_spec(fieldSpec, m_fieldSpecs);
  }
}
void TransferOptionHelper::add_field_spec(const std::string& fieldName)
{
  impl::add_field_spec(fieldName, m_fieldSpecs);
}
void TransferOptionHelper::add_field_spec(const std::string& fieldName, stk::mesh::FieldState fieldState)
{
  impl::add_field_spec(fieldName, fieldState, m_fieldSpecs);
}
void TransferOptionHelper::add_field_spec(const std::string& fieldName,  stk::mesh::FieldState fieldState, unsigned int fieldIndex)
{
  impl::add_field_spec(fieldName, fieldState, fieldIndex, m_fieldSpecs);
}
void TransferOptionHelper::add_part(const std::string& partName)
{
  stk::util::insert_keep_sorted_and_unique(partName, m_partNames);
}
void TransferOptionHelper::add_part(const std::vector<std::string>& partNames)
{
  for(const auto& partName : partNames) {
    stk::util::insert_keep_sorted_and_unique(partName, m_partNames);
  }
}

void TransferOptionHelper::set_interpolation_type(InterpolationType interpolationType)
{
  m_interpolationType = interpolationType;
  set_option(TransferOptions::INTERPOLATION_TYPE);
}
InterpolationType TransferOptionHelper::get_interpolation_type() const { return m_interpolationType; }
bool TransferOptionHelper::check_interpolation_type() const { return has_option(TransferOptions::INTERPOLATION_TYPE); }
bool TransferOptionHelper::make_interpolation_type()
{
  return check_prerequisite(TransferOptions::INTERPOLATION_TYPE, "make_interpolation_type()");
}

void TransferOptionHelper::set_extrapolate_option(stk::search::ObjectOutsideDomainPolicy extrapolateOption)
{
  m_extrapolateOption = extrapolateOption;
  set_option(TransferOptions::EXTRAPOLATION_TYPE);
}
stk::search::ObjectOutsideDomainPolicy TransferOptionHelper::get_extrapolate_option() const { return m_extrapolateOption; }
bool TransferOptionHelper::check_extrapolate_option() const { return has_option(TransferOptions::EXTRAPOLATION_TYPE); }
bool TransferOptionHelper::make_extrapolate_option()
{
  return check_prerequisite(TransferOptions::EXTRAPOLATION_TYPE, "make_extrapolate_option()");
}

void TransferOptionHelper::set_patch_recovery_type(PatchRecoveryEvaluationType patchEvalType)
{
  m_patchEvalType = patchEvalType;
  set_option(TransferOptions::PATCH_RECOVERY_TYPE);
}
PatchRecoveryEvaluationType TransferOptionHelper::get_patch_recovery_type() const { return m_patchEvalType; }
bool TransferOptionHelper::check_patch_recovery_type() const { return has_option(TransferOptions::PATCH_RECOVERY_TYPE); }
bool TransferOptionHelper::make_patch_recovery_type()
{
  return check_prerequisite(TransferOptions::PATCH_RECOVERY_TYPE, "make_patch_recovery_type()");
}

void TransferOptionHelper::set_external_point_handler(std::shared_ptr<stk::search::HandleExternalPointInterface> externalPointHandler)
{
  m_externalPointHandler = externalPointHandler;
  set_option(TransferOptions::EXTERNAL_POINT_HANDLER);
}
std::shared_ptr<stk::search::HandleExternalPointInterface> TransferOptionHelper::get_external_point_handler() const { return m_externalPointHandler; }
bool TransferOptionHelper::check_external_point_handler() const { return has_option(TransferOptions::EXTERNAL_POINT_HANDLER); }
bool TransferOptionHelper::make_external_point_handler()
{
  static std::string preamble{"make_external_point_handler()"};
  if(has_option(TransferOptions::EXTERNAL_POINT_HANDLER)) return true;

  if(!make_prerequisite(TransferOptions::EXTRAPOLATION_TYPE, preamble)) { return false; }

  switch(m_extrapolateOption) {
  case stk::search::ObjectOutsideDomainPolicy::TRUNCATE:
  {
    if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }
    if(!make_prerequisite(TransferOptions::MASTER_ELEMENT_PROVIDER, preamble)) { return false; }
    if(!make_prerequisite(TransferOptions::PARAMETRIC_TOLERANCE, preamble)) { return false; }
    if(!make_prerequisite(TransferOptions::GEOMETRIC_TOLERANCE, preamble)) { return false; }

    set_external_point_handler(std::make_shared<stk::search::MasterElementExternalPointTruncation>(
        m_bulk, m_coordinateField, m_masterElemProvider, m_geometricTolerance));
    break;
  }
  case stk::search::ObjectOutsideDomainPolicy::PROJECT:
  {
    if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }
    if(!make_prerequisite(TransferOptions::MASTER_ELEMENT_PROVIDER, preamble)) { return false; }
    if(!make_prerequisite(TransferOptions::PARAMETRIC_TOLERANCE, preamble)) { return false; }
    if(!make_prerequisite(TransferOptions::GEOMETRIC_TOLERANCE, preamble)) { return false; }

    set_external_point_handler(std::make_shared<stk::search::MasterElementExternalPointProjection>(
        m_bulk, m_coordinateField, m_masterElemProvider, m_parametricTolerance, m_geometricTolerance));
    break;
  }
  default:
  {
    if(!make_prerequisite(TransferOptions::PARAMETRIC_TOLERANCE, preamble)) { return false; }

    set_external_point_handler(std::make_shared<stk::search::ExternalPointNoOpHandler>(m_parametricTolerance));
    break;
  }
  }

  return has_option(TransferOptions::EXTERNAL_POINT_HANDLER);
}

void TransferOptionHelper::set_master_element_provider(std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider)
{
  m_masterElemProvider = masterElemProvider;
  set_option(TransferOptions::MASTER_ELEMENT_PROVIDER);
}
std::shared_ptr<stk::search::MasterElementProviderInterface>
TransferOptionHelper::get_master_element_provider() const { return m_masterElemProvider; }
bool TransferOptionHelper::check_master_element_provider() const { return has_option(TransferOptions::MASTER_ELEMENT_PROVIDER); }
bool TransferOptionHelper::make_master_element_provider()
{
  if(has_option(TransferOptions::MASTER_ELEMENT_PROVIDER)) return true;

  // TODO: No default implementation .. add in Intrepid2 formulation later

  return has_option(TransferOptions::MASTER_ELEMENT_PROVIDER);
}

void TransferOptionHelper::set_parametric_coordinates_finder(std::shared_ptr<stk::search::FindParametricCoordsInterface> parametricCoordsFinder)
{
  m_parametricCoordsFinder = parametricCoordsFinder;
  set_option(TransferOptions::PARAMETRIC_COORDINATE_FINDER);
}
std::shared_ptr<stk::search::FindParametricCoordsInterface>
TransferOptionHelper::get_parametric_coordinates_finder() const { return m_parametricCoordsFinder; }
bool TransferOptionHelper::check_parametric_coordinates_finder() const { return has_option(TransferOptions::PARAMETRIC_COORDINATE_FINDER); }
bool TransferOptionHelper::make_parametric_coordinates_finder()
{
  static std::string preamble{"make_parametric_coordinates_finder()"};
  if(has_option(TransferOptions::PARAMETRIC_COORDINATE_FINDER)) return true;

  if(!make_prerequisite(TransferOptions::SEND_MESH_TYPE, preamble)) { return false; }

  switch(m_sendMeshType) {
  case SendMeshType::NODE:
  {
    if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }

    set_parametric_coordinates_finder(std::make_shared<stk::search::NodeParametricCoordsFinder>(m_bulk, m_coordinateField));
    break;
  }
  case SendMeshType::EDGE:
  case SendMeshType::FACE:
  case SendMeshType::ELEMENT:
  {
    if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }
    if(!make_prerequisite(TransferOptions::MASTER_ELEMENT_PROVIDER, preamble)) { return false; }
    if(!make_prerequisite(TransferOptions::PARAMETRIC_TOLERANCE, preamble)) { return false; }

    set_parametric_coordinates_finder(std::make_shared<stk::search::MasterElementParametricCoordsFinder>(
        m_bulk, m_coordinateField, m_masterElemProvider, m_parametricTolerance));
    break;
  }
  case SendMeshType::INVALID:
  default:
    break;
  }

  return has_option(TransferOptions::PARAMETRIC_COORDINATE_FINDER);
}

void TransferOptionHelper::set_field_interpolator(std::shared_ptr<stk::transfer::InterpolateFieldsInterface> fieldInterpolator)
{
  m_fieldInterpolator = fieldInterpolator;
  set_option(TransferOptions::FIELD_INTERPOLATOR);
}
std::shared_ptr<stk::transfer::InterpolateFieldsInterface>
TransferOptionHelper::get_field_interpolator() const { return m_fieldInterpolator; }
bool TransferOptionHelper::check_field_interpolator() const { return has_option(TransferOptions::FIELD_INTERPOLATOR); }
bool TransferOptionHelper::make_field_interpolator()
{
  static std::string preamble{"make_field_interpolator()"};
  if(has_option(TransferOptions::FIELD_INTERPOLATOR)) return true;

  if(!make_prerequisite(TransferOptions::INTERPOLATION_TYPE, preamble)) { return false; }

  switch(m_interpolationType) {
  case InterpolationType::COPY:
  {
    set_field_interpolator(std::make_shared<stk::transfer::CopyFieldInterpolator>(m_bulk));
    break;
  }
  case InterpolationType::SUM:
  {
    set_field_interpolator(std::make_shared<stk::transfer::SumFieldInterpolator>(m_bulk));
    break;
  }
  case InterpolationType::MASTER_ELEMENT:
  {
    if(!make_prerequisite(TransferOptions::MASTER_ELEMENT_PROVIDER, preamble)) { return false; }

    set_field_interpolator(std::make_shared<MasterElementFieldInterpolator>(m_bulk, m_masterElemProvider));
    break;
  }
  case InterpolationType::PATCH_RECOVERY:
  {
    if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }
    if(!make_prerequisite(TransferOptions::PATCH_RECOVERY_TYPE, preamble)) { return false; }

    const stk::mesh::Part *meshPart = m_defaultPart;
    stk::mesh::PartVector parts = get_parts(m_bulk.mesh_meta_data(), m_partNames, m_transferName);
    stk::mesh::Selector selector(stk::mesh::selectUnion(parts));
    stk::mesh::Selector* selectorPtr = (parts.size() > 0) ? &selector : nullptr;
    set_field_interpolator(std::make_shared<PatchRecoveryFieldInterpolator>(
        m_bulk, m_coordinateField, meshPart, selectorPtr, m_patchEvalType));
    break;
  }
  case InterpolationType::INVALID:
  default:
    break;
  }

  return has_option(TransferOptions::FIELD_INTERPOLATOR);
}

void TransferOptionHelper::set_point_evaluator(std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator)
{
  m_pointEvaluator = pointEvaluator;
  set_option(TransferOptions::POINT_EVALUATOR);
}
std::shared_ptr<stk::search::PointEvaluatorInterface> TransferOptionHelper::get_point_evaluator() const { return m_pointEvaluator; }
bool TransferOptionHelper::check_point_evaluator() const { return has_option(TransferOptions::POINT_EVALUATOR); }
bool TransferOptionHelper::make_point_evaluator()
{
  static std::string preamble{"make_point_evaluator()"};
  if(has_option(TransferOptions::POINT_EVALUATOR)) return true;

  if(!make_prerequisite(TransferOptions::RECV_MESH_TYPE, preamble)) { return false; }

  switch(m_recvMeshType) {
  case RecvMeshType::EDGE_CENTROID:
  case RecvMeshType::FACE_CENTROID:
  case RecvMeshType::ELEMENT_CENTROID:
  {
    if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }

    set_point_evaluator(std::make_shared<stk::search::CentroidEvaluator>(m_bulk, m_coordinateField));
    break;
  }
  case RecvMeshType::EDGE_GAUSS_POINT:
  case RecvMeshType::FACE_GAUSS_POINT:
  case RecvMeshType::ELEMENT_GAUSS_POINT:
  {
    if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }
    if(!make_prerequisite(TransferOptions::MASTER_ELEMENT_PROVIDER, preamble)) { return false; }

    set_point_evaluator(std::make_shared<stk::search::MasterElementGaussPointEvaluator>(
        m_bulk, m_coordinateField, m_masterElemProvider));
    break;
  }
  default:
    break;
  }

  return has_option(TransferOptions::POINT_EVALUATOR);
}

void TransferOptionHelper::set_parametric_tolerance(double parametricTolerance)
{
  m_parametricTolerance = parametricTolerance;
  set_option(TransferOptions::PARAMETRIC_TOLERANCE);
}
double TransferOptionHelper::get_parametric_tolerance() const { return m_parametricTolerance; }
bool TransferOptionHelper::check_parametric_tolerance() const { return has_option(TransferOptions::PARAMETRIC_TOLERANCE); }
bool TransferOptionHelper::make_parametric_tolerance()
{
  if(has_option(TransferOptions::PARAMETRIC_TOLERANCE)) return true;
  set_parametric_tolerance(impl::__default_parametric_tolerance);
  return has_option(TransferOptions::PARAMETRIC_TOLERANCE);
}

void TransferOptionHelper::set_geometric_tolerance(double geometricTolerance)
{
  m_geometricTolerance = geometricTolerance;
  set_option(TransferOptions::GEOMETRIC_TOLERANCE);
}
double TransferOptionHelper::get_geometric_tolerance() const { return m_geometricTolerance; }
bool TransferOptionHelper::check_geometric_tolerance() const { return has_option(TransferOptions::GEOMETRIC_TOLERANCE); }
bool TransferOptionHelper::make_geometric_tolerance()
{
  static std::string preamble{"get_geometric_tolerance()"};
  if(has_option(TransferOptions::GEOMETRIC_TOLERANCE)) return true;

  if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }

  double scaleFactor = impl::__default_mesh_bounding_box_scale_factor;
  set_geometric_tolerance(mesh_bounding_box_geometric_tolerance(m_bulk, m_coordinateField, scaleFactor));
  return has_option(TransferOptions::GEOMETRIC_TOLERANCE);
}

void TransferOptionHelper::set_search_expansion_factor(double searchExpansionFactor)
{
  m_searchExpansionFactor = searchExpansionFactor;
  set_option(TransferOptions::BOX_EXPANSION_FACTOR);
}
double TransferOptionHelper::get_search_expansion_factor() const { return m_searchExpansionFactor; }
bool TransferOptionHelper::check_search_expansion_factor() const { return has_option(TransferOptions::BOX_EXPANSION_FACTOR); }
bool TransferOptionHelper::make_search_expansion_factor()
{
  if(has_option(TransferOptions::BOX_EXPANSION_FACTOR)) return true;
  set_search_expansion_factor(impl::__default_search_expansion_factor);
  return has_option(TransferOptions::BOX_EXPANSION_FACTOR);
}

void TransferOptionHelper::set_search_expansion_padding(double searchExpansionPadding)
{
  m_searchExpansionPadding = searchExpansionPadding;
  set_option(TransferOptions::BOX_EXPANSION_PADDING);
}
double TransferOptionHelper::get_search_expansion_padding() const { return m_searchExpansionPadding; }
bool TransferOptionHelper::check_search_expansion_padding() const { return has_option(TransferOptions::BOX_EXPANSION_PADDING); }
bool TransferOptionHelper::make_search_expansion_padding()
{
  static std::string preamble{"make_search_expansion_padding()"};
  if(has_option(TransferOptions::BOX_EXPANSION_PADDING)) return true;

  if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }

  double scaleFactor = impl::__default_mesh_bounding_box_scale_factor;
  set_search_expansion_padding(mesh_bounding_box_geometric_tolerance(m_bulk, m_coordinateField, scaleFactor));
  return has_option(TransferOptions::BOX_EXPANSION_PADDING);
}

void TransferOptionHelper::set_active_selector(const stk::mesh::Selector& activeSelector)
{
  m_activeSelector = activeSelector;
  set_option(TransferOptions::ACTIVE_SELECTOR);
}
const stk::mesh::Selector& TransferOptionHelper::get_active_selector() const { return m_activeSelector; }
bool TransferOptionHelper::check_active_selector() const { return has_option(TransferOptions::ACTIVE_SELECTOR); }
bool TransferOptionHelper::make_active_selector()
{
  static std::string preamble{"make_active_selector()"};
  if(has_option(TransferOptions::ACTIVE_SELECTOR)) return true;

  set_active_selector(stk::mesh::Selector().complement());
  return has_option(TransferOptions::ACTIVE_SELECTOR);
}

void TransferOptionHelper::set_default_part(const stk::mesh::Part* defaultPart)
{
  m_defaultPart = defaultPart;
  set_option(TransferOptions::DEFAULT_PART);
}
const stk::mesh::Part* TransferOptionHelper::get_default_part() const { return m_defaultPart; }
bool TransferOptionHelper::check_default_part() const { return has_option(TransferOptions::DEFAULT_PART); }
bool TransferOptionHelper::make_default_part()
{
  return has_option(TransferOptions::DEFAULT_PART);
}

void TransferOptionHelper::set_node_send_mesh(std::shared_ptr<stk::transfer::spmd::NodeSendMesh> nodeSendMesh)
{
  m_nodeSendMesh = nodeSendMesh;
  set_option(TransferOptions::NODE_SEND_MESH);
}
std::shared_ptr<stk::transfer::spmd::NodeSendMesh> TransferOptionHelper::get_node_send_mesh() const { return m_nodeSendMesh; }
bool TransferOptionHelper::check_node_send_mesh() const { return has_option(TransferOptions::NODE_SEND_MESH); }
bool TransferOptionHelper::make_node_send_mesh()
{
  static std::string preamble{"make_node_send_mesh()"};
  if(has_option(TransferOptions::NODE_SEND_MESH)) return true;

  if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::ACTIVE_SELECTOR, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::FIELD_INTERPOLATOR, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::GEOMETRIC_TOLERANCE, preamble)) { return false; }

  const stk::mesh::MetaData& meta = m_bulk.mesh_meta_data();
  stk::mesh::PartVector sendParts = get_parts(meta, m_partNames, m_transferName, m_defaultPart);

  if(!m_fieldInterpolator->fields_are_set()) {
    m_fieldInterpolator->set_fields(m_fieldSpecs);
  }

  set_node_send_mesh(std::make_shared<stk::transfer::spmd::NodeSendMesh>(&m_bulk, m_coordinateField, m_fieldSpecs,
                                                                         sendParts, m_activeSelector, m_bulk.parallel(),
                                                                         m_fieldInterpolator, m_geometricTolerance));

  if(has_option(TransferOptions::EXTRAPOLATION_TYPE)) {
    m_nodeSendMesh->set_extrapolate_option(m_extrapolateOption);
  }

  return has_option(TransferOptions::NODE_SEND_MESH);
}

void TransferOptionHelper::set_element_send_mesh(std::shared_ptr<stk::transfer::spmd::ElementSendMesh> elemSendMesh)
{
  m_elemSendMesh = elemSendMesh;
  set_option(TransferOptions::ELEMENT_SEND_MESH);
}
std::shared_ptr<stk::transfer::spmd::ElementSendMesh> TransferOptionHelper::get_element_send_mesh() const { return m_elemSendMesh; }
bool TransferOptionHelper::check_element_send_mesh() const { return has_option(TransferOptions::ELEMENT_SEND_MESH); }
bool TransferOptionHelper::make_element_send_mesh()
{
  static std::string preamble{"make_element_send_mesh()"};
  if(has_option(TransferOptions::ELEMENT_SEND_MESH)) return true;

  if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::ACTIVE_SELECTOR, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::PARAMETRIC_COORDINATE_FINDER, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::EXTERNAL_POINT_HANDLER, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::FIELD_INTERPOLATOR, preamble)) { return false; }

  if(has_option(TransferOptions::MASTER_ELEMENT_PROVIDER)) {
    if(!make_prerequisite(TransferOptions::MASTER_ELEMENT_PROVIDER, preamble)) { return false; }
  }

  const stk::mesh::MetaData& meta = m_bulk.mesh_meta_data();

  stk::mesh::PartVector sendParts = get_parts(meta, m_partNames, m_transferName, m_defaultPart);
  stk::mesh::EntityRank sendEntityRank = has_option(TransferOptions::ELEMENT_MESH_ENTITY_RANK) ?
                                         m_elementMeshEntityRank : get_transfer_rank(meta, m_partNames, m_transferName);

  if(!m_fieldInterpolator->fields_are_set()) {
    m_fieldInterpolator->set_fields(m_fieldSpecs);
  }

  if(has_option(TransferOptions::MASTER_ELEMENT_PROVIDER)) {
    set_element_send_mesh(std::make_shared<stk::transfer::spmd::ElementSendMesh>(&m_bulk, m_coordinateField, m_fieldSpecs,
                                                                                 sendEntityRank, sendParts, m_activeSelector,
                                                                                 m_bulk.parallel(), m_parametricCoordsFinder,
                                                                                 m_externalPointHandler, m_fieldInterpolator,
                                                                                 m_masterElemProvider));
  } else {
    set_element_send_mesh(std::make_shared<stk::transfer::spmd::ElementSendMesh>(&m_bulk, m_coordinateField, m_fieldSpecs,
                                                                                 sendEntityRank, sendParts, m_activeSelector,
                                                                                 m_bulk.parallel(), m_parametricCoordsFinder,
                                                                                 m_externalPointHandler, m_fieldInterpolator));
  }

  if(has_option(TransferOptions::EXTRAPOLATION_TYPE)) {
    m_elemSendMesh->set_extrapolate_option(m_extrapolateOption);
  }

  return has_option(TransferOptions::ELEMENT_SEND_MESH);
}

void TransferOptionHelper::set_node_recv_mesh(std::shared_ptr<stk::transfer::spmd::NodeRecvMesh> nodeRecvMesh)
{
  m_nodeRecvMesh = nodeRecvMesh;
  set_option(TransferOptions::NODE_RECV_MESH);
}
std::shared_ptr<stk::transfer::spmd::NodeRecvMesh> TransferOptionHelper::get_node_recv_mesh() const { return m_nodeRecvMesh; }
bool TransferOptionHelper::check_node_recv_mesh() const { return has_option(TransferOptions::NODE_RECV_MESH); }
bool TransferOptionHelper::make_node_recv_mesh()
{
  static std::string preamble{"make_node_recv_mesh()"};
  if(has_option(TransferOptions::NODE_RECV_MESH)) return true;

  if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::ACTIVE_SELECTOR, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::PARAMETRIC_TOLERANCE, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::GEOMETRIC_TOLERANCE, preamble)) { return false; }

  const stk::mesh::MetaData& meta = m_bulk.mesh_meta_data();
  stk::mesh::PartVector recvParts = get_parts(meta, m_partNames, m_transferName, m_defaultPart);

  set_node_recv_mesh(std::make_shared<stk::transfer::spmd::NodeRecvMesh>(&m_bulk, m_coordinateField, m_fieldSpecs,
                                                                         recvParts, m_activeSelector, m_bulk.parallel(),
                                                                         m_parametricTolerance, m_geometricTolerance));

  return has_option(TransferOptions::NODE_RECV_MESH);
}

void TransferOptionHelper::set_element_recv_mesh(std::shared_ptr<stk::transfer::spmd::ElementRecvMesh> elemRecvMesh)
{
  m_elemRecvMesh = elemRecvMesh;
  set_option(TransferOptions::ELEMENT_RECV_MESH);
}
std::shared_ptr<stk::transfer::spmd::ElementRecvMesh> TransferOptionHelper::get_element_recv_mesh() const { return m_elemRecvMesh; }
bool TransferOptionHelper::check_element_recv_mesh() const { return has_option(TransferOptions::ELEMENT_RECV_MESH); }
bool TransferOptionHelper::make_element_recv_mesh()
{
  static std::string preamble{"make_element_recv_mesh()"};
  if(has_option(TransferOptions::ELEMENT_RECV_MESH)) return true;

  if(!make_prerequisite(TransferOptions::COORDINATE_FIELD, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::ACTIVE_SELECTOR, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::POINT_EVALUATOR, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::PARAMETRIC_TOLERANCE, preamble)) { return false; }
  if(!make_prerequisite(TransferOptions::GEOMETRIC_TOLERANCE, preamble)) { return false; }

  const stk::mesh::MetaData& meta = m_bulk.mesh_meta_data();

  stk::mesh::PartVector recvParts = get_parts(meta, m_partNames, m_transferName, m_defaultPart);
  stk::mesh::EntityRank recvEntityRank = has_option(TransferOptions::ELEMENT_MESH_ENTITY_RANK) ?
                                         m_elementMeshEntityRank : get_transfer_rank(meta, m_partNames, m_transferName);

  set_element_recv_mesh(std::make_shared<stk::transfer::spmd::ElementRecvMesh>(&m_bulk, m_coordinateField, m_fieldSpecs,
                                                                               recvEntityRank, recvParts, m_activeSelector,
                                                                               m_bulk.parallel(), m_pointEvaluator,
                                                                               m_parametricTolerance, m_geometricTolerance));

  return has_option(TransferOptions::ELEMENT_RECV_MESH);
}

void TransferOptionHelper::set_element_mesh_entity_rank(stk::mesh::EntityRank elementMeshEntityRank)
{
  m_elementMeshEntityRank = elementMeshEntityRank;
  set_option(TransferOptions::ELEMENT_MESH_ENTITY_RANK);
}
stk::mesh::EntityRank TransferOptionHelper::get_element_mesh_entity_rank() const { return m_elementMeshEntityRank; }
bool TransferOptionHelper::check_element_mesh_entity_rank() const { return has_option(TransferOptions::ELEMENT_MESH_ENTITY_RANK); }
bool TransferOptionHelper::make_element_mesh_entity_rank()
{
  return check_prerequisite(TransferOptions::ELEMENT_MESH_ENTITY_RANK, "make_element_mesh_entity_rank()");
}

void TransferOptionHelper::set_parallel_machine(stk::ParallelMachine parallelMachine)
{
  m_parallelMachine = parallelMachine;
  set_option(TransferOptions::PARALLEL_MACHINE);
}
stk::ParallelMachine TransferOptionHelper::get_parallel_machine() const { return m_parallelMachine; }
bool TransferOptionHelper::check_parallel_machine() const { return has_option(TransferOptions::PARALLEL_MACHINE); }
bool TransferOptionHelper::make_parallel_machine()
{
  return has_option(TransferOptions::PARALLEL_MACHINE);
}

void GeometricTransferOptions::check_field_validity() const
{
  impl::PartChecker partChecker(*this, m_sendBulk.parallel());

  const stk::transfer::spmd::GeometricSendTransferOptions& sendOptions = get_send_options();
  const stk::transfer::spmd::GeometricRecvTransferOptions& recvOptions = get_recv_options();

  const std::vector<stk::transfer::FieldSpec>& sendFieldSpecs = m_sendOptions.get_fields();
  const std::vector<stk::transfer::FieldSpec>& recvFieldSpecs = m_recvOptions.get_fields();

  const stk::mesh::MetaData& sendMeta = sendOptions.get_bulk().mesh_meta_data();
  const stk::mesh::MetaData& recvMeta = recvOptions.get_bulk().mesh_meta_data();

  STK_ThrowRequireMsg(sendFieldSpecs.size() == recvFieldSpecs.size(),
      "Number of source fields does not match number of destination fields in transfer named: " << name());

  for(size_t i = 0; i < sendFieldSpecs.size(); i++) {
    const stk::mesh::FieldBase* sendField = stk::mesh::get_field_by_name(sendFieldSpecs[i].name, sendMeta);
    const stk::mesh::FieldBase* recvField = stk::mesh::get_field_by_name(recvFieldSpecs[i].name, recvMeta);

    STK_ThrowRequireMsg(nullptr != sendField,
          "Could not find send field named: " << sendFieldSpecs[i].name << " in transfer named: " << name());
    STK_ThrowRequireMsg(nullptr != recvField,
          "Could not find recv field named: " << recvFieldSpecs[i].name << " in transfer named: " << name());

    partChecker.process_fields(*sendField, *recvField);
  }

  partChecker.verify_field_part_definitions();
}

TransferType get_transfer_type(GeometricTransferOptions& transferOptions, bool throwOnError)
{
  TransferType invalidReturnValue{TransferType::INVALID};

  stk::transfer::spmd::GeometricSendTransferOptions& sendOptions = transferOptions.get_send_options();
  stk::transfer::spmd::GeometricRecvTransferOptions& recvOptions = transferOptions.get_recv_options();

  bool sendMeshTypeSet = sendOptions.has_option(TransferOptions::SEND_MESH_TYPE);
  bool isSendNode = impl::is_send_node_mesh(sendOptions);
  bool isSendElem = impl::is_send_element_mesh(sendOptions);

  ThrowOrReturnOnError(throwOnError, sendMeshTypeSet, "get_transfer_type(): SEND_MESH_TYPE option not set", invalidReturnValue);

  bool recvMeshTypeSet = recvOptions.has_option(TransferOptions::RECV_MESH_TYPE);
  bool isRecvNode      = impl::is_recv_node_mesh(recvOptions);
  bool isRecvGauss     = impl::is_recv_gauss_point_mesh(recvOptions);
  bool isRecvCentroid  = impl::is_recv_centroid_mesh(recvOptions);

  ThrowOrReturnOnError(throwOnError, recvMeshTypeSet, "get_transfer_type(): RECV_MESH_TYPE option not set", invalidReturnValue);

  if(isSendNode && isRecvNode) {
    // Node-node transfer
    return TransferType::NODE_TO_NODE;
  }

  if(isSendNode && (isRecvGauss || isRecvCentroid)) {
    // Node-element transfer
    return TransferType::NODE_TO_ELEMENT;
  }

  if(isSendElem && isRecvNode) {
    // Element-node transfer
    return TransferType::ELEMENT_TO_NODE;
  }

  if(isSendElem && (isRecvGauss || isRecvCentroid)) {
    // Element-element transfer
    return TransferType::ELEMENT_TO_ELEMENT;
  }

  return invalidReturnValue;
}

GeometricTransferPair create_node_to_node_transfer(GeometricTransferOptions& transferOptions, bool throwOnError)
{
  // Node-node transfer
  GeometricTransferPair emptyReturnValue(nullptr, nullptr);

  stk::transfer::spmd::GeometricSendTransferOptions& sendOptions = transferOptions.get_send_options();
  stk::transfer::spmd::GeometricRecvTransferOptions& recvOptions = transferOptions.get_recv_options();

  double boxExpansionFactor = transferOptions.get_search_expansion_factor();
  double boxExpansionPadding = transferOptions.get_search_expansion_padding();

  std::shared_ptr<stk::transfer::spmd::GeometricTransferDispatchBase> dispatch;
  std::shared_ptr<stk::transfer::TransferBase>                        transfer;

  ThrowOrReturnOnError(throwOnError, sendOptions.make_node_send_mesh(),
      "create_node_to_node_transfer(): could not create NODE send mesh", emptyReturnValue);
  ThrowOrReturnOnError(throwOnError, recvOptions.make_node_recv_mesh(),
      "create_node_to_node_transfer(): could not create NODE recv mesh", emptyReturnValue);

  using SPMDInterp   = stk::transfer::spmd::GeometricInterp<stk::transfer::spmd::NodeSendMesh, stk::transfer::spmd::NodeRecvMesh>;
  using SPMDTransfer = stk::transfer::spmd::GeometricTransfer< SPMDInterp >;

  if(transferOptions.has_option(TransferOptions::PARALLEL_MACHINE)) {
    transfer = std::make_shared< SPMDTransfer >(sendOptions.get_node_send_mesh(), recvOptions.get_node_recv_mesh(),
                                                transferOptions.name(), transferOptions.get_parallel_machine(),
                                                boxExpansionFactor, boxExpansionPadding);
  } else {
    transfer = std::make_shared< SPMDTransfer >(sendOptions.get_node_send_mesh(), recvOptions.get_node_recv_mesh(),
                                                transferOptions.name(), boxExpansionFactor, boxExpansionPadding);
  }

  using SPMDDispatch = stk::transfer::spmd::GeometricTransferDispatch<SPMDTransfer, SPMDInterp>;
  dispatch = std::make_shared< SPMDDispatch >(transfer);

  return std::make_pair(transfer, dispatch);
}

GeometricTransferPair create_node_to_element_transfer(GeometricTransferOptions& transferOptions, bool throwOnError)
{
  // Node-element transfer
  GeometricTransferPair emptyReturnValue(nullptr, nullptr);

  stk::transfer::spmd::GeometricSendTransferOptions& sendOptions = transferOptions.get_send_options();
  stk::transfer::spmd::GeometricRecvTransferOptions& recvOptions = transferOptions.get_recv_options();

  double boxExpansionFactor = transferOptions.get_search_expansion_factor();
  double boxExpansionPadding = transferOptions.get_search_expansion_padding();

  std::shared_ptr<stk::transfer::spmd::GeometricTransferDispatchBase> dispatch;
  std::shared_ptr<stk::transfer::TransferBase>                        transfer;

  ThrowOrReturnOnError(throwOnError, sendOptions.make_node_send_mesh(),
      "create_node_to_element_transfer(): could not create NODE send mesh", emptyReturnValue);
  ThrowOrReturnOnError(throwOnError, recvOptions.make_element_recv_mesh(),
      "create_node_to_element_transfer(): could not create ELEMENT recv mesh", emptyReturnValue);

  using SPMDInterp   = stk::transfer::spmd::GeometricInterp<stk::transfer::spmd::NodeSendMesh, stk::transfer::spmd::ElementRecvMesh>;
  using SPMDTransfer = stk::transfer::spmd::GeometricTransfer< SPMDInterp >;

  if(transferOptions.has_option(TransferOptions::PARALLEL_MACHINE)) {
    transfer = std::make_shared< SPMDTransfer >(sendOptions.get_node_send_mesh(), recvOptions.get_element_recv_mesh(),
                                                transferOptions.name(), transferOptions.get_parallel_machine(),
                                                boxExpansionFactor, boxExpansionPadding);
  } else {
    transfer = std::make_shared< SPMDTransfer >(sendOptions.get_node_send_mesh(), recvOptions.get_element_recv_mesh(),
                                                transferOptions.name(), boxExpansionFactor, boxExpansionPadding);
  }

  using SPMDDispatch = stk::transfer::spmd::GeometricTransferDispatch<SPMDTransfer, SPMDInterp>;
  dispatch = std::make_shared< SPMDDispatch >(transfer);

  return std::make_pair(transfer, dispatch);
}

GeometricTransferPair create_element_to_node_transfer(GeometricTransferOptions& transferOptions, bool throwOnError)
{
  // Element-node transfer
  GeometricTransferPair emptyReturnValue(nullptr, nullptr);

  stk::transfer::spmd::GeometricSendTransferOptions& sendOptions = transferOptions.get_send_options();
  stk::transfer::spmd::GeometricRecvTransferOptions& recvOptions = transferOptions.get_recv_options();

  double boxExpansionFactor = transferOptions.get_search_expansion_factor();
  double boxExpansionPadding = transferOptions.get_search_expansion_padding();

  std::shared_ptr<stk::transfer::spmd::GeometricTransferDispatchBase> dispatch;
  std::shared_ptr<stk::transfer::TransferBase>                        transfer;

  ThrowOrReturnOnError(throwOnError, sendOptions.make_element_send_mesh(),
      "create_element_to_node_transfer(): could not create ELEMENT send mesh", emptyReturnValue);
  ThrowOrReturnOnError(throwOnError, recvOptions.make_node_recv_mesh(),
      "create_element_to_node_transfer(): could not create NODE recv mesh", emptyReturnValue);

  using SPMDInterp   = stk::transfer::spmd::GeometricInterp<stk::transfer::spmd::ElementSendMesh, stk::transfer::spmd::NodeRecvMesh>;
  using SPMDTransfer = stk::transfer::spmd::GeometricTransfer< SPMDInterp >;

  if(transferOptions.has_option(TransferOptions::PARALLEL_MACHINE)) {
    transfer = std::make_shared< SPMDTransfer >(sendOptions.get_element_send_mesh(), recvOptions.get_node_recv_mesh(),
                                                transferOptions.name(), transferOptions.get_parallel_machine(),
                                                boxExpansionFactor, boxExpansionPadding);
  } else {
    transfer = std::make_shared< SPMDTransfer >(sendOptions.get_element_send_mesh(), recvOptions.get_node_recv_mesh(),
                                                transferOptions.name(), boxExpansionFactor, boxExpansionPadding);
  }

  using SPMDDispatch = stk::transfer::spmd::GeometricTransferDispatch<SPMDTransfer, SPMDInterp>;
  dispatch = std::make_shared< SPMDDispatch >(transfer);

  return std::make_pair(transfer, dispatch);
}

GeometricTransferPair create_element_to_element_transfer(GeometricTransferOptions& transferOptions, bool throwOnError)
{
  // Element-element transfer
  GeometricTransferPair emptyReturnValue(nullptr, nullptr);

  stk::transfer::spmd::GeometricSendTransferOptions& sendOptions = transferOptions.get_send_options();
  stk::transfer::spmd::GeometricRecvTransferOptions& recvOptions = transferOptions.get_recv_options();

  double boxExpansionFactor = transferOptions.get_search_expansion_factor();
  double boxExpansionPadding = transferOptions.get_search_expansion_padding();

  std::shared_ptr<stk::transfer::spmd::GeometricTransferDispatchBase> dispatch;
  std::shared_ptr<stk::transfer::TransferBase>                        transfer;

  ThrowOrReturnOnError(throwOnError, sendOptions.make_element_send_mesh(),
      "create_element_to_element_transfer(): could not create ELEMENT send mesh", emptyReturnValue);
  ThrowOrReturnOnError(throwOnError, recvOptions.make_element_recv_mesh(),
      "create_element_to_element_transfer(): could not create ELEMENT recv mesh", emptyReturnValue);

  using SPMDInterp   = stk::transfer::spmd::GeometricInterp<stk::transfer::spmd::ElementSendMesh, stk::transfer::spmd::ElementRecvMesh>;
  using SPMDTransfer = stk::transfer::spmd::GeometricTransfer< SPMDInterp >;

  if(transferOptions.has_option(TransferOptions::PARALLEL_MACHINE)) {
    transfer = std::make_shared< SPMDTransfer >(sendOptions.get_element_send_mesh(), recvOptions.get_element_recv_mesh(),
                                                transferOptions.name(), transferOptions.get_parallel_machine(),
                                                boxExpansionFactor, boxExpansionPadding);
  } else {
    transfer = std::make_shared< SPMDTransfer >(sendOptions.get_element_send_mesh(), recvOptions.get_element_recv_mesh(),
                                                transferOptions.name(), boxExpansionFactor, boxExpansionPadding);
  }

  using SPMDDispatch = stk::transfer::spmd::GeometricTransferDispatch<SPMDTransfer, SPMDInterp>;
  dispatch = std::make_shared< SPMDDispatch >(transfer);

  return std::make_pair(transfer, dispatch);
}

GeometricTransferPair create_transfer(GeometricTransferOptions& transferOptions, bool throwOnError)
{
  GeometricTransferPair emptyReturnValue(nullptr, nullptr);

  // Some preliminary checks
  bool hasExpansionFactor = transferOptions.make_search_expansion_factor();
  ThrowOrReturnOnError(throwOnError, hasExpansionFactor,
      "create_transfer(): could not obtain search expansion factor", emptyReturnValue);

  bool hasExpansionPadding = transferOptions.make_search_expansion_padding();
  ThrowOrReturnOnError(throwOnError, hasExpansionPadding,
      "create_transfer(): could not obtain search expansion padding", emptyReturnValue);

  TransferType transferType = get_transfer_type(transferOptions, throwOnError);

  if(transferType == TransferType::NODE_TO_NODE) {
    // Node-node transfer
   return create_node_to_node_transfer(transferOptions, throwOnError);
  }

  if(transferType == TransferType::NODE_TO_ELEMENT) {
    // Node-element transfer
    return create_node_to_element_transfer(transferOptions, throwOnError);
  }

  if(transferType == TransferType::ELEMENT_TO_NODE) {
    // Element-node transfer
    return create_element_to_node_transfer(transferOptions, throwOnError);
  }

  if(transferType == TransferType::ELEMENT_TO_ELEMENT) {
    // Element-element transfer
    return create_element_to_element_transfer(transferOptions, throwOnError);
  }

  return emptyReturnValue;
}

std::shared_ptr<stk::transfer::spmd::GeometricTransferDispatchBase>
create_dispatch(std::shared_ptr<stk::transfer::TransferBase> inputTransfer)
{
  {
    // Check Node-Node
    using SPMDInterp   = stk::transfer::spmd::GeometricInterp<stk::transfer::spmd::NodeSendMesh, stk::transfer::spmd::NodeRecvMesh>;
    using SPMDTransfer = stk::transfer::spmd::GeometricTransfer< SPMDInterp >;
    using SPMDDispatch = stk::transfer::spmd::GeometricTransferDispatch<SPMDTransfer, SPMDInterp>;

    std::shared_ptr<SPMDTransfer> transfer = std::dynamic_pointer_cast<SPMDTransfer>(inputTransfer);

    if(transfer) {
      return std::make_shared< SPMDDispatch >(transfer);
    }
  }

  {
    // Check Node-Element
    using SPMDInterp   = stk::transfer::spmd::GeometricInterp<stk::transfer::spmd::NodeSendMesh, stk::transfer::spmd::ElementRecvMesh>;
    using SPMDTransfer = stk::transfer::spmd::GeometricTransfer< SPMDInterp >;
    using SPMDDispatch = stk::transfer::spmd::GeometricTransferDispatch<SPMDTransfer, SPMDInterp>;

    std::shared_ptr<SPMDTransfer> transfer = std::dynamic_pointer_cast<SPMDTransfer>(inputTransfer);

    if(transfer) {
      return std::make_shared< SPMDDispatch >(transfer);
    }
  }

  {
    // Check Element-Node
    using SPMDInterp   = stk::transfer::spmd::GeometricInterp<stk::transfer::spmd::ElementSendMesh, stk::transfer::spmd::NodeRecvMesh>;
    using SPMDTransfer = stk::transfer::spmd::GeometricTransfer< SPMDInterp >;
    using SPMDDispatch = stk::transfer::spmd::GeometricTransferDispatch<SPMDTransfer, SPMDInterp>;

    std::shared_ptr<SPMDTransfer> transfer = std::dynamic_pointer_cast<SPMDTransfer>(inputTransfer);

    if(transfer) {
      return std::make_shared< SPMDDispatch >(transfer);
    }
  }

  {
    // Check Element-Element
    using SPMDInterp   = stk::transfer::spmd::GeometricInterp<stk::transfer::spmd::ElementSendMesh, stk::transfer::spmd::ElementRecvMesh>;
    using SPMDTransfer = stk::transfer::spmd::GeometricTransfer< SPMDInterp >;
    using SPMDDispatch = stk::transfer::spmd::GeometricTransferDispatch<SPMDTransfer, SPMDInterp>;

    std::shared_ptr<SPMDTransfer> transfer = std::dynamic_pointer_cast<SPMDTransfer>(inputTransfer);

    if(transfer) {
      return std::make_shared< SPMDDispatch >(transfer);
    }
  }

  return std::shared_ptr<stk::transfer::spmd::GeometricTransferDispatchBase>(nullptr);
}

} // namespace spmd
} // namespace transfer
} // namespace stk

