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

#ifndef STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICTRANSFEROPTIONS_HPP_
#define STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICTRANSFEROPTIONS_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_mesh/base/FieldState.hpp"                 // for FieldState
#include "stk_mesh/base/FieldBase.hpp"                  // for FieldBase
#include "stk_mesh/base/Selector.hpp"
#include <stk_search/ObjectOutsideDomainPolicy.hpp>
#include "stk_search_util/ExternalPointHandler.hpp"   // for ExternalPoin...
#include "stk_search_util/FindParametricCoordinates.hpp"  // for FindParametr...
#include "stk_search_util/MasterElementProvider.hpp"  // for ProvideMaste...
#include "stk_search_util/MeshUtility.hpp"
#include "stk_search_util/PointEvaluator.hpp"   // for EvaluatePointsInte...
#include "stk_transfer_util/InterpolateFields.hpp"
#include <stk_transfer_util/PointInterpolation.hpp>
#include "stk_util/parallel/Parallel.hpp"
#include <ostream>
#include <memory>                                       // for shared_ptr
#include <string>                                       // for string
#include <utility>                                      // for pair
#include <vector>                                       // for vector
#include <map>
#include <cstdint>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk { namespace transfer { namespace spmd { class NodeSendMesh; } } }
namespace stk { namespace transfer { namespace spmd { class NodeRecvMesh; } } }
namespace stk { namespace transfer { namespace spmd { class ElementSendMesh; } } }
namespace stk { namespace transfer { namespace spmd { class ElementRecvMesh; } } }
namespace stk { namespace transfer { namespace spmd { class GeometricTransferDispatchBase; } } }

namespace stk { namespace transfer { class TransferBase; } }

namespace stk {
namespace transfer {
namespace spmd {

enum class SendMeshType {
  NODE,
  EDGE,
  FACE,
  ELEMENT,
  INVALID
};

std::ostream & operator<<(std::ostream &out, SendMeshType type);

enum class RecvMeshType {
  NODE,
  EDGE_CENTROID,
  EDGE_GAUSS_POINT,
  FACE_CENTROID,
  FACE_GAUSS_POINT,
  ELEMENT_CENTROID,
  ELEMENT_GAUSS_POINT,
  INVALID
};

std::ostream & operator<<(std::ostream &out, RecvMeshType type);

enum class InterpolationType {
  COPY,
  SUM,
  MASTER_ELEMENT,
  PATCH_RECOVERY,
  INVALID
};

std::ostream & operator<<(std::ostream &out, InterpolationType type);

enum class TransferType {
  NODE_TO_NODE,

  NODE_TO_ELEMENT,
  NODE_TO_ELEM = NODE_TO_ELEMENT,

  ELEMENT_TO_NODE,
  ELEM_TO_NODE = ELEMENT_TO_NODE,

  ELEMENT_TO_ELEMENT,
  ELEM_TO_ELEM = ELEMENT_TO_ELEMENT,

  INVALID
};

std::ostream & operator<<(std::ostream &out, TransferType type);

enum TransferOptions {
  NONE                         = 0,
  SEND_MESH_TYPE               = (1L <<  0),
  RECV_MESH_TYPE               = (1L <<  1),
  COORDINATE_FIELD_NAME        = (1L <<  2),
  COORDINATE_FIELD_STATE       = (1L <<  3),
  COORDINATE_FIELD             = (1L <<  4),
  INTERPOLATION_TYPE           = (1L <<  5),
  EXTRAPOLATION_TYPE           = (1L <<  6),
  PATCH_RECOVERY_TYPE          = (1L <<  7),
  POINT_EVALUATOR              = (1L <<  8),
  EXTERNAL_POINT_HANDLER       = (1L <<  9),
  MASTER_ELEMENT_PROVIDER      = (1L << 10),
  PARAMETRIC_COORDINATE_FINDER = (1L << 11),
  FIELD_INTERPOLATOR           = (1L << 12),
  PARAMETRIC_TOLERANCE         = (1L << 13),
  GEOMETRIC_TOLERANCE          = (1L << 14),
  BOX_EXPANSION_FACTOR         = (1L << 15),
  BOX_EXPANSION_PADDING        = (1L << 16),
  ACTIVE_SELECTOR              = (1L << 17),
  DEFAULT_PART                 = (1L << 18),
  NODE_SEND_MESH               = (1L << 19),
  ELEMENT_SEND_MESH            = (1L << 20),
  NODE_RECV_MESH               = (1L << 21),
  ELEMENT_RECV_MESH            = (1L << 22),
  ELEMENT_MESH_ENTITY_RANK     = (1L << 23),
  PARALLEL_MACHINE             = (1L << 24)
};

std::ostream & operator<<(std::ostream &out, TransferOptions option);

using ExternalPointHandlerSharedPtr = std::shared_ptr<stk::search::HandleExternalPointInterface>;
using MasterElemProviderSharedPtr   = std::shared_ptr<stk::search::MasterElementProviderInterface>;
using FindParametricCoordsSharedPtr = std::shared_ptr<stk::search::FindParametricCoordsInterface>;
using InterpolateFieldsSharePtr     = std::shared_ptr<stk::transfer::InterpolateFieldsInterface>;
using PointEvaluatorSharedPtr       = std::shared_ptr<stk::search::PointEvaluatorInterface>;

using NodeSendMeshSharedPtr = std::shared_ptr<stk::transfer::spmd::NodeSendMesh>;
using NodeRecvMeshSharedPtr = std::shared_ptr<stk::transfer::spmd::NodeRecvMesh>;
using ElemSendMeshSharedPtr = std::shared_ptr<stk::transfer::spmd::ElementSendMesh>;
using ElemRecvMeshSharedPtr = std::shared_ptr<stk::transfer::spmd::ElementRecvMesh>;

class TransferOptionHelper;

class TransferOptionHelper {
 public:
  TransferOptionHelper(const std::string& transferName, stk::mesh::BulkData& bulk, const std::string& prefix)
    : m_transferName(transferName)
    , m_bulk(bulk)
    , m_prefix(prefix)
  {
    initialize_prerequisite_map();
  }

  TransferOptionHelper() = delete;
  ~TransferOptionHelper() = default;

  const std::string& name() const { return m_transferName; }
  const stk::mesh::BulkData& get_bulk() const { return m_bulk; }

  void set_option(TransferOptions flag) { m_assignedOptions |= flag; }
  bool has_option(TransferOptions flag) const { return (m_assignedOptions & flag); }

  void set_send_mesh_type(SendMeshType meshType);
  SendMeshType get_send_mesh_type() const;
  bool check_send_mesh_type() const;
  bool make_send_mesh_type();

  void set_recv_mesh_type(RecvMeshType meshType);
  RecvMeshType get_recv_mesh_type() const;
  bool check_recv_mesh_type() const;
  bool make_recv_mesh_type();

  void set_coordinate_field_name(const std::string& coordinateFieldName);
  const std::string& get_coordinate_field_name() const;
  bool check_coordinate_field_name() const;
  bool make_coordinate_field_name();

  void set_coordinate_field_state(stk::mesh::FieldState coordinateFieldState);
  stk::mesh::FieldState get_coordinate_field_state() const;
  bool check_coordinate_field_state() const;
  bool make_coordinate_field_state();

  void set_coordinate_field(const stk::mesh::FieldBase* coordinateField);
  const stk::mesh::FieldBase* get_coordinate_field() const;
  bool check_coordinate_field() const;
  bool make_coordinate_field();

  void add_field_spec(const FieldSpec& fieldSpec);
  void add_field_spec(const std::vector<FieldSpec>& fieldSpecs);
  void add_field_spec(const std::string& fieldName);
  void add_field_spec(const std::string& fieldName, stk::mesh::FieldState fieldState);
  void add_field_spec(const std::string& fieldName,  stk::mesh::FieldState fieldState, unsigned int fieldIndex);

  void add_part(const std::string& partName);
  void add_part(const std::vector<std::string>& partNames);

  void set_interpolation_type(InterpolationType interpolationType);
  InterpolationType get_interpolation_type() const;
  bool check_interpolation_type() const;
  bool make_interpolation_type();

  void set_extrapolate_option(stk::search::ObjectOutsideDomainPolicy extrapolateOption);
  stk::search::ObjectOutsideDomainPolicy get_extrapolate_option() const;
  bool check_extrapolate_option() const;
  bool make_extrapolate_option();

  void set_patch_recovery_type(PatchRecoveryEvaluationType patchEvalType);
  PatchRecoveryEvaluationType get_patch_recovery_type() const;
  bool check_patch_recovery_type() const;
  bool make_patch_recovery_type();

  void set_external_point_handler(std::shared_ptr<stk::search::HandleExternalPointInterface> externalPointHandler);
  std::shared_ptr<stk::search::HandleExternalPointInterface> get_external_point_handler() const;
  bool check_external_point_handler() const;
  bool make_external_point_handler();

  void set_master_element_provider(std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider);
  std::shared_ptr<stk::search::MasterElementProviderInterface> get_master_element_provider() const;
  bool check_master_element_provider() const;
  bool make_master_element_provider();

  void set_parametric_coordinates_finder(std::shared_ptr<stk::search::FindParametricCoordsInterface> parametricCoordsFinder);
  std::shared_ptr<stk::search::FindParametricCoordsInterface> get_parametric_coordinates_finder() const;
  bool check_parametric_coordinates_finder() const;
  bool make_parametric_coordinates_finder();

  void set_point_evaluator(std::shared_ptr<stk::search::PointEvaluatorInterface> pointEvaluator);
  std::shared_ptr<stk::search::PointEvaluatorInterface> get_point_evaluator() const;
  bool check_point_evaluator() const;
  bool make_point_evaluator();

  void set_field_interpolator(std::shared_ptr<stk::transfer::InterpolateFieldsInterface> fieldInterpolator);
  std::shared_ptr<stk::transfer::InterpolateFieldsInterface> get_field_interpolator() const;
  bool check_field_interpolator() const;
  bool make_field_interpolator();

  void set_parametric_tolerance(double parametricTolerance);
  double get_parametric_tolerance() const;
  bool check_parametric_tolerance() const;
  bool make_parametric_tolerance();

  void set_geometric_tolerance(double geometricTolerance);
  double get_geometric_tolerance() const;
  bool check_geometric_tolerance() const;
  bool make_geometric_tolerance();

  void set_search_expansion_factor(double searchExpansionFactor);
  double get_search_expansion_factor() const;
  bool check_search_expansion_factor() const;
  bool make_search_expansion_factor();

  void set_search_expansion_padding(double searchExpansionPadding);
  double get_search_expansion_padding() const;
  bool check_search_expansion_padding() const;
  bool make_search_expansion_padding();

  void set_active_selector(const stk::mesh::Selector& activeSelector);
  const stk::mesh::Selector& get_active_selector() const;
  bool check_active_selector() const;
  bool make_active_selector();

  void set_default_part(const stk::mesh::Part* defaultPart);
  const stk::mesh::Part* get_default_part() const;
  bool check_default_part() const;
  bool make_default_part();

  void set_node_send_mesh(std::shared_ptr<stk::transfer::spmd::NodeSendMesh> nodeSendMesh);
  std::shared_ptr<stk::transfer::spmd::NodeSendMesh> get_node_send_mesh() const;
  bool check_node_send_mesh() const;
  bool make_node_send_mesh();

  void set_element_send_mesh(std::shared_ptr<stk::transfer::spmd::ElementSendMesh> elemSendMesh);
  std::shared_ptr<stk::transfer::spmd::ElementSendMesh> get_element_send_mesh() const;
  bool check_element_send_mesh() const;
  bool make_element_send_mesh();

  void set_node_recv_mesh(std::shared_ptr<stk::transfer::spmd::NodeRecvMesh> nodeRecvMesh);
  std::shared_ptr<stk::transfer::spmd::NodeRecvMesh> get_node_recv_mesh() const;
  bool check_node_recv_mesh() const;
  bool make_node_recv_mesh();

  void set_element_recv_mesh(std::shared_ptr<stk::transfer::spmd::ElementRecvMesh> elemRecvMesh);
  std::shared_ptr<stk::transfer::spmd::ElementRecvMesh> get_element_recv_mesh() const;
  bool check_element_recv_mesh() const;
  bool make_element_recv_mesh();

  void set_element_mesh_entity_rank(stk::mesh::EntityRank elementMeshEntityRank);
  stk::mesh::EntityRank get_element_mesh_entity_rank() const;
  bool check_element_mesh_entity_rank() const;
  bool make_element_mesh_entity_rank();

  void set_parallel_machine(stk::ParallelMachine parallelMachine);
  stk::ParallelMachine get_parallel_machine() const;
  bool check_parallel_machine() const;
  bool make_parallel_machine();

  void set_verbose(bool flag) { m_verboseError = flag; }
  bool get_verbose() const { return m_verboseError; }

  bool make_prerequisite(TransferOptions option, const std::string& preamble);
  bool check_prerequisite(TransferOptions option, const std::string& preamble);

  void initialize_prerequisite_map();

  const std::vector<stk::transfer::FieldSpec>& get_fields() const { return m_fieldSpecs; }
  const std::vector<std::string>& get_part_names() const { return m_partNames; }

 private:
  std::string m_transferName{"UN-INITIALIZED TRANSFER NAME"};

  stk::mesh::BulkData& m_bulk;

  std::string m_prefix;

  // Must be big enough to hold all the bit shifted options
  uint32_t m_assignedOptions{TransferOptions::NONE};

  SendMeshType m_sendMeshType{SendMeshType::INVALID};
  RecvMeshType m_recvMeshType{RecvMeshType::INVALID};

  std::string                 m_coordinateFieldName;
  stk::mesh::FieldState       m_coordinateFieldState{stk::mesh::FieldState::StateNone};
  const stk::mesh::FieldBase* m_coordinateField{nullptr};

  std::vector<stk::transfer::FieldSpec> m_fieldSpecs;
  std::vector<std::string> m_partNames;

  InterpolationType                      m_interpolationType{InterpolationType::INVALID};
  stk::search::ObjectOutsideDomainPolicy m_extrapolateOption{stk::search::ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG};
  PatchRecoveryEvaluationType            m_patchEvalType{PatchRecoveryEvaluationType::UNDEFINED_PATCH_RECOVERY_EVALUATION_TYPE};

  std::shared_ptr<stk::search::HandleExternalPointInterface>   m_externalPointHandler;
  std::shared_ptr<stk::search::MasterElementProviderInterface> m_masterElemProvider;
  std::shared_ptr<stk::search::FindParametricCoordsInterface>  m_parametricCoordsFinder;
  std::shared_ptr<stk::transfer::InterpolateFieldsInterface>   m_fieldInterpolator;
  std::shared_ptr<stk::search::PointEvaluatorInterface>        m_pointEvaluator;

  std::shared_ptr<stk::transfer::spmd::NodeSendMesh>    m_nodeSendMesh;
  std::shared_ptr<stk::transfer::spmd::NodeRecvMesh>    m_nodeRecvMesh;
  std::shared_ptr<stk::transfer::spmd::ElementSendMesh> m_elemSendMesh;
  std::shared_ptr<stk::transfer::spmd::ElementRecvMesh> m_elemRecvMesh;

  double m_parametricTolerance{0.0};
  double m_geometricTolerance{0.0};

  double m_searchExpansionFactor{0.0};
  double m_searchExpansionPadding{0.0};

  stk::mesh::Selector m_activeSelector;

  const stk::mesh::Part* m_defaultPart{nullptr};

  stk::mesh::EntityRank m_elementMeshEntityRank{stk::topology::INVALID_RANK};

  stk::ParallelMachine m_parallelMachine{stk::parallel_machine_world()};

  bool m_verboseError{false};

  using PrerequisiteMap = std::map<TransferOptions, std::pair<bool(*)(TransferOptionHelper*), std::string>>;
  PrerequisiteMap m_prerequisiteMap;

  bool m_mapInitialized{false};

  using DispatchFunction = bool (*)(stk::transfer::spmd::TransferOptionHelper *);
  void register_prerequisite_map_entry(TransferOptions option, DispatchFunction dispatchFunction, const std::string& optionString);
};

class GeometricSendTransferOptions {
 public:
  GeometricSendTransferOptions(const std::string& transferName, stk::mesh::BulkData& sendBulk)
    : m_transferName(transferName)
    , m_bulk(sendBulk)
    , m_helper(transferName, sendBulk, "Send")
  {  }

  GeometricSendTransferOptions() = delete;
  ~GeometricSendTransferOptions() = default;

  const std::string& name() const { return m_transferName; }
  const stk::mesh::BulkData& get_bulk() const { return m_bulk; }

  bool has_option(TransferOptions flag) const { return m_helper.has_option(flag); }

  // Allowable implementations for source mesh
  void set_mesh_type(SendMeshType meshType) { m_helper.set_send_mesh_type(meshType); }
  SendMeshType get_mesh_type() const        { return m_helper.get_send_mesh_type(); }
  bool check_mesh_type() const              { return m_helper.check_send_mesh_type(); }
  bool make_mesh_type()                     { return m_helper.make_send_mesh_type(); }

  void set_coordinate_field_name(const std::string& coordinateFieldName) { m_helper.set_coordinate_field_name(coordinateFieldName); }
  const std::string& get_coordinate_field_name() const                   { return m_helper.get_coordinate_field_name(); }
  bool check_coordinate_field_name() const                               { return m_helper.check_coordinate_field_name(); }
  bool make_coordinate_field_name()                                      { return m_helper.make_coordinate_field_name(); }

  void set_coordinate_field_state(stk::mesh::FieldState coordinateFieldState) { m_helper.set_coordinate_field_state(coordinateFieldState); }
  stk::mesh::FieldState get_coordinate_field_state() const                    { return m_helper.get_coordinate_field_state(); }
  bool check_coordinate_field_state() const                                   { return m_helper.check_coordinate_field_state(); }
  bool make_coordinate_field_state()                                          { return m_helper.make_coordinate_field_state(); }

  void set_coordinate_field(const stk::mesh::FieldBase* coordinateField)  { m_helper.set_coordinate_field(coordinateField); }
  const stk::mesh::FieldBase* get_coordinate_field() const                { return m_helper.get_coordinate_field(); }
  bool check_coordinate_field() const                                     { return m_helper.check_coordinate_field(); }
  bool make_coordinate_field()                                            { return m_helper.make_coordinate_field(); }

  void add_field_spec(FieldSpec& fieldSpec)                                                        { m_helper.add_field_spec(fieldSpec); }
  void add_field_spec(const std::vector<FieldSpec>& fieldSpecs)                                    { m_helper.add_field_spec(fieldSpecs); }
  void add_field_spec(const std::string& fName)                                                    { m_helper.add_field_spec(fName); }
  void add_field_spec(const std::string& fName, stk::mesh::FieldState fState)                      { m_helper.add_field_spec(fName, fState); }
  void add_field_spec(const std::string& fName, stk::mesh::FieldState fState, unsigned int fIndex) { m_helper.add_field_spec(fName, fState, fIndex); }

  void add_part(const std::string& partName)                { m_helper.add_part(partName); }
  void add_part(const std::vector<std::string>& partNames)  { m_helper.add_part(partNames); }

  void set_interpolation_type(InterpolationType interpolationType) { m_helper.set_interpolation_type(interpolationType); }
  InterpolationType get_interpolation_type() const                 { return m_helper.get_interpolation_type(); }
  bool check_interpolation_type() const                            { return m_helper.check_interpolation_type(); }
  bool make_interpolation_type()                                   { return m_helper.make_interpolation_type(); }

  void set_extrapolate_option(stk::search::ObjectOutsideDomainPolicy extrapolateOption) { m_helper.set_extrapolate_option(extrapolateOption); }
  stk::search::ObjectOutsideDomainPolicy get_extrapolate_option() const                 { return m_helper.get_extrapolate_option(); }
  bool check_extrapolate_option() const                                                 { return m_helper.check_extrapolate_option(); }
  bool make_extrapolate_option()                                                        { return m_helper.make_extrapolate_option(); }

  void set_patch_recovery_type(PatchRecoveryEvaluationType patchEvalType) { m_helper.set_patch_recovery_type(patchEvalType); }
  PatchRecoveryEvaluationType get_patch_recovery_type() const             { return m_helper.get_patch_recovery_type(); }
  bool check_patch_recovery_type() const                                  { return m_helper.check_patch_recovery_type(); }
  bool make_patch_recovery_type()                                         { return m_helper.make_patch_recovery_type(); }

  void set_external_point_handler(ExternalPointHandlerSharedPtr externalPointHandler) { m_helper.set_external_point_handler(externalPointHandler); }
  ExternalPointHandlerSharedPtr get_external_point_handler() const                    { return m_helper.get_external_point_handler(); }
  bool check_external_point_handler() const                                           { return m_helper.check_external_point_handler(); }
  bool make_external_point_handler()                                                  { return m_helper.make_external_point_handler(); }

  void set_master_element_provider(MasterElemProviderSharedPtr masterElemProvider) { m_helper.set_master_element_provider(masterElemProvider); }
  MasterElemProviderSharedPtr get_master_element_provider() const                  { return m_helper.get_master_element_provider(); }
  bool check_master_element_provider() const                                       { return m_helper.check_master_element_provider(); }
  bool make_master_element_provider()                                              { return m_helper.make_master_element_provider(); }

  void set_parametric_coordinates_finder(FindParametricCoordsSharedPtr parametricCoordsFinder) { m_helper.set_parametric_coordinates_finder(parametricCoordsFinder); }
  FindParametricCoordsSharedPtr get_parametric_coordinates_finder() const                      { return m_helper.get_parametric_coordinates_finder(); }
  bool check_parametric_coordinates_finder() const                                             { return m_helper.check_parametric_coordinates_finder(); }
  bool make_parametric_coordinates_finder()                                                    { return m_helper.make_parametric_coordinates_finder(); }

  void set_field_interpolator(InterpolateFieldsSharePtr fieldInterpolator) { m_helper.set_field_interpolator(fieldInterpolator); }
  InterpolateFieldsSharePtr get_field_interpolator() const                 { return m_helper.get_field_interpolator(); }
  bool check_field_interpolator() const                                    { return m_helper.check_field_interpolator(); }
  bool make_field_interpolator()                                           { return m_helper.make_field_interpolator(); }

  void set_parametric_tolerance(double parametricTolerance) { m_helper.set_parametric_tolerance(parametricTolerance); }
  double get_parametric_tolerance() const                   { return m_helper.get_parametric_tolerance(); }
  bool check_parametric_tolerance() const                   { return m_helper.check_parametric_tolerance(); }
  bool make_parametric_tolerance()                          { return m_helper.make_parametric_tolerance(); }

  void set_geometric_tolerance(double geometricTolerance) { m_helper.set_geometric_tolerance(geometricTolerance); }
  double get_geometric_tolerance() const                  { return m_helper.get_geometric_tolerance(); }
  bool check_geometric_tolerance() const                  { return m_helper.check_geometric_tolerance(); }
  bool make_geometric_tolerance()                         { return m_helper.make_geometric_tolerance(); }

  void set_active_selector(const stk::mesh::Selector& activeSelector) { m_helper.set_active_selector(activeSelector); }
  const stk::mesh::Selector& get_active_selector() const              { return m_helper.get_active_selector(); }
  bool check_active_selector() const                                  { return m_helper.check_active_selector(); }
  bool make_active_selector()                                         { return m_helper.make_active_selector(); }

  void set_default_part(const stk::mesh::Part* defaultPart) { m_helper.set_default_part(defaultPart); }
  const stk::mesh::Part* get_default_part() const           { return m_helper.get_default_part(); }
  bool check_default_part() const                           { return m_helper.check_default_part(); }
  bool make_default_part()                                  { return m_helper.make_default_part(); }

  void set_node_send_mesh(NodeSendMeshSharedPtr nodeSendMesh) { m_helper.set_node_send_mesh(nodeSendMesh); }
  NodeSendMeshSharedPtr get_node_send_mesh() const            { return m_helper.get_node_send_mesh(); }
  bool check_node_send_mesh() const                           { return m_helper.check_node_send_mesh(); }
  bool make_node_send_mesh()                                  { return m_helper.make_node_send_mesh(); }

  void set_element_send_mesh(ElemSendMeshSharedPtr elemSendMesh) { m_helper.set_element_send_mesh(elemSendMesh); }
  ElemSendMeshSharedPtr get_element_send_mesh() const            { return m_helper.get_element_send_mesh(); }
  bool check_element_send_mesh() const                           { return m_helper.check_element_send_mesh(); }
  bool make_element_send_mesh()                                  { return m_helper.make_element_send_mesh(); }

  void set_element_mesh_entity_rank(stk::mesh::EntityRank elementMeshEntityRank) { m_helper.set_element_mesh_entity_rank(elementMeshEntityRank); }
  stk::mesh::EntityRank get_element_mesh_entity_rank() const                     { return m_helper.get_element_mesh_entity_rank(); }
  bool check_element_mesh_entity_rank() const                                    { return m_helper.check_element_mesh_entity_rank(); }
  bool make_element_mesh_entity_rank()                                           { return m_helper.make_element_mesh_entity_rank(); }

  void set_verbose(bool flag) { m_helper.set_verbose(flag); }
  bool get_verbose() const    { return m_helper.get_verbose(); }

  bool make_prerequisite(TransferOptions option, const std::string& preamble)  { return m_helper.make_prerequisite(option, preamble); }
  bool check_prerequisite(TransferOptions option, const std::string& preamble) { return m_helper.check_prerequisite(option, preamble); }

  const std::vector<stk::transfer::FieldSpec>& get_fields() const { return m_helper.get_fields(); }
  const std::vector<std::string>& get_part_names() const { return m_helper.get_part_names(); }

 private:
  std::string m_transferName{"UN-INITIALIZED TRANSFER NAME"};
  stk::mesh::BulkData& m_bulk;
  TransferOptionHelper m_helper;
};

class GeometricRecvTransferOptions {
 public:
  GeometricRecvTransferOptions(const std::string& transferName, stk::mesh::BulkData& recvBulk)
    : m_transferName(transferName)
    , m_bulk(recvBulk)
    , m_helper(transferName, recvBulk, "Recv")
  { }

  GeometricRecvTransferOptions() = delete;
  ~GeometricRecvTransferOptions() = default;

  const std::string& name() const { return m_transferName; }
  const stk::mesh::BulkData& get_bulk() const { return m_bulk; }

  bool has_option(TransferOptions flag) const { return m_helper.has_option(flag); }

  // Allowable implementations for destination mesh
  void set_mesh_type(RecvMeshType meshType) { m_helper.set_recv_mesh_type(meshType); }
  RecvMeshType get_mesh_type() const        { return m_helper.get_recv_mesh_type(); }
  bool check_mesh_type() const              { return m_helper.check_recv_mesh_type(); }
  bool make_mesh_type()                     { return m_helper.make_recv_mesh_type(); }

  void set_coordinate_field_name(const std::string& coordinateFieldName) { m_helper.set_coordinate_field_name(coordinateFieldName); }
  const std::string& get_coordinate_field_name() const                   { return m_helper.get_coordinate_field_name(); }
  bool check_coordinate_field_name() const                               { return m_helper.check_coordinate_field_name(); }
  bool make_coordinate_field_name()                                      { return m_helper.make_coordinate_field_name(); }

  void set_coordinate_field_state(stk::mesh::FieldState coordinateFieldState) { m_helper.set_coordinate_field_state(coordinateFieldState); }
  stk::mesh::FieldState get_coordinate_field_state() const                    { return m_helper.get_coordinate_field_state(); }
  bool check_coordinate_field_state() const                                   { return m_helper.check_coordinate_field_state(); }
  bool make_coordinate_field_state()                                          { return m_helper.make_coordinate_field_state(); }

  void set_coordinate_field(const stk::mesh::FieldBase* coordinateField)  { m_helper.set_coordinate_field(coordinateField); }
  const stk::mesh::FieldBase* get_coordinate_field() const                { return m_helper.get_coordinate_field(); }
  bool check_coordinate_field() const                                     { return m_helper.check_coordinate_field(); }
  bool make_coordinate_field()                                            { return m_helper.make_coordinate_field(); }

  void add_field_spec(FieldSpec& fieldSpec)                                                        { m_helper.add_field_spec(fieldSpec); }
  void add_field_spec(const std::vector<FieldSpec>& fieldSpecs)                                    { m_helper.add_field_spec(fieldSpecs); }
  void add_field_spec(const std::string& fName)                                                    { m_helper.add_field_spec(fName); }
  void add_field_spec(const std::string& fName, stk::mesh::FieldState fState)                      { m_helper.add_field_spec(fName, fState); }
  void add_field_spec(const std::string& fName, stk::mesh::FieldState fState, unsigned int fIndex) { m_helper.add_field_spec(fName, fState, fIndex); }

  void add_part(const std::string& partName)                { m_helper.add_part(partName); }
  void add_part(const std::vector<std::string>& partNames)  { m_helper.add_part(partNames); }

  void set_point_evaluator(PointEvaluatorSharedPtr pointEvaluator) { m_helper.set_point_evaluator(pointEvaluator); }
  PointEvaluatorSharedPtr get_point_evaluator() const              { return m_helper.get_point_evaluator(); }
  bool check_point_evaluator() const                               { return m_helper.check_point_evaluator(); }
  bool make_point_evaluator()                                      { return m_helper.make_point_evaluator(); }

  void set_master_element_provider(MasterElemProviderSharedPtr masterElemProvider) { m_helper.set_master_element_provider(masterElemProvider); }
  MasterElemProviderSharedPtr get_master_element_provider() const                  { return m_helper.get_master_element_provider(); }
  bool check_master_element_provider() const                                       { return m_helper.check_master_element_provider(); }
  bool make_master_element_provider()                                              { return m_helper.make_master_element_provider(); }

  void set_parametric_tolerance(double parametricTolerance) { m_helper.set_parametric_tolerance(parametricTolerance); }
  double get_parametric_tolerance() const                   { return m_helper.get_parametric_tolerance(); }
  bool check_parametric_tolerance() const                   { return m_helper.check_parametric_tolerance(); }
  bool make_parametric_tolerance()                          { return m_helper.make_parametric_tolerance(); }

  void set_geometric_tolerance(double geometricTolerance) { m_helper.set_geometric_tolerance(geometricTolerance); }
  double get_geometric_tolerance() const                  { return m_helper.get_geometric_tolerance(); }
  bool check_geometric_tolerance() const                  { return m_helper.check_geometric_tolerance(); }
  bool make_geometric_tolerance()                         { return m_helper.make_geometric_tolerance(); }

  void set_active_selector(const stk::mesh::Selector& activeSelector) { m_helper.set_active_selector(activeSelector); }
  const stk::mesh::Selector& get_active_selector() const              { return m_helper.get_active_selector(); }
  bool check_active_selector() const                                  { return m_helper.check_active_selector(); }
  bool make_active_selector()                                         { return m_helper.make_active_selector(); }

  void set_default_part(const stk::mesh::Part* defaultPart) { m_helper.set_default_part(defaultPart); }
  const stk::mesh::Part* get_default_part() const           { return m_helper.get_default_part(); }
  bool check_default_part() const                           { return m_helper.check_default_part(); }
  bool make_default_part()                                  { return m_helper.make_default_part(); }

  void set_node_recv_mesh(NodeRecvMeshSharedPtr nodeRecvMesh) { m_helper.set_node_recv_mesh(nodeRecvMesh); }
  NodeRecvMeshSharedPtr get_node_recv_mesh() const            { return m_helper.get_node_recv_mesh(); }
  bool check_node_recv_mesh() const                           { return m_helper.check_node_recv_mesh(); }
  bool make_node_recv_mesh()                                  { return m_helper.make_node_recv_mesh(); }

  void set_element_recv_mesh(ElemRecvMeshSharedPtr elemRecvMesh) { m_helper.set_element_recv_mesh(elemRecvMesh); }
  ElemRecvMeshSharedPtr get_element_recv_mesh() const            { return m_helper.get_element_recv_mesh(); }
  bool check_element_recv_mesh() const                           { return m_helper.check_element_recv_mesh(); }
  bool make_element_recv_mesh()                                  { return m_helper.make_element_recv_mesh(); }

  void set_element_mesh_entity_rank(stk::mesh::EntityRank elementMeshEntityRank) { m_helper.set_element_mesh_entity_rank(elementMeshEntityRank); }
  stk::mesh::EntityRank get_element_mesh_entity_rank() const                     { return m_helper.get_element_mesh_entity_rank(); }
  bool check_element_mesh_entity_rank() const                                    { return m_helper.check_element_mesh_entity_rank(); }
  bool make_element_mesh_entity_rank()                                           { return m_helper.make_element_mesh_entity_rank(); }

  void set_verbose(bool flag) { m_helper.set_verbose(flag); }
  bool get_verbose() const    { return m_helper.get_verbose(); }

  bool make_prerequisite(TransferOptions option, const std::string& preamble)  { return m_helper.make_prerequisite(option, preamble); }
  bool check_prerequisite(TransferOptions option, const std::string& preamble) { return m_helper.check_prerequisite(option, preamble); }

  const std::vector<stk::transfer::FieldSpec>& get_fields() const { return m_helper.get_fields(); }
  const std::vector<std::string>& get_part_names() const { return m_helper.get_part_names(); }

 private:
  std::string m_transferName{"UN-INITIALIZED TRANSFER NAME"};
  stk::mesh::BulkData& m_bulk;
  TransferOptionHelper m_helper;
};

class GeometricTransferOptions {
 public:
  GeometricTransferOptions(const std::string& transferName, stk::mesh::BulkData& sendBulk, stk::mesh::BulkData& recvBulk)
  : m_transferName(transferName)
  , m_sendBulk(sendBulk)
  , m_recvBulk(recvBulk)
  , m_helper(transferName, sendBulk, "Transfer")
  , m_sendOptions(transferName, sendBulk)
  , m_recvOptions(transferName, recvBulk)
  { }

  GeometricTransferOptions() = delete;
  ~GeometricTransferOptions() = default;

  GeometricSendTransferOptions& get_send_options() { return m_sendOptions; }
  GeometricRecvTransferOptions& get_recv_options() { return m_recvOptions; }

  const GeometricSendTransferOptions& get_send_options() const { return m_sendOptions; }
  const GeometricRecvTransferOptions& get_recv_options() const { return m_recvOptions; }

  const std::string& name() const { return m_transferName; }

  const stk::mesh::BulkData& get_send_bulk() const { return m_sendBulk; }
  const stk::mesh::BulkData& get_recv_bulk() const { return m_recvBulk; }

  bool has_option(TransferOptions flag) const { return m_helper.has_option(flag); }

  void set_interpolation_type(InterpolationType interpolationType) { m_helper.set_interpolation_type(interpolationType);
                                                                     m_sendOptions.set_interpolation_type(interpolationType); }
  InterpolationType get_interpolation_type() const                 { return m_helper.get_interpolation_type(); }
  bool check_interpolation_type() const                            { return m_helper.check_interpolation_type(); }
  bool make_interpolation_type()                                   { return m_helper.make_interpolation_type(); }

  void set_search_expansion_factor(double searchExpansionFactor) { m_helper.set_search_expansion_factor(searchExpansionFactor); }
  double get_search_expansion_factor() const                     { return m_helper.get_search_expansion_factor(); }
  bool check_search_expansion_factor() const                     { return m_helper.check_search_expansion_factor(); }
  bool make_search_expansion_factor()                            { return m_helper.make_search_expansion_factor(); }

  void set_search_expansion_padding(double searchExpansionPadding) { m_helper.set_search_expansion_padding(searchExpansionPadding); }
  double get_search_expansion_padding() const                      { return m_helper.get_search_expansion_padding(); }
  bool check_search_expansion_padding() const                      { return m_helper.check_search_expansion_padding(); }
  bool make_search_expansion_padding()                             { return m_helper.make_search_expansion_padding(); }

  void set_parallel_machine(stk::ParallelMachine parallelMachine) { m_helper.set_parallel_machine(parallelMachine); }
  stk::ParallelMachine get_parallel_machine() const               { return m_helper.get_parallel_machine(); }
  bool check_parallel_machine() const                             { return m_helper.check_parallel_machine(); }
  bool make_parallel_machine()                                    { return m_helper.make_parallel_machine(); }

  void set_verbose(bool flag) { m_helper.set_verbose(flag); m_sendOptions.set_verbose(flag); m_recvOptions.set_verbose(flag); }
  bool get_verbose() const { return m_helper.get_verbose(); }

  void check_field_validity() const;

 private:
  std::string m_transferName{"UN-INITIALIZED TRANSFER NAME"};

  stk::mesh::BulkData& m_sendBulk;
  stk::mesh::BulkData& m_recvBulk;

  TransferOptionHelper m_helper;

  GeometricSendTransferOptions m_sendOptions;
  GeometricRecvTransferOptions m_recvOptions;
};

using GeometricTransferPair =
    std::pair<std::shared_ptr<stk::transfer::TransferBase>, std::shared_ptr<stk::transfer::spmd::GeometricTransferDispatchBase>>;

TransferType get_transfer_type(GeometricTransferOptions& transferOptions, bool throwOnError = false);

GeometricTransferPair create_node_to_node_transfer(GeometricTransferOptions& transferOptions, bool throwOnError = true);
GeometricTransferPair create_node_to_element_transfer(GeometricTransferOptions& transferOptions, bool throwOnError = true);
GeometricTransferPair create_element_to_node_transfer(GeometricTransferOptions& transferOptions, bool throwOnError = true);
GeometricTransferPair create_element_to_element_transfer(GeometricTransferOptions& transferOptions, bool throwOnError = true);

GeometricTransferPair create_transfer(GeometricTransferOptions& transferOptions, bool throwOnError = true);

std::shared_ptr<stk::transfer::spmd::GeometricTransferDispatchBase>
create_dispatch(std::shared_ptr<stk::transfer::TransferBase> inputTransfer);

} // namespace spmd
} // namespace transfer
} // namespace stk

#endif /* STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICTRANSFEROPTIONS_HPP_ */
