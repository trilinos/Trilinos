/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <string>
#include <ostream>
#include <stk_search/ObjectOutsideDomainPolicy.hpp>
#include <stk_transfer_util/TransferMainSettings.hpp>

namespace stk {
namespace transfer_util {

stk::transfer::spmd::RecvMeshType convert_to_recv_mesh_type(const std::string& id)
{
  if(id == "NODE") return stk::transfer::spmd::RecvMeshType::NODE;
  if(id == "EDGE_CENTROID") return stk::transfer::spmd::RecvMeshType::EDGE_CENTROID;
  if(id == "EDGE_GAUSS_POINT") return stk::transfer::spmd::RecvMeshType::EDGE_GAUSS_POINT;
  if(id == "FACE_CENTROID") return stk::transfer::spmd::RecvMeshType::FACE_CENTROID;
  if(id == "FACE_GAUSS_POINT") return stk::transfer::spmd::RecvMeshType::FACE_GAUSS_POINT;
  if(id == "ELEMENT_CENTROID") return stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID;
  if(id == "ELEMENT_GAUSS_POINT") return stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT;

  return stk::transfer::spmd::RecvMeshType::INVALID;
}

TransferMainSettings::TransferMainSettings()
  : m_numInputProcessors(1),
    m_numOutputProcessors(1),
    m_recvTypeString("NODE"),
    m_transferType("INTERP"),
    m_coordTransformZExpr("z=0")
{
  m_OODP = stk::search::get_object_outside_domain_policy("IGNORE");
  m_recvType = convert_to_recv_mesh_type(m_recvTypeString);
}

void TransferMainSettings::set_num_input_processors(unsigned numInputProcs)
{
  m_numInputProcessors = numInputProcs;
}

void TransferMainSettings::set_num_output_processors(unsigned numOutputProcs)
{
  m_numOutputProcessors = numOutputProcs;
}

void TransferMainSettings::set_sendMesh_filename(const std::string& sendMesh)
{
 m_sendMesh = sendMesh;
}

void TransferMainSettings::set_recvMesh_filename(const std::string& recvMesh)
{
  m_recvMesh = recvMesh;
}

void TransferMainSettings::set_outputMesh_filename(const std::string& outputMesh)
{
  m_outputMesh = outputMesh;
}

void TransferMainSettings::set_transfer_field(const std::pair<std::string, std::string>& fieldNamePair)
{
  m_transferFields.push_back(fieldNamePair);
}

void TransferMainSettings::set_transfer_send_parts(const std::vector<std::string>& partName)
{
  m_transferSendParts = partName;
}

void TransferMainSettings::set_transfer_recv_parts(const std::vector<std::string>& partName)
{
  m_transferRecvParts = partName;
}

bool TransferMainSettings::set_extrapolate_option(const std::string& policy)
{
  m_OODP = stk::search::get_object_outside_domain_policy(policy);
  return (m_OODP != stk::search::ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG);
}

void TransferMainSettings::set_master_elements_name(const std::string& masterElementsName)
{
  m_masterElementsName = masterElementsName;
}

bool TransferMainSettings::set_recv_type(const std::string& recvType)
{
  m_recvTypeString = recvType;
  m_recvType = convert_to_recv_mesh_type(m_recvTypeString);
  return (m_recvType != stk::transfer::spmd::RecvMeshType::INVALID);
}

void TransferMainSettings::set_time_steps_spec(const std::string& timeSteps)
{
  m_timeSteps = timeSteps;
}

void TransferMainSettings::set_coord_transf_z_expr(const std::string& coordTransformZExpression)
{
  m_coordTransformZExpr = coordTransformZExpression;
}

bool TransferMainSettings::set_transfer_type(const std::string& xferType)
{
  m_transferType = xferType;
  return (m_transferType == "INTERP" || m_transferType == "COPY" || m_transferType == "PATCH");
}

bool TransferMainSettings::set_2d_to_3d_mapping_type(const std::string& mappingType)
{
  m_2d_to_3d_mappingType = mappingType;
  return (m_2d_to_3d_mappingType == "ZPLANE" || m_2d_to_3d_mappingType == "EXTRUDE");
}

void TransferMainSettings::set_2d_to_3d_axisymmetric_transfer_params(const TwoDTo3DAxisymmetricParams& params)
{
  m_enable2Dto3DAxisymmetric = true;
  m_twoDToThreeDAxisymmetricParams = params;
}

void TransferMainSettings::set_3d_to_3d_axisymmetric_transfer_params(const ThreeDTo3DAxisymmetricParams& params)
{
  m_enable3Dto3DAxisymmetric = true;
  m_threeDToThreeDAxisymmetricParams = params;
}

unsigned TransferMainSettings::get_num_input_processors() const
{
  return m_numInputProcessors;
}

unsigned TransferMainSettings::get_num_output_processors() const
{
  return m_numOutputProcessors;
}

const std::string& TransferMainSettings::get_sendMesh_filename() const
{
  return m_sendMesh;
}
const std::string& TransferMainSettings::get_recvMesh_filename() const
{
  return m_recvMesh;
}

const std::string& TransferMainSettings::get_outputMesh_filename() const
{
  return m_outputMesh;
}

const std::vector<std::pair<std::string, std::string>>& TransferMainSettings::get_transfer_fields() const
{
  return m_transferFields;
}

const std::vector<std::string>& TransferMainSettings::get_transfer_send_parts() const {
  return m_transferSendParts;
}

const std::vector<std::string>& TransferMainSettings::get_transfer_recv_parts() const {
  return m_transferRecvParts;
}

stk::search::ObjectOutsideDomainPolicy TransferMainSettings::get_extrapolate_option() const
{
  return m_OODP;
}

std::string TransferMainSettings::get_extrapolate_option_string() const
{
  return stk::search::get_object_outside_domain_policy(m_OODP);
}

std::string TransferMainSettings::get_field_list_string() const
{
  std::string fieldListString;

  for (unsigned i = 0; i < m_transferFields.size(); i++) {
    if (i != 0) { fieldListString += ", "; }
    fieldListString += m_transferFields[i].first + ":" + m_transferFields[i].second;
  }

  return fieldListString;
}

std::string TransferMainSettings::get_send_part_list_string() const
{
  std::string sendPartString;

  for (unsigned i = 0; i < m_transferSendParts.size(); i++) {
    if (i != 0) { sendPartString += ", "; }
    sendPartString += m_transferSendParts[i];
  }

  return sendPartString;
}

std::string TransferMainSettings::get_recv_part_list_string() const
{
  std::string recvPartString;

  for (unsigned i = 0; i < m_transferRecvParts.size(); i++) {
    if (i != 0) { recvPartString += ", "; }
    recvPartString += m_transferRecvParts[i];
  }

  return recvPartString;
}

std::string TransferMainSettings::get_part_list_string() const
{
  std::string partListString;

  partListString += get_send_part_list_string() + ":" + get_recv_part_list_string();

  return partListString;
}

const std::string& TransferMainSettings::get_master_elements_name() const
{
  return m_masterElementsName;
}

stk::transfer::spmd::RecvMeshType TransferMainSettings::get_recv_type() const
{
  return m_recvType;
}

std::string TransferMainSettings::get_recv_type_string() const
{
  return m_recvTypeString;
}

stk::mesh::EntityRank TransferMainSettings::get_default_recv_part_rank() const
{
  if(m_recvType == stk::transfer::spmd::RecvMeshType::NODE) return stk::topology::ELEMENT_RANK;
  if(m_recvType == stk::transfer::spmd::RecvMeshType::EDGE_CENTROID) return stk::topology::EDGE_RANK;
  if(m_recvType == stk::transfer::spmd::RecvMeshType::EDGE_GAUSS_POINT) return stk::topology::EDGE_RANK;
  if(m_recvType == stk::transfer::spmd::RecvMeshType::FACE_CENTROID) return stk::topology::FACE_RANK;
  if(m_recvType == stk::transfer::spmd::RecvMeshType::FACE_GAUSS_POINT) return stk::topology::FACE_RANK;
  if(m_recvType == stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID) return stk::topology::ELEMENT_RANK;
  if(m_recvType == stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT) return stk::topology::ELEMENT_RANK;

  return stk::topology::INVALID_RANK;
}

stk::mesh::EntityRank TransferMainSettings::get_recv_field_rank() const
{
  if(m_recvType == stk::transfer::spmd::RecvMeshType::NODE) return stk::topology::NODE_RANK;
  if(m_recvType == stk::transfer::spmd::RecvMeshType::EDGE_CENTROID) return stk::topology::EDGE_RANK;
  if(m_recvType == stk::transfer::spmd::RecvMeshType::EDGE_GAUSS_POINT) return stk::topology::EDGE_RANK;
  if(m_recvType == stk::transfer::spmd::RecvMeshType::FACE_CENTROID) return stk::topology::FACE_RANK;
  if(m_recvType == stk::transfer::spmd::RecvMeshType::FACE_GAUSS_POINT) return stk::topology::FACE_RANK;
  if(m_recvType == stk::transfer::spmd::RecvMeshType::ELEMENT_CENTROID) return stk::topology::ELEMENT_RANK;
  if(m_recvType == stk::transfer::spmd::RecvMeshType::ELEMENT_GAUSS_POINT) return stk::topology::ELEMENT_RANK;

  return stk::topology::INVALID_RANK;
}

std::string TransferMainSettings::get_time_steps_spec() const
{
  return m_timeSteps;
}

std::string TransferMainSettings::get_coord_transf_z_expr() const
{
  return m_coordTransformZExpr;
}

std::string TransferMainSettings::get_transfer_type() const
{
  return m_transferType;
}

std::string TransferMainSettings::get_2d_to_3d_mapping_type() const
{
  return m_2d_to_3d_mappingType;
}

bool TransferMainSettings::get_enable_2d_to_3d_axisymmetric_transfer() const
{
  return m_enable2Dto3DAxisymmetric;
}

const TwoDTo3DAxisymmetricParams& TransferMainSettings::get_2d_to_3d_axisymmetric_transfer_params() const
{
  return m_twoDToThreeDAxisymmetricParams;
}

bool TransferMainSettings::get_enable_3d_to_3d_axisymmetric_transfer() const
{
  return m_enable3Dto3DAxisymmetric;
}

const ThreeDTo3DAxisymmetricParams& TransferMainSettings::get_3d_to_3d_axisymmetric_transfer_params() const
{
  return m_threeDToThreeDAxisymmetricParams;
}

} // namespace transfer_util
} // namespace stk
