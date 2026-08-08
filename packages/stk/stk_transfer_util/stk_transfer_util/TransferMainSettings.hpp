/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_TRANSFER_UTIL_TRANSFERMAINSETTINGS_HPP
#define STK_TRANSFER_UTIL_TRANSFERMAINSETTINGS_HPP

#include <utility>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <string>
#include <ostream>
#include <stk_search/ObjectOutsideDomainPolicy.hpp>
#include <stk_transfer_util/spmd/GeometricTransferOptions.hpp>

namespace stk {
namespace transfer_util {

struct TwoDTo3DAxisymmetricParams
{
  std::array<double, 3> axis = {0, 0, 1};
  std::array<double, 3> trans = {0, 0, 0};
};

struct ThreeDTo3DAxisymmetricParams
{
  double theta_min = 0;
  double theta_max = 0;
  std::array<double, 3> axis = {0, 0, 1};
  std::array<double, 3> trans = {0, 0, 0};
};

class TransferMainSettings
{
public:

  TransferMainSettings();

  void set_num_input_processors(unsigned numInputProcs);
  void set_num_output_processors(unsigned numOutputProcs);
  void set_sendMesh_filename(const std::string& sendMesh);
  void set_recvMesh_filename(const std::string& recvMesh);
  void set_outputMesh_filename(const std::string& recvMesh);
  void set_transfer_field(const std::pair<std::string, std::string>& fieldName);
  void set_transfer_send_parts(const std::vector<std::string>& partName);
  void set_transfer_recv_parts(const std::vector<std::string>& partName);
  bool set_extrapolate_option(const std::string& policy);
  void set_master_elements_name(const std::string& masterElementsName);
  bool set_recv_type(const std::string& recvType);
  void set_time_steps_spec(const std::string& timeSteps);
  void set_coord_transf_z_expr(const std::string& coordTransformExpression);
  bool set_transfer_type(const std::string& xferType);
  bool set_2d_to_3d_mapping_type(const std::string& mappingType);
  void set_2d_to_3d_axisymmetric_transfer_params(const TwoDTo3DAxisymmetricParams& params);
  void set_3d_to_3d_axisymmetric_transfer_params(const ThreeDTo3DAxisymmetricParams& params);

  unsigned get_num_input_processors() const;
  unsigned get_num_output_processors() const;
  const std::string& get_sendMesh_filename() const;
  const std::string& get_recvMesh_filename() const;
  const std::string& get_outputMesh_filename() const;
  const std::vector<std::pair<std::string, std::string>>& get_transfer_fields() const;
  const std::vector<std::string>& get_transfer_send_parts() const;
  const std::vector<std::string>& get_transfer_recv_parts() const;
  stk::search::ObjectOutsideDomainPolicy get_extrapolate_option() const;
  std::string get_extrapolate_option_string() const;
  std::string get_field_list_string() const;
  std::string get_send_part_list_string() const;
  std::string get_recv_part_list_string() const;
  std::string get_part_list_string() const;
  const std::string& get_master_elements_name() const;
  stk::transfer::spmd::RecvMeshType get_recv_type() const;
  std::string get_recv_type_string() const;
  std::string get_transfer_type() const;
  std::string get_2d_to_3d_mapping_type() const;
  std::string get_time_steps_spec() const;
  std::string get_coord_transf_z_expr() const;
  stk::mesh::EntityRank get_default_recv_part_rank() const;
  stk::mesh::EntityRank get_recv_field_rank() const;
  bool get_enable_2d_to_3d_axisymmetric_transfer() const;
  const TwoDTo3DAxisymmetricParams& get_2d_to_3d_axisymmetric_transfer_params() const;
  bool get_enable_3d_to_3d_axisymmetric_transfer() const;
  const ThreeDTo3DAxisymmetricParams& get_3d_to_3d_axisymmetric_transfer_params() const;

private:
  unsigned m_numInputProcessors;
  unsigned m_numOutputProcessors;
  std::string m_sendMesh;
  std::string m_recvMesh;
  std::string m_outputMesh;
  std::vector<std::pair<std::string, std::string>> m_transferFields;
  std::vector<std::string> m_transferSendParts;
  std::vector<std::string> m_transferRecvParts;
  stk::search::ObjectOutsideDomainPolicy m_OODP;
  std::string m_masterElementsName;
  stk::transfer::spmd::RecvMeshType m_recvType;
  std::string m_recvTypeString;
  std::string m_transferType;
  std::string m_2d_to_3d_mappingType;
  std::string m_timeSteps;
  bool m_enable2Dto3DAxisymmetric = false;
  TwoDTo3DAxisymmetricParams m_twoDToThreeDAxisymmetricParams;
  bool m_enable3Dto3DAxisymmetric = false;
  ThreeDTo3DAxisymmetricParams m_threeDToThreeDAxisymmetricParams;
  std::string m_coordTransformZExpr;
};

} // namespace transfer_util
} // namespace stk

#endif // STK_TRANSFER_UTIL_TRANSFERMAINSETTINGS_HPP
