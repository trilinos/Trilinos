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

namespace stk {
namespace transfer_util {

class TransferMainSettings
{
public:

  TransferMainSettings();

  void set_num_input_processors(unsigned numInputProcs);
  void set_num_output_processors(unsigned numOutputProcs);
  void set_sendMesh_filename(const std::string& sendMesh);
  void set_recvMesh_filename(const std::string& recvMesh);
  void set_transfer_field(const std::pair<std::string, std::string>& fieldName);
  bool set_extrapolate_option(const std::string& policy);

  unsigned get_num_input_processors() const;
  unsigned get_num_output_processors() const;
  const std::string& get_sendMesh_filename() const;
  const std::string& get_recvMesh_filename() const;
  const std::vector<std::pair<std::string, std::string>>& get_transfer_fields() const;
  stk::search::ObjectOutsideDomainPolicy get_extrapolate_option() const;
  std::string get_extrapolate_option_string() const;
  std::string get_field_list_string() const;

private:
  unsigned m_numInputProcessors;
  unsigned m_numOutputProcessors;
  std::string m_sendMesh;
  std::string m_recvMesh;
  std::vector<std::pair<std::string, std::string>> m_transferFields;
  stk::search::ObjectOutsideDomainPolicy m_OODP;

};

} // namespace transfer_util
} // namespace stk

#endif // STK_TRANSFER_UTIL_TRANSFERMAINSETTINGS_HPP
