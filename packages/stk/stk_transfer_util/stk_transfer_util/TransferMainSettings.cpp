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

TransferMainSettings::TransferMainSettings()
  : m_numInputProcessors(1),
    m_numOutputProcessors(1)
{
  m_OODP = stk::search::get_object_outside_domain_policy("IGNORE");
}

void TransferMainSettings::set_num_input_processors(unsigned numInputProcs)
{
  m_numInputProcessors = numInputProcs;
}

void TransferMainSettings::set_num_output_processors(unsigned numOutputProcs)
{
  m_numOutputProcessors = numOutputProcs;
}

void TransferMainSettings::set_fromMesh_filename(const std::string& fromMesh)
{
 m_fromMesh = fromMesh; 
}

void TransferMainSettings::set_toMesh_filename(const std::string& toMesh)
{
  m_toMesh = toMesh;
}
void TransferMainSettings::set_transfer_field(const std::string& fieldName)
{
  m_transferFields.push_back(fieldName); 
}

bool TransferMainSettings::set_extrapolate_option(const std::string& policy)
{
  m_OODP = stk::search::get_object_outside_domain_policy(policy);
  return (m_OODP != stk::search::ObjectOutsideDomainPolicy::UNDEFINED_OBJFLAG);
}

unsigned TransferMainSettings::get_num_input_processors() const
{
  return m_numInputProcessors;
}

unsigned TransferMainSettings::get_num_output_processors() const
{
  return m_numOutputProcessors;
}

std::string TransferMainSettings::get_fromMesh_filename() const
{
  return m_fromMesh;
}
std::string TransferMainSettings::get_toMesh_filename() const
{
  return m_toMesh;
}

std::vector<std::string> TransferMainSettings::get_transfer_fields() const
{
  return m_transferFields;
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
    fieldListString += m_transferFields[i];
  }

  return fieldListString;
}

} // namespace transfer_util
} // namespace stk
