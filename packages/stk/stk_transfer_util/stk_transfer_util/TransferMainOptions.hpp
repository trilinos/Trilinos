/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_TRANSFER_UTIL_TRANSFERMAINOPTIONS_HPP
#define STK_TRANSFER_UTIL_TRANSFERMAINOPTIONS_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <string>
#include <ostream>

namespace stk {
namespace transfer_util {

struct OptionNames
{
  const std::string fromMesh = "from-mesh";
  const std::string toMesh = "to-mesh";
};

class TransferMainOptions
{
public:
  TransferMainOptions(MPI_Comm comm, int argc, const char** argv);

  std::string get_fromMesh_filename() const;
  std::string get_toMesh_filename() const;

  stk::CommandLineParser::ParseState get_parse_state() const;

  void print_usage_help(std::ostream& os) const;

private:
  OptionNames m_optionNames;
  stk::CommandLineParserParallel m_cmdLineParser;
  stk::CommandLineParser::ParseState m_parseState;
};

} // namespace transfer_util
} // namespace stk

#endif // STK_TRANSFER_UTIL_TRANSFERMAINOPTIONS_HPP
