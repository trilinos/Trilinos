/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef MOCK_CONFIGURATION_INFO_HPP
#define MOCK_CONFIGURATION_INFO_HPP

#include <map>
#include <string>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp> // for CommBuffer

namespace mock_coupling
{

class ConfigurationInfo
{
public:
  std::map<std::string, double> dvals;
  std::map<std::string, int> ivals;
  std::map<std::string, std::string> svals;

  ConfigurationInfo exchange(stk::ParallelMachine global, stk::ParallelMachine local);

private:
  void pack(stk::CommBuffer & b) const;
  void unpack(stk::CommBuffer & b);
};

}

#endif /* MOCK_CONFIGURATION_INFO_HPP */
