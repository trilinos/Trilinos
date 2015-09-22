/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_search_use_cases_Performance_hpp
#define stk_search_use_cases_Performance_hpp

#include <string>
#include <stk_util/parallel/Parallel.hpp>

namespace stk_classic {
namespace search {
class Options
{
public:
  Options()
    : mesh_filename(""), mesh_type("exodusii"), entity(""),
      offset(0.0), scale(1.0) {}
  std::string mesh_filename;
  std::string mesh_type;
  std::string entity;
  double      offset;
  double      scale;
private:
  Options(const Options&); // Do not implement
  Options& operator=(const Options&); // Do not implement
};
}
}

void performance_driver(stk_classic::ParallelMachine  comm,
                        const std::string &working_directory,
                        const std::string &search_type,
                        const stk_classic::search::Options &range,
                        const stk_classic::search::Options &domain,
                        bool               performance);

#endif
