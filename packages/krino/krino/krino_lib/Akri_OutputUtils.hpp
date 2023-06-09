/*
 * Akri_OutputUtils.hpp
 *
 *  Created on: Nov 21, 2022
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_HPP_
#include <string>

#include <stk_io/DatabasePurpose.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace stk { namespace mesh{ class BulkData; } }
namespace stk { namespace mesh{ class Part; } }
namespace krino { class Faceted_Surface; }

namespace krino {

void output_field_with_mesh(const stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart, const std::string & fileName, int step, double time, stk::io::DatabasePurpose purpose = stk::io::WRITE_RESULTS);
std::string create_file_name(const std::string & fileBaseName, const int fileIndex);
void write_facets( const int dim, const Faceted_Surface & facetedSurface, const std::string & fileBaseName, const int fileIndex, const stk::ParallelMachine comm);

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_HPP_ */
