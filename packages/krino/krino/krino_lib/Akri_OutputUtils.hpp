/*
 * Akri_OutputUtils.hpp
 *
 *  Created on: Nov 21, 2022
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_HPP_
#include <string>
#include <vector>

#include <stk_io/DatabasePurpose.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>

namespace krino {

template<class FACET>
class Faceted_Surface;

void output_mesh_with_fields(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const std::string & fileName, int step, double time, stk::io::DatabasePurpose purpose = stk::io::WRITE_RESULTS);
void output_composed_mesh_with_fields(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const std::string & fileName, int step, double time, stk::io::DatabasePurpose purpose = stk::io::WRITE_RESULTS);
std::string create_file_name(const std::string & fileBaseName, const int fileIndex);
std::string create_filename_from_base_filename(const std::string & baseFileName, const int numFileRevisions);
stk::mesh::PartVector turn_off_output_for_empty_io_parts(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector);

template<class FACET>
void write_facets( const std::vector<FACET> & facetedSurface, const std::string & fileBaseName, const int fileIndex, const stk::ParallelMachine comm);

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_HPP_ */
