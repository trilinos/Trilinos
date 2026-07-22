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
#include <stk_math/StkVector.hpp>
#include <Akri_FieldRef.hpp>

namespace krino {

class FacetedSurfaceBase;
class Facet3d;

bool is_parallel_io_enabled();
bool does_file_exist(const std::string& filename);

void output_mesh_with_fields(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const std::string & fileName, int step, double time, stk::io::DatabasePurpose purpose = stk::io::WRITE_RESULTS);
void output_composed_mesh_with_fields(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const std::string & fileName, int step, double time, stk::io::DatabasePurpose purpose = stk::io::WRITE_RESULTS);
void fix_ownership_and_output_composed_mesh_with_fields(stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector, const std::string & fileName, int step, double time, stk::io::DatabasePurpose purpose = stk::io::WRITE_RESULTS);
std::string create_file_name(const std::string & fileBaseName, const int fileIndex);
stk::mesh::PartVector turn_off_output_for_empty_io_parts(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & outputSelector);

void write_facets(const stk::mesh::BulkData & mesh, const FacetedSurfaceBase & facets, const std::string & fileBaseName, const int fileIndex);
void write_shells_for_surfaces(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const std::vector<stk::mesh::Part*> & surfaces, const stk::mesh::Part & activePart, const std::string & fileBaseName, const int fileIndex);

std::string stl_parallel_file_name(const std::string &fileBaseName, const stk::ParallelMachine comm);
void write_stl(const std::string &filename, const std::vector<const Facet3d *> & facets);
void write_stl(const std::string &filename, const std::vector<Facet3d> & facets);
void write_stl(const std::string &filename,
  const std::vector<stk::math::Vector3d> & vertices,
  const std::vector<std::array<unsigned,3>> & facetConnectivity);
}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_OUTPUTUTILS_HPP_ */
