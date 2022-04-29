// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Region_h
#define Akri_Region_h

#include <stk_util/diag/Timer.hpp>
#include <string>
#include <memory>

namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace io { class StkMeshIoBroker; } }
namespace krino { class BoundingBoxMesh; }
namespace krino { class RegionForwarder; }
namespace krino { class ResultsOutputOptions; }
namespace krino { class Simulation; }
namespace Ioss { class Region; }

namespace krino{

class Region {
public:
  Region(Simulation & owning_simulation, const std::string & regionName);
  virtual ~Region();

  virtual void commit();
  virtual void initialize();
  virtual void execute();

  double time_step() const;
  const std::string & name() const { return my_name; }
  unsigned spatial_dimension() const;
  const stk::mesh::BulkData& get_stk_mesh_bulk_data() const;
  stk::mesh::BulkData& get_stk_mesh_bulk_data();
  const stk::mesh::MetaData& get_stk_mesh_meta_data() const;
  stk::mesh::MetaData& get_stk_mesh_meta_data();
  stk::diag::Timer & getRegionTimer() const { return my_timerRegion; }
  stk::diag::Timer & getMeshInputTimer() const { return my_timerMeshInput; }
  stk::diag::Timer & getMeshOutputTimer() const { return my_timerMeshOutput; }

  stk::io::StkMeshIoBroker & stk_IO();
  std::string name_of_input_mesh() const;
  Ioss::Region * get_input_io_region();
  void associate_input_mesh(const std::string & model_name, bool assert_32bit_ids, bool force_64bit_ids);
  void set_generated_mesh_domain();
  void create_output_mesh();
  void declare_output_variables(size_t result_output_index);
  void process_output(bool forceOutput);
  ResultsOutputOptions * get_results_options() { return my_results_options.get(); }
  void set_initial_refinement_levels(int levels) { my_initial_refinement_levels = levels; }

private:
  Simulation & my_simulation;
  stk::mesh::MetaData * my_meta;
  stk::mesh::BulkData * my_bulk;
  std::unique_ptr<stk::io::StkMeshIoBroker> myIOBroker;
  std::unique_ptr<BoundingBoxMesh> my_generated_mesh;
  std::unique_ptr<ResultsOutputOptions> my_results_options;
  std::string my_name;
  std::string my_input_model_name;

  mutable stk::diag::Timer my_timerRegion;           ///< Region's root timer
  mutable stk::diag::Timer my_timerInitialize;       ///< Initialize timer
  mutable stk::diag::Timer my_timerExecute;          ///< Execute timer
  mutable stk::diag::Timer my_timerMeshInput;        ///< Mesh input timer
  mutable stk::diag::Timer my_timerMeshOutput;       ///< Mesh output timer
  size_t my_output_file_index;
  bool my_output_file_created;
  int my_initial_refinement_levels;
};

} // namespace krino

#endif // Akri_Region_h
