// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.



#define IN_MESH_ADAPT 1
#define MESH_ADAPT_CPP11_MVI 0

#include <sys/unistd.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>
#include <stdint.h>
#include <map>
#include <list>
#include <string>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>
#include <percept/RunEnvironment.hpp>
#include <percept/ProgressMeter.hpp>
#include <percept/Histograms.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/memory_util.hpp>

#include <adapt/RefinerUtil.hpp>

#include <adapt/main/RunAdaptRun.hpp>
#include <adapt/TransitionElementAdapter.hpp>
#include <adapt/ElementRefinePredicate.hpp>

#include <adapt/UniformRefiner.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <Ioss_Utils.h>
#include <Ioss_SerializeIO.h>

#include <adapt/AdaptedMeshVerifier.hpp>

#include <percept/mesh/mod/smoother/SpacingFieldUtil.hpp>
#include <percept/mesh/geometry/stk_geom/3D/FitGregoryPatches.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <percept/DihedralAngleCheck.hpp>
#include <percept/mesh/geometry/stk_geom/LocalCubicSplineFit.hpp>
#endif
#if HAVE_OPENNURBS
#include <percept/mesh/geometry/recovery/GeometryRecoverySplineFit.hpp>
#endif

#define ALLOW_MEM_TEST 1

#include "AdaptMain.hpp"
#include "RunAdaptRun.hpp"
#include <stk_mesh/base/MeshUtils.hpp>

class PGeom;

namespace stk { 
namespace diag {
  class TimerSet;
  class Timer;
}
}

namespace percept {

  void print_timer_table(stk::diag::Timer &timer);

  class MemoryInfo;

  class MeshAdapt {

  public:

    enum GeometryType {
      GEOM_NONE,
      OPENNURBS,
      PGEOM_OPENNURBS,
      MESH_BASED_GREGORY_PATCH,
      PGEOM_ACIS,
      N_GeometryType
    };

    // member variable definition and initialization (c++11) 
#include "MeshAdaptMemberVarInit.hpp"    

    // constructor
    MeshAdapt();

    ~MeshAdapt();

    bool has_suffix(const std::string &str, const std::string &suffix) {
      return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }
    // true if s is a well-formed integer, and set its value in i
    bool is_number(std::string &s, int &i);

    // Determine the basename and refinement level from a given an ExodusII file name.
    // This is the last integer before the final '.e' or '.exo' in the file name.
    // If no level was found, returns 0 as the level.
    //
    // Return:
    //    string - the basename which is everything before the number
    //    int    - the refinement level
    std::pair<std::string,int> get_basename_and_level(const std::string& input_mesh);
    void to_vec_char( vector<char>& v, const std::string &s );
    void to_vec_char( vector<char>& v, const char *s );
    // return true and set found_h to true if input string s is the help string h, 
    // followed by an = sign or nothing
    bool check_optionstring( const char *s, const char *h, bool &found_h, std::string &value_h, std::string default_value = category_name[COMMON] );
    // return true if any help option was found
    bool check_help_options(int argc, char **argv, 
                           bool &found_print_version, bool &found_help_list, bool &found_help, bool &found_help_full); 
    // for printing words lined up in columns
    void print_in_columns(const std::vector<std::string> &words) const;
    void print_help(std::vector< std::string > &help_strings, HelpType help_type) const;
    int do_help(ParserSystem &ps, HelpType help_type);
    // set the default operation and output files, if needed
    void set_operation_defaults(int argc, char **argv); 

    // non-trivial initializations that aren't safe to put in the constructor
    void init(int argc, char **argv);
    int checkInput(const std::string &option, const std::string &value, const std::string &allowed_values);
    void fill_help_strings(int argc, char **argv);
    // determine default output_mesh, previous and next meshes based on input mesh name
    void default_output_files(std::string input_mesh_s, std::string operation, int number_refines,
                              std::string &output_mesh_s, std::string &previous_adapted_mesh_s, std::string &next_adapted_mesh_s);
    // true if there is an error
    bool set_input_geometry_type();
    // true if there is an error
    bool check_parsing_results();
    // returns index to the "operation" keyword, e.g. "adapt". Returns 0 if it wasn't simple syntax
    int check_for_simple_options(int argc, char **argv);
    // true if simple options were specified
    int adapt_main_simple_options(int argc, char **argv, bool &parsing_ok, bool &do_adapt);

    // what to print when a parsing error occurred, generic syntax help
    void print_syntax_error_instructions(bool generic_syntax_error = true);
    // clean up (finalize) libraries and exit
    void exit_safely(int exit_code); 

    struct EntitySelectorUCF {
      PerceptMesh& m_eMesh;
      int m_dim;
      EntitySelectorUCF(PerceptMesh& eMesh) : m_eMesh(eMesh), m_dim(m_eMesh.get_spatial_dim()) {}

      bool operator()(stk::mesh::Entity e_from, stk::mesh::Entity e_to)
      {
        VERIFY_OP_ON(m_eMesh.m_unprojected_coordinates, !=, 0, "bad m_unprojected_coordinates");
        double *c_from = static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , e_from ));
        double *uc_to = stk::mesh::field_data( *m_eMesh.m_unprojected_coordinates , e_to );
        if (std::fabs(uc_to[3]) < 1.e-6)
          {
            for (int j=0; j < m_dim; j++)
              {
                uc_to[j] = c_from[j];
              }
            if (m_dim != 3)
              uc_to[2] = 0.0;
            uc_to[3] = 1.0;
          }
        // don't do further processing
        return false;
      }

    };

    // helpers (could be made virtual and specialized by derived,
    //  for now, we use if-stmts for special cases)
    void create_refine_pattern();
    void create_refiner();
    void pre_refiner_tasks(int iBreak);
    int  bucket_in_block_names(stk::mesh::Bucket& bucket, std::string& block);

    // utilities
    void setup_m2g_parts(std::string input_geometry);
    void initialize_m2g_geometry(std::string input_geometry);

    // sub-algorithms of do_run algorithms
    BlockNamesType process_block_names();
    void mesh_based_geometry_setup();
    void mesh_based_geometry_fitting();

    void verify_mesh_util(bool isInit);
    void verify_mesh_init();
    void verify_mesh_final();
    void do_histograms_init();
    void do_histograms_final();
    void compute_hmesh_sizes_init();
    void compute_hmesh_sizes_final();
    void run_adapt_run();
    void do_dihedral_angle_check();

    // do_run algorithms
    void pre_open();
    int do_run_pre_commit();
    int do_run_post_commit();
    int do_post_proc();

    // version - return true if built in Sierra
    bool get_version(std::string* v=0);

    void log_usage( bool status = true );

    // main routines
    int adapt_main(int argc, char **argv); // parse then adapt
    int adapt_main_full_options(int argc, char **argv, bool &parsing_ok, bool &do_adapt); // parsing
    int adapt_main_do_adapt(); // adaptation, after all parsing is done and info has been transfered to class variables
    int main(int argc, char **argv);

    std::shared_ptr<FitGregoryPatches> fitter; // 3D
#if HAVE_OPENNURBS
    std::shared_ptr<GeometryRecoverySplineFit> grsf; // 2D
#endif

    std::shared_ptr<UniformRefinerPatternBase> pattern;
    std::shared_ptr<Refiner> refiner;

    stk::mesh::FieldBase* proc_rank_field_ptr = 0;

    std::shared_ptr<PerceptMesh> eMeshP;
    std::shared_ptr<AdaptedMeshVerifier> adaptedMeshVerifier;
    PGeom * m_PGeomPntr = NULL;

    BlockNamesType m_block_names;
    stk::mesh::Selector block_selector;

    std::shared_ptr<ElementRefinePredicate> element_refine_predicate;
    std::shared_ptr<stk::mesh::Selector> univ_selector;

    typedef std::map<std::string, int> StringIntMap;
    StringIntMap block_names_x_map;
    std::shared_ptr<DihedralAngleCheck> m_dihedral_angle_check;

  private:
    stk::ParallelMachine m_comm;
    MemoryInfo m_meminfo;
 
    stk::diag::TimerSet m_timerSet;
    stk::diag::Timer m_timer;
 };


}

