// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


// Note: this file is included in MeshAdapt.hpp and .cpp to provide a path to using c++11 member var init

// options start...
// ====================================================================================================================================
std::string options_description_desc = "adapt options";

// NOTE: Options --directory --output-log --runtest are handled/defined in RunEnvironment
std::string input_mesh  = "";
std::string input_geometry  = "";
#if IN_MESH_ADAPT
GeometryType input_geometry_type = GEOM_NONE;
#endif
std::string fit_geometry_file  = "";
double fit_angle = 115.0;
std::string fit_3d_file = "";
std::string output_mesh  = "";
std::string block_name_inc = "";
std::string block_name_exc = "";
int use_transition_elements = 0;
std::string previous_adapted_mesh  = "";
std::string next_adapted_mesh  = "";
std::string memory_logfile_name  = "";

// for Salinas
#if defined(STK_BUILT_FOR_SIERRA)
std::string rbar_blocks = ""; 
#endif
// for Salinas and other codes
//std::string ignore_blocks = "";
// just use block_name_inc to exclude....

#if defined(STK_BUILT_FOR_SIERRA)
std::string version_prefix = "Sierra_";
#else
std::string version_prefix = "NonSierra_";
#endif

std::string version = "1.0";
int print_version = 0;

// Default values are reasonable values to use if the corresponding operation is specified.
// A default value for "enrich" does not mean that "enrich" will be done by default.
std::string convert  = ""; // depends on element type
std::string refine  = "DEFAULT";
std::string enrich  = "DEFAULT";
int generated_mesh = 1;

int dihedral_angle_check = 0;
int skin_mesh = 0;
int output_active_elements_only = 0;
std::string verify_meshes = "0";
int number_refines = 1;
int progress_meter = 0;
int print_timers = 0;
int smooth_geometry = 0;
int smoother_niter = 1000;
double smoother_tol = 1.e-4;
std::string adapt_input = "";
int smooth_use_reference_mesh = 1;
int fix_all_block_boundaries = 0;
std::string ioss_write_options = "";
std::string ioss_read_options = "";
#if !defined(NO_GEOM_SUPPORT)
int respect_spacing = 1;
#endif
int smooth_surfaces = 0;
//double min_spacing_factor = 0.25; // range [0,0.5]
int remove_geometry_blocks = 0;
int dump_geometry_file = 0;
int sync_io_regions = 1;
std::string property_map = "";
std::string compute_hmesh = "";
int print_hmesh_surface_normal = 0;
int m_excludeWedgesNotConnectedToPyramids = 0;

double hmesh_factor = 0.0;
double hmesh_min_max_ave_factor[3] = {0,  0, 0};
std::string histograms_root = "cout";
//std::string histogram_options = "mesh: [edge_length, quality_edge, quality_vol_edge_ratio, volume]";
std::string histogram_options = "";

#if IN_MESH_ADAPT
std::string convert_options = UniformRefinerPatternBase::s_convert_options;
std::string refine_options  = UniformRefinerPatternBase::s_refine_options;
std::string enrich_options  = UniformRefinerPatternBase::s_enrich_options;
#endif

// : if not specified, use input mesh name appended with _{converted,refined_#refines,enriched}");

std::string help_list = "COMMON";
std::string help = "COMMON";
std::string help_full = "COMMON";
// ====================================================================================================================================
// options end
#ifndef NDEBUG
int debug = 1;
#else
int debug = 0;
#endif

int s_spatialDim = 0;

double t0   = 0.0;
double t1   = 0.0;
double cpu0 = 0.0;
double cpu1 = 0.0;

int i_pass = 0;
int m_M = 1;
int m_W = 1;
int m_iW = 0;
int m_M_0 = 0;
int m_M_1 = 0;
int m_iM = 0;

std::string histogram_basic_options  = "{file_root: " + histograms_root + ", mesh: [edge_length, quality_edge, quality_vol_edge_ratio, volume] }";
int result = 0;
unsigned failed_proc_rank = 0u;

unsigned p_rank = 0;
unsigned p_size = 0;

double mMeshInputTime       = 0.0;
double mMeshOutputTime      = 0.0;
double mAdaptTimeOverall    = 0.0;
double mAdaptCPUTimeOverall = 0.0;

// general help stings, that are not specific to a particular command option
const std::string syntax_error = "\nCommand line syntax error: ";
const std::string operation_error = "exactly one operation and its value must be specified: --refine, --enrich, --convert (or --generated_mesh).";
const std::string default_operation = "if no operation or value is specified, the defaults are --refine=DEFAULT";
const std::string syntax_example = "e.g., 'mesh_adapt --input_mesh=myfile.g' is the same as 'mesh_adapt --refine=DEFAULT --input_mesh=myfile.g'";
const std::string simple_operation_error = "exactly one simple operation must be specified : adapt, refine, enrich, or convert.";
// non-trivial, filled in by init
std::string help_instructions; // "Run 'mesh_adapt --help' for help.";
std::string simple_std_help;     // standard syntax requirements
std::string simple_usage_string; // simple syntax
const std::string simple_usage_string_one_line = "mesh_adapt  {convert, refine, enrich, adapt}  input_mesh_filename  output_mesh_filename";
const std::string help_hint = "For a list of all options, try '--help_list=ALL'";
