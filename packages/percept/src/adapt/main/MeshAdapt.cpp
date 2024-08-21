// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <numeric>
#include <sys/stat.h>

#include <adapt/main/MeshAdapt.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <adapt/AdaptHelperFunctions.hpp>
#include <adapt/FixSideSets.hpp>
#include <adapt/ExcludeWedgesNotConnectedToPyramids.hpp>

#include <percept/PerceptUtils.hpp>

#include <stk_util/registry/ProductRegistry.hpp>
#include <stk_util/util/string_case_compare.hpp>

#include <sys/ioctl.h>

#if defined(STK_BUILT_FOR_SIERRA)
#include <AuditLog.h>
#endif

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <percept/mesh/geometry/kernel/GeometryKernelGregoryPatch.hpp>
#include <percept/mesh/geometry/kernel/GeometryKernelOpenNURBS.hpp>

#ifdef HAVE_CUBIT
#include "PGeom.hpp"
#include "PGeomAssoc.hpp"
#ifdef HAVE_ACIS
#include "PGeomACIS.hpp"
#endif

#include <percept/mesh/geometry/kernel/GeometryKernelPGEOM.hpp>
#endif

#include <percept/mesh/geometry/kernel/MeshGeometry.hpp>
#include <percept/mesh/geometry/kernel/GeometryFactory.hpp>

#endif

namespace percept {

  MeshAdapt::MeshAdapt() :
    m_timerSet(sierra::Diag::TIMER_ALL),
    m_timer(stk::diag::createRootTimer("MeshAdapt", m_timerSet))
  {
    m_timer.start();
  }

  void print_timer_table(stk::diag::Timer &timer)
  {
    std::ostringstream str;
    timer.stop();
    stk::diag::printTimersTable(str, timer,
                                stk::diag::METRICS_ALL, false,
                                MPI_COMM_WORLD);
    if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
      {
        std::cout << str.str() << std::endl;
      }
  }

  // Determine the basename and refinement level from a given an ExodusII file name.
  // This is the last integer before the final '.e' or '.exo' in the file name.
  // If no level was found, returns 0 as the level.
  //
  // Return:
  //    string - the basename which is everything before the number
  //    int    - the refinement level
  std::pair<std::string,int> MeshAdapt::get_basename_and_level(const std::string& input_mesh)
  {
    // find the last number before the first dot character
    size_t dot_index = input_mesh.find_last_of(".");
    std::string prefix(input_mesh, 0, dot_index);
    size_t last_index = prefix.find_last_not_of("0123456789");

    int level = 0; // level will be one if we don't find a level
    if (last_index + 1 < dot_index) {
      std::string level_str(prefix, last_index + 1, std::string::npos);
      level = std::stoi(level_str);
    }
    std::string basename(prefix, 0, last_index + 1);

    return std::make_pair(basename, level);
  }

  int MeshAdapt::checkInput(const std::string &option, const std::string &value, const std::string &allowed_values)
  {
    std::vector<std::string> vals = Util::split(allowed_values, ", ");
    for (auto &s : vals)
    {
      if (s == value)
        return 0;
    }
    
    if (p_rank==0) std::cout << "\n" << syntax_error
                             << "bad value for option, " << option << "=" << value 
                             << "\n  allowed values = " << allowed_values << std::endl;
    return 1;
  }

  void MeshAdapt::exit_safely(int exit_code)
  {
#if defined( STK_HAS_MPI )
    Kokkos::finalize();
    stk::parallel_machine_finalize();
#endif
  }      

  void MeshAdapt::fill_help_strings(int argc, char **argv)
  {
    std::stringstream ss;
    ss << 
      "Simple syntax:\n" << 
      "       " << argv[0] << " [refine|enrich] input_mesh_filename  [output_mesh_filename] [number_levels]\n"
      "       " << argv[0] << " convert input_mesh_filename [output_mesh_filename]\n"
      "       " << argv[0] << " adapt input_mesh_filename [analysis_results_filename]\n";
    simple_usage_string = ss.str();

    simple_std_help = "\n" + simple_usage_string + "\n" + "Standard syntax requirements:\n    " + operation_error + 
                      "\n    " + default_operation +
                      "\n    " + syntax_example;

    help_instructions = "Run '" + std::string(argv[0]) + " --help' for help.";
  }

  int MeshAdapt::check_for_simple_options(int argc, char **argv)
  {
    for (int i = 1; i < argc; i++)
      {
        char *s = argv[i];
        if (strcmp( s, "adapt")  == 0 ||
            strcmp( s, "refine") == 0 ||
            strcmp( s, "enrich") == 0 ||
            strcmp( s, "convert")== 0 )
        {
          return i;
        }
      }
    return 0;
  }

  void MeshAdapt::to_vec_char( vector<char>& v, const char *s )
  {
    v.clear();
    while ( *s != '\0' )
      v.push_back( *s++ );
    v.push_back( '\0' );
  }
  void MeshAdapt::to_vec_char( vector<char>& v, const std::string &s )
  {
    v.clear(); 
    v.reserve(s.size());
    std::copy( s.begin(), s.end(), std::back_inserter(v) );
    v.push_back('\0');
  }

  void MeshAdapt::default_output_files(std::string input_mesh_s, std::string operation, int number_refines,
    std::string &output_mesh_s, std::string &previous_adapted_mesh_s, std::string &next_adapted_mesh_s)
  {

    output_mesh_s = input_mesh_s;
    previous_adapted_mesh_s.clear();
    next_adapted_mesh_s.clear();

    // auto-generate next mesh file name, based on input mesh file name
    size_t dot_index = input_mesh_s.find_last_of(".");
    std::string suffix = input_mesh_s.substr(dot_index);

    // set previous and next adapted mesh, even if operation==refine, because maybe rar-info was given
    std::pair<std::string,int> basename_and_level = get_basename_and_level(input_mesh_s);
    std::string basename(basename_and_level.first);
    int input_level = basename_and_level.second;
    if (input_level > 0) 
      previous_adapted_mesh_s = basename + std::to_string(input_level) + "_ft" + suffix;
    // else leave previous_adapted_mesh="";
    std::string next_level(std::to_string(input_level + 1));
    next_adapted_mesh_s = basename + next_level + "_ft" + suffix;

    if (operation == "adapt") 
    {      
      output_mesh_s = basename + next_level + suffix;
    }
    else 
    {
      // add correct-English suffix adapt->_adapted, refine->_refined3, etc.
      bool ends_in_e = !operation.empty() && operation[operation.length()-1] == 'e';
      std::string new_suffix = "_" + operation + (ends_in_e ? "d" : "ed" );
      if ( operation == "refine" )
        new_suffix += "_" + std::to_string(number_refines);
      new_suffix += suffix;

      Util::replace(output_mesh_s, suffix, new_suffix);
    }

  }

  bool MeshAdapt::is_number(std::string &s, int &i)
  {
    bool isInt = false;
    {
      try {
        i = std::stoi(s);
        (void)i;
        isInt = true;
      }
      catch( ... )
      {
      }
      if (!isInt)
        i = 1;
    }
    return isInt;
  }

  void find_next_nonoption( int argc, char **argv, int &i )
  {
    while ( i < argc && Util::startsWith(argv[i], "--") )
      ++i;
  }

  int MeshAdapt::adapt_main_simple_options(int argc, char **argv, bool &parsing_ok, bool &do_adapt)
  {    
    parsing_ok = false;
    do_adapt = false;

    int simple_options_index = check_for_simple_options(argc, argv);
    if (simple_options_index == 0)
      return 0; // simple options were not specified, try standard syntax

    // simple options were specified, return 1 later

    // check_for_simple_options above ensures option is one of "adapt", "refine", "enrich" or "convert"
    std::string option = argv[simple_options_index];

    // too few or too many options?
    // we return 1 because we found the simple syntax
    // we set parsing_ok = false to trigger an error message, and set do_adapt = false

    // input file
    int next_i = simple_options_index+1;
    find_next_nonoption( argc, argv, next_i );
    std::string input_mesh_s    ( argc > next_i ? argv[next_i++] : std::string() );
    // next option could be output_file, analysis_results_filename, or number_levels
    // we could have both output_file and number_levels, or only one of them
    find_next_nonoption( argc, argv, next_i );
    std::string output_mesh_s   ( argc > next_i ? argv[next_i++] : std::string() );
    find_next_nonoption( argc, argv, next_i );
    std::string number_refines_s( argc > next_i ? argv[next_i++] : std::string() );

    if (input_mesh_s.empty())
    {
      if (p_rank==0) std::cout << "Error: simple syntax \"" << option << "\" was specified, but the input mesh was missing." << std::endl;
      return 1; 
    }

    // if only one option was specified, was it output_mesh or number_refines?
    int nref(1); // integer version of number_refines
    {
      bool number_refines_is_number = is_number(number_refines_s,nref);
      // if number_refines_is_number, then leave it as is
      if (!number_refines_is_number)
      {
        int num2;
        bool output_mesh_is_number = is_number(output_mesh_s,num2);
        // the wrong one was a number, swap them
        if (output_mesh_is_number)
        {
          nref = num2;
          output_mesh_s.swap(number_refines_s);
        }
        // neither was a number, default number_refines ="1"
        else 
        {
          nref = 1;
          number_refines_s = "1";
        }
      }
    }

    // generate default names
    std::string def_output_mesh_s;
    default_output_files(input_mesh_s, option, nref, def_output_mesh_s, previous_adapted_mesh, next_adapted_mesh);
    if (option == "adapt") {

      // set input mesh name to analysis mesh, if it was specified
      // i.e. interpret output_mesh as the analysis mesh, not the output mesh
      if (!output_mesh_s.empty())
        input_mesh_s = output_mesh_s; // ??? should we clear the output_mesh_s or replace it with the default??

      // use the default output mesh name, based on the input mesh, not the analysis mesh
      output_mesh_s = def_output_mesh_s;
      if (p_rank==0) std::cout << "defaulting to output_mesh_filename \"" << output_mesh_s << "\"" << std::endl;

      number_refines_s = "1";
      nref = 1;
    }
    else 
    {
      if (output_mesh_s.empty())
      {
        output_mesh_s = def_output_mesh_s;
        if (p_rank==0) std::cout << "defaulting to output_mesh_filename \"" << output_mesh_s << "\"" << std::endl;
      }
    }


    // build a new set of normal options, based on interpretation of the simple options

    // create non-const char* for the converted-sytax options
    // ok to use these until the vectors go out of scope
    std::vector<char> in_mesh;  to_vec_char( in_mesh, "--input_mesh=" + input_mesh_s );

    std::vector<char> out_mesh; to_vec_char( out_mesh, "--output_mesh=" + output_mesh_s );

    std::vector<char> num_refines; to_vec_char( num_refines, "--number_refines=" + number_refines_s );

    std::vector<char> ref;        to_vec_char(ref, "--refine=DEFAULT");
    std::vector<char> spacing;    to_vec_char(spacing, "--respect_spacing=0");
    std::vector<char> rar;        to_vec_char(rar,  "--RAR_info=\"adapt.yaml\"");
    std::vector<char> rich;       to_vec_char(rich, "--enrich=DEFAULT");
    std::vector<char> cvrt;       to_vec_char(cvrt, "--convert=Hex8_Tet4_24");
    std::vector<char> ioss_read;  to_vec_char(ioss_read,   "--ioss_read_options=\"large\"");
    std::vector<char> ioss_write; to_vec_char(ioss_write, "--ioss_write_options=\"large\"");

    std::vector<char *> argv_new;
    argv_new.reserve(argc + 11); 
    argv_new.push_back( argv[0] );

    // ok to use until the above vectors go out of scope
    argv_new.push_back( &in_mesh[0] ); 
    argv_new.push_back( &out_mesh[0] ); 
    if (option == "refine")
      argv_new.push_back( &ref[0] );
    else if (option == "adapt") {
      argv_new.push_back( &ref[0] );
      argv_new.push_back( &spacing[0] );
      argv_new.push_back( &rar[0] );
    }
    else if (option == "enrich")
      argv_new.push_back( &rich[0] );
    else
      argv_new.push_back( &cvrt[0] );
    argv_new.push_back( &num_refines[0] );
    argv_new.push_back( &ioss_read[0] );
    argv_new.push_back( &ioss_write[0] );

    // push all standard (non-simple) options *after* so they take precedence over the auto-derived ones
    for (int i = 1; i < argc; ++i)
    {
      if (Util::startsWith(argv[i], "--"))
        argv_new.push_back(argv[i]);
    }


    // if we got to here, then we found simple options and should return 1,
    adapt_main_full_options(argv_new.size(), &argv_new[0], parsing_ok, do_adapt);

    if (!parsing_ok)
      if (p_rank==0) std::cout << syntax_error << "malformed simple syntax for \"" << option << "\"" << std::endl;

    return 1; 
  }

  void MeshAdapt::print_syntax_error_instructions(bool generic_syntax_error)
  {
    if (generic_syntax_error) 
      if (p_rank==0) std::cout << "\n" << syntax_error << "\n";
    if (p_rank==0) std::cout << help_instructions << "\n" << std::endl;
  }

  int MeshAdapt::adapt_main(int argc, char **argv)
  {
    bool parsing_ok = false;
    bool do_adapt = false;
    try
    {
      init(argc, argv);
      if (!adapt_main_simple_options(argc, argv, parsing_ok, do_adapt))
        adapt_main_full_options(argc, argv, parsing_ok, do_adapt);
    }
    catch(...)
    {
      print_syntax_error_instructions();
      return 1;        
    }
    
    if (!parsing_ok)
    {
      print_syntax_error_instructions();
      return 1;                
    }
    if (do_adapt)
    {
      // call actual mesh adaptation
      return adapt_main_do_adapt();
    }
    return 0; // all is well, we printed some sort of help
  }

  bool MeshAdapt::check_optionstring( const char *s, const char *h, bool &found_h, std::string &value_h, std::string default_value )
  {
    // found_h = false; don't set it to false if it is already true
    bool found_this_time = false;
    const int h_len = strlen(h);
    if (strncmp(s,h,h_len) == 0)
    {
      if ( s[h_len] == '\0' || s[h_len] == '=' )
      {
        found_h = true;
        found_this_time = true;
        if (s[h_len] == '=')        
          value_h = &s[h_len+1]; // copy string from character after = to end
        // else, leave value_h at its default value, whatever that is

        // if value_h is empty, regardless of why it is empty, replace it with COMMON
        if (value_h.empty())
          value_h = default_value;
      }
    }

    return found_this_time;
  }

  bool match_prefix( const std::string &s, const std::string &pre )
  {
    const auto len = pre.length();
    return ( s.length()>=len && s.substr(0,len) == pre );
  }

  bool MeshAdapt::check_help_options(int argc, char **argv, 
    bool &found_print_version, bool &found_help_list, bool &found_help, bool &found_help_full) 
  {

    // if no options were specified, treat it as --help
    if (argc<=1)
    {
      found_help = true;
      return true;
    }


    for ( int i = 1; i < argc; ++i ) // skip 0
    {
      char *s = argv[i];
      if (strcmp(s,"--version") == 0 || strcmp(s,"-v") == 0 )   found_print_version = true;
      else if ( check_optionstring( s, "--help_full", found_help_full, help_full ) ) {}
      else if ( check_optionstring( s, "--help_list", found_help_list, help_list ) ) {}
      else if ( check_optionstring( s, "--help", found_help, help ) ) {}
      // malformed help option?, check this last
      else 
      {
      // copy the part of s beore the '=', if any
        std::string ss(s);
        size_t equal_pos = ss.find( '=' );
        if (equal_pos != std::string::npos)
          ss.resize(equal_pos);
        if (!ss.empty())
        {
      // make ss lower case
          transform(ss.begin(), ss.end(), ss.begin(), ::tolower);
          if (ss == "h" || ss == "-h" || ss == "--h" || match_prefix(ss, "help") || match_prefix( ss, "-help" ) || match_prefix(ss, "--help") )
          {
            found_help = true;
            if (help.empty()) help = category_name[COMMON];
            if (p_rank==0) std::cout << "\nInterpreting malformed input option '" << s << "' as '--help=" << help << "'" << std::endl;
            if (p_rank==0) std::cout << help_hint << "\n" << std::endl;
          }
        }
      }
    }
    return found_print_version || found_help_list   || found_help || found_help_full;
  }

  void MeshAdapt::set_operation_defaults(int argc, char **argv)
  {
    // check if convert, enrich, or refine was explicitly mentioned in input
    // convert, enrich, refine, generated_mesh may have default values, so we can't just check their string length after parsing
    // but don't bother with their default values, we already set those
    bool found_convert(false), found_enrich(false), found_refine(false), found_generate(false);
    std::string junk;
    for (int i = 1; i < argc; ++i)
    {
      char *s = argv[i];
      check_optionstring(s,"--convert",        found_convert,  junk);
      check_optionstring(s,"--enrich",         found_enrich,   junk);
      check_optionstring(s,"--refine",         found_refine,   junk);
      check_optionstring(s,"--generated_mesh", found_generate, junk);
    }

    // if no operation, then default to refine
    if ( !(found_convert || found_enrich || found_refine || found_generate) )
    {
      if (p_rank==0) std::cout << "defaulting to --refine" << std::endl;
      found_refine = true;
    }

    // if more than one operation found, remove all but one
    if (found_refine)   { found_enrich=false; found_convert=false; found_generate=false;}
    if (found_enrich)   { found_refine=false; found_convert=false; found_generate=false;}
    if (found_convert)  { found_enrich=false; found_refine=false; found_generate=false;}
    if (found_generate) { found_enrich=false; found_convert=false; found_refine=false;}

    // remove the default values of the non-selected operations
    if (!found_convert ) convert.clear();
    if (!found_enrich  ) enrich .clear();
    if (!found_refine  ) refine .clear();
    if (!found_generate) generated_mesh =0;

    // reset the default values if Teuchos command line processor cleared them.
    // e.g. --refine with no equals or value clears the default value!
    if (found_refine   && refine  .empty()) {refine   = "DEFAULT";} 
    if (found_enrich   && enrich  .empty()) {enrich   = "DEFAULT";} 
    // if (found_convert  && convert .empty()) {convert  = "DEFAULT";}  // no reasonable default value
    if (found_generate && generated_mesh <= 0) {generated_mesh = 0;}

    // set default output file names, provided there was an input file
    const bool err_input_mesh = input_mesh.empty();

    // determine default outputs based on inputs
    if (!err_input_mesh)
    {
      // set the operation string, for default file names
      std::string operation;
      if (found_convert ) operation = "convert";
      if (found_enrich  ) operation = "enrich";
      if (found_refine  ) operation = "refine";
      if (found_generate) operation = "generate";

      std::string def_output_mesh, def_previous_adapted_mesh, def_next_adapted_mesh;
      default_output_files(input_mesh, operation, number_refines,
                           def_output_mesh, def_previous_adapted_mesh, def_next_adapted_mesh);

      // If prev and next were set by the simple options, use those values, 
      // because in the adapt case with an analysis file, they are different.
      // Otherwise use the default values generated here.
      if (previous_adapted_mesh.empty()) previous_adapted_mesh = def_previous_adapted_mesh;
      if (next_adapted_mesh.empty()) next_adapted_mesh = def_next_adapted_mesh;
      // replace unspecified strings with defaults, and let the user know!
      if (output_mesh.empty() && !def_output_mesh.empty()) 
      {
        output_mesh = def_output_mesh;
        if (p_rank==0) std::cout << "defaulting to --output_mesh=" << output_mesh << std::endl;
      }
    }
  }

  bool MeshAdapt::set_input_geometry_type()
  {

    std::string error_string;

    // check geometry file extensions
    if (input_geometry.length())
    {

#if defined( STK_HAS_MPI )
      MPI_Barrier( MPI_COMM_WORLD );
#endif

      // 3dm == opennurbs
      if (input_geometry.find(".3dm") != std::string::npos )
      {
        // could be pgeom or regular opennurbs
        // replace extension with "m2g" and see if file exists
        std::string m2gFile = input_geometry.substr(0,input_geometry.length()-3) + "m2g";
        struct stat s;
        if (0 == stat(m2gFile.c_str(), &s))
        {
          input_geometry_type = PGEOM_OPENNURBS;
#ifndef HAVE_CUBIT
          error_string = "\nERROR: CUBIT is needed for PGeom OpenNURBS, but this mesh_adapt was not built with CUBIT.";
#endif              
        }
        else
        {
          input_geometry_type = OPENNURBS;
        }
      }
      else if (has_suffix( input_geometry, ".g"  ) || 
               has_suffix( input_geometry, ".e"  ) || 
               has_suffix( input_geometry, ".exo") )
      {
        input_geometry_type = MESH_BASED_GREGORY_PATCH;
      }
      else if (has_suffix( input_geometry, ".sat"))
      {
        input_geometry_type = PGEOM_ACIS;
#ifndef HAVE_ACIS
        error_string = "\nERROR: ACIS is needed for PGeom ACIS files, but this mesh_adapt was not built with ACIS.";
#endif
      }
      else
      {
        error_string = "\n ERROR: Unrecognized geometry type.\n" 
        "Recognized file extensions are "
        ".3dm (OpenNURBS); "
        ".3dm.m2g (PGeom OpenNURBS); "
        "{.g .e .exo} (Exodus mesh based Gregory patch); and "
        ".sat (ACIS)";
      }

      // print error message, if any
      if (p_rank == 0)
      {
        // debug / diagnostic
        if (p_rank==0) std::cout << "mesh_adapt: --input_geometry=" << input_geometry << " type is ";
        switch (input_geometry_type)
        {
          case PGEOM_OPENNURBS:
          if (p_rank==0) std::cout << "OpenNURBS (PGeom)";
          break;
          case OPENNURBS:
          if (p_rank==0) std::cout << "OpenNURBS";
          break;
          case MESH_BASED_GREGORY_PATCH:
          if (p_rank==0) std::cout << " mesh-based (Gregory patch)";
          break;
          case PGEOM_ACIS:
          if (p_rank==0) std::cout << "ACIS";
          break;
          default:
          case GEOM_NONE:
          case N_GeometryType:
          if (p_rank==0) std::cout << "undefined";
          break;
        }
        // error, if any
        if (p_rank==0) std::cout << error_string << std::endl;
      }
    }
    return !error_string.empty();
  }

  bool MeshAdapt::check_parsing_results()
  {
    const bool no_operation = (convert.length()+enrich.length()+refine.length()+generated_mesh == 0);

    const bool bad_operation = 
      (convert.length() && checkInput("convert", convert, convert_options)) ||
      ( enrich.length() && checkInput("enrich",  enrich,  enrich_options)) ||
      ( refine.length() && checkInput("refine",  refine,  refine_options)) ; 
      // todo, check input_mesh string when generated_mesh=1
      // todo, check if file exists when generated_mesh=0
    const bool err_input_mesh  =  input_mesh.empty(); 
    const bool err_output_mesh = output_mesh.empty();

    const std::string verify_options("0, 1, 2, FV");
    bool err_verify = (!verify_meshes.empty() && checkInput("verify_meshes", verify_meshes, verify_options));

    if ( no_operation || bad_operation || err_input_mesh || err_output_mesh || err_verify )
    {
      // bad operation already prints syntax_error, otherwise we want to print it now
      if (!bad_operation) 
        if (p_rank==0) std::cout << "\n" << syntax_error;
      if ( err_input_mesh )
        if (p_rank==0) std::cout << "--input_mesh must be specified. ";
      if ( err_output_mesh )
        if (p_rank==0) std::cout << "--output_mesh was not specified or derived. ";
      if ( no_operation )
      {
        if (p_rank==0) std::cout << "\nEither " << operation_error << "\nOr " << simple_operation_error << "\n";
      }
      if (p_rank==0) std::cout << std::endl;
      return true;
    }
    return false;

  }
  
  int MeshAdapt::adapt_main_full_options(int argc, char **argv, bool &parsing_ok, bool &do_adapt)
  {
    ParserSystem ps;
    parsing_ok = false;
    do_adapt = false;

    // check for help ( or version )
    bool found_print_version(false), found_help_list(false), found_help(false), found_help_full(false);
    bool found_any_help = check_help_options(argc, argv, found_print_version, found_help_list, found_help, found_help_full);
    if (found_any_help)
    {
      parsing_ok = true;
      do_adapt = false;
    }

    if (found_print_version)
    {
      std::string v;
      bool vb = get_version(&v);
      if (p_rank == 0) std::cout << "Version: " << v << std::endl;
      return ( vb ? 0 : 1 ); // true means everything was fine, a version was given
    }

    // create the parser parameter list
    //================================== basic options ===========================================================
    // help
    // the first three options take a keyword value, the others don't
    const std::string keyword_options = ", for {SIMPLE, COMMON, ALL}, option type (e.g. IO), or command (e.g. input_mesh)";
    std::string keyword_options_ext="option types are {";
    for ( int i = 0; i < num_operations-1; ++i)
    {
      keyword_options_ext += operation_name[i] + ", ";
    }
    keyword_options_ext += operation_name[num_operations-1] +"}";

    // always give the full help for "help," to help the user get started
    set_option(ps, "help"                       , &help              , HELP, "short one-line help" + keyword_options + 
                              "\n                                                 "  // indent, exceptional use of two lines for ONELINE help
                              + keyword_options_ext, "", "", true, false);
    set_option(ps, "help_list"                  , &help_list         , HELP, "list command options", "", "", false, false);
    set_option(ps, "help_full"                  , &help_full         , HELP, "long help, including syntax options", "", "", false, false);
    set_option(ps, "version"                    , &print_version     , HELP,  "print version and exit", "", "", false, true);

    // files
    set_option(ps, "input_mesh"               , &input_mesh               , IO, "input mesh name", "", "", false, false);
    set_option(ps, "output_mesh"              , &output_mesh              , IO, "output mesh name", "", "", false, false);

    // operation type
    set_option(ps, "convert"                  , &convert                  , OPERATION, 
      "split each element into several of a new type", "", "\n\tOptions are " + convert_options);

    set_option(ps, "enrich"                   , &enrich                   , OPERATION, 
      "add nodes to make elements higher order", ". E.g. DEFAULT",  "\n\tOptions are " + enrich_options);

    set_option(ps, "refine"                   , &refine                   ,  OPERATION, "split elements", "", "\n\tOptions are " + refine_options, true, false);

    set_option(ps, "number_refines"           , &number_refines           , PARAMETER, "number of refinement passes", "",  ".  Must be >= N in --blocks=b:Nx", true, false);
    // subsetting
    
    set_option(ps, "blocks"                   , &block_name_inc           , PARAMETER,  
      "blocks to operate on, and block-specific refinement levels", "",
      "\n"
      "\t(0) by default, all blocks are included\n"
      "\t(1) [+]block_<n> will include that block, e.g., block_2\n"
      "\t(2) -block_<n> will exclude that block, e.g., -block_3\n"
      "\t(3) file:<filename> reads values from the file, e.g., file:blocklist.dat\n"
      "\t(4) block_1..block_10 specifies the range of blocks 1 to 10\n"
      "\t(5) :Nx (or :NX) suffix, where N specifies the number of times to refine\n"
      "\t    e.g. --blocks=1..3:2x,5,6:3x refines blocks 1,2 and 3 twice, block 5 'number_refines' times,\n"
      "\t         and block 6 three times. Block 4 would not be refined at all.\n"
      "\tNote any combination of + and - and ranges can be specified.\n"
      "\tNote `block_` is optional, so 'block_1' is equivalent to '1'. The + sign is optional.");


    // spacing
#if !defined(NO_GEOM_SUPPORT)
    set_option(ps, "respect_spacing"          , &respect_spacing          , PARAMETER, "preserve the input-mesh size-gradation during refinement", "", "");
#endif

    // ==== end basic options

    //================================== advanced options ===========================================================
    // subsetting
    set_option(ps, "use_transition_elements"  , &use_transition_elements  ,
                                  PARAMETER, "create a conformal mesh", "", 
                                  "\n\tUse transition elements between adjacent more-refined and less-refined blocks");

    // ioss options -  we'd like to change this to set automatically, based on other criteria
    set_option(ps, "ioss_read_options"        , &ioss_read_options        , IO, "input mesh IOSS/Exodus file layout", "",
      "\n\tioss_read_options=\"{large, auto-join:<yes | no>\" e.g. \"auto-join:no\""); 
    set_option(ps, "ioss_write_options"       , &ioss_write_options       , IO, "output mesh IOSS/Exodus file layout", "",
      "\n\tioss_write_options=\"{large, auto-decomp:<yes | no>}\" e.g. \"large,auto-decomp:yes\"");

    // run-adapt-run
    set_option(ps, "RAR_info"                 , &adapt_input              , IO, "name of input file for advanced usage");

    // query/control
    set_option(ps, "verify_meshes"            , &verify_meshes            , INFO, "verify positive volumes for input and output meshes", "",
      "\n\tset to 1 for finite element volume checks"
      "\n\tset to FV for finite volume checks");
    set_option(ps, "progress_meter"           , &progress_meter           , INFO, "progress meter on or off");
    set_option(ps, "print_timers"             , &print_timers             , INFO, "print more detailed timing info");
    set_option(ps, "skin_mesh"                , &skin_mesh                , PARAMETER, "produce a boundary sideset for any exposed boundaries");
    set_option(ps, "dihedral_angle_check"     , &dihedral_angle_check     , INFO, "list obtuse dihedral angles", "", 
      "\n\tOnly tetrahedral meshes are supported."
      "\n\tHigher values print more info, value must be >= 0");
    set_option(ps, "memory_logfile_name"      , &memory_logfile_name      , INFO, "write memory usage to file");

    // properties
    set_option(ps, "property_map"             , &property_map             , PARAMETER, "YAML-based list of string pairs", "",
      "\n\te.g. {smoother_type: Newton, smoother_niter: 200, smoother_tol: 1.e-3}");

    // geometry
    set_option(ps, "input_geometry"           , &input_geometry           , IO, "input geometry name");
    set_option(ps, "smooth_geometry"          , &smooth_geometry          , PARAMETER, "smooth nodes that lie on geometry", "", 
      "\n\tmoves nodes after projecting them to the geometry to try to avoid bad-quality meshes");
    set_option(ps, "smoother_niter"           , &smoother_niter           , PARAMETER, "mesh smoother number of iterations");
    set_option(ps, "smoother_tol"             , &smoother_tol             , PARAMETER, "mesh smoother convergence tolerance");
    set_option(ps, "smooth_surfaces"          , &smooth_surfaces          , PARAMETER, "allow nodes to move on surfaces when smoothing");
    set_option(ps, "dump_geometry_file"       , &dump_geometry_file       , INFO, "debug print geometry (OpenNURBS 3dm) file contents");
    set_option(ps, "fit_geometry_file"        , &fit_geometry_file        , IO, "create a smooth geometry from the mesh boundary", "", 
          "\n\t2D meshes only."
          "\n\tcreate an OpenNURBS 3dm file fitting cubic splines to all mesh boundaries");
    set_option(ps, "fit_angle"                , &fit_angle                ,
                                  PARAMETER, "dihedral angle criterion for determining corners of MBG", "",
                                  "\n\tSpecified in degrees. 2D meshes only."
                                  "\n\tIf the dihedral is less than fit_angle, the node is considered a corner."
                                  "\n\tNote a flat surface has 180 degree angles everywhere.");
    set_option(ps, "fit_3d_file"              , &fit_3d_file              ,
                                  PARAMETER, "input YAML file. 3D meshes only.", "",
                                  "\n\tFit bi-cubics to mesh surface geometry, and store in fields."
                                  "\n\tFile contains control information on which surfaces"
                                  "to fit and angle criterion for surface seams."
                                  "\n\tA sample YAML file will be printed if --fit_3d_file=sample.yaml.");

    // smoothing
    set_option(ps, "smooth_use_reference_mesh", &smooth_use_reference_mesh, PARAMETER, "use the reference mesh when smoothing", "",
      "\n\tSet to 0 to smooth without a reference mesh."
      "\n\tFor most cases, set to 1 (default).");
    set_option(ps, "fix_all_block_boundaries" , &fix_all_block_boundaries , PARAMETER, "nodes on block boundaries are not smoothed", "",
      "\n\tEspecially useful when smoothing without geometry.");

    // mesh query
    set_option(ps, "compute_hmesh"            , &compute_hmesh            , PARAMETER, "compute mesh parameter using method eigens|edges");
    set_option(ps, "print_hmesh_surface_normal"  , &print_hmesh_surface_normal            , IO, "prints a table of normal mesh spacing at each surface");

    // histograms
    set_option(ps, "histogram"  , &histogram_options  , INFO, "print a histogram, of fields and mesh quality at a timestep.", "",
                                  "\n\tValue may be a filename or a string."
                                  "\n\tFilename: file contains parameters in YAML format (yaml.org)"
                                  "\n\tString: parameter format is \"{ fields: [field_1,...,field_n], file_root: my_histograms,"
                                  "\n\t  mesh: [edge_length, quality_edge, quality_vol_edge_ratio, volume], time: 0.1, step: 2 }\"");

    set_option(ps, "histogram_file_root"      , &histograms_root          , INFO, "root name of histogram files, or cout (screen).");

    set_option(ps, "generated_mesh"           , &generated_mesh           , OPERATION, "generate an NxMxP hexahedral mesh and write to --output_mesh", "",
                  "\n\tSpecify the mesh size using --input_mesh=NxMxP|bbox:xmin,xmax,ymin,ymax,zmin,zmax|sideset:xXyYzZ"
                  "\n\tBbox and sideset are optional. sideset args are where to put sidesets"
                  "\n\tE.g: '--generated_mesh --input_mesh=\"10x10x10|bbox:-1,-1,-1,2,2,1|sideset:xXyZ\"'"
                  "\n\tgenerates a 10x10x10 mesh inside sidesets on minX, maxX, minY and maxZ sides."
                  );

    set_option(ps, "remove_geometry_blocks"   , &remove_geometry_blocks   , IO, "remove geometry blocks from output Exodus file", "",
      "\n\tRemoval is done after any refinement and geometry projection.");

    // internal
    set_option(ps, "sync_io_regions"          , &sync_io_regions          , IO, "synchronize input & output regions' Exodus ids", "", 
                                  "\n\tEnsures output mesh has the same block ids and names as the input mesh."
                                  "\n\tThis only makes sense for refine (not enrich or convert).", true);

#if defined(STK_BUILT_FOR_SIERRA)
    // Salinas
    set_option(ps, "rbar_blocks"              , &rbar_blocks              , PARAMETER, "blocks to treat as special Salinas RBARs.", "", 
      "\n\tRBARs will connect new nodes between two surfaces."
      "\n\tUse help_full=--blocks for help on the option format.");
#endif

    // ad simple usage help to the command line processor 
    // ps.clp.setDocString( simple_usage_string.c_str() );  // zzyk delete me

    // help
    bool err = false;
    if (found_help_list)
    {
      err = do_help(ps, LIST) || err;
    }
    if (found_help)
    {
      err = do_help(ps, ONELINE) || err;
    }
    if (found_help_full)
    {
      err = do_help(ps, FULL) || err;
    }

    if ( found_any_help )
      return (int) err;

    // try to parse the option to run, not just print help
    parsing_ok = false;
    do_adapt = false;

    int err_clp(0);
    int bad_option(0);
    try 
    {
      err_clp = processCommandLine(ps.clp, argc, argv, bad_option);
    }
    catch (const std::exception &exc) 
    {
      if (p_rank==0) std::cout << "\n" << syntax_error << " unknown parameters or unallowed values.\n";
      print_syntax_error_instructions(false);
      exit_safely(1);
    }
    if (err_clp) 
    {
      if (bad_option)
      {
        help_full = argv[bad_option];
        do_help(ps, FULL);
      }
      return err_clp;
    }

    set_operation_defaults(argc, argv);

    if ( set_input_geometry_type() )
      return 1;

    if ( check_parsing_results() )
      return 1;

    // all is well if got to here
    parsing_ok = true;
    do_adapt = true;
    return 0;
  }


  void MeshAdapt::print_in_columns(const std::vector<std::string> &words) const
  {
    if (words.empty()) return;

    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    int screen_char_width = (int) w.ws_col;

    const int num_words = (int) words.size();
    // see how many columns can fit in screen_char_width
    std::vector<int> col_width, test_width; // start with zero columns
    while(true)
    {
      int cols = col_width.size()+1;
      int rows = (int) ceil ( (double) num_words / (double) cols );
      // set initial width to all zeros
      test_width.clear();
      test_width.resize(cols, 0);
      for(int r=0; r<rows; r++)
      {
        for (int c = 0; c < cols; ++c)
        {
          auto i = r + c * rows;
          if (i<num_words)
          {
            const int padding = 4;
            test_width[c] = std::max( test_width[c], (int) words[i].length() + padding );
          }
        }
      }
      auto width_sum = std::accumulate(test_width.begin(), test_width.end(), 0);
      if (width_sum<screen_char_width)
      {
        col_width.swap(test_width);
      }
      else
        break;
      if ( rows <= 1 )
        break;
    }

    // now print it
    std::stringstream ss;
    const int cols = (int) col_width.size();
    int rows = (int) ceil ( (double) num_words / (double) cols );
    for (int r = 0; r < rows; ++r)
    {
      for (int c = 0; c < cols; ++c)
      {
        auto i = r + c * rows;
        if (i<num_words)
        {
          ss << std::left << std::setw(col_width[c]) << words[i];
        }
      }
      ss << "\n"; // end of row
    }
    ss << "\n"; // extra
    //PRINT_INFO("%s", ss.str().c_str());
  }


  int MeshAdapt::do_help(ParserSystem &ps, HelpType help_type)
  {
    std::string help_prefix;
    std::string option;
    switch (help_type)
    {
      case LIST:
      help_prefix = "\nListing options";
      option = help_list;
      break;
      case ONELINE:
      default:
      help_prefix = "\nHelp";
      option = help;
      break;
      case FULL: 
      help_prefix = "\nFull help";
      option = help_full;
      break;
    }
    // strip any preceeding "-" characters
    while (!option.empty() && option[0]=='-')
      option.erase(0,1);
    // strip anything after '=', including the '='
    auto eq_pos = option.find('=');
    if (eq_pos != std::string::npos)
      option.erase(eq_pos, option.length());
    help_prefix += " for `" + option + "':\n";

    std::vector< std::string > help_strings;

    // is the option an operation name?
    for (int op = 0; op < num_operations; ++op)
      if (option == operation_name[op])
      {
        if (p_rank==0) std::cout << help_prefix;
        // search both basic and advanced, regardless of operation
        for (auto & odm : ps.basic_option_help )
          if ( odm.second.op == op )
            help_strings.push_back( odm.second.get_help( help_type ) );
        for (auto & odm : ps.advanced_option_help)
          if ( odm.second.op == op )
            help_strings.push_back( odm.second.get_help( help_type ) );
        if (p_rank==0) print_help(help_strings, help_type);
        if (p_rank==0) std::cout << std::endl;
        return 0;
      }

    // is it one of the categories? 
    if (option == category_name[SIMPLE] ||
        option == category_name[COMMON] ||
        option == category_name[ALL] )
    {
      if (help_type == LIST)
        std::cout << "\n" << simple_usage_string_one_line << std::endl;
      else
        std::cout << simple_std_help << std::endl;
      if (option == category_name[SIMPLE])
        return 0;
      if (option == category_name[COMMON] ||
          option == category_name[ALL] )
      {
        // common 
        std::cout << help_prefix;
        for ( auto &odm : ps.basic_option_help )
          help_strings.push_back( odm.second.get_help( help_type ) );
      }
      if ( option == category_name[ALL] )
      {
        // advanced
        for ( auto &odm : ps.advanced_option_help )
          help_strings.push_back( odm.second.get_help( help_type ) );
      }
      if (p_rank==0) print_help(help_strings, help_type);
      if (p_rank==0) std::cout << std::endl;
      return 0;
    }
 
    // is it a specific keyword?
    std::string help_string;
    if ( 0 == get_help(ps, option, help_string, help_type ) )
    {
      std::cout << help_prefix << help_string << std::endl;
      return 0;
    }

    // was it none of the above, an unrecognized options?
    std::cout << "\nHelp for '" << option << "'' not found\n" << std::endl;
    std::cout << help_hint << "\n" << std::endl;

    return 1;

  }

  void MeshAdapt::print_help(std::vector< std::string > &help_strings, HelpType help_type) const
  {
    if (help_type == LIST)
      print_in_columns(help_strings);
    else
    {
      for ( auto &w : help_strings )
      {
        std::cout << w;
        // extra blank line between multi-line help
        if (help_type == FULL && std::count(w.begin(), w.end(), '\n') > 1)
          std::cout << "\n"; 
      }
      std::cout << std::endl;
    }
  }


  void MeshAdapt::mesh_based_geometry_setup()
  {
    FitGregoryPatches::SurfaceSets surfaceSets;
    FitGregoryPatches::AngleMap angleMap;
    fitter.reset( new FitGregoryPatches(*eMeshP, surfaceSets, angleMap, 15.0));
    if (fit_3d_file.length())
      {
        fitter->parse(fit_3d_file);
        bool doRegister=true;
        fitter->register_or_set_fields(doRegister);
      }
    else if (input_geometry_type == MESH_BASED_GREGORY_PATCH)
      {
        bool doRegister=true;
        fitter->register_or_set_fields(doRegister);
        eMeshP->add_input_field(eMeshP->m_gregory_control_points_field);
        eMeshP->add_input_field(eMeshP->m_gregory_control_points_field_shell);
        eMeshP->add_input_field(eMeshP->m_node_normals);
      }

#if HAVE_OPENNURBS
    grsf.reset(new GeometryRecoverySplineFit(*eMeshP, fit_angle));
#if defined(STK_PERCEPT_HAS_GEOMETRY)
    if (fit_geometry_file != "")
      {
        grsf->fit_geometry_create_parts_meta();
      }
#elif defined(NO_GEOM_SUPPORT)
    std::ostringstream oss;
    oss << "\nERROR: Geometry and/or smoothing is not currently supported on this platform. Try running with geometry turned off.";
    throw std::runtime_error(oss.str());
#endif
#endif
  }
  void MeshAdapt::mesh_based_geometry_fitting()
  {
    if (fit_3d_file.length())
      {
        if (!eMeshP->get_rank()) std::cout << "Fitting 3D geometry using control data from file: " << fit_3d_file
                                           << " QA.activate= " << fitter->m_QA.m_activate << " QA.file= " << fitter->m_QA.m_file
                                           << std::endl;
        fitter->computeControlPoints();
        if (fitter->m_QA.m_activate && fitter->m_QA.m_file.length())
          {
            fitter->createQA(fitter->m_QA.m_file+"_"+input_mesh);
          }
      }

#if defined(STK_PERCEPT_HAS_GEOMETRY)
#if HAVE_OPENNURBS
    if (fit_geometry_file != "")
      {
        grsf->fit_geometry(fit_geometry_file);
      }
#endif
#endif
  }

  void MeshAdapt::verify_mesh_util(bool isInit)
  {
    if (verify_meshes == "1" || verify_meshes == "2" || verify_meshes == "FV")
      {
        bool print_table=true;
        double badJac=0.0;
        bool use_finite_volume = false;
        if (verify_meshes == "FV")
          use_finite_volume = true;
        int dump_all_elements = 0;
        if (verify_meshes == "2")
          dump_all_elements = 1;
        std::string type = (isInit ? "input":"refined");
        if (!eMeshP->get_rank()) std::cout << "Verify " << type << " mesh..." << std::endl;
        if (eMeshP->check_mesh_volumes(print_table, badJac, dump_all_elements, use_finite_volume))
          {
            throw std::runtime_error("ERROR: verify_meshes shows a bad "+type+" mesh");
          }
        if (isInit) adaptedMeshVerifier.reset( new percept::AdaptedMeshVerifier(true));
        if (!adaptedMeshVerifier->isValid(*eMeshP, isInit))
          throw std::runtime_error("ERROR: AdaptedMeshVerifier found invalid "+type+" mesh");
      }
    else
      {
        if (verify_meshes != "0")
          VERIFY_MSG("--verify_meshes option unrecognized, use 0, 1, 2 or FV, your option: "+verify_meshes);
      }
  }

  void MeshAdapt::do_histograms_init()
  {
    if (histogram_options.size() != 0)
      {
        Histograms<double> histograms(histograms_root);
        HistogramsParser<double> hparser(histogram_options);
        hparser.create(histograms);

        double hopt_time = histograms.m_database_time;
        int hopt_step = histograms.m_database_step;
        int current_step = eMeshP->get_current_database_step();
        if (hopt_time >= 0.0)
          {
            eMeshP->read_database_at_time(hopt_time);
          }
        else if (hopt_step >= 0)
          {
            eMeshP->read_database_at_step(hopt_step);
          }
        else
          {
            // get last step
            int step = eMeshP->get_database_time_step_count();
            eMeshP->read_database_at_step(step);
            //std::cout << "step= " << step << " current_step= " << current_step << std::endl;
            //eMeshP->read_database_at_step(step?step:1);
          }

        eMeshP->mesh_field_stats(&histograms);
        histograms.compute_uniform_bins(10);

        if (!p_rank) {
          std::cout << "Before refine, user-requested histograms= " << std::endl;
          histograms.print(true);
        }
        // reset
        eMeshP->read_database_at_step(current_step);
      }
  }

  void MeshAdapt::do_histograms_final()
  {
    if (histogram_options.size() != 0)
      {
        Histograms<double> histograms(histograms_root);
        HistogramsParser<double> hparser(histogram_options);
        hparser.create(histograms);

        eMeshP->mesh_field_stats(&histograms);
        histograms.compute_uniform_bins(10);

        if (!p_rank) {
          std::cout << "After refine, user-requested histograms= " << std::endl;
          histograms.print(true);
        }
      }
  }

  void MeshAdapt::compute_hmesh_sizes_init()
  {
    if (compute_hmesh.size() != 0)
      {
        double hmesh=0.0;
        Histograms<double> histograms(histograms_root);
        HistogramsParser<double> hparser(histogram_basic_options);
        hparser.create(histograms);
        if (compute_hmesh == "eigens")
          {
            hmesh = eMeshP->hmesh_stretch_eigens(hmesh_min_max_ave_factor, &histograms["stretch_eigens"].m_data, &histograms["quality_edge_eigens"].m_data);
            histograms["stretch_eigens"].set_titles("Stretch Eigens Histogram");
            histograms["quality_edge_eigens"].set_titles("Stretch Eigens Max/Min Quality Histogram");
          }
        else if (compute_hmesh == "edges")
          {
            hmesh = eMeshP->hmesh_edge_lengths(hmesh_min_max_ave_factor, &histograms["edge_length"].m_data, &histograms["quality_edge"].m_data);
            histograms["edge_length"].set_titles("Edge Length Histogram");
            histograms["quality_edge_eigens"].set_titles("Edge Max/Min Quality Histogram");
          }
        else
          {
            throw std::runtime_error("unknown option for compute_hmesh: "+compute_hmesh);
          }
        histograms.compute_uniform_bins(10);

        if (!p_rank) {
          std::cout << "Before refine, Mesh size (h-parameter) = " << hmesh
                    << " (min = " << hmesh_min_max_ave_factor[0]
                    << " max = " << hmesh_min_max_ave_factor[1]
                    << " ave = " << hmesh_min_max_ave_factor[2]
                    << ") "
                    << std::endl;
          histograms.print(true);
        }
        hmesh_factor = hmesh;
      }
  }

  void MeshAdapt::compute_hmesh_sizes_final()
  {
    if (compute_hmesh.size() != 0)
      {
        double hmesh=0.0;
        double min_max_ave[3];
        Histograms<double> histograms(histograms_root);
        HistogramsParser<double> hparser(histogram_basic_options);
        hparser.create(histograms);

        if (compute_hmesh == "eigens")
          {
            hmesh = eMeshP->hmesh_stretch_eigens(min_max_ave, &histograms["stretch_eigens"].m_data, &histograms["quality_edge_eigens"].m_data);
            histograms["stretch_eigens"].set_titles("Stretch Eigens Histogram");
            histograms["quality_edge_eigens"].set_titles("Stretch Eigens Max/Min Quality Histogram");
          }
        else if (compute_hmesh == "edges")
          {
            hmesh = eMeshP->hmesh_edge_lengths(min_max_ave, &histograms["edge_length"].m_data, &histograms["quality_edge"].m_data);
            histograms["edge_length"].set_titles("Edge Length Histogram");
            histograms["quality_edge_eigens"].set_titles("Edge Max/Min Quality Histogram");
          }
        else
          {
            throw std::runtime_error("unknown option for compute_hmesh: "+compute_hmesh);
          }
        histograms.compute_uniform_bins(10);

        hmesh_factor /= hmesh;
        hmesh_min_max_ave_factor[0] /= min_max_ave[0];
        hmesh_min_max_ave_factor[1] /= min_max_ave[1];
        hmesh_min_max_ave_factor[2] /= min_max_ave[2];
        if (!p_rank) {
          std::cout << "After refine, Mesh size (h-parameter) = " << hmesh << " oldH/newH factor= " << hmesh_factor
                    << "\n (new min = " << min_max_ave[0]
                    << " max = " << min_max_ave[1]
                    << " ave = " << min_max_ave[2]
                    << ") "
                    << "\n (old/new min = " << hmesh_min_max_ave_factor[0]
                    << " max = " << hmesh_min_max_ave_factor[1]
                    << " ave = " << hmesh_min_max_ave_factor[2]
                    << ") "
                    << std::endl;
          histograms.print(true);
        }
      }
  }

  void MeshAdapt::create_refine_pattern()
  {
    if (use_transition_elements)
      {
        eMeshP->register_and_set_refine_fields();
        std::set<stk::mesh::Part *> pl;
        if (1)
          {
            for (unsigned ii=0; ii < m_block_names[eMeshP->element_rank()].size(); ++ii)
              {
                std::string bn = m_block_names[eMeshP->element_rank()][ii];
                VERIFY_OP_ON ((bn[0] == '+' || bn[0] == '-'), ==, true, "bad block name: "+bn);
                std::string bname = bn.substr(1);
                stk::mesh::Part *part = eMeshP->get_fem_meta_data()->get_part(bname);

                if (part && !stk::equal_case(part->name(), bname)) {
                  const std::string alias = part->name();

                    if (debug && !eMeshP->get_rank())
                      std::cout << "block= " << bname << " replaced with alias=" << alias << std::endl;

                    const int mult = block_names_x_map[bname];
                    block_names_x_map[alias] = mult;
                }

                VERIFY_OP_ON(part, !=, 0, "couldn't find part: "+bname);

                if (bn[0] == '+')
                  {
                    pl.insert(part);
                  }
              }
            for (unsigned ii=0; ii < m_block_names[eMeshP->element_rank()].size(); ++ii)
              {
                std::string bn = m_block_names[eMeshP->element_rank()][ii];
                std::string bname = bn.substr(1);
                stk::mesh::Part *part = eMeshP->get_fem_meta_data()->get_part(bname);

                if (bn[0] == '-')
                  {
                    if (pl.find(part) != pl.end())
                      pl.erase(part);
                  }
              }
          }
        stk::mesh::PartVector pv(pl.begin(), pl.end());
        block_selector = stk::mesh::selectUnion( pv );

        if (debug && !eMeshP->get_rank())
          std::cout << "block_names= " << m_block_names << "\nfrom parts = " << eMeshP->print_part_vector_string(pv)
                    << "\nblock_selector= " << block_selector << std::endl;

        pattern = Teuchos::get_shared_ptr( make_local_break_pattern(*eMeshP) );
      }
    else
      {
        if (convert == "Wedge6_Pyramid5_Tet4_Exclude_Unconnected_Wedges") {
          convert = "Hex8_Wedge6_Pyramid5_Tet4";
          m_excludeWedgesNotConnectedToPyramids = 1;
        }

        pattern = Teuchos::get_shared_ptr(UniformRefinerPatternBase::createPattern(refine, enrich, convert, *eMeshP, m_block_names));
      }
  }

  void MeshAdapt::create_refiner()
  {
    if (use_transition_elements)
      {
        univ_selector.reset(new stk::mesh::Selector(eMeshP->get_fem_meta_data()->universal_part()));
        element_refine_predicate.reset(new ElementRefinePredicate(*eMeshP, univ_selector.get(), eMeshP->m_refine_field, 0.0));
        refiner.reset(new TransitionElementAdapter<ElementRefinePredicate>(*element_refine_predicate, *eMeshP, *pattern, 0));
        refiner->setRemoveOldElements(false);
        refiner->setAlwaysInitializeNodeRegistry(false);
      }
    else
      {
        refiner.reset(new UniformRefiner(*eMeshP, *pattern, proc_rank_field_ptr));
        refiner->setRemoveOldElements(true);
        if (block_name_inc.size())
          {
            refiner->setRemoveOldElements(false);
            refiner->setAlwaysInitializeNodeRegistry(false);
            eMeshP->output_active_children_only(true);
          }
        if (m_excludeWedgesNotConnectedToPyramids) {
          refiner->setRemoveOldElements(false);
        }
      }
    refiner->setAlternateRootTimer(&m_timer);
  }

  int MeshAdapt::bucket_in_block_names(stk::mesh::Bucket& bucket, std::string& block)
  {
    stk::mesh::PartVector pv = bucket.supersets();
    for (unsigned ii=0; ii < pv.size(); ++ii)
      {
        std::string partName = pv[ii]->name();

        stk::mesh::Part& part = *pv[ii];
        bool auto_part = 0 != part.attribute<AutoPart>();
        if (stk::mesh::is_auto_declared_part(part) || auto_part)
          continue;

        if (block_names_x_map.find(partName) != block_names_x_map.end())
          {
            int mult = block_names_x_map[partName];
            block = partName;
            return mult;
          }
      }
    block = "";
    return 0;
  }

  void MeshAdapt::pre_refiner_tasks(int iBreak)
  {
    if (m_excludeWedgesNotConnectedToPyramids) {
      stk::mesh::Part* exclude_part = eMeshP->get_fem_meta_data()->get_part("exclude_part");
      stk::mesh::PartVector pv(1, exclude_part);
      excludeWedgesNotConnectedToPyramids(*eMeshP, pv);
      
      refiner->setExcludeParts(pv);
    }

    if (use_transition_elements)
      {
        const stk::mesh::BucketVector & buckets = eMeshP->get_bulk_data()->buckets( eMeshP->element_rank() );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            RefineFieldType::value_type val = static_cast<RefineFieldType::value_type>(0);
            if (block_selector(bucket))
              {
                std::string pname;
                int mult = bucket_in_block_names(bucket, pname);
                if (iBreak < mult)
                  {
                    val = static_cast<RefineFieldType::value_type>(1);
                  }
              }
            const unsigned num_elements_in_bucket = bucket.size();
            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                stk::mesh::Entity element = bucket[iElement];
                RefineFieldType::value_type *fdata = stk::mesh::field_data( *eMeshP->m_refine_field , element );
                fdata[0] = static_cast<RefineFieldType::value_type>(val);
              }
          }
      }
    else if (block_name_inc.size())
      {
        stk::mesh::PartVector pv = eMeshP->get_fem_meta_data()->get_mesh_parts();
        typedef std::set<stk::mesh::Part *> PartSet;
        PartSet parts_all,  parts_in, parts_not_in;
        for (unsigned ii=0; ii < pv.size(); ++ii)
          {
            bool auto_part = (0 != pv[ii]->attribute<AutoPart>());
            if (auto_part)
              continue;

            if (pv[ii]->primary_entity_rank() == eMeshP->element_rank())
              {
                parts_all.insert(pv[ii]);
              }
          }
        for (StringIntMap::iterator it = block_names_x_map.begin(); it != block_names_x_map.end(); ++it)
          {
            std::string partName = it->first;
            int mult = it->second;
            stk::mesh::Part *part = eMeshP->get_fem_meta_data()->get_part(partName);
            if (debug) std::cout << "part= " << part->name() << " iBreak= " << iBreak << " mult= " << mult << std::endl;
            if (iBreak < mult)
              {
                parts_in.insert(part);
              }
          }
        std::set_difference(parts_all.begin(), parts_all.end(), parts_in.begin(), parts_in.end(), std::inserter(parts_not_in, parts_not_in.begin()));
        pv.assign(parts_not_in.begin(), parts_not_in.end());
        if (debug && !eMeshP->get_rank())
          {
            std::cout << " without use_transition_elements excluded parts= " << eMeshP->print_part_vector_string(pv,"\n") << std::endl;
          }
        if (pv.size())
          refiner->setExcludeParts(pv);
      }
  }

  // define static vars from RunAdaptRunInfo
  double RunAdaptRunInfo::mMeshInputTime = 0.0;
  double RunAdaptRunInfo::mMeshOutputTime = 0.0;
  double RunAdaptRunInfo::mAdaptTimeOverall = 0.0;
  double RunAdaptRunInfo::mAdaptCPUTimeOverall = 0.0;

  void MeshAdapt::run_adapt_run()
  {

    // This is not currently used in our test suite
    // to support legacy uses, added an extension to the output_mesh name if there is none
    if ( !(has_suffix(output_mesh, ".e") || has_suffix(output_mesh, ".exo"))) {
      output_mesh = output_mesh + ".e";
    }

    RunAdaptRunInfo::run_adapt_run(previous_adapted_mesh, input_mesh,
                                   next_adapted_mesh, output_mesh,
                                   adapt_input, input_geometry, smooth_geometry,
                                   ioss_read_options,
                                   ioss_write_options,
                                   property_map,
                                   verify_meshes,
                                   memory_logfile_name,
                                   m_timer,
                                   false);

    mMeshInputTime = RunAdaptRunInfo::mMeshInputTime;
    mMeshOutputTime = RunAdaptRunInfo::mMeshOutputTime;
    mAdaptTimeOverall = RunAdaptRunInfo::mAdaptTimeOverall;
    mAdaptCPUTimeOverall = RunAdaptRunInfo::mAdaptCPUTimeOverall;
  }

  void MeshAdapt::do_dihedral_angle_check()
  {
    m_dihedral_angle_check->find_simplex_elements_with_obtuse_angles();
  }

  BlockNamesType MeshAdapt::process_block_names()
  {
    // FIXME move this next block of code to a method on UniformRefiner
    BlockNamesType block_names(percept::EntityRankEnd+1u);
    BlockNamesType block_names_rbar(percept::EntityRankEnd+1u);

#if defined(STK_BUILT_FOR_SIERRA)
    if (rbar_blocks.length())
      {
        BlockNamesType rbar_names(percept::EntityRankEnd+1u);
        std::string block_name_inc_orig = block_name_inc;
        block_names_rbar = RefinerUtil::getBlockNames(block_name_inc, eMeshP->get_rank(), *eMeshP);
        if (rbar_blocks.length())
          {
            rbar_names = RefinerUtil::getBlockNames(rbar_blocks, eMeshP->get_rank(), *eMeshP);
            if (!eMeshP->get_rank())
              std::cout << "rbar_blocks= " << rbar_blocks << " rbar_names= " << rbar_names << std::endl;
          }
        for (unsigned ii=0; ii < rbar_names[eMeshP->element_rank()].size(); ii++)
          {
            std::string srb = rbar_names[eMeshP->element_rank()][ii];
            Util::replace(srb, "+", "-");
            block_names_rbar[eMeshP->element_rank()].push_back(srb);
          }
        block_name_inc = "";
        for (unsigned ii=0; ii < block_names_rbar[eMeshP->element_rank()].size(); ii++)
          {
            block_name_inc += (ii ? "," : "") + block_names_rbar[eMeshP->element_rank()][ii];
          }
        if (!eMeshP->get_rank())
          std::cout << "rbar: original block_name option = " << block_name_inc_orig << " new = " << block_name_inc << std::endl;
      }
#endif

    if (block_name_inc.length())
      {
        if (debug && !eMeshP->get_rank()) std::cout << "block_names: original (or after rbar) block_name option = " << block_name_inc << std::endl;
        block_names = RefinerUtil::getBlockNames(block_name_inc, eMeshP->get_rank(), *eMeshP);

        if (1)
          {
            int number_refines_max = 0;
            for (unsigned ii=0; ii < block_names[eMeshP->element_rank()].size(); ++ii)
              {
                std::string& bn = block_names[eMeshP->element_rank()][ii];
                //std::cout << "bn= " << bn << std::endl;
                size_t pos = bn.find(":");
                if (pos != std::string::npos)
                  {
                    std::string xp = bn.substr(pos+1);
                    size_t xpos = xp.find("x");
                    if (xpos == std::string::npos)
                      xpos = xp.find("X");
                    VERIFY_OP_ON(xpos, !=, std::string::npos, "syntax for multiplier for --blocks option is missing 'x' or 'X'");
                    xp = xp.substr(0, xpos);
                    std::string bnNew = bn.substr(0,pos);
                    if (!eMeshP->get_rank())
                      std::cout << "Note: found multiplier in block name, bn = " << bn << " mult= " << xp << " bnNew= " << bnNew << std::endl;
                    int mult = toInt(xp);
                    size_t posP = bnNew.find("+");
                    size_t posM = bnNew.find("-");
                    if (posP == std::string::npos && posM == std::string::npos)
                      {
                        VERIFY_MSG("ERROR: block name after processing not preceded by '+' or '-', bnNew= " +bnNew);
                      }
                    if (posP != std::string::npos)
                      {
                        bn = bnNew;
                        bnNew = bnNew.substr(1);
                        block_names_x_map[bnNew] = mult;
                        number_refines_max = std::max(number_refines_max, mult);
                        if (debug && !eMeshP->get_rank())
                          std::cout << "Note: bnNew= " << bnNew << std::endl;
                      }
                  }
                else
                  {
                    size_t posP = bn.find("+");
                    size_t posM = bn.find("-");
                    if (posP == std::string::npos && posM == std::string::npos)
                      {
                        VERIFY_MSG("ERROR: block name after processing not preceded by '+' or '-', bn= " +bn);
                      }
                    if (posP != std::string::npos)
                      {
                        std::string bnNew = bn;
                        bnNew = bnNew.substr(1);
                        block_names_x_map[bnNew] = number_refines;
                        if (debug && !eMeshP->get_rank())
                          std::cout << "Note: bnNew= " << bnNew << std::endl;
                      }
                  }
              }
            if (number_refines < number_refines_max)
              {
                if (!eMeshP->get_rank())
                  std::cout << "WARNING: number_refines must be set >= to the max found in --blocks option with the optional :Nx specification = " << number_refines_max
                            << "\n  NOTE: resetting number_refines to the max value." << std::endl;
              }
            if (debug && !eMeshP->get_rank())
              {
                std::cout << "block_names_x_map=\n" << block_names_x_map << std::endl;
              }
          }

        if (1)
          {
            eMeshP->commit();
            block_names = RefinerUtil::correctBlockNamesForPartPartConsistency(*eMeshP, block_names, input_geometry);

            eMeshP->close();
            pre_open();
            eMeshP->open(input_mesh);

            if (smooth_geometry) {
                const bool output = smooth_geometry > 1;
                eMeshP->add_coordinate_state_fields(output);
            }
#if !defined(NO_GEOM_SUPPORT)
            if (respect_spacing >= 1) {
              const bool output = respect_spacing > 1;
              eMeshP->set_respect_spacing(true);
              eMeshP->add_spacing_fields(output);
            }
#endif
            if (smooth_surfaces == 1) eMeshP->set_smooth_surfaces(true);

          }
        if (debug && !eMeshP->get_rank()) std::cout << "block_names after processing: " << block_names << std::endl;
      }
    return block_names;
  }


  void MeshAdapt::pre_open()
  {
    eMeshP->set_avoid_add_all_mesh_fields_as_input_fields(true);
  }

  int MeshAdapt::do_run_pre_commit()
  {
    if (ioss_read_options.length() || ioss_write_options.length())
      {
        if (!eMeshP->get_rank())
          {
            std::cout << "INFO: ioss_read_options=" << ioss_read_options << " ioss_write_options=" << ioss_write_options << std::endl;
          }
      }

    if (smooth_geometry)
      {
        if (smooth_surfaces == 1) eMeshP->set_smooth_surfaces(true);
        eMeshP->setProperty("smoother_niter", toString(smoother_niter));
        eMeshP->setProperty("smoother_tol", toString(smoother_tol));
        eMeshP->setProperty("smooth_use_reference_mesh", (smooth_use_reference_mesh?"1":"0"));
      }
    if (ioss_read_options.length())  eMeshP->set_ioss_read_options(ioss_read_options);
    if (ioss_write_options.length()) eMeshP->set_ioss_write_options(ioss_write_options);

    if (adapt_input.length()) {
      run_adapt_run();
      return -1;
    }

    pre_open();
    {
      stk::diag::Timer timerReadMesh_("ReadMesh", m_timer);
      stk::diag::TimeBlock tbReadMeshg_(timerReadMesh_);
      
      double t0 = stk::wall_time();
      if (generated_mesh)
        eMeshP->new_mesh(GMeshSpec(input_mesh));
      else
        eMeshP->open(input_mesh);
      double t1 = stk::wall_time();
      mMeshInputTime = t1 - t0;
    }
    if (smooth_surfaces == 1) eMeshP->set_smooth_surfaces(true);

    if (smooth_geometry) {
        const bool output = smooth_geometry > 1;
        eMeshP->add_coordinate_state_fields(output);
    }
#if !defined(NO_GEOM_SUPPORT)
    if (respect_spacing >= 1) {
      const bool output = respect_spacing > 1;
      eMeshP->set_respect_spacing(true);
      eMeshP->add_spacing_fields(output);
    }
#endif
    if (sync_io_regions)
      {
        if (convert.length() || enrich.length())
          {
            if (!eMeshP->get_rank())
              {
                std::cout       << "WARNING: removing property original_topology_type from input (and output mesh) since topology is changed by convert or enrich" << std::endl;
              }
            eMeshP->set_remove_io_orig_topo_type(true);
          }
      }
    eMeshP->set_sync_io_regions(sync_io_regions);
    if (!s_spatialDim) s_spatialDim = eMeshP->get_spatial_dim();

    Util::setRank(eMeshP->get_rank());

    if ((input_geometry_type != PGEOM_ACIS) && (input_geometry_type != PGEOM_OPENNURBS))
    {
      m_block_names = process_block_names();
      create_refine_pattern();
    }

    if (fix_all_block_boundaries)
      {
        bool make_part_io_part=true;
        eMeshP->add_part("inner_skin_part", make_part_io_part);
      }

    // Mesh-based geometry Fitting - setup parts
    mesh_based_geometry_setup();

    if (input_geometry.length()) {
      eMeshP->register_and_set_smoothing_fields(); 

      setup_m2g_parts(input_geometry);
    }

    eMeshP->add_registered_refine_fields_as_input_fields();

    eMeshP->get_ioss_mesh_data()->add_all_mesh_fields_as_input_fields();

    if (dihedral_angle_check != 0)
      {
        m_dihedral_angle_check.reset(new DihedralAngleCheck(eMeshP.get(), (dihedral_angle_check > 0 ? dihedral_angle_check : 1)));
      }

    if (m_excludeWedgesNotConnectedToPyramids)
    {
      eMeshP->get_fem_meta_data()->declare_part("exclude_part", eMeshP->element_rank());
    }

    return 0;
  }

  int MeshAdapt::do_run_post_commit()
  {

    if (DO_MEMORY) {
      std::string hwm = print_memory_both(eMeshP->parallel());
      if (!eMeshP->get_rank()) std::cout << "MEM: " << hwm << " initial memory after opening input mesh."  << std::endl;
    }

    if (fix_all_block_boundaries)
      {
        eMeshP->get_skin_part("inner_skin_part", true);
      }

    if (dihedral_angle_check != 0)
      {
        do_dihedral_angle_check();
      }

    if (skin_mesh)
      {
        eMeshP->skin_mesh();
      }

    if (input_geometry_type == MESH_BASED_GREGORY_PATCH)
      {
        stk::mesh::FieldVector from(1, eMeshP->get_coordinates_field()), to(1, eMeshP->m_unprojected_coordinates);
        PerceptMesh::copy_fields(*eMeshP, *eMeshP, from, to, eMeshP->node_rank(), EntitySelectorUCF(*eMeshP));
      }

    // Mesh-based-geometry fitting
    mesh_based_geometry_fitting();

    verify_mesh_util(true);

    do_histograms_init();

    compute_hmesh_sizes_init();

    if (print_hmesh_surface_normal)
      {
        std::string msg="before refine";
        eMeshP->print_hmesh_surface_normal(msg, std::cout);
      }

#if !defined(NO_GEOM_SUPPORT)
    if (respect_spacing)
      {
        SpacingFieldUtil sfu(*eMeshP);
        sfu.compute_spacing_field();
      }
#endif

      write_memory_logfile(m_comm, READ_MESH, memory_logfile_name);

    // doRefineMesh
      {
        t0 =  stk::wall_time();
        cpu0 = stk::cpu_time();

        create_refiner();

#if defined(STK_PERCEPT_LITE) && STK_PERCEPT_LITE == 0
        ProgressMeter pm(*refiner);
#endif

        if (input_geometry != "") {
          refiner->setGeometryFile(input_geometry);
          refiner->setSmoothGeometry(smooth_geometry);
          
          initialize_m2g_geometry(input_geometry);
        }

        refiner->setFixAllBlockBoundaries(fix_all_block_boundaries);
        refiner->setDoProgressMeter(progress_meter == 1 && 0 == p_rank);
#if defined(STK_BUILT_FOR_SIERRA)
        if (rbar_blocks.length())
          {
            BlockNamesType rbar_names(percept::EntityRankEnd+1u);
            if (rbar_blocks.length())
              {
                if (!eMeshP->get_rank())
                  std::cout << "For RBAR treatment: rbar_blocks= " << rbar_blocks << std::endl;

                rbar_names = RefinerUtil::getBlockNames(rbar_blocks, eMeshP->get_rank(), *eMeshP);
                if (!eMeshP->get_rank())
                  std::cout << "For RBAR treatment: rbar_names (after getBlockNames)= " << rbar_names << std::endl;
              }
            refiner->set_rbar_special_treatment(rbar_names);
          }
#endif

        for (int iBreak = 0; iBreak < number_refines; iBreak++)
          {
            if (!eMeshP->get_rank())
              {
                std::cout << "Refinement pass # " << (iBreak+1) << " start..." << std::endl;
              }

            pre_refiner_tasks(iBreak);

            refiner->doBreak();

            if (!eMeshP->get_rank())
              {
                std::cout << std::endl;
                bool printAllTopologies = false;
                refiner->getRefinementInfo().printTable(std::cout, iBreak, printAllTopologies);
                std::cout << std::endl;
              }

            write_memory_logfile(m_comm, REFINE_MESH, memory_logfile_name);

            if (DO_MEMORY) {
              std::string hwm = print_memory_both(eMeshP->parallel());
              if (!eMeshP->get_rank()) std::cout << "MEM: " << hwm << " memory after refining mesh for pass # " << iBreak  << std::endl;
            }


          } // iBreak


          {
            if (DO_MEMORY) {
              std::string hwm = print_memory_both(eMeshP->parallel());
              if (!eMeshP->get_rank()) std::cout << "MEM: " << hwm << " memory before deleteParentElements" << std::endl;
            }

            refiner->deleteParentElements();

            if (DO_MEMORY) {
              std::string hwm = print_memory_both(eMeshP->parallel());
              if (!eMeshP->get_rank()) std::cout << "MEM: " << hwm << " memory after deleteParentElements"  << std::endl;
            }
          }


        if (number_refines == 0 && smooth_geometry)
          {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
            refiner->setSmoothGeometry(smooth_geometry);
            const bool snap_geometry = false;
            refiner->snapAndSmooth(snap_geometry, input_geometry, (smooth_use_reference_mesh == 1) );
#elif defined(NO_GEOM_SUPPORT)
            std::ostringstream oss;
            oss << "\nERROR: Geometry and/or smoothing is not currently supported on this platform. Try running with geometry turned off.";
            throw std::runtime_error(oss.str());
#endif
          }

        refiner->deleteParentElements();

        t1 =  stk::wall_time();
        cpu1 = stk::cpu_time();

        mAdaptTimeOverall = t1 - t0;
        mAdaptCPUTimeOverall = cpu1 - cpu0;

        verify_mesh_util(false);

        do_histograms_final();

        compute_hmesh_sizes_final();

        if (print_hmesh_surface_normal)
          {
            std::string msg="after refine";
            eMeshP->print_hmesh_surface_normal(msg, std::cout);
          }

        if (remove_geometry_blocks) eMeshP->remove_geometry_blocks_on_output(input_geometry);

        {
          stk::diag::Timer timerWriteMesh_("WriteMesh", m_timer);
          stk::diag::TimeBlock tbWriteMeshg_(timerWriteMesh_);

          double t0 = stk::wall_time();
          eMeshP->save_as(output_mesh);
          double t1 = stk::wall_time();
          mMeshOutputTime = t1 - t0;
        }

        write_memory_logfile(m_comm, WRITE_MESH, memory_logfile_name);

      } // doRefineMesh
    return 0;
  }

  int MeshAdapt::do_post_proc()
  {
    if (DO_MEMORY) {
      std::string hwm = print_memory_both(eMeshP->parallel());
      if (!eMeshP->get_rank()) std::cout << "MEM: " << hwm << " final memory after refining mesh."  << std::endl;
    }

    double cpuMax = (cpu1-cpu0);
    double wallMax = (t1-t0);
    double wallMin = wallMax;
    double cpuSum = (cpu1-cpu0);
    double cpuMin = cpuMax;

    stk::all_reduce( eMeshP->get_bulk_data()->parallel(), stk::ReduceSum<1>( &cpuSum ) );
    stk::all_reduce( eMeshP->get_bulk_data()->parallel(), stk::ReduceMax<1>( &cpuMax ) );
    stk::all_reduce( eMeshP->get_bulk_data()->parallel(), stk::ReduceMin<1>( &cpuMin ) );
    stk::all_reduce( eMeshP->get_bulk_data()->parallel(), stk::ReduceMax<1>( &wallMax ) );
    stk::all_reduce( eMeshP->get_bulk_data()->parallel(), stk::ReduceMin<1>( &wallMin ) );

    if (0 == p_rank)
      {
        std::cout << "P[" << p_rank << ", " << p_size << "]  min wall clock time = " << wallMin << " (sec)" << std::endl;
        std::cout << "P[" << p_rank << ", " << p_size << "]  max wall clock time = " << wallMax << " (sec)" << std::endl;
        std::cout << "P[" << p_rank << ", " << p_size << "]  min cpu  clock time = " << cpuMin << " (sec)" << std::endl;
        std::cout << "P[" << p_rank << ", " << p_size << "]  max cpu  clock time = " << cpuMax << " (sec)" << std::endl;
        std::cout << "P[" << p_rank << ", " << p_size << "]  sum cpu  clock time = " << cpuSum << " (sec)" << std::endl;
      }

    return 0;
  }

  // version - return true if build in Sierra
  bool MeshAdapt::get_version(std::string* v)
  {
    if (v) *v = version_prefix+version;
#if defined(STK_BUILT_FOR_SIERRA)
    return true;
#else
    return false;
#endif
  }

  void MeshAdapt::log_usage( bool status )
  {
#if defined(STK_BUILT_FOR_SIERRA)
    const bool disable_audit = !sierra::Env::get_param("noaudit").empty() || std::getenv("SIERRA_USAGE_METRICS_OFF") != NULL;

    size_t hwm_max = 0, hwm_min = 0, hwm_avg = 0;

    stk::get_memory_high_water_mark_across_processors(eMeshP->parallel(), hwm_max, hwm_min, hwm_avg);

    stk::all_reduce( eMeshP->parallel(), stk::ReduceSum<1>( &mAdaptCPUTimeOverall ) );

    if (eMeshP->get_parallel_rank() == 0) {
      // Audit log
      bool runtest = sierra::Env::get_param("runtest").empty() ? false : true;
      if ((!disable_audit) && (!runtest)) {
        const double bytes_in_MB = 1024*1024;
        auditdata data;
        AuditLogDefaults(&data, "mesh_adapt", sierra::ProductRegistry::version(), eMeshP->get_parallel_size());
        data.starttime = sierra::format_time(sierra::Env::start_time(), "%Y%m%d%H%M%S");
        data.purpose = "meshing";

        data.num_proc       = eMeshP->get_parallel_size();

        data.time_elapsed   = mAdaptTimeOverall;
        data.cputimesum     = mAdaptCPUTimeOverall;
        data.meshinputtime  = mMeshInputTime;
        data.meshoutptutime = mMeshOutputTime;
        data.hwm_min        = hwm_min / bytes_in_MB;
        data.hwm_max        = hwm_max / bytes_in_MB;
        data.hwm_avg        = hwm_avg / bytes_in_MB;
        data.status = status ? "success" : "fail";

        OutputAuditLog(&data);
      }
    }
#endif    
  }


  void MeshAdapt::init(int argc, char **argv)
  {
    fill_help_strings(argc, argv);

#if defined( STK_HAS_MPI )
    m_comm = stk::ParallelMachine(stk::parallel_machine_init(&argc, &argv));
    Kokkos::initialize(argc,argv);
#endif

    EXCEPTWATCH;

    p_rank = stk::parallel_machine_rank(m_comm);
    p_size = stk::parallel_machine_size(m_comm);
  }

  // the actual mesh adaptation! no parsing.
  int MeshAdapt::adapt_main_do_adapt()
  {
      write_memory_logfile(m_comm, INITIALIZE, memory_logfile_name);

    t0   = 0.0;
    t1   = 0.0;
    cpu0 = 0.0;
    cpu1 = 0.0;

#if defined( STK_PERCEPT_HAS_GEOMETRY )
    if (dump_geometry_file && input_geometry.find(".3dm") != std::string::npos)
      {
#if HAVE_OPENNURBS
        GeometryKernelOpenNURBS gko;
        gko.debug_dump_file(input_geometry);
#endif
      }
#endif

//    std::string input_mesh_save = input_mesh;
//    std::string output_mesh_save = output_mesh;

    eMeshP.reset(new percept::PerceptMesh);

    if (output_active_elements_only)
      eMeshP->output_active_children_only(true);

    if (1)
      {
        if (property_map.length())
          {
            eMeshP->parse_property_map_string(property_map);
          }
        const char * env_val = std::getenv("Percept_property_map");
        if (env_val)
          {
            std::string vv(env_val);
            if (eMeshP->get_rank() == 0) std::cout << "found Percept_property_map = " << vv << std::endl;
            eMeshP->parse_property_map_string(vv);
          }
        if (eMeshP->getProperty("MeshAdapt.debug") == "true")
          debug = 1;
      }

    if (progress_meter && eMeshP->get_rank() == 0)
      {
        std::cout << "Stage: Open mesh..." << " cpu: " << eMeshP->cpu_time() << " [sec]" << std::endl;
      }

    int res1 = do_run_pre_commit();

    // trap for RunAdaptRun, etc.
    if (res1 < 0)
      {
        if (print_timers) print_timer_table(m_timer);

        log_usage();

        eMeshP.reset();

        write_memory_logfile(m_comm, FINALIZE, memory_logfile_name);

        return 0;
      }

    if (progress_meter && eMeshP->get_rank() == 0)
      {
        std::cout << "Stage: Commit mesh..." << " cpu: " << eMeshP->cpu_time() << " [sec]" << std::endl;
      }

    eMeshP->commit();
    if (progress_meter && eMeshP->get_rank() == 0)
      {
        std::cout << "Stage: Commit mesh...done" << " cpu: " << eMeshP->cpu_time() << " [sec]" << std::endl;
      }

    do_run_post_commit();

    do_post_proc();

    if (print_timers) print_timer_table(m_timer);

    log_usage();

    eMeshP.reset();

    write_memory_logfile(m_comm, FINALIZE, memory_logfile_name);

    return result;
  }


void MeshAdapt::setup_m2g_parts(std::string input_geometry)
{
  if ( (input_geometry_type != PGEOM_ACIS) && (input_geometry_type != PGEOM_OPENNURBS) )	  return;

#ifdef HAVE_CUBIT        
#ifdef HAVE_ACIS
  //Madison Brewer: what I don't like about this is that we're not really utilizing the geometry interface layer/kernel to its full extent.
  //It essentially just gets called for snapping and disappears after that.
  //BIG QUESTION: Do we want to try and refactor percept to truly interact with its geometry through these kernels? How much refactoring would this incur?

  m_PGA = new PGeomACIS;
  m_PGeomPntr = m_PGA;
#else
  m_PGeomPntr = new PGeom;
#endif

  if (input_geometry_type == PGEOM_ACIS) {
    m_PGeomPntr->initialize(ACIS_GEOMETRY_ENGINE);
    m_PGeomPntr->import_acis_file(input_geometry.c_str());
  }
  else if (input_geometry_type == PGEOM_OPENNURBS) {
    m_PGeomPntr->initialize(OPENNURBS_GEOMETRY_ENGINE);
    m_PGeomPntr->import_open_nurbs_file(input_geometry.c_str());
  }

  stk::mesh::MetaData * md = eMeshP->get_fem_meta_data();

  std::vector<int> surfIDs;
  m_PGeomPntr->get_surfaces(surfIDs);
  std::vector<std::string> quadNames(surfIDs.size());
  std::vector<std::string> triNames(surfIDs.size());

  //make parts that we'll store new mesh entities on
  for(unsigned i = 0; i < surfIDs.size(); i++) { 
    std::string name = "geom_surface_quad_";
    name = name + std::to_string(surfIDs[i]);
    stk::mesh::Part& part = md->declare_part_with_topology(name,
                                                           stk::topology::SHELL_QUAD_4);
    if (dump_geometry_file) stk::io::put_io_part_attribute(part);
    quadNames[i] = name;
  } 

  for(unsigned i = 0; i<surfIDs.size();i++){
    std::string name = "geom_surface_tri_";
    name = name + std::to_string(surfIDs[i]);
    stk::mesh::Part& part = md->declare_part_with_topology(name,
                                                           stk::topology::SHELL_TRI_3);
    if (dump_geometry_file) stk::io::put_io_part_attribute(part);
    triNames[i] = name;
  }

  std::vector<int> curveIDs;
  m_PGeomPntr->get_curves(curveIDs);
  std::vector<std::string> curveNames(curveIDs.size());
  for (unsigned i = 0; i < curveIDs.size(); i++) { 
    std::string name = "geom_curve_";
    name = name + std::to_string(curveIDs[i]);
    stk::mesh::Part& part = md->declare_part_with_topology(name,
                                                           stk::topology::BEAM_2);
    if (dump_geometry_file) stk::io::put_io_part_attribute(part);
    curveNames[i] = name;
  }

  //setup refinement: bcarnes: why is this here?
  m_block_names = process_block_names();
  create_refine_pattern();
#endif
}
  
void MeshAdapt::initialize_m2g_geometry(std::string input_geometry)
{
  if( (input_geometry_type != PGEOM_ACIS) && (input_geometry_type != PGEOM_OPENNURBS) ) return;
#ifdef HAVE_CUBIT  

  std::string m2gFile = input_geometry.substr(0,input_geometry.length()-3) + "m2g";

  int THIS_PROC_NUM = stk::parallel_machine_rank( MPI_COMM_WORLD);

  stk::mesh::MetaData* md = eMeshP->get_fem_meta_data();
  stk::mesh::BulkData* bd = eMeshP->get_bulk_data();

  std::vector<int> curveIDs;
  m_PGeomPntr->get_curves(curveIDs);

  std::vector<int> surfIDs;
  m_PGeomPntr->get_surfaces(surfIDs);

  eMeshP->initializeIdServer();

  PGeomAssoc<stk::mesh::BulkData, stk::mesh::Entity, stk::mesh::Entity,
    stk::mesh::Entity, stk::mesh::Entity> geom_assoc(m_PGeomPntr);
  geom_assoc.set_node_callback(get_node_from_id);
  geom_assoc.set_edge_callback(get_beam_from_ids);
  geom_assoc.set_face_callback(get_shell_from_ids);
  geom_assoc.set_elem_callback(get_hex_from_id);
  geom_assoc.set_validate_nodes_callback(validate_node_ownership);
  geom_assoc.set_validate_element_callback(validate_element_ownership);
  geom_assoc.set_mesh(bd);
  geom_assoc.set_fill_curve_and_surface_maps_during_import(false);
  const bool geometry_exists = true;
  geom_assoc.import_m2g_file(m2gFile.c_str(), geometry_exists);

  bd->modification_begin();

  std::vector<stk::mesh::Entity> nodesToCheck;
  std::vector<int> procsSharedTo;
  for (unsigned i = 0; i < curveIDs.size(); i++) { //create beams and put them into corresponding curve parts

    std::vector<stk::mesh::Part *> add_parts_beams(1, static_cast<stk::mesh::Part*>(0));

    add_parts_beams[0] = md->get_part("geom_curve_" + std::to_string(curveIDs[i]));

    std::vector<std::vector<int>> edge_node_ids;
    geom_assoc.get_curve_edge_nodes(curveIDs[i], edge_node_ids);

    for (unsigned ii = 0; ii < edge_node_ids.size(); ii++) {

      bool toDeclare = true;

      std::vector<stk::mesh::EntityId> beamNodeIDs;
      int lowestRank = std::numeric_limits<int>::max();

      for (unsigned j = 0; j < edge_node_ids[ii].size(); j++)
        beamNodeIDs.push_back((stk::mesh::EntityId) edge_node_ids[ii][j]);

      nodesToCheck.clear();
      for (unsigned j = 0; j < edge_node_ids[ii].size(); j++) {

        stk::mesh::Entity cur_node = bd->get_entity(stk::topology::NODE_RANK, beamNodeIDs[j]);

        nodesToCheck.push_back(cur_node);
      }

      bd->shared_procs_intersection(nodesToCheck, procsSharedTo);
      procsSharedTo.push_back(THIS_PROC_NUM); //find all processes that own or have this node shared to it
      for (size_t iii = 0; iii < procsSharedTo.size(); iii++) {
        if (procsSharedTo[iii] < lowestRank)
          lowestRank = procsSharedTo[iii]; //lowest ranking process is responsible for creation of this entity
        //		QUESTION: does this create a significant load imbalance for creation?
      }

      if (lowestRank != THIS_PROC_NUM)
        toDeclare = false;

      if (toDeclare) {

        stk::mesh::EntityId id2 = eMeshP->getNextId(stk::topology::ELEMENT_RANK);
        stk::mesh::declare_element(*eMeshP->get_bulk_data(),
                                   add_parts_beams, id2, beamNodeIDs);
      }
    }
  }

  std::vector<std::pair<std::vector<int>, int>> uncreatedEnts;
  for (unsigned i = 0; i < surfIDs.size(); i++) {
    std::vector<stk::mesh::Part *> add_parts_shells(1,static_cast<stk::mesh::Part*>(0));

    std::vector<std::vector<int>> face_node_ids;
    geom_assoc.get_surface_face_nodes(surfIDs[i], face_node_ids);

    for (unsigned ii = 0; ii < face_node_ids.size(); ii++) {

      std::vector<stk::mesh::EntityId> shellNodeIDs;
      for (unsigned j = 0; j < face_node_ids[ii].size(); j++)
        shellNodeIDs.push_back((stk::mesh::EntityId) face_node_ids[ii][j]);

      if (shellNodeIDs.size() == 3)
        add_parts_shells[0] = md->get_part("geom_surface_tri_"
                                           + std::to_string(surfIDs[i])); //store these parts in a partvector for FASTER access
      else {
        add_parts_shells[0] = md->get_part("geom_surface_quad_"
                                           + std::to_string(surfIDs[i]));
      }

      bool toDeclare = true;

      int lowestRank = std::numeric_limits<int>::max();
      std::vector<stk::mesh::Entity> entitiesToCheck;
     
      procsSharedTo.clear(); //std::vector<int> procsSharedTo;

      stk::mesh::Entity cur_node;
      for (unsigned j = 0; j < face_node_ids[ii].size(); j++) {

        cur_node = bd->get_entity(stk::topology::NODE_RANK, shellNodeIDs[j]);

        entitiesToCheck.push_back(cur_node);
      }

      bd->shared_procs_intersection(entitiesToCheck, procsSharedTo);
      procsSharedTo.push_back(THIS_PROC_NUM);//find all processes that either own or have these nodes shared to them
      for (size_t iii = 0; iii < procsSharedTo.size(); iii++) {
        if (procsSharedTo[iii] < lowestRank)
          lowestRank = procsSharedTo[iii]; //lowest ranking process is responsible for creation
        //		QUESTION: does this create a significant load imbalance for creation?
      }

      if (lowestRank != THIS_PROC_NUM) {
        toDeclare = false;
      }

      if (toDeclare) {
        stk::mesh::EntityId id2 = eMeshP->getNextId(stk::topology::ELEMENT_RANK);

        stk::mesh::declare_element(*bd, add_parts_shells,
                                   id2, shellNodeIDs);
      }
    }
  }
  bd->modification_end();

  geom_assoc.fill_curve_and_surface_maps();
  delete m_PGeomPntr; //bad things happen if you don't explicitly reallocate the memory here
#endif
}

  int MeshAdapt::main(int argc, char **argv) 
  {
    int res = adapt_main(argc, argv);
    exit_safely(res);
    return res;
  }

}

