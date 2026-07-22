// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <math.h>
#include "uns_inline_decomp.h"
#include "element_density_function.h"
#include "geometry_transform.h"
#include "inline_mesh_desc.h"
#include "legacy_inline_mesh_desc.h"
#include "radial_inline_mesh_desc.h"
#include "radial_trisection_inline_mesh_desc.h"
#include "brick_inline_mesh_desc.h"
#include <sstream>
#include <fstream>
#include "../parser/parser.h"
#include "../parser/parse_table.h"
#include "../parser/InputBlock.h"
#include "../parser/token_stream.h"
#include "parse_routines.h"
#include <stdlib.h>

namespace PAMGEN_NEVADA {


  enum ParamType { 
    P_STARTING_BLOCK_NUMBER = 0,
    P_OFFSET,
    P_GMIN,
    P_GMAX ,
    P_NTHETA,
    P_NUMZ,
    P_NUMR,
    P_NUMA,
    P_INITIAL_RADIUS,
    P_TRISECTION_BLOCKS,
    P_FIRST_SIZE,
    P_LAST_SIZE,
    P_INTERVAL,
    P_ZBLOCK,
    P_RBLOCK,
    P_ABLOCK,
    P_NPHI,
    P_NR,
    P_NX,
    P_NY,
    P_NZ,
    P_BX,
    P_BY,
    P_BZ,
    P_BTHETA,
    P_BPHI,
    P_BR,
    P_RI,
    P_RO,
    P_THETA,
    P_MINTHETA,
    P_PHI,
    P_XMIN,
    P_YMIN,
    P_ZMIN,
    P_ZMAX,
    P_AMR,
    P_PERIODIC,
    P_NODESET,
    P_SIDESET,
    P_BLOCK_NODESET,
    P_BLOCK_SIDESET,
    P_DEBUG,
    P_SQUARED,
    P_ENFORCE_PERIODIC,
    P_NUMELARRAY,
    P_TRANSITION_RADIUS,
    P_NOCUTSI,
    P_NOCUTSJ,
    P_NOCUTSK,
    MAX_PARAM };


  static long long sdimension;
  static long long * parse_error_count = NULL;
  static bool has_specified_geometry_type=false;

  void Allow_New_Mesh_Specification(){ has_specified_geometry_type=false;}

  //This takes an integer and returns it in an English representation
  std::string EnglishNumber(long long the_number) {
    // set up our core numbers from which other numbers can be constructed
    const long long numbers[] = {0,1,2,3,4,5,6,7,8,9,
      10,11,12,13,14,15,16,17,18,19,20,30,40,50,60,70,80,90,
      100,1000,1000000LL,1000000000LL,1000000000000LL};

    const std::string english[] = {"zero", "one", "two", "three", "four", "five", "six",
      "seven", "eight", "nine", "ten", "eleven", "twelve", "thirteen", "fourteen",
      "fifteen", "sixteen", "seventeen", "eighteen", "nineteen", "twenty",
      "thirty", "fourty", "fifty", "sixty", "seventy", "eighty", "ninety",
      "hundred", "thousand", "million","billion","quatrillion"};

    // count of core numbers (size of above arrays)
    const int count = 33;

    // two strings, one holding our number so far, and another temp.
    std::string returnString = "";
    std::string temp;

    if(the_number < 0) {
      the_number = -the_number;
      returnString = "negative";
    }

    // simply return the number if it already exists in the array above
    for(int i = 0; i < count; i++) {
      if(the_number == numbers[i]) {
        if(returnString.size()) // this is a negative number, so add qualifier
          return returnString + " " + english[i];

        return english[i]; // positive number, just return it
      }
    }

    // otherwise we have to contruct it
    for(int i = count - 1; i > 0; i--) {
      if(the_number >= numbers[i]) {
        temp = english[i];

        // e.g. if n = 5000, then we have "thousand" + n/1000 which equals 5 so,
        //      "five" + " " + "thousand"
        if(the_number >= 100)
          temp = EnglishNumber(the_number / numbers[i]) + " " + temp;

        if(returnString.size() > 0) returnString = returnString + " ";
        returnString = returnString + temp; // add this number to our existing number

        // remove this number from n, to complete construction of the number
        the_number = the_number - (long long)(the_number / numbers[i]) * numbers[i];
      }
    }

    return returnString;
  }

  /****************************************************************************/
  Inline_Mesh_Desc * Parse_Inline_Mesh(std::string & file_name,
      std::stringstream & input_stream,
      long long & sparse_error_count,
      long long incoming_dim,
      long long max_int)
    /****************************************************************************/
  {
    // makes this a noop if an Inline_Mesh_Desc has alread been created.
    if(Inline_Mesh_Desc::im_static_storage) return Inline_Mesh_Desc::im_static_storage;
    Allow_New_Mesh_Specification();

    sdimension = incoming_dim;
    parse_error_count = & sparse_error_count;

    Parse_Table * parse_table = new Parse_Table;

    Keyword keyword_table[] = {
      {"MESH", 0, PAMGEN_NEVADA::Parse_Inline_Mesh_Tok},
      {"TWO D MESH", 0, PAMGEN_NEVADA::Parse_Inline_Mesh_2D_Tok},
      {"THREE D MESH", 0, PAMGEN_NEVADA::Parse_Inline_Mesh_3D_Tok},
    };

    parse_table->merge(keyword_table,3);

    InputBlock * input_tree = new InputBlock;
    input_tree->set( file_name.c_str(), 0 );

    Token_Stream token_stream(input_stream, Inline_Mesh_Desc::echo_stream, input_tree,0,0,true);

    PAMGEN_NEVADA::Parse(&token_stream,parse_table,TK_EXIT);

    if(!Inline_Mesh_Desc::im_static_storage) {
      delete parse_table;
      return Inline_Mesh_Desc::im_static_storage;
    }

    if(token_stream.Error_Count() != 0) {
      delete parse_table;
      return NULL;
    }

    PAMGEN_NEVADA::Inline_Mesh_Desc * ims = Inline_Mesh_Desc::first_im_static_storage;
    while(ims){

      /*flip to sequential if any topology is changed*/
      ims->brokerDecompositionStrategy();

      long long error_code = ims->Check_Blocks();
      if(error_code){
        (*parse_error_count) ++;
        Inline_Mesh_Desc::echo_stream << "SETUP ERROR IN CHECK_BLOCKS " << ims->getErrorString() << std::endl;
        return NULL;
      }

      error_code = ims->Check_Block_BC_Sets();
      if(error_code){
        (*parse_error_count) ++;
        Inline_Mesh_Desc::echo_stream << "SETUP ERROR IN CHECK_BLOCK_BC_SETS " << ims->getErrorString() << std::endl;
        return NULL;
      }

      error_code = ims->Set_Up();
      if(error_code){
        (*parse_error_count) ++;
        Inline_Mesh_Desc::echo_stream << "SETUP ERROR " << ims->getErrorString() << std::endl;
        return NULL;
      }

      error_code = ims->Check_Spans();
      if(error_code){
        (*parse_error_count) ++;
        Inline_Mesh_Desc::echo_stream << "SPAN ERROR " << ims->getErrorString() << std::endl;
        return NULL;
      }
      ims = ims->next;
    }


    ims = Inline_Mesh_Desc::first_im_static_storage;
    while (ims)
    {
      long long total_el_count;
      long long total_node_count;
      long long total_edge_count;

      ims->calculateSize(total_el_count,total_node_count,total_edge_count);

      long long error_code = ims->reportSize(total_el_count,
          total_node_count,
          total_edge_count,
          ims->error_stream,
          max_int);
      if(error_code){
        (*parse_error_count)++;
        Inline_Mesh_Desc::echo_stream << "SIZING ERROR " << ims->getErrorString();
        return NULL;
      }
      {
        ims->info_stream << "Inline mesh specification requested: \n";
        ims->info_stream << "\t" << total_el_count;
        //       ims->info_stream << "\t" << EnglishNumber(total_el_count);
        ims->info_stream << " Elements \n";
        ims->info_stream << "\t" << total_node_count;
        //       ims->info_stream << "\t" << EnglishNumber(total_node_count);
        ims->info_stream << " Nodes and\n ";
        ims->info_stream << "\t" << total_edge_count;
        //       ims->info_stream << "\t" << EnglishNumber(total_edge_count);
        ims->info_stream << " Edges.\n";
      }
      ims = ims->next;
    }


    if(parse_table)delete parse_table;
    if(input_tree) delete input_tree;
    return Inline_Mesh_Desc::first_im_static_storage;
  }

  /****************************************************************************/
  Token Parse_Inline_Mesh_3D_Tok(Token_Stream *token_stream, int value)
    /****************************************************************************/
  {
    sdimension = 3;
    return Parse_Inline_Mesh_Tok(token_stream,value);
  }

  /****************************************************************************/
  Token Parse_Inline_Mesh_2D_Tok(Token_Stream *token_stream, int value)
    /****************************************************************************/
  {
    sdimension = 2;
    return Parse_Inline_Mesh_Tok(token_stream,value);
  }

  /****************************************************************************/
  Token Parse_Inline_Mesh_Tok(Token_Stream *token_stream, int value)
    /****************************************************************************/
  {
    // if this function is used to parse a region_io block of input
    // then only a default set of keywords is parsed:

    Token token;
    long long block_index = 0;
    long long cubit_axis = 0;
    bool geometry_spec_error = false;
    {
      token = token_stream->Shift();
      std::string mesh_type = token.As_String();

      if(mesh_type == "RECTILINEAR" || mesh_type == "2DR" || mesh_type == "3DR" || mesh_type == "2DC"){
        Inline_Mesh_Desc::addDisc( new Cartesian_Inline_Mesh_Desc(sdimension ));
        Inline_Mesh_Desc::im_static_storage->inline_geometry_type = INLINE_CARTESIAN;
        if(has_specified_geometry_type) geometry_spec_error = true; else has_specified_geometry_type = true;
      }
      else if(mesh_type == "SPHERICAL"){
        Inline_Mesh_Desc::addDisc( new Spherical_Inline_Mesh_Desc(sdimension ));
        Inline_Mesh_Desc::im_static_storage->inline_geometry_type = INLINE_SPHERICAL;
        if(has_specified_geometry_type) geometry_spec_error = true; else has_specified_geometry_type = true;
      }
      else if(mesh_type == "CYLINDRICAL"){
        Inline_Mesh_Desc::addDisc( new Cylindrical_Inline_Mesh_Desc(sdimension ));
        Inline_Mesh_Desc::im_static_storage->inline_geometry_type = INLINE_CYLINDRICAL;
        if(has_specified_geometry_type) geometry_spec_error = true; else has_specified_geometry_type = true;
      }
      else if(mesh_type == "CUBIT RADIAL"){
        Inline_Mesh_Desc::addDisc( new Radial_Inline_Mesh_Desc(sdimension ));
        Inline_Mesh_Desc::im_static_storage->inline_geometry_type = RADIAL;
        if(has_specified_geometry_type) geometry_spec_error = true; else has_specified_geometry_type = true;
      }
      else if(mesh_type == "RADIAL"){
        Inline_Mesh_Desc::addDisc( new Radial_Inline_Mesh_Desc(sdimension ));
        Inline_Mesh_Desc::im_static_storage->inline_geometry_type = RADIAL;
        if(has_specified_geometry_type) geometry_spec_error = true; else has_specified_geometry_type = true;
      }
      else if(mesh_type == "BRICK"){
        Inline_Mesh_Desc::addDisc( new Brick_Inline_Mesh_Desc(sdimension ));
        Inline_Mesh_Desc::im_static_storage->inline_geometry_type = RADIAL;
        if(has_specified_geometry_type) geometry_spec_error = true; else has_specified_geometry_type = true;
      }
      else if(mesh_type == "CUBIT RADIAL TRISECTION"){
        Inline_Mesh_Desc::addDisc( new Radial_Trisection_Inline_Mesh_Desc(sdimension ));
        Inline_Mesh_Desc::im_static_storage->inline_geometry_type = RADIAL_TRISECTION;
        if(has_specified_geometry_type) geometry_spec_error = true; else has_specified_geometry_type = true;
      }
      else if(mesh_type == "RADIAL TRISECTION"){
        Inline_Mesh_Desc::addDisc( new Radial_Trisection_Inline_Mesh_Desc(sdimension ));
        Inline_Mesh_Desc::im_static_storage->inline_geometry_type = RADIAL_TRISECTION;
        if(has_specified_geometry_type) geometry_spec_error = true; else has_specified_geometry_type = true;
      }
      else if(mesh_type == "SET ASSIGN"){
      }
      else if(mesh_type == "USER DEFINED ELEMENT DENSITY"){
        Parse_User_Defined_Element_Density(token_stream,0);
      }
      else if(mesh_type == "USER DEFINED GEOMETRY TRANSFORMATION"){
        Parse_User_Defined_Geometry_Transformation(token_stream,0);
      }
      else if(mesh_type == "DECOMPOSITION STRATEGY"){
        Parse_Decomposition_Strategy(token_stream,0);
      }
      else if(mesh_type == "TOPOLOGY MODIFICATION"){
        Parse_Topology_Modification(token_stream,0);
      }
      else{
        //DMH trouble here Inline_Mesh_Desc::im_static_storage not yet in existence
        token_stream->Parse_Error("incorrect keyword ",
            "expected a geometry type (CYLINDRICAL,SPHERICAL,RECTILINEAR, RADIAL, RADIAL TRISECTION, USER DEFINED ELEMENT DENSITY, USER DEFINED GEOMETRY TRANSFORMATION, DECOMPOSITION STRATEGY) or SET ASSIGN");
      }
    }

    if(geometry_spec_error)  {
      token_stream->Parse_Error("incorrect keyword ",
                                "cannot specify more than one geometry type in a single mesh object");
    }


    Keyword parameter_table[] = {
      {"STARTING BLOCK ID OFFSET", P_STARTING_BLOCK_NUMBER, Get_Real_Token},
      {"OFFSET", P_OFFSET, Get_Real_Token},
      {"GMIN", P_GMIN, Get_Real_Token},
      {"GMAX", P_GMAX, Get_Real_Token},
      {"NTHETA", P_NTHETA, Get_Real_Token},
      {"NUMZ", P_NUMZ, Get_Real_Token},
      {"NUMR", P_NUMR, Get_Real_Token},
      {"NUMX", P_NUMR, Get_Real_Token},
      {"NUMA", P_NUMA, Get_Real_Token},
      {"NUMY", P_NUMA, Get_Real_Token},
      {"INITIAL RADIUS", P_INITIAL_RADIUS, Get_Real_Token},
      {"TRISECTION BLOCKS", P_TRISECTION_BLOCKS, Get_Real_Token},
      {"FIRST SIZE", P_FIRST_SIZE, Get_Real_Token},
      {"LAST SIZE", P_LAST_SIZE, Get_Real_Token},
      {"INTERVAL", P_INTERVAL, Get_Real_Token},
      {"ZBLOCK", P_ZBLOCK, Get_Real_Token},
      {"RBLOCK", P_RBLOCK, Get_Real_Token},
      {"XBLOCK", P_RBLOCK, Get_Real_Token},
      {"ABLOCK", P_ABLOCK, Get_Real_Token},
      {"YBLOCK", P_ABLOCK, Get_Real_Token},
      {"NPHI", P_NPHI, Get_Real_Token},
      {"NR", P_NR, Get_Real_Token},
      {"NX", P_NX, Get_Real_Token},
      {"NY", P_NY, Get_Real_Token},
      {"NZ", P_NZ, Get_Real_Token},
      {"BX", P_BX, Get_Real_Token},
      {"BY", P_BY, Get_Real_Token},
      {"BZ", P_BZ, Get_Real_Token},
      {"BTHETA", P_BTHETA, Get_Real_Token},
      {"BPHI", P_BPHI, Get_Real_Token},
      {"BR", P_BR, Get_Real_Token},
      {"RI", P_RI, Get_Real_Token},
      {"RO", P_RO, Get_Real_Token},
      {"THETA", P_THETA, Get_Real_Token},
      {"MINTHETA", P_MINTHETA, Get_Real_Token},
      {"PHI", P_PHI, Get_Real_Token},
      {"XMIN", P_XMIN, Get_Real_Token},
      {"YMIN", P_YMIN, Get_Real_Token},
      {"ZMIN", P_ZMIN, Get_Real_Token},
      {"ZMAX", P_ZMAX, Get_Real_Token},
      {"AMR", P_AMR, Get_Real_Token},
      {"NODESET", P_NODESET, Get_Real_Token},
      {"SIDESET", P_SIDESET, Get_Real_Token},
      {"BLOCK NODESET", P_BLOCK_NODESET, Get_Real_Token},
      {"BLOCK SIDESET", P_BLOCK_SIDESET, Get_Real_Token},
      {"DEBUG", P_DEBUG, Get_Real_Token},
      {"SQUARED", P_SQUARED, Get_Real_Token},
      {"ENFORCE PERIODIC", P_ENFORCE_PERIODIC, Get_Real_Token},
      {"NUMELARRAY", P_NUMELARRAY, Get_Real_Token},
      {"TRANSITION RADIUS", P_TRANSITION_RADIUS, Get_Real_Token},
      {"NOCUTSI", P_NOCUTSI, Get_Real_Token},
      {"NOCUTSJ", P_NOCUTSJ, Get_Real_Token},
      {"NOCUTSK", P_NOCUTSK, Get_Real_Token},
    };

    long long N = sizeof(parameter_table)/sizeof(Keyword);

    qsort(parameter_table,
        N,
        sizeof(Keyword),
        PAMGEN_Keyword_Compare);

    token = token_stream->Shift();

    while(token.Type() != TK_END){
      // check token for end of material model input


      // evaluate model parameter and value

      if (token.Type() != TK_IDENTIFIER) {
        token_stream->Parse_Error("an inline mesh keyword was expected",
            PAMGEN_NEVADA::Concatenate_Legal_Commands(parameter_table,N));
      }

      char *name = (char*) token.As_String();
      const Keyword *match = (Keyword*)bsearch( name,
          (char*)parameter_table,
          N,
          sizeof(Keyword),
          PAMGEN_Cstring_Keyword_Compare );
      if(!match){
        token_stream->Parse_Error(std::string("Unknown identifier in this context: ") + name);
      }

      PAMGEN_NEVADA::Check_for_Ambiguity( token_stream,
          name,
          match,
          (Keyword const *) parameter_table,
          N );
      long long param_id = match->argument;

      switch(param_id){
        case P_STARTING_BLOCK_NUMBER:{
			 int start_block_id = token_stream->Parse_Integer();

                         if(start_block_id <=0){
                           std::stringstream ss;
                           ss << "The  starting block number " << start_block_id << " must be positive.";
                           token_stream->Semantics_Error(ss.str());
                         }
			 Inline_Mesh_Desc::im_static_storage->inline_block_start = start_block_id;

                        break;}
        case P_OFFSET:{
                        Inline_Mesh_Desc::im_static_storage->inline_offset[0] = token_stream->Parse_Real();
                        Inline_Mesh_Desc::im_static_storage->inline_offset[1] = token_stream->Parse_Real();
                        if(sdimension == 3){
                          Inline_Mesh_Desc::im_static_storage->inline_offset[2] = token_stream->Parse_Real();
                        }
                        break;}
        case P_GMIN:{
                      Inline_Mesh_Desc::im_static_storage->inline_gmin[0] = token_stream->Parse_Real();
                      Inline_Mesh_Desc::im_static_storage->inline_gmin[1] = token_stream->Parse_Real();
                      if(sdimension == 3){
                        Inline_Mesh_Desc::im_static_storage->inline_gmin[2] = token_stream->Parse_Real();
                      }
                      break;}
        case P_GMAX:{
                      Inline_Mesh_Desc::im_static_storage->inline_gmax[0] = token_stream->Parse_Real();
                      Inline_Mesh_Desc::im_static_storage->inline_gmax[1] = token_stream->Parse_Real();
                      if(sdimension == 3){
                        Inline_Mesh_Desc::im_static_storage->inline_gmax[2] = token_stream->Parse_Real();
                      }
                      break;}
        case P_RI:{
                    Inline_Mesh_Desc::im_static_storage->inline_gmin[0] = token_stream->Parse_Real();
                    break;}
        case P_RO:{
                    Inline_Mesh_Desc::im_static_storage->inline_gmax[0] = token_stream->Parse_Real();
                    break;}
        case P_THETA:{
                       Inline_Mesh_Desc::im_static_storage->inline_gmax[1] = token_stream->Parse_Real();
                       if(sdimension == 2){
                         if(Inline_Mesh_Desc::im_static_storage->inline_gmax[1] >= 360.0){
                           Inline_Mesh_Desc::im_static_storage->inline_gmax[1] = 360.0;
                           Inline_Mesh_Desc::im_static_storage->periodic_j = true;
                         }
                       }
                       else{
                         if(Inline_Mesh_Desc::im_static_storage->inline_gmax[1] >= 180.0){
                           Inline_Mesh_Desc::im_static_storage->inline_gmax[1] = 180.0;
                           // An in this case there should be a degenerate type
                           // bc applied along the x-axis
                         }
                       }
                       break;}
        case P_PHI:{
                     Inline_Mesh_Desc::im_static_storage->inline_gmax[2] = token_stream->Parse_Real();
                     if(Inline_Mesh_Desc::im_static_storage->inline_gmax[2] >= 360.0){
                       Inline_Mesh_Desc::im_static_storage->inline_gmax[2] = 360.0;
                       Inline_Mesh_Desc::im_static_storage->periodic_k = true;
                     }
                     break;}
        case P_XMIN:{
                      Inline_Mesh_Desc::im_static_storage->inline_gmin[0] = token_stream->Parse_Real();
                      break;}
        case P_YMIN:{
                      Inline_Mesh_Desc::im_static_storage->inline_gmin[1] = token_stream->Parse_Real();
                      break;}
        case P_MINTHETA:{
                          Inline_Mesh_Desc::im_static_storage->inline_gmin[1] = token_stream->Parse_Real();
                          break;}
        case P_ZMIN:{
                      Inline_Mesh_Desc::im_static_storage->inline_gmin[2] = token_stream->Parse_Real();
                      break;}
        case P_ZMAX:{
                      Inline_Mesh_Desc::im_static_storage->inline_gmax[2] = token_stream->Parse_Real();
                      break;}
        case P_NX:{
                    Inline_Mesh_Desc::im_static_storage->inline_n[0] = token_stream->Parse_Integer();
                    break;}
        case P_NY:{
                    Inline_Mesh_Desc::im_static_storage->inline_n[1] = token_stream->Parse_Integer();
                    break;}
        case P_NZ:{
                    if(sdimension == 2){
                      token_stream->Semantics_Error(std::string("NZ may not be specified in 2D simulations."));
                    }
                    else{
                      Inline_Mesh_Desc::im_static_storage->inline_n[2] = token_stream->Parse_Integer();
                    }
                    break;}
        case P_NTHETA:{
                        Inline_Mesh_Desc::im_static_storage->inline_n[1] = token_stream->Parse_Integer();
                        break;}


        case P_NUMZ:{
                      cubit_axis = 2;
                      Inline_Mesh_Desc::im_static_storage->inline_b[2] = token_stream->Parse_Integer();

                      if(Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis] = new Real[Inline_Mesh_Desc::im_static_storage->inline_b[2]];

                      if(Inline_Mesh_Desc::im_static_storage->c_block_dist[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->c_block_dist[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->c_block_dist[cubit_axis] = new Real[Inline_Mesh_Desc::im_static_storage->inline_b[2]+1];

                      if(Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis] = new Real[Inline_Mesh_Desc::im_static_storage->inline_b[2]];

                      if(Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis] = new Real[Inline_Mesh_Desc::im_static_storage->inline_b[2]];

                      if(Inline_Mesh_Desc::im_static_storage->interval[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->interval[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->interval[cubit_axis] = new long long[Inline_Mesh_Desc::im_static_storage->inline_b[2]];

                      for(long long i = 0; i < Inline_Mesh_Desc::im_static_storage->inline_b[2]; i ++){
                        Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis][i] = 0.;
                        Inline_Mesh_Desc::im_static_storage->c_block_dist[cubit_axis][i] = 0.;
                        Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis][i] = 0.;
                        Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis][i] = 0.;
                        Inline_Mesh_Desc::im_static_storage->interval[cubit_axis][i] = 0;
                      }

                      break;}
        case P_NUMR:{
                      cubit_axis = 0;
                      Inline_Mesh_Desc::im_static_storage->inline_b[0] = token_stream->Parse_Integer();

                      if(Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis] = new Real[Inline_Mesh_Desc::im_static_storage->inline_b[0]];

                      if(Inline_Mesh_Desc::im_static_storage->c_block_dist[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->c_block_dist[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->c_block_dist[cubit_axis] = new Real[Inline_Mesh_Desc::im_static_storage->inline_b[0]+1];

                      if(Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis] = new Real[Inline_Mesh_Desc::im_static_storage->inline_b[0]];

                      if(Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis] = new Real[Inline_Mesh_Desc::im_static_storage->inline_b[0]];

                      if(Inline_Mesh_Desc::im_static_storage->interval[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->interval[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->interval[cubit_axis] = new long long[Inline_Mesh_Desc::im_static_storage->inline_b[0]];

                      for(long long i = 0; i < Inline_Mesh_Desc::im_static_storage->inline_b[0]; i ++){
                        Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis][i] = 0.;
                        Inline_Mesh_Desc::im_static_storage->c_block_dist[cubit_axis][i] = 0.;
                        Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis][i] = 0.;
                        Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis][i] = 0.;
                        Inline_Mesh_Desc::im_static_storage->interval[cubit_axis][i] = 0;
                      }

                      break;}
        case P_NUMA:{
                      cubit_axis = 1;
                      Inline_Mesh_Desc::im_static_storage->inline_b[1] = token_stream->Parse_Integer();

                      if(Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis] = new Real[Inline_Mesh_Desc::im_static_storage->inline_b[1]];

                      if(Inline_Mesh_Desc::im_static_storage->c_block_dist[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->c_block_dist[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->c_block_dist[cubit_axis] = new Real[Inline_Mesh_Desc::im_static_storage->inline_b[1]+1];

                      if(Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis] = new Real[Inline_Mesh_Desc::im_static_storage->inline_b[1]];

                      if(Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis] = new Real[Inline_Mesh_Desc::im_static_storage->inline_b[1]];

                      if(Inline_Mesh_Desc::im_static_storage->interval[cubit_axis])delete Inline_Mesh_Desc::im_static_storage->interval[cubit_axis];
                      Inline_Mesh_Desc::im_static_storage->interval[cubit_axis] = new long long[Inline_Mesh_Desc::im_static_storage->inline_b[1]];
                      for(long long i = 0; i < Inline_Mesh_Desc::im_static_storage->inline_b[1]; i ++){
                        Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis][i] = 0.;
                        Inline_Mesh_Desc::im_static_storage->c_block_dist[cubit_axis][i] = 0.;
                        Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis][i] = 0.;
                        Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis][i] = 0.;
                        Inline_Mesh_Desc::im_static_storage->interval[cubit_axis][i] = 0;
                      }

                      break;}

        case P_INITIAL_RADIUS:{
                                Inline_Mesh_Desc::im_static_storage->inline_gmin[0] = token_stream->Parse_Real();
                                break;}

        case P_TRISECTION_BLOCKS:{
                                   Inline_Mesh_Desc::im_static_storage->trisection_blocks = token_stream->Parse_Integer();
                                   break;}


        case P_FIRST_SIZE:{
                            Inline_Mesh_Desc::im_static_storage->first_size[cubit_axis][block_index] = token_stream->Parse_Real();
                            break;}
        case P_LAST_SIZE:{
                           Inline_Mesh_Desc::im_static_storage->last_size[cubit_axis][block_index] = token_stream->Parse_Real();
                           break;}
        case P_INTERVAL:{
                          Inline_Mesh_Desc::im_static_storage->interval[cubit_axis][block_index] = token_stream->Parse_Integer();
                          break;}


        case P_ZBLOCK:{
                        block_index = token_stream->Parse_Integer()-1;
                        if(block_index+1 > Inline_Mesh_Desc::im_static_storage->inline_b[2]){
                          std::stringstream ss;
                          ss << "Zblock index out of range, should run from 1 to " << Inline_Mesh_Desc::im_static_storage->inline_b[2] << ".";
                          token_stream->Semantics_Error(ss.str());
                        }
                        Real val = token_stream->Parse_Real();
                        if(val <= 0.){
                          std::stringstream ss;
                          ss << "Block size must be positive." << Inline_Mesh_Desc::im_static_storage->inline_b[2] << ".";
                          token_stream->Semantics_Error(ss.str());
                        } 
                        Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis][block_index] = val;
                        break;}
        case P_RBLOCK:{
                        block_index = token_stream->Parse_Integer()-1;
                        if(block_index+1 > Inline_Mesh_Desc::im_static_storage->inline_b[0]){
                          std::stringstream ss;
                          ss << "Rblock/Xblock index out of range, should run from 1 to " << Inline_Mesh_Desc::im_static_storage->inline_b[0] << ".";
                          token_stream->Semantics_Error(ss.str());
                        }
                        Real val = token_stream->Parse_Real();
                        if(val <= 0.){
                          std::stringstream ss;
                          ss << "Block size must be positive." << Inline_Mesh_Desc::im_static_storage->inline_b[2] << ".";
                          token_stream->Semantics_Error(ss.str());
                        } 
                        Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis][block_index] = val;
                        break;}
        case P_ABLOCK:{
                        block_index = token_stream->Parse_Integer()-1;
                        if(block_index+1 > Inline_Mesh_Desc::im_static_storage->inline_b[1]){
                          std::stringstream ss;
                          ss << "Ablock/Yblock index out of range, should run from 1 to " << Inline_Mesh_Desc::im_static_storage->inline_b[1] << ".";
                          token_stream->Semantics_Error(ss.str());
                        }
                        Real val = token_stream->Parse_Real();
                        if(val <= 0.){
                          std::stringstream ss;
                          ss << "Block size must be positive." << Inline_Mesh_Desc::im_static_storage->inline_b[2] << ".";
                          token_stream->Semantics_Error(ss.str());
                        } 
                        Inline_Mesh_Desc::im_static_storage->block_dist[cubit_axis][block_index] = val;
                        break;}


        case P_NPHI:{
                      Inline_Mesh_Desc::im_static_storage->inline_n[2] = token_stream->Parse_Integer();
                      break;}
        case P_NR:{
                    Inline_Mesh_Desc::im_static_storage->inline_n[0] = token_stream->Parse_Integer();
                    break;}
        case P_BX:{
                    Inline_Mesh_Desc::im_static_storage->inline_b[0] = token_stream->Parse_Integer();
                    break;}
        case P_BY:{
                    Inline_Mesh_Desc::im_static_storage->inline_b[1] = token_stream->Parse_Integer();
                    break;}
        case P_BZ:{
                    if(sdimension == 2){
                      token_stream->Semantics_Error(std::string("BZ may not be specified in 2D simulations."));
                    }
                    else{
                      Inline_Mesh_Desc::im_static_storage->inline_b[2] = token_stream->Parse_Integer();
                    }
                    break;}
        case P_BTHETA:{
                        Inline_Mesh_Desc::im_static_storage->inline_b[1] = token_stream->Parse_Integer();
                        break;}
        case P_BPHI:{
                      Inline_Mesh_Desc::im_static_storage->inline_b[2] = token_stream->Parse_Integer();
                      break;}
        case P_BR:{
                    Inline_Mesh_Desc::im_static_storage->inline_b[0] = token_stream->Parse_Integer();
                    break;}
        case P_AMR:{

                     break;}
        case P_DEBUG:{

                       break;}
        case P_SQUARED:{
                         Inline_Mesh_Desc::im_static_storage->try_squared = true;
                         break;}
        case P_ENFORCE_PERIODIC:{
                                  Inline_Mesh_Desc::im_static_storage->enforce_periodic = true;
                                  break;}
        case P_NUMELARRAY:{
                            Inline_Mesh_Desc::im_static_storage->inc_nels[0] = token_stream->Parse_Integer();
                            Inline_Mesh_Desc::im_static_storage->inc_nels[1] = token_stream->Parse_Integer();
                            Inline_Mesh_Desc::im_static_storage->inc_nels[2] = token_stream->Parse_Integer();
                            break;}
        case P_TRANSITION_RADIUS:{
                                   Inline_Mesh_Desc::im_static_storage->transition_radius = token_stream->Parse_Real();
                                   break;}
        case P_NOCUTSI:{
                         Inline_Mesh_Desc::im_static_storage->inc_nocuts[0] = true;
                         break;}
        case P_NOCUTSJ:{
                         Inline_Mesh_Desc::im_static_storage->inc_nocuts[1] = true;
                         break;}
        case P_NOCUTSK:{
                         Inline_Mesh_Desc::im_static_storage->inc_nocuts[2] = true;
                         break;}
        case P_NODESET:{
                         const char* type_token;

                         token = token_stream->Shift();
                         type_token  = token.As_String();
                         long long the_id = (long long)token_stream->Parse_Integer();
                         if(the_id <=0){
                           std::stringstream ss;
                           ss << "Nodeset ids must be positive " << the_id;
                           token_stream->Semantics_Error(ss.str());
                         }
                         Topo_Loc the_loc = NUM_TOPO_CONNECTIONS;
                         // case insensitive string compare
                         if(!strncasecmp(type_token,"ihi",3))the_loc = PLUS_I;
                         else if(!strncasecmp(type_token,"all",3))the_loc = ALL_NODES;
                         else if(!strncasecmp(type_token,"xhi",3))the_loc = PLUS_I;
                         else if(!strncasecmp(type_token,"rhi",3))the_loc = PLUS_I;
                         else if(!strncasecmp(type_token,"jhi",3))the_loc = PLUS_J;
                         else if(!strncasecmp(type_token,"yhi",3))the_loc = PLUS_J;
                         else if(!strncasecmp(type_token,"thetahi",7))the_loc = PLUS_J;
                         else if(!strncasecmp(type_token,"khi",3))the_loc = PLUS_K;
                         else if(!strncasecmp(type_token,"zhi",3))the_loc = PLUS_K;
                         else if(!strncasecmp(type_token,"ilo",3) ||
                             !strncasecmp(type_token,"xlo",3) ||
                             !strncasecmp(type_token,"rlo",3)){
                           the_loc = MINUS_I;
                           if(Inline_Mesh_Desc::im_static_storage->inline_geometry_type == RADIAL_TRISECTION){
                             the_loc = Z_AXIS;
                           }
                         }
                         else if(!strncasecmp(type_token,"jlo",3))the_loc = MINUS_J;
                         else if(!strncasecmp(type_token,"ylo",3))the_loc = MINUS_J;
                         else if(!strncasecmp(type_token,"thetalo",7))the_loc = MINUS_J;
                         else if(!strncasecmp(type_token,"klo",3))the_loc = MINUS_K;
                         else if(!strncasecmp(type_token,"zlo",3))the_loc = MINUS_K;
                         else if(!strncasecmp(type_token,"zaxis",5))the_loc = Z_AXIS;
                         else if(!strncasecmp(type_token,"e00",3))the_loc = EDGE0;
                         else if(!strncasecmp(type_token,"e01",3))the_loc = EDGE1;
                         else if(!strncasecmp(type_token,"e02",3))the_loc = EDGE2;
                         else if(!strncasecmp(type_token,"e03",3))the_loc = EDGE3;
                         else if(!strncasecmp(type_token,"e04",3))the_loc = EDGE4;
                         else if(!strncasecmp(type_token,"e05",3))the_loc = EDGE5;
                         else if(!strncasecmp(type_token,"e06",3))the_loc = EDGE6;
                         else if(!strncasecmp(type_token,"e07",3))the_loc = EDGE7;
                         else if(!strncasecmp(type_token,"e08",3))the_loc = EDGE8;
                         else if(!strncasecmp(type_token,"e09",3))the_loc = EDGE9;
                         else if(!strncasecmp(type_token,"e10",3))the_loc = EDGE10;
                         else if(!strncasecmp(type_token,"e11",3))the_loc = EDGE11;
                         else if(!strncasecmp(type_token,"v00",3))the_loc = VERTEX0;
                         else if(!strncasecmp(type_token,"v01",3))the_loc = VERTEX1;
                         else if(!strncasecmp(type_token,"v02",3))the_loc = VERTEX2;
                         else if(!strncasecmp(type_token,"v03",3))the_loc = VERTEX3;
                         else if(!strncasecmp(type_token,"v04",3))the_loc = VERTEX4;
                         else if(!strncasecmp(type_token,"v05",3))the_loc = VERTEX5;
                         else if(!strncasecmp(type_token,"v06",3))the_loc = VERTEX6;
                         else if(!strncasecmp(type_token,"v07",3))the_loc = VERTEX7;
                         else {
                           std::stringstream ss;
                           ss << "Inappropriate token in nodesets:" << type_token;
                           token_stream->Semantics_Error(ss.str());
                           break;
                         }

                         if(sdimension == 2 && (the_loc == MINUS_K || the_loc == PLUS_K))break;

                         PG_BC_Specification * the_ns = Inline_Mesh_Desc::im_static_storage->getNodeset_by_id(the_id);
                         if(the_ns){
                           the_ns->addEntry(the_loc,false,0);
                         }
                         else{
                           PG_BC_Specification * bcs = new PG_BC_Specification(the_id,the_loc,false,0);
                           Inline_Mesh_Desc::im_static_storage->nodeset_list.push_back(bcs);
                         }

                         break;}
        case P_BLOCK_NODESET:{
                               const char* type_token;

                               token = token_stream->Shift();
                               type_token  = token.As_String();
                               long long the_id = (long long)token_stream->Parse_Integer();
                               if(the_id <=0){
                                 std::stringstream ss;
                                 ss << "Nodeset ids must be positive " << the_id;
                                 token_stream->Semantics_Error(ss.str());
                               }
                               long long the_block = (long long) token_stream->Parse_Integer();

                               Topo_Loc the_loc = NUM_TOPO_CONNECTIONS;
                               // case insensitive string compare
                               if(!strncasecmp(type_token,"ihi",3))the_loc = PLUS_I;
                               else if(!strncasecmp(type_token,"xhi",3))the_loc = PLUS_I;
                               else if(!strncasecmp(type_token,"rhi",3))the_loc = PLUS_I;
                               else if(!strncasecmp(type_token,"jhi",3))the_loc = PLUS_J;
                               else if(!strncasecmp(type_token,"yhi",3))the_loc = PLUS_J;
                               else if(!strncasecmp(type_token,"thetahi",7))the_loc = PLUS_J;
                               else if(!strncasecmp(type_token,"khi",3))the_loc = PLUS_K;
                               else if(!strncasecmp(type_token,"zhi",3))the_loc = PLUS_K;
                               else if(!strncasecmp(type_token,"ilo",3))the_loc = MINUS_I;
                               else if(!strncasecmp(type_token,"xlo",3))the_loc = MINUS_I;
                               else if(!strncasecmp(type_token,"rlo",3))the_loc = MINUS_I;
                               else if(!strncasecmp(type_token,"jlo",3))the_loc = MINUS_J;
                               else if(!strncasecmp(type_token,"ylo",3))the_loc = MINUS_J;
                               else if(!strncasecmp(type_token,"thetalo",7))the_loc = MINUS_J;
                               else if(!strncasecmp(type_token,"klo",3))the_loc = MINUS_K;
                               else if(!strncasecmp(type_token,"zlo",3))the_loc = MINUS_K;
                               else if(!strncasecmp(type_token,"e00",3))the_loc = EDGE0;
                               else if(!strncasecmp(type_token,"e01",3))the_loc = EDGE1;
                               else if(!strncasecmp(type_token,"e02",3))the_loc = EDGE2;
                               else if(!strncasecmp(type_token,"e03",3))the_loc = EDGE3;
                               else if(!strncasecmp(type_token,"e04",3))the_loc = EDGE4;
                               else if(!strncasecmp(type_token,"e05",3))the_loc = EDGE5;
                               else if(!strncasecmp(type_token,"e06",3))the_loc = EDGE6;
                               else if(!strncasecmp(type_token,"e07",3))the_loc = EDGE7;
                               else if(!strncasecmp(type_token,"e08",3))the_loc = EDGE8;
                               else if(!strncasecmp(type_token,"e09",3))the_loc = EDGE9;
                               else if(!strncasecmp(type_token,"e10",3))the_loc = EDGE10;
                               else if(!strncasecmp(type_token,"e11",3))the_loc = EDGE11;
                               else if(!strncasecmp(type_token,"v00",3))the_loc = VERTEX0;
                               else if(!strncasecmp(type_token,"v01",3))the_loc = VERTEX1;
                               else if(!strncasecmp(type_token,"v02",3))the_loc = VERTEX2;
                               else if(!strncasecmp(type_token,"v03",3))the_loc = VERTEX3;
                               else if(!strncasecmp(type_token,"v04",3))the_loc = VERTEX4;
                               else if(!strncasecmp(type_token,"v05",3))the_loc = VERTEX5;
                               else if(!strncasecmp(type_token,"v06",3))the_loc = VERTEX6;
                               else if(!strncasecmp(type_token,"v07",3))the_loc = VERTEX7;
                               else {
                                 std::stringstream ss;
                                 ss << "Inappropriate token in block nodesets:" << type_token;
                                 token_stream->Semantics_Error(ss.str());
                                 break;
                               }

                               if(sdimension == 2 && (the_loc == MINUS_K || the_loc == PLUS_K))break;

                               PG_BC_Specification * the_ns = Inline_Mesh_Desc::im_static_storage->getNodeset_by_id(the_id);
                               if(the_ns){
                                 the_ns->addEntry(the_loc,true,the_block);
                               }
                               else{
                                 PG_BC_Specification * bcs = new PG_BC_Specification(the_id,the_loc,true,the_block);
                                 Inline_Mesh_Desc::im_static_storage->nodeset_list.push_back(bcs);
                               }
                               break;}
        case P_SIDESET:{
                         const char* type_token;
                         token = token_stream->Shift();
                         type_token  = token.As_String();
                         long long the_id = (long long)token_stream->Parse_Integer();
                         if(the_id <=0){
                           std::stringstream ss;
                           ss << "Sideset ids must be positive " << the_id;
                           token_stream->Semantics_Error(ss.str());
                         }
                         Topo_Loc the_loc = NUM_TOPO_CONNECTIONS;
                         // case insensitive string compare
                         if(!strncasecmp(type_token,"ihi",3))the_loc = PLUS_I;
                         else if(!strncasecmp(type_token,"xhi",3))the_loc = PLUS_I;
                         else if(!strncasecmp(type_token,"rhi",3))the_loc = PLUS_I;
                         else if(!strncasecmp(type_token,"jhi",3))the_loc = PLUS_J;
                         else if(!strncasecmp(type_token,"yhi",3))the_loc = PLUS_J;
                         else if(!strncasecmp(type_token,"thetahi",7))the_loc = PLUS_J;
                         else if(!strncasecmp(type_token,"khi",3))the_loc = PLUS_K;
                         else if(!strncasecmp(type_token,"zhi",3))the_loc = PLUS_K;
                         else if(!strncasecmp(type_token,"ilo",3) ||
                             !strncasecmp(type_token,"xlo",3) ||
                             !strncasecmp(type_token,"rlo",3)){
                           the_loc = MINUS_I;
                           if(Inline_Mesh_Desc::im_static_storage->inline_geometry_type == RADIAL_TRISECTION){
                             std::stringstream ss;
                             ss << "Sideset location disallowed for Trisection Meshes:" << type_token;
                             token_stream->Semantics_Error(ss.str());
                           }
                         }
                         else if(!strncasecmp(type_token,"jlo",3))the_loc = MINUS_J;
                         else if(!strncasecmp(type_token,"ylo",3))the_loc = MINUS_J;
                         else if(!strncasecmp(type_token,"thetalo",7))the_loc = MINUS_J;
                         else if(!strncasecmp(type_token,"klo",3))the_loc = MINUS_K;
                         else if(!strncasecmp(type_token,"zlo",3))the_loc = MINUS_K;
                         else {
                           std::stringstream ss;
                           ss << "Inappropriate token in sidesets:" << type_token;
                           token_stream->Semantics_Error(ss.str());
                           break;
                         }



                         if(sdimension == 2 && (the_loc == MINUS_K || the_loc == PLUS_K))break;

                         PG_BC_Specification * the_ss = Inline_Mesh_Desc::im_static_storage->getSideset_by_id(the_id);
                         if(the_ss){
                           the_ss->addEntry(the_loc,false,0);
                         }
                         else{
                           PG_BC_Specification * bcs = new PG_BC_Specification(the_id,the_loc,false,0);
                           Inline_Mesh_Desc::im_static_storage->sideset_list.push_back(bcs);
                         }


                         break;}
        case P_BLOCK_SIDESET:{
                               const char* type_token;
                               token = token_stream->Shift();
                               type_token  = token.As_String();
                               long long the_id = (long long)token_stream->Parse_Integer();
                               if(the_id <=0){
                                 std::stringstream ss;
                                 ss << "Sideset ids must be positive " << the_id;
                                 token_stream->Semantics_Error(ss.str());
                               }
                               long long the_block = (long long) token_stream->Parse_Integer();

                               Topo_Loc the_loc = NUM_TOPO_CONNECTIONS;
                               // case insensitive string compare
                               if(!strncasecmp(type_token,"ihi",3))the_loc = PLUS_I;
                               else if(!strncasecmp(type_token,"xhi",3))the_loc = PLUS_I;
                               else if(!strncasecmp(type_token,"rhi",3))the_loc = PLUS_I;
                               else if(!strncasecmp(type_token,"jhi",3))the_loc = PLUS_J;
                               else if(!strncasecmp(type_token,"yhi",3))the_loc = PLUS_J;
                               else if(!strncasecmp(type_token,"thetahi",7))the_loc = PLUS_J;
                               else if(!strncasecmp(type_token,"khi",3))the_loc = PLUS_K;
                               else if(!strncasecmp(type_token,"zhi",3))the_loc = PLUS_K;
                               else if(!strncasecmp(type_token,"ilo",3))the_loc = MINUS_I;
                               else if(!strncasecmp(type_token,"xlo",3))the_loc = MINUS_I;
                               else if(!strncasecmp(type_token,"rlo",3))the_loc = MINUS_I;
                               else if(!strncasecmp(type_token,"jlo",3))the_loc = MINUS_J;
                               else if(!strncasecmp(type_token,"ylo",3))the_loc = MINUS_J;
                               else if(!strncasecmp(type_token,"thetalo",7))the_loc = MINUS_J;
                               else if(!strncasecmp(type_token,"klo",3))the_loc = MINUS_K;
                               else if(!strncasecmp(type_token,"zlo",3))the_loc = MINUS_K;
                               else {
                                 std::stringstream ss;
                                 ss << "Inappropriate token in block sidesets:" << type_token;
                                 token_stream->Semantics_Error(ss.str());
                                 break;
                               }
                               if(sdimension == 2 &&(the_loc == MINUS_K || the_loc == PLUS_K))break;

                               if(the_loc == MINUS_I && the_block == 1){
                                 if(Inline_Mesh_Desc::im_static_storage->inline_geometry_type == RADIAL_TRISECTION){
                                   std::stringstream ss;
                                   ss << "Block sideset location disallowed on this block for Trisection Meshes:" << type_token;
                                   token_stream->Semantics_Error(ss.str());
                                 }
                               }

                               PG_BC_Specification * the_ss = Inline_Mesh_Desc::im_static_storage->getSideset_by_id(the_id);

                               if(the_ss){
                                 the_ss->addEntry(the_loc,true,the_block);
                               }
                               else{
                                 PG_BC_Specification * bcs = new PG_BC_Specification(the_id,the_loc,true,the_block);
                                 Inline_Mesh_Desc::im_static_storage->sideset_list.push_back(bcs);
                               }
                               break;}
        default:{
                  std::stringstream ss;
                  ss << "Unknown identifier in this context: " << name;
                  token_stream->Parse_Error(ss.str());
                }
      }
      token = token_stream->Shift();
    }

    Token next = token_stream->Lookahead();
    if(next.Type() != TK_END){
      Parse_Inline_Mesh_Tok(token_stream,value);
      return token;
    }

    std::string intervalResults = Inline_Mesh_Desc::im_static_storage->Calc_Intervals();
    if(!intervalResults.empty()){
      token_stream->Parse_Error(intervalResults);
    }

    token = token_stream->Shift();
    return token;
  }

  /*****************************************************************************/
  Token Parse_User_Defined_Element_Density(Token_Stream *token_stream, int)
    /*****************************************************************************/
  {
    assert(token_stream);

    // USER DEFINED ELEMENT DENSITY, I|J|K
    //   "
    //   function using coordinates[0,1,2] and assigning values
    //   to return_value[0,1,2]
    //   "
    // END

    //get the variable that the function will affect
    long long function_index = 0;
    Token token = token_stream->Shift();
    std::string direction = token.As_String();
    if(token == "I" || token == "X" || token == "R"){
      function_index = 0;
    }
    else if(token == "J" || token == "Y" || token == "THETA"){
      function_index = 1;
    }
    else if(token == "K" || token == "Z" || token == "PHI"){
      function_index = 2;
    }
    else{
      token_stream->Parse_Error("expected I|J|K");
    }

    //get function body
    Token token2 = token_stream->Shift();
    if (token2.Type() != TK_STRING){
      token_stream->Parse_Error("expected a quoted C-language function body");
    }
    std::string funcBody = token2.As_Stripped_String();



    Token next = token_stream->Lookahead();
    if(next.Type() != TK_END){
      token_stream->Parse_Error("missing END after User Defined Element Density Function");
    }

    // may be replacing previously defined function.
    Element_Density_Function *EDF = Inline_Mesh_Desc::im_static_storage->Element_Density_Functions[function_index];
    if(EDF)delete EDF;
    std::stringstream estring;
    EDF = new Element_Density_Function(funcBody,direction,estring);
    if(!estring.str().empty()){
      token_stream->Parse_Error(estring.str());
    }
    Inline_Mesh_Desc::im_static_storage->Element_Density_Functions[function_index] = EDF;
    return Token();
  }

  /*****************************************************************************/
  Token Parse_User_Defined_Geometry_Transformation(Token_Stream *token_stream, int)
    /*****************************************************************************/
  {
    assert(token_stream);

    // USER DEFINED GEOMETRY TRANSFORMATION
    //   "
    //   function using coordinates[0,1,2] and assigning values
    //   to return_value[0,1,2]
    //   "
    // END

    //get the variable that the function will affect

    //get function body
    Token token2 = token_stream->Shift();
    if (token2.Type() != TK_STRING){
      token_stream->Parse_Error("expected a quoted C-language function body");
    }
    std::string funcBody = token2.As_Stripped_String();



    Token next = token_stream->Lookahead();
    if(next.Type() != TK_END){
      token_stream->Parse_Error("missing END after User Defined Geometry Transformation");
    }

    // may be replacing previously defined function.


    Geometry_Transform *EDF = Inline_Mesh_Desc::im_static_storage->Geometry_Transform_Function;
    if(EDF)delete EDF;
    std::stringstream estring;
    EDF = new Geometry_Transform(funcBody,estring);
    if(!estring.str().empty()){
      token_stream->Parse_Error(estring.str());
    }
    Inline_Mesh_Desc::im_static_storage->Geometry_Transform_Function = EDF;

    return Token();
  }

  /*****************************************************************************/
  Token Parse_Topology_Modification(Token_Stream *token_stream, int)
    /*****************************************************************************/
  {
    assert(token_stream);

    // TOPOLOGY MODIFICATION
    //  SUPPRESS BLOCK 1
    //  SUPPRESS BLOCK 2
    //
    //
    //
    // END


    Token next = token_stream->Lookahead();
    Token token = next;

    while(next.Type() != TK_END){
      token = token_stream->Shift();

      if(token == "SUPPRESS BLOCK"){
        int the_block = token_stream->Parse_Integer();
        Inline_Mesh_Desc::im_static_storage->addSuppressedBlock(the_block);
      }
      else {
        token_stream->Parse_Error("Expecting SUPPRESS BLOCK int");
      }
      next = token_stream->Lookahead();
    }
    return token;
  }
  /*****************************************************************************/
  Token Parse_Decomposition_Strategy(Token_Stream *token_stream, int)
    /*****************************************************************************/
  {
    assert(token_stream);

    // DECOMPOSITION STRATEGY
    //   BISECT | PROCESSOR LAYOUT | SEQUENTIAL
    //   I,3
    //   J,2
    //   k,4
    // END


    Token next = token_stream->Lookahead();
    Token token = next;

    while(next.Type() != TK_END){
      token = token_stream->Shift();

      if(token == "BISECT"){
        Inline_Mesh_Desc::im_static_storage->inline_decomposition_type = BISECTION;
      }
      if(token == "SEQUENTIAL"){
        Inline_Mesh_Desc::im_static_storage->inline_decomposition_type = SEQUENTIAL;
      }
      if(token == "RANDOM"){
        Inline_Mesh_Desc::im_static_storage->inline_decomposition_type = RANDOM;
      }
      else if(token == "PROCESSOR LAYOUT"){
        //This keyord is deprecated and should eventually be removed
        Inline_Mesh_Desc::im_static_storage->inline_decomposition_type = PROCESSOR_LAYOUT;
      }
      else if (token == "NUMPROCS I" || token == "NUMPROCS X" || token == "NUMPROCS R"){
        Inline_Mesh_Desc::im_static_storage->inline_decomposition_type = PROCESSOR_LAYOUT;
        Inline_Mesh_Desc::im_static_storage->inline_nprocs[0] = token_stream->Parse_Integer();
      }
      else if (token == "NUMPROCS J" || token == "NUMPROCS Y" || token == "NUMPROCS THETA"){
        Inline_Mesh_Desc::im_static_storage->inline_decomposition_type = PROCESSOR_LAYOUT;
        Inline_Mesh_Desc::im_static_storage->inline_nprocs[1] = token_stream->Parse_Integer();
      }
      else if (token == "NUMPROCS K" || token == "NUMPROCS Z" || token == "NUMPROCS PHI"){
        Inline_Mesh_Desc::im_static_storage->inline_decomposition_type = PROCESSOR_LAYOUT;
        Inline_Mesh_Desc::im_static_storage->inline_nprocs[2] = token_stream->Parse_Integer();
      }
      next = token_stream->Lookahead();
    }
    return token;
  }


  /*****************************************************************************/
  void Token_Stream::Parse_Error(const std::string &s, const std::string &v)
    /*****************************************************************************/
  {
    Semantics_Error(s, v);

    // Never returns
    longjmp(recovery_context, 1);

  }

  /*****************************************************************************/
  void Token_Stream::Semantics_Error(const std::string &s, const std::string &v)
    /*****************************************************************************/
  {
    if (!recovery_flag){
      error_count++;
      (*PAMGEN_NEVADA::parse_error_count) ++;
      output << "\n***** ERROR Parsing: " << s << std::endl;
      output << v << std::endl;
      output << std::flush;
    }
  }

}// end namespace PAMGEN_NEVADA
