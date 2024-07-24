// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "create_inline_mesh.h"
#include "inline_mesh_desc.h"
#include <iostream>
#include <strings.h>
#include <cstring>

#include "../mesh_spec_lt/pamgen_mesh_specification.h"
#include "inline_mesh_desc.h"
ms_lt::Mesh_Specification * buildMeshSpecification_LT(
    PAMGEN_NEVADA::Inline_Mesh_Desc* imd,
    long long rank,
    long long num_procs
    );
ms_lt::Mesh_Specification * consolidateMeshSpecification_LT(
    ms_lt::Mesh_Specification * bms
    );

/*****************************************************************************/
long long Delete_Pamgen_Mesh()
{
  if(PAMGEN_NEVADA::Inline_Mesh_Desc::first_im_static_storage){
    delete PAMGEN_NEVADA::Inline_Mesh_Desc::first_im_static_storage;
  }
  PAMGEN_NEVADA::Inline_Mesh_Desc::im_static_storage = NULL;
  PAMGEN_NEVADA::Inline_Mesh_Desc::first_im_static_storage = NULL;

  if(ms_lt::Mesh_Specification::first_ms_static_storage){
    delete ms_lt::Mesh_Specification::first_ms_static_storage;
  }
  ms_lt::Mesh_Specification::first_ms_static_storage = NULL;

  return 0;
}

/*****************************************************************************/
long long Create_Pamgen_Mesh(
    const char * file_char_array,
    long long dimension,
    long long rank,
    long long num_procs,
    long long max_int
    )
{
  PAMGEN_NEVADA::Inline_Mesh_Desc * imd = NULL;
  PAMGEN_NEVADA::Inline_Mesh_Desc * fimd = NULL;
  std::string fn("PAMGEN LIBRARY");

  PAMGEN_NEVADA::Partition::partition_count = 0;

  // copy input into stream, no file operations in library.
  std::stringstream input_stream;
  long long sfca = strlen(file_char_array);
  input_stream.write(file_char_array,sfca);

  long long pec = 0;

  fimd = PAMGEN_NEVADA::Parse_Inline_Mesh(fn,
      input_stream,
      pec,
      dimension,
      max_int);

  if(pec > 0)return ERROR_PARSING_DEFINITION;

  if(!fimd)return ERROR_CREATING_IMD;
  imd = fimd;
  ms_lt::Mesh_Specification * ams = NULL;
  while(imd){

    ams = buildMeshSpecification_LT(imd,
        rank,
        num_procs);

    if(!ams)return ERROR_CREATING_MS;

    ms_lt::Mesh_Specification::Add_MS(ams);

    imd = imd->next;
  }

  ms_lt::Mesh_Specification * nms =  ms_lt::Mesh_Specification::first_ms_static_storage->consolidateMS();

  ms_lt::Mesh_Specification::Replace_MS(nms);

  return ERROR_FREE_CREATION;
}


/*****************************************************************************/
char * getPamgenEchoStream(char * car)
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::echo_stream.str();
  const char * cst = st.c_str();
  strcpy(car,cst);
  return car;
}

/*****************************************************************************/
long long getPamgenEchoStreamSize()
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::echo_stream.str();
  const char * cst = st.c_str();
  long long stsz = strlen(cst);
  return stsz;
}

/*****************************************************************************/
long long getPamgenErrorStreamSize()
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::first_im_static_storage->getErrorString();
  const char * cst = st.c_str();
  long long stsz = strlen(cst);
  return stsz;
}

/*****************************************************************************/
long long getPamgenWarningStreamSize()
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::first_im_static_storage->getWarningString();
  const char * cst = st.c_str();
  long long stsz = strlen(cst);
  return stsz;
}


/*****************************************************************************/
long long getPamgenInfoStreamSize()
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::first_im_static_storage->getInfoString();
  const char * cst = st.c_str();
  long long stsz = strlen(cst);
  return stsz;
}

/*****************************************************************************/
char * getPamgenErrorStream(char * car)
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::first_im_static_storage->getErrorString();
  const char * cst = st.c_str();
  strcpy(car,cst);
  return car;
}

/*****************************************************************************/
char * getPamgenWarningStream(char * car)
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::first_im_static_storage->getWarningString();
  const char * cst = st.c_str();
  strcpy(car,cst);
  return car;
}

/*****************************************************************************/
char * getPamgenInfoStream(char * car)
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::first_im_static_storage->getInfoString();
  const char * cst = st.c_str();
  strcpy(car,cst);
  return car;
}
