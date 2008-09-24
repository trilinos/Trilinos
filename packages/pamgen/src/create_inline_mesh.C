#include "create_inline_mesh.h"
#include "inline_mesh_desc.h"
#include <iostream>
#include <strings.h>
#include <cstring>

#include "../mesh_spec_lt/mesh_specification.h" 
#include "inline_mesh_desc.h"
ms_lt::Mesh_Specification * buildMeshSpecification_LT(PAMGEN_NEVADA::Inline_Mesh_Desc* imd,int rank, int num_procs);

/*****************************************************************************/
int Delete_Pamgen_Mesh()
/*****************************************************************************/
{
  if(PAMGEN_NEVADA::Inline_Mesh_Desc::static_storage){
    delete PAMGEN_NEVADA::Inline_Mesh_Desc::static_storage;
    PAMGEN_NEVADA::Inline_Mesh_Desc::static_storage = NULL;
  }
  if(ms_lt::Mesh_Specification::static_storage){
    delete ms_lt::Mesh_Specification::static_storage;
    ms_lt::Mesh_Specification::static_storage = NULL;
  }
  return 0;
}

/*****************************************************************************/
int Create_Pamgen_Mesh(const char * file_char_array, 
		       int dimension,
		       int rank,
		       int num_procs)
/*****************************************************************************/
{
  PAMGEN_NEVADA::Inline_Mesh_Desc * imd = NULL;
  std::string fn("PAMGEN LIBRARY");

  // copy input into stream, no file operations in library.
  std::stringstream input_stream;
  int sfca = strlen(file_char_array);
  input_stream.write(file_char_array,sfca);


  int pec = 0;

  imd = PAMGEN_NEVADA::Parse_Inline_Mesh(fn,
			  input_stream,
			  pec,
			  dimension);

  if(pec > 0)return ERROR_PARSING_DEFINITION;
  
  if(!imd)return ERROR_CREATING_IMD;


  ms_lt::Mesh_Specification * ams = buildMeshSpecification_LT(imd,
							      rank, 
							      num_procs);
  if(!ams)return ERROR_CREATING_MS;

  return ERROR_FREE_CREATION;
}


/*****************************************************************************/
char * getPamgenEchoStream(char * car)
/*****************************************************************************/
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::echo_stream.str();
  const char * cst = st.c_str();
  strcpy(car,cst);
  return car;
}

/*****************************************************************************/
int getPamgenEchoStreamSize()
/*****************************************************************************/
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::echo_stream.str();
  const char * cst = st.c_str();
  int stsz = strlen(cst);
  return stsz;
}

/*****************************************************************************/
int getPamgenErrorStreamSize()
/*****************************************************************************/
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::static_storage->getErrorString();
  const char * cst = st.c_str();
  int stsz = strlen(cst);
  return stsz;
}

/*****************************************************************************/
int getPamgenWarningStreamSize()
/*****************************************************************************/
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::static_storage->getWarningString();
  const char * cst = st.c_str();
  int stsz = strlen(cst);
  return stsz;
}


/*****************************************************************************/
int getPamgenInfoStreamSize()
/*****************************************************************************/
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::static_storage->getInfoString();
  const char * cst = st.c_str();
  int stsz = strlen(cst);
  return stsz;
}

/*****************************************************************************/
char * getPamgenErrorStream(char * car)
/*****************************************************************************/
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::static_storage->getErrorString();
  const char * cst = st.c_str();
  strcpy(car,cst);
  return car;
}

/*****************************************************************************/
char * getPamgenWarningStream(char * car)
/*****************************************************************************/
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::static_storage->getWarningString();
  const char * cst = st.c_str();
  strcpy(car,cst);
  return car;
}

/*****************************************************************************/
char * getPamgenInfoStream(char * car)
/*****************************************************************************/
{
  std::string st = PAMGEN_NEVADA::Inline_Mesh_Desc::static_storage->getInfoString();
  const char * cst = st.c_str();
  strcpy(car,cst);
  return car;
}
