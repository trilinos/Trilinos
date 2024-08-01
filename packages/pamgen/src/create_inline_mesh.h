// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef create_inline_meshH
#define create_inline_meshH

#ifdef __cplusplus
extern "C"
{
#endif

  long long Delete_Pamgen_Mesh();

  long long Create_Pamgen_Mesh(const char * file_char_array, 
			       long long dimension,
			       long long rank,
			       long long num_procs,
			       long long int_max);

  char * getPamgenEchoStream(char *);
  long long getPamgenEchoStreamSize();

  char * getPamgenErrorStream(char *);
  long long getPamgenErrorStreamSize();

  char * getPamgenWarningStream(char *);
  long long getPamgenWarningStreamSize();

  char * getPamgenInfoStream(char *);
  long long getPamgenInfoStreamSize();
    
#define ERROR_FREE_CREATION 0
#define ERROR_CREATING_IMD 1
#define ERROR_CREATING_MS 2
#define ERROR_PARSING_DEFINITION 3


#ifdef __cplusplus
}
#endif
#endif
