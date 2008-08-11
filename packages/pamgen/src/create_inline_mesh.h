#ifndef create_inline_meshH
#define create_inline_meshH

#ifdef __cplusplus
extern "C"
{
#endif

  int Delete_Pamgen_Mesh();

  int Create_Pamgen_Mesh(const char * file_char_array, 
			 int dimension,
			 int rank,
			 int num_procs);

  char * getPamgenEchoStream(char *);
  int getPamgenEchoStreamSize();

  char * getPamgenErrorStream(char *);
  int getPamgenErrorStreamSize();

  char * getPamgenWarningStream(char *);
  int getPamgenWarningStreamSize();

  char * getPamgenInfoStream(char *);
  int getPamgenInfoStreamSize();
    
#define ERROR_FREE_CREATION 0
#define ERROR_CREATING_IMD 1
#define ERROR_CREATING_MS 2
#define ERROR_PARSING_DEFINITION 3


#ifdef __cplusplus
}
#endif
#endif
