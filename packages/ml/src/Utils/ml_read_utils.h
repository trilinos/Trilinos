
#ifndef _MLREADER_
#define _MLREADER_

#include <stdio.h>

struct reader_context {
   int id;
   int N_levels, nsmooth, maxcoarsesize, coarse_its;
   int N_dofPerNode;
   double agg_thresh;
   char *smoother, *agg_coarsen_scheme, *coarse_solve, *krylov;
   char *partition_file;
   int  output;
   double tol, agg_damping;
   char *agg_spectral_norm;
};

#define MAX_INPUT_STR_LN 101
#define MAX_TOKENS 50


#ifdef __cplusplus
extern "C"
{
#endif

extern void ML_Reader_GetGeneralSpecs(FILE *ifp, 
                                      struct reader_context *context);
extern int  ML_Reader_LookFor(FILE *ifp, char *string, char input[], 
                              int ch_term);
extern int  ML_Reader_ReadString(FILE *ifp, char string[], char ch);
extern void ML_Reader_Strip(char string[]);
extern void ML_Reader_GetSolutionSpecs(FILE *ifp, 
                            struct reader_context *context);
extern void ML_Reader_GetAggregationSpecs(FILE *ifp, 
                            struct reader_context *context);
extern void ML_Reader_InitContext(struct reader_context *context);
extern void ML_Reader_ReadInput(char *cmd_file_name, 
                                 struct reader_context **context);

extern int ML_strcmp(char *input, char *string);


#ifdef __cplusplus
}
#endif


#endif
