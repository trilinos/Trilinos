
#ifndef _MLREADER_
#define _MLREADER_

#include <stdio.h>
#include "ml_common.h"
#include "ml_memory.h"

struct reader_context {
   int id;
   int N_levels, nsmooth, maxcoarsesize, coarse_its;
   int max_outer_its;
   int N_dofPerNode;
   double agg_thresh;
   char smoother[80], agg_coarsen_scheme[80], coarse_solve[80], krylov[80];
   char partition_file[80];
   int  output, output_level;
   double tol, agg_damping;
   char agg_spectral_norm[80];
   char subsmoother[80];
   char cycle[80];
};

#define MAX_INPUT_STR_LN 101
#define MAX_TOKENS 50


#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{
#endif
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


#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif


#endif
