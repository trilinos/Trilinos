/*****************************************************************************/
/*                                                                           *
 *  This file contains routines for reading input to drive ML from a file.   *
 *  These routines work with the example 'ml_readfile.c'.                    * 
 *  The current input file has the following form:                           *
 *                                                                           *
#
#  Test input file ML 
#       
# Notes: 
#   1) Captilization should not matter
#   2) Lines starting with # are comments
#   3) Parallel Partitioning File is not used
#      when only 1 processor is present.
#   4) comments in [] indicate available options.
#   5) The matrix must be stored in a file according
#      to Aztec's AZ_read_msr() format.
#   6) Including the file 'rhsfile' will cause this 
#      data (stored in Aztec's AZ_read_msr() format)
#      to be used as righthand side.
#   7) Including the file 'initguessfile' will cause this 
#      data (stored in Aztec's AZ_read_msr() format)
#      to be used as righthand side.
#   8) rigid body mode information can be input by
#      keeping files 'rigid_body_mode%d' (stored
#      in Aztec's AZ_read_msr() format) where %d
#      starts from 0 and increases. Each file
#      should contain 1 rigid body mode.
#   9) The utility ml/util/az_capt2read.c can be used
#      to convert matlab type '*.dat' data into 
#      AZ_read_msr() format.
#
-----------------------------------------------
      General Problem Specifications
-----------------------------------------------
Number of DOF per node       = 1
Parallel Partitioning File   = myfile  
Output Frequency             = 2       
Tolerance                    = 1.0e-11
Print Level                  = 1
#                              [0,1,...]

-----------------------------------------------
      Solution Specifications
-----------------------------------------------
Max Number of Levels         = 4
Type of Smoother             = SymGaussSeidel 
#                              [Parasails, GaussSeidel, SymGaussSeidel, Poly,
#                               BlockGaussSeidel, Aggregate, Jacobi, Metis]
Smoother steps per level     = 7
Coarse grid solver           = SuperLU
#                              [Parasails, GaussSeidel, SymGaussSeidel, Poly,
#                               BlockGaussSeidel, Aggregate, Jacobi, Metis,
#                               SuperLU]
Coarse Grid iterations       = 1
Outer Iteration              = Cg
#                              [Cg, Bicgstab, Tfqmr, Gmres] 

-----------------------------------------------
      Aggregation Specifications
-----------------------------------------------
Type of Aggregation          = Mis
#                              [Mis, Uncoupled, Coupled]
Aggregate threshold          = 0.0
Max coarse size              = 30
Smoothed aggregation damping = 1.5 
Spectral norm calculation    = Anorm
#                              [Anorm, Calc]
# end of inputfile

 *                                                                           *
 *  In the example ml_readfile.c, this information should be put in a        *
 *  a file entitled 'ml_inputfile'.                                          *
 *****************************************************************************/
/*****************************************************************************/
#include "ml_read_utils.h"
#include "ml_struct.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void ML_Reader_ReadInput(char *cmd_file_name, struct reader_context **context)

/* Read the input file for multigrid driver according to the tumi/shadi
 * format (not written down anywhere) and set algorithm and problem dependent
 * parameters accordingly.
 *
 *    Authors:  John N. Shadid, 1421
 *              Harry K. Moffat 1421
 *              Scott A. Hutchinson 1421
 *              Gary L. Hennigan 1421
 */

{
  FILE *ifp;


  /*********************** BEGIN EXECUTION ***********************************/
  
  *context =(struct reader_context *) ML_allocate(sizeof(struct reader_context));
  ML_Reader_InitContext(*context);

  /* Open the command input file */

  if ( (ifp = fopen(cmd_file_name, "r")) != NULL) {
      ;
/*    (void) fprintf(stderr, "ML input mode\n"); */
  } else {
    (void) fprintf(stderr, "read_input_file: Can't open input file, %s,",
                   cmd_file_name);
    (void) fprintf(stderr, " for reading\n");
    exit(-1);
  }

  ML_Reader_GetGeneralSpecs(ifp, *context);
  ML_Reader_GetSolutionSpecs(ifp, *context);
  ML_Reader_GetAggregationSpecs(ifp, *context);


}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void ML_Reader_GetGeneralSpecs(FILE *ifp, struct reader_context *context)

{
/* Read the General Problem specification section of the input file. */

  char        input[MAX_INPUT_STR_LN], *c_srch;
  static char yo[] = "get_general_specs";
  int         output_level;

  c_srch = "general problem specifications";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '=')) {
    c_srch = "general specifications";
    if (!ML_Reader_LookFor(ifp, c_srch, input, '=')) {
      fprintf(stderr, "%s: ERROR, couldn't find \"%s\" or\n \"%s\"!\n", yo,
                       "General Problem Specifications", c_srch);
      exit(-1);
    }
  }

  /* Read the number of degrees of freedom per node information */

  c_srch = "number of dof per node";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '='))
    context->N_dofPerNode = 1;       /* Defaults to 1 */
  else {
    ML_Reader_ReadString(ifp, input, '\n');
    if (sscanf(input, "%d", &(context->N_dofPerNode)) != 1) {
      fprintf(stderr, "%s ERROR: can\'t interp int while looking for \"%s\"\n",
              yo, c_srch);
      exit(-1);
    }
  }

  /* Select a file containing partioning information */

  c_srch = "parallel partitioning file";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '=')) {
    fprintf(stderr, "%s: ERROR, couldn't find \"%s\"!\n", yo, c_srch);
    exit(-1);
  }
  ML_Reader_ReadString(ifp,input, '\n');
  ML_Reader_Strip(input);
  strcpy(context->partition_file,(const char *) input);

  /* How often should the residual be printed */

  c_srch = "output frequency";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '='))
    context->output = 1;       /* Defaults to 1 */
  else {
    ML_Reader_ReadString(ifp, input, '\n');
    if (sscanf(input, "%d", &(context->output)) != 1) {
      fprintf(stderr, "%s ERROR: can\'t interp int while looking for \"%s\"\n",
              yo, c_srch);
      exit(-1);
    }
  }

  /* Solve to what tolerance */

  c_srch = "tolerance";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '='))
    context->tol = 1.0e-6;       /* Default */
  else {
    ML_Reader_ReadString(ifp, input, '\n');
    if (sscanf(input, "%lf", &(context->tol)) != 1) {
      fprintf(stderr, "%s ERROR: can\'t interp double while looking for \"%s\"\n",
              yo, c_srch);
      exit(-1);
    }
  }

  /* Determine the amount of information that ML should print out. */

  c_srch = "print level";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '='))
    output_level = 0;       /* Defaults to 0 */
  else {
    ML_Reader_ReadString(ifp, input, '\n');
    if (sscanf(input, "%d", &output_level) != 1) {
      fprintf(stderr, "%s ERROR: can\'t interp int while looking for \"%s\"\n",
              yo, c_srch);
      exit(-1);
    }
  }
  context->output_level = output_level;


}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ML_Reader_LookFor(FILE *ifp, char *string, char input[], int ch_term)

/* Scan the input file (reading in strings according to 'ML_Reader_ReadString(
 * ifp)' specifications) until the character pattern in 'string' is matched.
 *
 * Author:                 Ray S. Tuminaro Div 1422
 * revised:                10/2/90 John N. Shadid
 *
 * Parameter list:
 *
 * ifp    == pointer to file "input"
 * string == contains string pattern to be matched.
 * input  == buffer array to hold characters that are read ing.
 * ch_term  == termination character. When scanning a line of input
 *             is read until either a newline, the 'ch' termination
 *             character is read, or the end-of-file is read.
 */

{

  long  file_pos;

  /*-----------------------------Execution Begins----------------------------*/

  /* Store the current position in the file */
  file_pos = ftell(ifp);

  /* Begin reading the file */
  if (ML_Reader_ReadString(ifp, input, ch_term) == -1) {

    /* Upon failure, reset the file position and return */
    fseek(ifp, file_pos, SEEK_SET);

    return 0;
  }

  ML_Reader_Strip(input);

  while (ML_strcmp(input, string) != 0 ) {
    if (ML_Reader_ReadString(ifp, input, ch_term) == -1) {

      /* Upon failure, reset the file position and return */
      fseek(ifp, file_pos, SEEK_SET);

      return 0;
    }
    ML_Reader_Strip(input);
  }

  return 1;

} /* ML_Reader_LookFor */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ML_Reader_ReadString(FILE *ifp, char string[], char ch)

/* This routine reads the standard input until encountering the end-of-file, a
 * newline, the character 'ch' or until MAX_INPUT_STR_LN characters are read.
 * The inputted characters are read into 'string'.  If an error occurs, -1 is
 * returned and an error message is written to standard error.  Upon successful
 * completion, the number of characters read plus 1 for the null character at
 * the end of the string is returned.
 *
 *      Author:                 Ray S. Tuminaro Div 1422
 *
 *      Parameter list:
 *
 *      ifp    == pointer to file "input"
 *      string == On output 'string' contains the characters read
 *                from the input stream.
 *      ch     == Termination character. That is, input function
 *                stops when 'ch' is read.
 */

{
  int i = 0;
  int new_ch = 0;

  while ( (i < MAX_INPUT_STR_LN) && ((new_ch = getc(ifp)) != ch)
          && (new_ch != '\n')
          && (new_ch != EOF) ) {
    if ((new_ch >= 'A') && (new_ch <= 'Z'))
       new_ch = 'a' + new_ch - 'A'; 
    string[i++] = new_ch;
  }

  if (new_ch == EOF) return -1;

  if (i == MAX_INPUT_STR_LN) {
    (void) fprintf (stderr, "ML_Reader_ReadString ERROR: scanned %d characters and could not find (%c)\n", MAX_INPUT_STR_LN, ch);
    return -1;
  }

  string[i] = '\0';

  return (i+1);

} /* ML_Reader_ReadString */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void ML_Reader_Strip(char string[])

/* This routine strips off blanks and tabs (only leading and trailing
   characters) in 'string'.

   Author:                 Ray S. Tuminaro Div 1422

   Parameter list:

   string == On output 'string' contains the same characters as on input except
             the leading and trailing blanks and tabs have been removed.
*/

{
  int  i, j;
  char ch;

  /* find first real character */

  i = 0;
  while ( ((ch = string[i]) != '\0' ) && ((ch == ' ') || (ch == '\t')) )
    i++;

  /* move real part of string to the front */

  j = 0;
  while ( (ch = string[j+i]) != '\0' ) {
    string[j] = string[j+i];
    j++;
  }
  string[j] = string[j+i];
  j--;

  /* remove trailing blanks */

  while ( (j != -1 ) && (((ch = string[j]) == ' ')||(ch == '\t') ||
                         (ch == '\n')) )
    j--;
  string[j+1] = '\0';

} /* ML_Reader_Strip */


void ML_Reader_GetSolutionSpecs(FILE *ifp, struct reader_context *context)

{
/* Read the Solution specification section of the input file. */

  char *yo = "get_solution_specs";

  char  input[MAX_INPUT_STR_LN], *c_srch;

  /* Read in Solution Specifications */

  c_srch = "solution specifications";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '=')) {
    fprintf(stderr, "%s: ERROR, couldn't find \"%s\"!\n", yo, c_srch);
    exit(-1);
  }

  /* select the maximum number of MG levels */

  c_srch = "max number of levels";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '='))
    context->N_levels = 1;       /* Defaults to 1 */
  else {
    ML_Reader_ReadString(ifp, input, '\n');
    if (sscanf(input, "%d", &(context->N_levels)) != 1) {
      fprintf(stderr, "%s ERROR: can\'t interp int while looking for \"%s\"\n",
              yo, c_srch);
      exit(-1);
    }
  }

  /* select the smoother on all levels */

  c_srch = "type of smoother";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '=')) {
    fprintf(stderr, "%s: ERROR, couldn't find \"%s\"!\n", yo, c_srch);
    exit(-1);
  }
  ML_Reader_ReadString(ifp,input, '\n');
  ML_Reader_Strip(input);
  strcpy(context->smoother,input);

  /* select the subsmoother (if applicable) */

  c_srch = "type of subsmoother";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '=')) {
    strcpy(context->subsmoother,"default");
  }
  else {
     ML_Reader_ReadString(ifp,input, '\n');
     ML_Reader_Strip(input);
     strcpy(context->subsmoother,input);
  }

  /* select the number of smoother iterations per level */

  c_srch = "smoother steps per level";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '='))
    context->nsmooth = 1;       /* Defaults to 1 */
  else {
    ML_Reader_ReadString(ifp, input, '\n');
    if (sscanf(input, "%d", &(context->nsmooth)) != 1) {
      fprintf(stderr, "%s ERROR: can\'t interp int while looking for \"%s\"\n",
              yo, c_srch);
      exit(-1);
    }
  }

  /* select the smoother on all levels */

  c_srch = "coarse grid solver";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '=')) {
    fprintf(stderr, "%s: ERROR, couldn't find \"%s\"!\n", yo, c_srch);
    exit(-1);
  }
  ML_Reader_ReadString(ifp,input, '\n');
  ML_Reader_Strip(input);
  strcpy(context->coarse_solve,input);

  /* select the number of coarse grid iterations */

  c_srch = "coarse grid iterations";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '='))
    context->coarse_its = 1;       /* Defaults to 1 */
  else {
    ML_Reader_ReadString(ifp, input, '\n');
    if (sscanf(input, "%d", &(context->coarse_its)) != 1) {
      fprintf(stderr, "%s ERROR: can\'t interp int while looking for \"%s\"\n",
              yo, c_srch);
      exit(-1);
    }
  }

  /* select the outer iterative method */

  c_srch = "outer iteration";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '=')) {
    fprintf(stderr, "%s: ERROR, couldn't find \"%s\"!\n", yo, c_srch);
    exit(-1);
  }
  ML_Reader_ReadString(ifp,input, '\n');
  ML_Reader_Strip(input);
  strcpy(context->krylov,input);

  /* select the maximum number of outer iterations */

  c_srch = "max number of outer iterations";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '='))
    context->max_outer_its = 500;       /* Defaults to 500 */
  else {
    ML_Reader_ReadString(ifp, input, '\n');
    if (sscanf(input, "%d", &(context->max_outer_its)) != 1) {
      fprintf(stderr, "%s ERROR: can\'t interp int while looking for \"%s\"\n",
              yo, c_srch);
      exit(-1);
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void ML_Reader_GetAggregationSpecs(FILE *ifp, struct reader_context *context)

{
/* Read the Aggregation specification section of the input file. */

  char *yo = "get_aggregation_specs";

  char  input[MAX_INPUT_STR_LN], *c_srch;

  /* Read in Aggregation Specifications */

  c_srch = "aggregation specifications";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '=')) {
    fprintf(stderr, "%s: ERROR, couldn't find \"%s\"!\n", yo, c_srch);
    exit(-1);
  }

  /* select aggregation method */

  c_srch = "type of aggregation";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '=')) {
    fprintf(stderr, "%s: ERROR, couldn't find \"%s\"!\n", yo, c_srch);
    exit(-1);
  }
  ML_Reader_ReadString(ifp,input, '\n');
  ML_Reader_Strip(input);
  strcpy(context->agg_coarsen_scheme,input);

  /* select aggregation drop tolerance */

  c_srch = "aggregate threshold";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '='))
    context->agg_thresh = 0.0;       /* Defaults to 0 */
  else {
    ML_Reader_ReadString(ifp, input, '\n');
    if (sscanf(input, "%lf", &(context->agg_thresh)) != 1) {
      fprintf(stderr,"%s ERROR: can\'t interp double while looking for \"%s\"\n",
              yo, c_srch);
      exit(-1);
    }
  }

  /* select max size of coarsest grid */

  c_srch = "max coarse size";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '='))
    context->maxcoarsesize = 100;       /* Defaults to 100 */
  else {
    ML_Reader_ReadString(ifp, input, '\n');
    if (sscanf(input, "%d", &(context->maxcoarsesize)) != 1) {
      fprintf(stderr,"%s ERROR: can\'t interp int while looking for \"%s\"\n",
              yo, c_srch);
      exit(-1);
    }
  }


  /* select aggregation damping */

  c_srch = "smoothed aggregation damping";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '='))
    context->agg_damping = 1.333333333333;       /* Default */
  else {
    ML_Reader_ReadString(ifp, input, '\n');
    if (sscanf(input, "%lf", &(context->agg_damping)) != 1) {
      fprintf(stderr,"%s ERROR: can\'t interp double while looking for \"%s\"\n",
              yo, c_srch);
      exit(-1);
    }
  }

  /* calculate spectral norm with CG or just use Anorm */

  c_srch = "spectral norm calculation";
  if (!ML_Reader_LookFor(ifp, c_srch, input, '=')) {
    fprintf(stderr, "%s: ERROR, couldn't find \"%s\"!\n", yo, c_srch);
    exit(-1);
  }
  ML_Reader_ReadString(ifp,input, '\n');
  ML_Reader_Strip(input);
  strcpy(context->agg_spectral_norm,input);

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void ML_Reader_InitContext(struct reader_context *context)
{
/* Initialize the context used for keeping track of multigrid input options */
   
   context->id            = 0;
   context->N_levels      = 0;
   context->nsmooth       = 0;
   context->maxcoarsesize = 0;
   context->coarse_its    = 0;
   context->N_dofPerNode  = 0;
   context->agg_thresh    = 0.0;
   context->smoother[0]      = '\0';
   context->agg_coarsen_scheme[0] = '\0';
   context->coarse_solve[0]  = '\0';
   context->krylov[0]        = '\0';
   context->partition_file[0] = '\0';
   context->output         = 1;
   context->tol            = 1.e-6;
   context->agg_damping    = 1.33333;
   context->agg_spectral_norm[0] = '\0';
   context->output_level = 0;


}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int ML_strcmp(char *input, char *string)
{
/* Similar to 'C' strcmp except this one converts everything to lower case. */
  int i;
  char *input_copy, *string_copy;

  input_copy = (char *) ML_allocate(sizeof(char)*(strlen(input)+1));
  string_copy = (char *) ML_allocate(sizeof(char)*(strlen(string)+1));
  strcpy(input_copy,input);
  strcpy(string_copy,string);
  i = 0;
  while ((input_copy[i] != '\0') && (input_copy[i] != '\n')) {
    if ((input_copy[i] >= 'A') && (input_copy[i] <= 'Z'))
       input_copy[i] = 'a' + input_copy[i] - 'A'; 
    i++;
  }
  i = 0;
  while ((string_copy[i] != '\0') && (string_copy[i] != '\n')) {
    if ((string_copy[i] >= 'A') && (string_copy[i] <= 'Z'))
       string_copy[i] = 'a' + string_copy[i] - 'A'; 
    i++;
  }

  i = strcmp(input_copy, string_copy);
  ML_free(input_copy);
  ML_free(string_copy);

  return(i);
}
