/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/* Function prototypes */

int LB_Get_Processor_Name(LB *, char *);
int LB_Build_Machine_Desc(LB *);
int LB_divide_machine(LB *, int, MPI_Comm, int *, int *, int *, int *,
                      double *);

/* Misc. constants */
#define MAX_PROC_NAME_LEN (MPI_MAX_PROCESSOR_NAME+1)
#define MACHINE_DESC_FILE_DEFAULT "/etc/local/Zoltan_Machine_Desc"

