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
 * $Name$
 *====================================================================*/
#ifndef lint
static char *cvs_dr_input_const = "$Id$";
#endif

#ifndef _DR_INPUT_CONST_H_
#define _DR_INPUT_CONST_H_

/* Define NEMESIS_IO to use Nemesis I/O.              */
/* Undefine NEMESIS_IO to avoid compiling in Nemesis. */
#define NEMESIS_IO

/* define the FEM input file types */
#define NEMESIS_FILE 0
#define CHACO_FILE   1

/* Structure used to store the information necessary for parallel I/O. */
struct Parallel_IO
{
  int     dsk_list_cnt;

  int    *dsk_list;
  int     rdisk;

  int     num_dsk_ctrlrs;            /* The number of disk controllers.     */
  int     pdsk_add_fact;             /* The offset from zero used by the    */
                                     /* the target machine.                 */

  int     zeros;        /* 1 - if the target machine uses leading zeros when */
                        /*     designating the disk number (eg - the paragon */
                        /*     uses /pfs/io_01)                              */
                        /* 0 - if it does not (eg - the tflop uses           */
                        /*     /pfs/tmp_1)                                   */

  int     file_type;    /* input file type */

  /* The root location of the parallel disks */
  char    pdsk_root[FILENAME_MAX+1];

  /* The subdirectory to write files to */
  char    pdsk_subdir[FILENAME_MAX+1];

  /* The base name of the parallel file. */
  char    pexo_fname[FILENAME_MAX+1];

};
typedef struct Parallel_IO  PARIO_INFO;
typedef struct Parallel_IO *PARIO_INFO_PTR;


/* Function prototypes */
extern
int read_cmd_file(
  char *filename,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info
);

int check_inp(
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info
);

void brdcst_cmd_info(
  int Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info
);

void gen_par_filename(
  char *scalar_fname,
  char *par_fname,
  PARIO_INFO_PTR pio_info,
  int proc_for,
  int nprocs
);


#endif /* _DR_INPUT_CONST_H_ */
