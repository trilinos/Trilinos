/*
 * Copyright(C) 1999-2021 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* This file contains the source code for the program used to test the
 * Nemesis distribution.
 *****************************************************************************
 * Written By: Gary L. Hennigan (SNL, 1421)
 *****************************************************************************
 * Functions contained in this file:
 *      main() - Entry point and main calling program.
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <exodusII.h>

#include "ne_nemesisI.h"

/* Constants for init_global functions */
#define NNG   100
#define NEG   50
#define NEBG  5
#define NEBCG 10
#define NNSG  2
#define NSSG  3

/* Constants for load-balance functions */
#define NPROC  10
#define NPROCF NPROC
#define NINTN  200
#define NBORN  10
#define NEXTN  5
#define NINTE  100
#define NBORE  20
#define NNCMAP 4
#define NECMAP 2

/* Constants for communication map */
#define NCNTCM 20
#define ECNTCM 17

int main(int argc, char *argv[])
{

  /* Local function calls */
  int ne_test_glbp(int);
  int ne_test_piinf(int);
  int ne_test_pinig(int);
  int ne_test_pelbid(int);
  int ne_test_pnsp(int);
  int ne_test_pssp(int);
  int ne_test_pnm(int);
  int ne_test_pem(int);
  int ne_test_pcmp(int);
  int ne_test_pncm(int);
  int ne_test_pecm(int);

  int ne_test_giinf(int);
  int ne_test_ginig(int);
  int ne_test_gelbid(int);
  int ne_test_gnsp(int);
  int ne_test_gssp(int);
  int ne_test_gnm(int);
  int ne_test_gem(int);
  int ne_test_gncm(int);
  int ne_test_gecm(int);

  int ne_test_plbpc(int);
  int ne_test_pcmpc(int);

  /* Uninitialized local variables */
  int   ne_file_id;
  char *file_name = "./ne_test.exoII";
  float version;

  /* Initialized local variables */
  int mode3 = EX_CLOBBER;
  int mode4 = EX_CLOBBER | EX_NETCDF4 | EX_NOCLASSIC;

  char *yo    = "main";
  int   io_ws = 0, cpu_ws = 0, t_pass = 0, t_fail = 0;

  /*---------------------------------------------------------------------------*/
  /*                      OUTPUT TEST SECTION                                  */
  /*---------------------------------------------------------------------------*/

  printf("*********************Output Tests***********************\n");

  /* Create the ExodusII/Nemesis file */
  printf("creating ExodusII file...");
  fflush(stdout);

  /* Attempt to create a netcdf4-format file; if it fails, then assume
     that the netcdf library does not support that mode and fall back
     to classic netcdf3 format.  If that fails, issue an error and
     die.
  */
  if ((ne_file_id = ex_create(file_name, mode4, &cpu_ws, &io_ws)) < 0) {
    /* netcdf4 create failed, try netcdf3 */
    if ((ne_file_id = ex_create(file_name, mode3, &cpu_ws, &io_ws)) < 0) {
      printf("FAILED\n");
      t_fail++;
      fprintf(stderr, "[%s]: ERROR, unable to create test file \"%s\"!\n", yo, file_name);
      exit(-1);
    }
    else {
      printf(" (netcdf3 format) ");
    }
  }
  else {
    printf(" (netcdf4 format) ");
  }
  printf("successful\n");
  fflush(stdout);
  t_pass++;

  /* Test the output of initial information */
  printf("testing init info output...");
  fflush(stdout);
  if (ne_test_piinf(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test the output of initial global information */
  printf("testing global init info output...");
  fflush(stdout);
  if (ne_test_pinig(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test the output of the global element block IDs */
  printf("testing global element block ID output...");
  fflush(stdout);
  if (ne_test_pelbid(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test the output of the global node-set info */
  printf("testing global node-set params output...");
  fflush(stdout);
  if (ne_test_pnsp(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test the output of the global side-set info */
  printf("testing global side-set params output...");
  fflush(stdout);
  if (ne_test_pssp(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test the output of the concatenated load-balance parameters */
  printf("testing concatenated load balance info output...");
  fflush(stdout);
  if (ne_test_plbpc(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
  }
  else {
    printf("successful\n");
    fflush(stdout);
  }

  /* Test the output of the node map */
  printf("testing node map output...");
  fflush(stdout);
  if (ne_test_pnm(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test the output of the element map */
  printf("testing element map output...");
  fflush(stdout);
  if (ne_test_pem(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test the output of the concatenated communication map params */
  printf("testing concatenated communication map params output...");
  fflush(stdout);
  if (ne_test_pcmpc(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test nodal communication map output */
  printf("testing nodal communication map output...");
  fflush(stdout);
  if (ne_test_pncm(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test elemental communication map output */
  printf("testing elemental communication map output...");
  fflush(stdout);
  if (ne_test_pecm(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Close the ExodusII/Nemesis test file */
  printf("closing ExodusII file...");
  fflush(stdout);
  if (ex_close(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
    fprintf(stderr, "[%s]: ERROR, unable to close test file \"%s\"!\n", yo, file_name);
    exit(-1);
  }
  printf("successful\n");
  fflush(stdout);
  t_pass++;

  /*---------------------------------------------------------------------------*/
  /*                       INPUT TEST SECTION                                  */
  /*---------------------------------------------------------------------------*/

  printf("**********************Input Tests***********************\n");

  /* Re-open the ExodusII/NemesisI file */
  printf("reopening ExodusII file...");
  fflush(stdout);
  if (ex_open(file_name, EX_READ, &cpu_ws, &io_ws, &version) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test read of of the initial information */
  printf("testing init info input...");
  fflush(stdout);
  if (ne_test_giinf(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test read of initial global information */
  printf("testing global init info input...");
  fflush(stdout);
  if (ne_test_ginig(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test read of global element block IDs */
  printf("testing global element block IDs input...");
  fflush(stdout);
  if (ne_test_gelbid(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test read of global node-set params */
  printf("testing global node-set params input...");
  fflush(stdout);
  if (ne_test_gnsp(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test read of global side-set params */
  printf("testing global side-set params input...");
  fflush(stdout);
  if (ne_test_gssp(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test read of load-balance params */
  printf("testing load-balance params input...");
  fflush(stdout);
  if (ne_test_glbp(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test read of the node map */
  printf("testing node map input...");
  fflush(stdout);
  if (ne_test_gnm(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test read of element map */
  printf("testing element map input...");
  fflush(stdout);
  if (ne_test_gem(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test read of nodal communication maps */
  printf("testing nodal communication map input...");
  fflush(stdout);
  if (ne_test_gncm(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Test read of elemental communication maps */
  printf("testing elemental communication map input...");
  fflush(stdout);
  if (ne_test_gecm(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
  }
  else {
    printf("successful\n");
    fflush(stdout);
    t_pass++;
  }

  /* Close the ExodusII/Nemesis test file */
  printf("closing ExodusII file...");
  fflush(stdout);
  if (ex_close(ne_file_id) < 0) {
    printf("FAILED\n");
    fflush(stdout);
    t_fail++;
    fprintf(stderr, "[%s]: ERROR, unable to close test file \"%s\"!\n", yo, file_name);
    exit(-1);
  }
  printf("successful\n");
  fflush(stdout);
  t_pass++;

  /* Output a test summary */
  printf("\n");
  printf("Tests Passed: %d\n", t_pass);
  printf("Tests Failed: %d\n", t_fail);

  return 0;
}

/*****************************************************************************/
int ne_test_piinf(int fileid)
{
  char ftype[3];

  ex_copy_string(ftype, "s", 3);

  return (ne_put_init_info(fileid, NPROC, NPROCF, ftype));
}

/*****************************************************************************/
int ne_test_pinig(int fileid)
{

  int nng = NNG, neg = NEG, nebg = NEBG, nnsg = NNSG, nssg = NSSG;

  return (ne_put_init_global(fileid, nng, neg, nebg, nnsg, nssg));
}

/*****************************************************************************/
int ne_test_pelbid(int fileid)
{
  int elblk_ids[NEBG], elblk_cnt[NEBG];

  for (int i = 0; i < NEBG; i++) {
    elblk_ids[i] = (i + 1);
    elblk_cnt[i] = NEBCG;
  }

  return (ne_put_eb_info_global(fileid, elblk_ids, elblk_cnt));
}

/*****************************************************************************/
int ne_test_pnsp(int fileid)
{
  int global_ids[NNSG], global_n_cnts[NNSG], global_df_cnts[NNSG];

  for (int i = 0; i < NNSG; i++) {
    global_ids[i]     = 2 * (i + 1);
    global_n_cnts[i]  = 3 * (i + 1);
    global_df_cnts[i] = 1;
  }

  return (ne_put_ns_param_global(fileid, global_ids, global_n_cnts, global_df_cnts));
}

/*****************************************************************************/
int ne_test_pssp(int fileid)
{
  int global_ids[NSSG], global_el_cnts[NSSG], global_df_cnts[NSSG];

  for (int i = 0; i < NSSG; i++) {
    global_ids[i]     = 3 * (i + 1);
    global_el_cnts[i] = 2 * (i + 1);
    global_df_cnts[i] = 1;
  }

  return (ne_put_ss_param_global(fileid, global_ids, global_el_cnts, global_df_cnts));
}

/*****************************************************************************/
int ne_test_pnm(int fileid)
{
  int j1 = 0;
  int node_mapi[NINTN], node_mapb[NBORN], node_mape[NEXTN];

  for (int iproc = 0; iproc < NPROCF; iproc++) {
    for (int j = 0; j < NINTN; node_mapi[j++] = j1++) {
      ;
    }
    for (int j = 0; j < NBORN; node_mapb[j++] = j1++) {
      ;
    }
    for (int j = 0; j < NEXTN; node_mape[j++] = j1++) {
      ;
    }
    j1        = 0;
    int error = ne_put_node_map(fileid, node_mapi, node_mapb, node_mape, iproc);
    if (error < 0) {
      return error;
    }
  }

  return 0;
}

/*****************************************************************************/
int ne_test_pem(int fileid)
{
  int j1 = 0;
  int elem_mapi[NINTE], elem_mapb[NBORE];

  for (int iproc = 0; iproc < NPROCF; iproc++) {
    for (int j = 0; j < NINTE; elem_mapi[j++] = j1++) {
      ;
    }
    for (int j = 0; j < NBORE; elem_mapb[j++] = j1++) {
      ;
    }
    j1        = 0;
    int error = ne_put_elem_map(fileid, elem_mapi, elem_mapb, iproc);
    if (error < 0) {
      return error;
    }
  }

  return 0;
}

/*****************************************************************************/
int ne_test_pcmp(int fileid)
{
  int node_map_ids[NNCMAP], node_map_node_cnts[NNCMAP];
  int elem_map_ids[NECMAP], elem_map_elem_cnts[NECMAP];

  for (int iproc = 0; iproc < NPROCF; iproc++) {
    for (int i = 0; i < NNCMAP; i++) {
      node_map_ids[i]       = (i + 1);
      node_map_node_cnts[i] = NCNTCM;
    }

    for (int i = 0; i < NECMAP; i++) {
      elem_map_ids[i]       = 2 * (i + 1);
      elem_map_elem_cnts[i] = ECNTCM;
    }

    int error = ne_put_cmap_params(fileid, node_map_node_cnts, node_map_ids, elem_map_elem_cnts,
                                   elem_map_ids, iproc);
    if (error < 0) {
      return error;
    }
  }

  return 0;
}

/*****************************************************************************/
int ne_test_pncm(int fileid)
{
  int node_map_ids[NNCMAP], node_ids[NCNTCM], proc_ids[NCNTCM];

  for (int iproc = 0; iproc < NPROCF; iproc++) {
    for (int i = 0; i < NNCMAP; i++) {
      node_map_ids[i] = (i + 1);
    }
    for (int i = 0; i < NCNTCM; i++) {
      node_ids[i] = 2 * (i + 1);
      proc_ids[i] = 3 * (i + 1);
    }

    for (int i = 0; i < NNCMAP; i++) {
      int error = ne_put_node_cmap(fileid, node_map_ids[i], node_ids, proc_ids, iproc);
      if (error < 0) {
        return error;
      }
    }
  }

  return 0;
}

/*****************************************************************************/
int ne_test_pecm(int fileid)
{
  int elem_map_ids[NECMAP], elem_ids[ECNTCM], side_ids[ECNTCM];
  int proc_ids[ECNTCM];

  for (int iproc = 0; iproc < NPROCF; iproc++) {
    for (int i = 0; i < NECMAP; i++) {
      elem_map_ids[i] = 2 * (i + 1);
    }
    for (int i = 0; i < ECNTCM; i++) {
      elem_ids[i] = 2 * (i + 1);
      side_ids[i] = 3 * (i + 1);
      proc_ids[i] = 4 * (i + 1);
    }

    for (int i = 0; i < NECMAP; i++) {
      int error = ne_put_elem_cmap(fileid, elem_map_ids[i], elem_ids, side_ids, proc_ids, iproc);
      if (error < 0) {
        return error;
      }
    }
  }

  return 0;
}

/*****************************************************************************/
int ne_test_giinf(int fileid)
{
  int  nproc, nprocf;
  char ftype[2];

  int error = ne_get_init_info(fileid, &nproc, &nprocf, ftype);

  if (error < 0) {
    return error;
  }

  if (nproc != NPROC) {
    return -1;
  }
  if (nprocf != NPROCF) {
    return -1;
  }
  if (strcmp(ftype, "s") != 0) {
    return -1;
  }

  return 0;
}

/*****************************************************************************/
int ne_test_ginig(int fileid)
{
  int num_nodes_g, num_elems_g, num_elem_blks_g, num_ns_g, num_ss_g;

  int error = ne_get_init_global(fileid, &num_nodes_g, &num_elems_g, &num_elem_blks_g, &num_ns_g,
                                 &num_ss_g);

  if (error < 0) {
    return error;
  }

  if (num_nodes_g != NNG) {
    return -1;
  }
  if (num_elems_g != NEG) {
    return -1;
  }
  if (num_elem_blks_g != NEBG) {
    return -1;
  }
  if (num_ns_g != NNSG) {
    return -1;
  }
  if (num_ss_g != NSSG) {
    return -1;
  }

  return 0;
}

/*****************************************************************************/
int ne_test_gelbid(int fileid)
{
  int el_blk_ids[NEBG], el_blk_cnt[NEBG];

  int error = ne_get_eb_info_global(fileid, el_blk_ids, el_blk_cnt);

  if (error < 0) {
    return error;
  }

  for (int i = 0; i < NEBG; i++) {
    if (el_blk_ids[i] != (i + 1)) {
      return -1;
    }
    if (el_blk_cnt[i] != NEBCG) {
      return -1;
    }
  }

  return 0;
}

/*****************************************************************************/
int ne_test_gnsp(int fileid)
{
  int global_ids[NNSG], global_n_cnts[NNSG], global_df_cnts[NNSG];

  int error = ne_get_ns_param_global(fileid, global_ids, global_n_cnts, global_df_cnts);

  if (error < 0) {
    return error;
  }

  for (int i = 0; i < NNSG; i++) {
    if (global_ids[i] != 2 * (i + 1)) {
      return -1;
    }
    if (global_n_cnts[i] != 3 * (i + 1)) {
      return -1;
    }
    if (global_df_cnts[i] != 1) {
      return -1;
    }
  }

  return 0;
}

/*****************************************************************************/
int ne_test_gssp(int fileid)
{
  int global_ids[NSSG], global_e_cnts[NSSG], global_df_cnts[NSSG];

  int error = ne_get_ss_param_global(fileid, global_ids, global_e_cnts, global_df_cnts);

  if (error < 0) {
    return error;
  }

  for (int i = 0; i < NSSG; i++) {
    if (global_ids[i] != 3 * (i + 1)) {
      return -1;
    }
    if (global_e_cnts[i] != 2 * (i + 1)) {
      return -1;
    }
    if (global_df_cnts[i] != 1) {
      return -1;
    }
  }

  return 0;
}

/*****************************************************************************/
int ne_test_glbp(int fileid)
{
  int nintn, nborn, nextn, ninte, nbore, nncmap, necmap;

  for (int iproc = 0; iproc < NPROCF; iproc++) {
    int error = ne_get_loadbal_param(fileid, &nintn, &nborn, &nextn, &ninte, &nbore, &nncmap,
                                     &necmap, iproc);

    if (error < 0) {
      return error;
    }

    if (nintn != NINTN) {
      return -1;
    }
    if (nborn != NBORN) {
      return -1;
    }
    if (nextn != NEXTN) {
      return -1;
    }
    if (ninte != NINTE) {
      return -1;
    }
    if (nbore != NBORE) {
      return -1;
    }
    if (nncmap != NNCMAP) {
      return -1;
    }
    if (necmap != NECMAP) {
      return -1;
    }
  }

  return 0;
}

/*****************************************************************************/
int ne_test_gnm(int fileid)
{
  int j1 = 0;
  int node_mapi[NINTN], node_mapb[NBORN], node_mape[NEXTN];

  for (int iproc = 0; iproc < NPROCF; iproc++) {
    int error = ne_get_node_map(fileid, node_mapi, node_mapb, node_mape, iproc);

    if (error < 0) {
      return error;
    }

    for (int j = 0; j < NINTN; j++) {
      if (node_mapi[j] != j1++) {
        return -1;
      }
    }
    for (int j = 0; j < NBORN; j++) {
      if (node_mapb[j] != j1++) {
        return -1;
      }
    }
    for (int j = 0; j < NEXTN; j++) {
      if (node_mape[j] != j1++) {
        return -1;
      }
    }

    j1 = 0;
  }

  return 0;
}

/*****************************************************************************/
int ne_test_gem(int fileid)
{
  int j1 = 0;
  int elem_mapi[NINTE], elem_mapb[NBORE];

  for (int iproc = 0; iproc < NPROCF; iproc++) {

    int error = ne_get_elem_map(fileid, elem_mapi, elem_mapb, iproc);

    if (error < 0) {
      return error;
    }

    for (int j = 0; j < NINTE; j++) {
      if (elem_mapi[j] != j1++) {
        return -1;
      }
    }
    for (int j = 0; j < NBORE; j++) {
      if (elem_mapb[j] != j1++) {
        return -1;
      }
    }
    j1 = 0;
  }

  return 0;
}

/*****************************************************************************/
int ne_test_gncm(int fileid)
{
  int node_map_ids[NNCMAP], node_map_cnts[NNCMAP];
  int node_ids[NCNTCM], proc_ids[NCNTCM];

  for (int iproc = 0; iproc < NPROCF; iproc++) {

    int error = ne_get_cmap_params(fileid, node_map_ids, node_map_cnts, NULL, NULL, iproc);

    if (error < 0) {
      return error;
    }

    for (int i = 0; i < NNCMAP; i++) {
      error = ne_get_node_cmap(fileid, node_map_ids[i], node_ids, proc_ids, iproc);

      if (error < 0) {
        return error;
      }

      for (int j = 0; j < NCNTCM; j++) {
        if (node_ids[j] != 2 * (j + 1)) {
          return -1;
        }
        if (proc_ids[j] != 3 * (j + 1)) {
          return -1;
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************/
int ne_test_gecm(int fileid)
{
  int elem_ids[ECNTCM], elem_map_cnts[NECMAP], proc_ids[ECNTCM];
  int side_ids[ECNTCM], elem_map_ids[NECMAP];

  for (int iproc = 0; iproc < NPROCF; iproc++) {

    int error = ne_get_cmap_params(fileid, NULL, NULL, elem_map_ids, elem_map_cnts, iproc);

    if (error < 0) {
      return error;
    }

    for (int i = 0; i < NECMAP; i++) {
      error = ne_get_elem_cmap(fileid, elem_map_ids[i], elem_ids, side_ids, proc_ids, iproc);

      if (error < 0) {
        return error;
      }

      for (int j = 0; j < ECNTCM; j++) {
        if (elem_ids[j] != 2 * (j + 1)) {
          return -1;
        }
        if (side_ids[j] != 3 * (j + 1)) {
          return -1;
        }
        if (proc_ids[j] != 4 * (j + 1)) {
          return -1;
        }
      }
    }
  }

  return 0;
}

/*****************************************************************************/
int ne_test_plbpc(int fileid)
{
  int num_int_nodes[NPROCF], num_bor_nodes[NPROCF], num_ext_nodes[NPROCF];
  int num_int_elems[NPROCF], num_bor_elems[NPROCF];
  int num_node_cmaps[NPROCF], num_elem_cmaps[NPROCF];

  /* Set up the vectors */
  for (int iproc = 0; iproc < NPROCF; iproc++) {
    num_int_nodes[iproc] = NINTN;
    num_bor_nodes[iproc] = NBORN;
    num_ext_nodes[iproc] = NEXTN;

    num_int_elems[iproc] = NINTE;
    num_bor_elems[iproc] = NBORE;

    num_node_cmaps[iproc] = NNCMAP;
    num_elem_cmaps[iproc] = NECMAP;
  }

  return (ne_put_loadbal_param_cc(fileid, num_int_nodes, num_bor_nodes, num_ext_nodes,
                                  num_int_elems, num_bor_elems, num_node_cmaps, num_elem_cmaps));
}

/*****************************************************************************/
int ne_test_pcmpc(int fileid)
{
  int nmap_ids[NNCMAP * NPROCF], nmap_n_cnts[NNCMAP * NPROCF];
  int nmap_proc_ptr[NPROCF + 1];
  int emap_ids[NECMAP * NPROCF], emap_e_cnts[NECMAP * NPROCF];
  int emap_proc_ptr[NPROCF + 1];

  nmap_proc_ptr[0] = 0;
  emap_proc_ptr[0] = 0;
  int n_cntr       = 0;
  int e_cntr       = 0;
  for (int iproc = 0; iproc < NPROCF; iproc++) {

    for (int j = 0; j < NNCMAP; j++) {
      nmap_ids[n_cntr]      = (j + 1);
      nmap_n_cnts[n_cntr++] = NCNTCM;
    }

    for (int j = 0; j < NECMAP; j++) {
      emap_ids[e_cntr]      = 2 * (j + 1);
      emap_e_cnts[e_cntr++] = ECNTCM;
    }

    nmap_proc_ptr[iproc + 1] = nmap_proc_ptr[iproc] + NNCMAP;
    emap_proc_ptr[iproc + 1] = emap_proc_ptr[iproc] + NECMAP;
  }

  return (ne_put_cmap_params_cc(fileid, nmap_ids, nmap_n_cnts, nmap_proc_ptr, emap_ids, emap_e_cnts,
                                emap_proc_ptr));
}
