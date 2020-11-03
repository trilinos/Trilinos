/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <assert.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "exodusII.h"

#define DEFAULT_FILE_NAME "mesh"
#define DEFAULT_MAP_ORIGIN 1
#define DEFAULT_NUM_DOMAINS 1
#define DEFAULT_NUM_ELEMENTS 1000000
#define DEFAULT_NUM_FIELDS 0
#define DEFAULT_NUM_TIMESTEPS 1

#define MAX_STRING_LEN 128
#define NUM_BYTES_PER_INT 4
#define NUM_NODES_PER_ELEM 8

#define EBLK_ID 100000
#define EXODUSII_FILE_TYPE ".e"

typedef double realtyp;

typedef int64_t INT;

INT StringToCount(char *size_str)
{
  INT  size = 0;
  char range;
  int  rc;

  rc = sscanf(size_str, "%" PRId64 "%c", &size, &range);
  if (rc == 2) {
    switch ((int)range) {
    case 'k':
    case 'K': size *= 1000; break;
    case 'm':
    case 'M': size *= 1000000; break;
    case 'g':
    case 'G': size *= 1000000000; break;
    }
  }
  else if (rc == 0) {
    size = -1;
  }
  return (size);
} /* StringToCount() */

void get_file_name(const char *base, const char *ext, int rank, int nprocs, const char *other,
                   char *output);

/* We need to do a cube-root to find the number of elements on each
 * side of the cube, but don't want to link in -lm just for
 * this. Use this routine which is better than
 * a brute force approach. Found at
 * http://www.hackersdelight.org/HDcode/icbrt.c
 */
INT icbrt(unsigned x)
{
  INT      s;
  unsigned y, b;

  s = 30;
  y = 0;
  while (s >= 0) { /* Do 11 times. */
    y = 2 * y;
    b = (3 * y * (y + 1) + 1) << s;
    s = s - 3;
    if (x >= b) {
      x = x - b;
      y = y + 1;
    }
  }
  return y;
}

/* Prototypes */
void create_rr_elem_map(INT loc_num_elements, INT *elem_map, INT map_origin, INT num_domains,
                        INT current_domain);

void create_elem_map(INT loc_num_elems, INT elem_num, INT *elem_map, INT map_origin);

void create_local_connect(INT *node_map, INT len_node_map, INT len_connect, INT *domain_connect,
                          INT *loc_connect, INT map_origin);

void extract_connect(INT element_offset, INT num_elem, INT *elem_map, INT *connect,
                     INT *domain_connect, INT map_origin);

void make_mesh(realtyp *x, realtyp *y, realtyp *z, INT *connect, INT map_origin,
               INT num_elements_1d);

void parse_input(int argc, char *argv[], bool *debug, INT *map_origin, INT *num_elements_1d,
                 INT *num_domains, INT *num_nodal_fields, INT *num_global_fields,
                 INT *num_element_fields, INT *num_timesteps, char *file_name, int *exodus,
                 int *compression_level, int *shuffle, int *int64bit);

void write_exo_mesh(int debug, char *file_name, INT map_origin, INT num_nodes, INT num_elements,
                    INT num_domains, INT num_nodal_fields, INT num_global_fields,
                    INT num_element_fields, INT num_timesteps, realtyp *x, realtyp *y, realtyp *z,
                    INT *connect, int compression_level, int shuffle, int int64bit);

void create_node_map(INT len_map, INT len_connect, INT *domain_connect, INT *node_map,
                     INT *loc_num_nodes, INT map_origin);

INT bin_search2(INT value, INT num, INT List[]);

/***********************************************************************
 *
 *  Main function
 *
 ***********************************************************************/

int main(int argc, char *argv[])
{
  INT *connect;
  bool debug = false; /* true, display debug information; false       */
  /* otherwise.                                 */
  static char file_name[MAX_STRING_LEN] = DEFAULT_FILE_NAME;
  int         exodus                    = true;
  INT         map_origin                = DEFAULT_MAP_ORIGIN;
  INT         num_domains               = DEFAULT_NUM_DOMAINS;
  INT         num_elements_1d;
  INT         num_elements       = DEFAULT_NUM_ELEMENTS;
  INT         num_nodal_fields   = DEFAULT_NUM_FIELDS;
  INT         num_global_fields  = DEFAULT_NUM_FIELDS;
  INT         num_element_fields = DEFAULT_NUM_FIELDS;
  INT         num_timesteps      = DEFAULT_NUM_TIMESTEPS;
  INT         num_nodes;
  int         compression_level = 0;
  int         shuffle           = 0;
  int         int64bit          = 0;
  size_t      size;

  realtyp *x;
  realtyp *y;
  realtyp *z;

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* Parse Input */
  parse_input(argc, argv, &debug, &map_origin, &num_elements, &num_domains, &num_nodal_fields,
              &num_global_fields, &num_element_fields, &num_timesteps, file_name, &exodus,
              &compression_level, &shuffle, &int64bit);

  /* Create Coordinates and Connectivity Array */
  num_elements_1d = icbrt(num_elements);
  num_nodes       = (num_elements_1d + 1) * (num_elements_1d + 1) * (num_elements_1d + 1);
  x               = malloc(num_nodes * sizeof(realtyp));
  y               = malloc(num_nodes * sizeof(realtyp));
  z               = malloc(num_nodes * sizeof(realtyp));
  assert(x != NULL && y != NULL && z != NULL);

  num_elements = num_elements_1d * num_elements_1d * num_elements_1d;
  size         = (size_t)NUM_NODES_PER_ELEM * num_elements * sizeof(INT);
  assert(size > 0);
  connect = malloc(size);
  assert(connect != NULL);

  fprintf(stderr, "Creating a 3D mesh of %" PRId64 " hex elements and %" PRId64 " nodes.\n",
          num_elements, num_nodes);

  make_mesh(x, y, z, connect, map_origin, num_elements_1d);
  fprintf(stderr, "\t...Mesh topology created.\n");

  /*
   *    Write Out Mesh
   */

  if (exodus) {
    write_exo_mesh(debug, file_name, map_origin, num_nodes, num_elements, num_domains,
                   num_nodal_fields, num_global_fields, num_element_fields, num_timesteps, x, y, z,
                   connect, compression_level, shuffle, int64bit);
  }

  free(x);
  free(y);
  free(z);
  free(connect);
  return 0;
} /* end of main() */

/***********************************************************************
 ***********************************************************************/
void parse_input(int argc, char *argv[], bool *debug, INT *map_origin, INT *num_elements_1d,
                 INT *num_domains, INT *num_nodal_fields, INT *num_global_fields,
                 INT *num_element_fields, INT *num_timesteps, char *file_name, int *exodus,
                 int *compression_level, int *shuffle, int *int64bit)
{
  int arg = 0; /* Argument index.      */

  while (++arg < argc) {
    if (strcmp("-c", argv[arg]) == 0) {
      if (++arg < argc) {
        *num_nodal_fields = atoi(argv[arg]);
      }
    }
    else if (strcmp("-compress", argv[arg]) == 0) {
      if (++arg < argc) {
        *compression_level = atoi(argv[arg]);
      }
    }
    else if (strcmp("-shuffle", argv[arg]) == 0) {
      *shuffle = 1;
    }
    else if (strcmp("-64", argv[arg]) == 0) {
      *int64bit = 1;
    }
    else if (strcmp("-nv", argv[arg]) == 0) {
      if (++arg < argc) {
        *num_nodal_fields = atoi(argv[arg]);
      }
    }
    else if (strcmp("-gv", argv[arg]) == 0) {
      if (++arg < argc) {
        *num_global_fields = atoi(argv[arg]);
      }
    }
    else if (strcmp("-ev", argv[arg]) == 0) {
      if (++arg < argc) {
        *num_element_fields = atoi(argv[arg]);
      }
    }
    else if (strcmp("-t", argv[arg]) == 0) {
      if (++arg < argc) {
        *num_timesteps = atoi(argv[arg]);
      }
    }
    else if (strcmp("-d", argv[arg]) == 0) {
      *debug = true;
    }
    else if (strcmp("-f", argv[arg]) == 0) {
      if (++arg < argc) {
        ex_copy_string(file_name, argv[arg], MAX_STRING_LEN);
      }
    }
    else if (strcmp("-m", argv[arg]) == 0) {
      if (++arg < argc) {
        *map_origin = atoi(argv[arg]);
      }
    }
    else if (strcmp("-n", argv[arg]) == 0) {
      if (++arg < argc) {
        *num_elements_1d = StringToCount(argv[arg]);
      }
    }
    else if (strcmp("-p", argv[arg]) == 0) {
      if (++arg < argc) {
        *num_domains = atoi(argv[arg]);
      }
    }
    else if (strcmp("-x", argv[arg]) == 0) {
      *exodus = true;
    }
    else if ((strcmp("-h", argv[arg]) == 0) || (strcmp("-u", argv[arg]) == 0)) {
      printf("                                                                \n");
      printf("NAME                                                            \n");
      printf("                                                                \n");
      printf("create_mesh - creates a mesh file for performance benchmarking. \n");
      printf("                                                                \n");
      printf("SYNOPSIS                                                        \n");
      printf("                                                                \n");
      printf("create_mesh [-c fields] [-t timesteps] [-d] [-f file_name] \n");
      printf("            [-m map_origin] [-n elements] [-p domains]          \n");
      printf("            [-nv number] [-ev number] [-gv number] ");
      printf("            [-r] [-u] [-h] ");
      printf("[-x]");
      printf("                 \n");
      printf("                                                                \n");
      printf("DESCRIPTION                                                     \n");
      printf("                                                                \n");
      printf("This program creates a 3-D mesh for performance benchmarking.   \n");
      printf("The EXODUSII II database file(s) created by this       \n");
      printf("prrogram is/are read by the rd_wt_mesh program to perform the   \n");
      printf("actual benchmark.                                               \n");
      printf("                                                                \n");
      printf("OPTIONS                                                         \n");
      printf("                                                                \n");
      printf("-c  fields     number of nodal   fields. Default: %d            \n",
             DEFAULT_NUM_FIELDS);
      printf("-nv fields     number of nodal   fields. Default: %d            \n",
             DEFAULT_NUM_FIELDS);
      printf("-ev fields     number of element fields. Default: %d           \n",
             DEFAULT_NUM_FIELDS);
      printf("-gv fields     number of global  fields. Default: %d            \n",
             DEFAULT_NUM_FIELDS);
      printf("-t timesteps   number of timesteps. Default: %d                 \n",
             DEFAULT_NUM_TIMESTEPS);
      printf("-d             display debug information.                         \n");
      printf("-f file_name   file name prefix for all created files:          \n");
      printf("                                                                \n");
      printf("                  'file_name'_n%s [EXODUSII II file]              \n",
             EXODUSII_FILE_TYPE);
      printf("                                                                \n");
      printf("               where n varies from 0 to number of domains-1.  \n");
      printf("               Default: %s                                      \n",
             DEFAULT_FILE_NAME);
      printf("-h             display help/usage information.                  \n");
      printf("-m map_origin  element map origin. Default: %d                  \n",
             DEFAULT_MAP_ORIGIN);
      printf("-n elements    number of elements in mesh   \n");
      printf("               Can suffix with 'k', 'm', 'g' for thousand, million, billion\n");
      printf("               elements/file = elements/number_of_domains.      \n");
      printf("               Default: %d                                      \n",
             DEFAULT_NUM_ELEMENTS);
      printf("-p domains     number of domains. Default: %d                   \n",
             DEFAULT_NUM_DOMAINS);
      printf("-compress val  set compression to level 'val' [0..9]            \n");
      printf("-shuffle       enable hdf5-shuffle                              \n");
      printf("-64            enable 64-bit integers                           \n");
      printf("-u             display help/usage information.                  \n");

      exit(0);
    }
    else {
      fprintf(stderr, "Unknown option: %s\n", argv[arg]);
      fprintf(stderr, "Enter create_mesh -h for description of valid options.\n");

      exit(0);
    }
  }
}

/***********************************************************************
 *
 *  Create the coordinates and connectivity array for the mesh
 *
 ***********************************************************************/

void make_mesh(realtyp *x, realtyp *y, realtyp *z, INT *connect, INT map_origin,
               INT num_elements_1d)
{
  size_t i, j, k, m, base, cnt;
  size_t elp1sq = (num_elements_1d + 1) * (num_elements_1d + 1);

  /* create global coordinates */

  for (m = 0, k = 0; m < (num_elements_1d + 1); m++) {
    for (i = 0; i < (num_elements_1d + 1); i++) {
      for (j = 0; j < (num_elements_1d + 1); j++, k++) {
        x[k] = (realtyp)j;
        y[k] = (realtyp)i;
        z[k] = (realtyp)m;
      }
    }
  }

  /* build connectivity array (node list) for mesh */

  for (m = 0, k = 0, cnt = 0; m < num_elements_1d; m++) {
    for (i = 0, k = 0; i < num_elements_1d; i++) {
      for (j = 0; j < num_elements_1d; j++, k++) {
        base           = (m * elp1sq) + k + i + map_origin;
        connect[cnt++] = base;
        connect[cnt++] = base + 1;
        connect[cnt++] = base + num_elements_1d + 2;
        connect[cnt++] = base + num_elements_1d + 1;

        connect[cnt++] = elp1sq + base;
        connect[cnt++] = elp1sq + base + 1;
        connect[cnt++] = elp1sq + base + num_elements_1d + 2;
        connect[cnt++] = elp1sq + base + num_elements_1d + 1;
      }
    }
  }
} /* end of make_mesh() */

/***********************************************************************
 ***********************************************************************/
void write_exo_mesh(int debug, char *file_name, INT map_origin, INT num_nodes, INT num_elements,
                    INT num_domains, INT num_nodal_fields, INT num_global_fields,
                    INT num_element_fields, INT num_timesteps, realtyp *x, realtyp *y, realtyp *z,
                    INT *connect, int compression_level, int shuffle, int int64bit)
{
  int  CPU_word_size = sizeof(realtyp);
  int  IO_word_size  = sizeof(realtyp);
  int  exoid, err, num_dim, num_elem_blk, num_node_sets, num_side_sets;
  INT  i, j, t, index, loc_num_elements, loc_num_nodes, len_connect;
  INT *elem_map = NULL, *node_map = NULL, *domain_connect = NULL, *loc_connect = NULL;
  int *elem_var_tab;
  INT  accum_num_elements = 0;
  INT  loc_node_size      = -1;

  realtyp *loc_xcoords = NULL;
  realtyp *loc_ycoords = NULL;
  realtyp *loc_zcoords = NULL;
  realtyp *globals     = NULL;

  char   temporary_name[MAX_STRING_LEN];
  char **var_name;

  accum_num_elements = 0;
  for (i = 0; i < num_domains; i++) {
    int mymode = EX_MAPS_INT64_API | EX_BULK_INT64_API | EX_IDS_INT64_API;
    if (int64bit) {
      mymode |= EX_MAPS_INT64_DB | EX_BULK_INT64_DB | EX_IDS_INT64_DB;
    }

    /* create the EXODUSII file */
    get_file_name(file_name, "e", i, num_domains, NULL, temporary_name);

    exoid = ex_create(temporary_name, EX_CLOBBER | mymode, &CPU_word_size, &IO_word_size);

    if (exoid < 0) {
      fprintf(stderr, "after ex_create, error = %d\n", exoid);
      exit(-1);
    }

    ex_set_option(exoid, EX_OPT_COMPRESSION_LEVEL, compression_level);
    ex_set_option(exoid, EX_OPT_COMPRESSION_SHUFFLE, shuffle);

    if (num_domains > 1) {
      /* Determine local number of elements */
      if (num_elements < num_domains) {
        fprintf(stderr, "number of elements is less than number of domains.\n");
        if (i < num_elements) {
          loc_num_elements = 1;
        }
        else {
          loc_num_elements = 0;
        }
      }
      else {
        loc_num_elements = num_elements / num_domains;
        if (i < (num_elements % num_domains)) {
          loc_num_elements++;
        }
      }

      len_connect = NUM_NODES_PER_ELEM * loc_num_elements;

      /* malloc things we need */

      if (i == 0) { /* first time through; max size arrays occur on
                       first iteration */
        elem_map       = malloc(loc_num_elements * sizeof(INT));
        domain_connect = malloc(len_connect * sizeof(INT));
        loc_connect    = malloc(len_connect * sizeof(INT));
        node_map       = malloc(num_nodes * sizeof(INT));
      }

      /* Create element local/global map */
      create_elem_map(loc_num_elements, accum_num_elements, elem_map, map_origin);

      /* Extract current domain's connectivity, referencing global node ids */
      extract_connect(accum_num_elements, loc_num_elements, elem_map, connect, domain_connect,
                      map_origin);

      accum_num_elements += loc_num_elements;

      /* The local/global node map is just the current domain's connectivity,
         sorted with duplicate entries removed */
      create_node_map(num_nodes, len_connect, domain_connect, node_map, &loc_num_nodes, map_origin);

      /* Using local/global node map, convert the domain connectivity
         (referencing global node ids) to local connectivity (referencing
         local node ids) */

      create_local_connect(node_map, loc_num_nodes, len_connect, domain_connect, loc_connect,
                           map_origin);
    }
    else {
      loc_num_elements = num_elements;
      loc_num_nodes    = num_nodes;
    }

    if (debug) {
      fprintf(stderr, "\n\n\n");

      fprintf(stderr, "\n domain: %" PRId64 "\n", i);
      fprintf(stderr, "\n loc_num_elements: %" PRId64 "\n", loc_num_elements);
      fprintf(stderr, "\n loc_num_nodes: %" PRId64 "\n", loc_num_nodes);
    }

    num_dim       = 3;
    num_elem_blk  = 1;
    num_node_sets = 0;
    num_side_sets = 0;

    err = ex_put_init(exoid, "This is an EXODUSII performance test.", num_dim, loc_num_nodes,
                      loc_num_elements, num_elem_blk, num_node_sets, num_side_sets);

    if (err) {
      fprintf(stderr, "after ex_put_init, error = %d\n", err);
      ex_close(exoid);
      exit(-1);
    }

    /* Extract the local x and y coordinates */
    if (num_domains > 1) {
      if (loc_num_nodes > loc_node_size) {
        loc_xcoords   = realloc(loc_xcoords, loc_num_nodes * sizeof(realtyp));
        loc_ycoords   = realloc(loc_ycoords, loc_num_nodes * sizeof(realtyp));
        loc_zcoords   = realloc(loc_zcoords, loc_num_nodes * sizeof(realtyp));
        loc_node_size = loc_num_nodes;
      }

      for (j = 0; j < loc_num_nodes; j++) {
        index          = node_map[j] - map_origin;
        loc_xcoords[j] = x[index];
        loc_ycoords[j] = y[index];
        loc_zcoords[j] = z[index];
      }

      err = ex_put_coord(exoid, loc_xcoords, loc_ycoords, loc_zcoords);
    }
    else {
      err = ex_put_coord(exoid, x, y, z);
    }
    if (err) {
      fprintf(stderr, "after ex_put_coord, error = %d\n", err);
      ex_close(exoid);
      exit(-1);
    }
    if (debug) {
      fprintf(stderr, "\tCoordinates output.\n");
    }
#if 1
    {
      INT   ids[1] = {EBLK_ID};
      INT   num_elem_per_block[1];
      char *names[1] = {"hex"};
      INT   num_node_per_elem[1];
      INT   num_attr_per_block[1];
      bool  write_map       = num_domains > 1 ? true : false;
      num_elem_per_block[0] = loc_num_elements;
      num_node_per_elem[0]  = NUM_NODES_PER_ELEM;
      num_attr_per_block[0] = 0;
      err = ex_put_concat_elem_block(exoid, ids, names, num_elem_per_block, num_node_per_elem,
                                     num_attr_per_block, write_map);
    }
#else
    err = ex_put_block(exoid, EX_ELEM_BLOCK, 10000000000, "hex", loc_num_elements,
                       NUM_NODES_PER_ELEM, 0);
#endif

    if (err) {
      fprintf(stderr, "after ex_put_elem_block, error = %d\n", err);
      ex_close(exoid);
      exit(-1);
    }

    if (num_domains > 1) {
      err = ex_put_conn(exoid, EX_ELEM_BLOCK, EBLK_ID, loc_connect, NULL, NULL);
    }
    else {
      err = ex_put_conn(exoid, EX_ELEM_BLOCK, EBLK_ID, connect, NULL, NULL);
    }

    if (err) {
      fprintf(stderr, "after ex_put_elem_conn, error = %d\n", err);
      ex_close(exoid);
      exit(-1);
    }

    if (debug) {
      fprintf(stderr, "\tConnectivity output.\n");
    }
    /* write out element and node maps */

    if (num_domains > 1) {
      err = ex_put_id_map(exoid, EX_NODE_MAP, node_map);

      if (err) {
        fprintf(stderr, "after ex_put_id_map, error = %d\n", err);
        ex_close(exoid);
        exit(-1);
      }

      err = ex_put_id_map(exoid, EX_ELEM_MAP, elem_map);

      if (err) {
        fprintf(stderr, "after ex_put_id_map, error = %d\n", err);
        ex_close(exoid);
        exit(-1);
      }

      if (debug) {
        fprintf(stderr, "\tMaps output.\n");
      }
    }

    /* write out simulated results fields;
       we'll just write out the x coordinate field 'num_nodal_fields' times */
    if (loc_num_nodes < loc_num_elements) {
      fprintf(stderr, "INTERNAL ERROR: Programmer assumed number of nodes > number of elements, "
                      "but that is not true.\n");
      ex_close(exoid);
      exit(-1);
    }

    if (num_element_fields > 0) {
      elem_var_tab = malloc(num_element_fields * sizeof(int));
      for (j = 0; j < num_element_fields; j++) {
        elem_var_tab[j] = 1;
      }
    }
    else {
      elem_var_tab = 0;
    }
    err = ex_put_all_var_param(exoid, num_global_fields, num_nodal_fields, num_element_fields,
                               elem_var_tab, 0, 0, 0, 0);
    if (err) {
      fprintf(stderr, "after ex_put_all_var_param, error = %d\n", err);
      ex_close(exoid);
      exit(-1);
    }

    if (num_nodal_fields > 0) {

      var_name = malloc(num_nodal_fields * sizeof(char *));
      for (j = 0; j < num_nodal_fields; j++) {
        var_name[j] = malloc((MAX_STRING_LEN + 1) * sizeof(char));
        sprintf(var_name[j], "node_field_%" PRId64, j + 1);
      }
      err = ex_put_variable_names(exoid, EX_NODAL, num_nodal_fields, var_name);
      for (j = 0; j < num_nodal_fields; j++) {
        free(var_name[j]);
      }
      free(var_name);
    }

    if (num_global_fields > 0) {
      globals  = malloc(num_global_fields * sizeof(realtyp));
      var_name = malloc(num_global_fields * sizeof(char *));
      for (j = 0; j < num_global_fields; j++) {
        var_name[j] = malloc((MAX_STRING_LEN + 1) * sizeof(char));
        sprintf(var_name[j], "global_field_%" PRId64, j + 1);
        globals[j] = j;
      }
      err = ex_put_variable_names(exoid, EX_GLOBAL, num_global_fields, var_name);
      for (j = 0; j < num_global_fields; j++) {
        free(var_name[j]);
      }
      free(var_name);
    }

    if (num_element_fields > 0) {
      free(elem_var_tab);
      var_name = malloc(num_element_fields * sizeof(char *));
      for (j = 0; j < num_element_fields; j++) {
        var_name[j] = malloc((MAX_STRING_LEN + 1) * sizeof(char));
        sprintf(var_name[j], "element_field_%" PRId64, j + 1);
      }
      err = ex_put_variable_names(exoid, EX_ELEM_BLOCK, num_element_fields, var_name);
      for (j = 0; j < num_element_fields; j++) {
        free(var_name[j]);
      }
      free(var_name);
    }

    if (num_nodal_fields + num_global_fields + num_element_fields > 0) {
      fprintf(stderr, "Domain %" PRId64 "/%" PRId64 ", Writing Timestep: ", i + 1, num_domains);
      for (t = 0; t < num_timesteps; t++) {
        realtyp time = t;
        ex_put_time(exoid, t + 1, &time);
        fprintf(stderr, "%" PRId64 ", ", t + 1);
        if (num_global_fields > 0) {
          err = ex_put_var(exoid, t + 1, EX_GLOBAL, 1, 0, num_global_fields, globals);
          if (err) {
            fprintf(stderr, "after ex_put_global_var, error = %d\n", err);
            ex_close(exoid);
            exit(-1);
          }
        }
        for (j = 0; j < num_nodal_fields; j++) {
          err = ex_put_var(exoid, t + 1, EX_NODAL, j + 1, 0, loc_num_nodes, x);
          if (err) {
            fprintf(stderr, "after ex_put_nodal_var, error = %d\n", err);
            ex_close(exoid);
            exit(-1);
          }
        }
        for (j = 0; j < num_element_fields; j++) {
          err = ex_put_var(exoid, t + 1, EX_ELEM_BLOCK, j + 1, EBLK_ID, loc_num_elements, x);
          if (err) {
            fprintf(stderr, "after ex_put_element_var, error = %d\n", err);
            ex_close(exoid);
            exit(-1);
          }
        }
      }
      fprintf(stderr, "\n");
    }

    err = ex_close(exoid);

    if (err) {
      fprintf(stderr, "after ex_close, error = %d\n", err);
      exit(-1);
    }
    if (debug) {
      fprintf(stderr, "\tFile written.\n");
    }
  }

  /*
   * Free Memory
   */

  if (num_domains > 1) {
    free(domain_connect);
    free(elem_map);
    free(loc_connect);
    free(loc_xcoords);
    free(loc_ycoords);
    free(loc_zcoords);
    free(node_map);
  }
  if (num_global_fields > 0) {
    free(globals);
  }
}

/***********************************************************************
 *
 * Create element local/global map
 *
 * This puts contiguous groups of elements in each domain.  This is
 * a somewhat reasonable map for a realistic application.
 *
 ***********************************************************************/
void create_elem_map(INT loc_num_elems, INT elem_num, INT *elem_map, INT map_origin)
{
  INT i;

  for (i = 0; i < loc_num_elems; i++) {
    elem_map[i] = map_origin + elem_num++;
  }
}

/***********************************************************************
 *
 * Extract current domain's connectivity, referencing global node ids
 *
 * This extracts the "domain connectivity," that is, the connectivity
 * of the elements in the current domain.  The node ids in the domain
 * connectivity reference global node ids.
 *
 ***********************************************************************/

void extract_connect(INT element_offset, INT num_elem, INT *elem_map, INT *connect,
                     INT *domain_connect, INT map_origin)
{
  INT i, j, k, m, offset;

  for (i = element_offset, j = 0, m = 0; j < num_elem; j++) {
    if (elem_map[j] == i + map_origin) { /* extract this element */
      offset = (i * NUM_NODES_PER_ELEM);
      for (k = offset; k < offset + NUM_NODES_PER_ELEM; k++) {
        domain_connect[m++] = connect[k];
      }
      i++;
    }
  }
}

/***********************************************************************
 *
 * The local/global node map is just the current domain's connectivity,
 * sorted, with duplicate entries removed.  This isn't obvious, but
 * trust me.
 *
 ***********************************************************************/
void create_node_map(INT len_map, INT len_connect, INT *domain_connect, INT *node_map,
                     INT *loc_num_nodes, INT map_origin)
{
  INT cnt, i;
  for (i = 0; i < len_map; i++) {
    node_map[i] = 0;
  }

  for (i = 0; i < len_connect; i++) {
    node_map[domain_connect[i] - map_origin] = 1;
  }

  cnt = 0;
  for (i = 0; i < len_map; i++) {
    if (node_map[i] > 0) {
      node_map[cnt++] = i + map_origin;
    }
  }
  *loc_num_nodes = cnt;
}

/***********************************************************************
 *
 * Using local/global node map, convert the domain connectivity
 * (referencing global node ids) to local connectivity (referencing
 * local node ids).
 *
 * This requires inverting the local/global map, a relatively expensive
 * operation.  The procedure is:
 *
 *   for every entry in the domain connectivity
 *     search the node map until found
 *     set the value of the entry in the local connectivity to
 *       the index of the located value in the node map
 *
 ***********************************************************************/
void create_local_connect(INT *node_map, INT len_node_map, INT len_connect, INT *domain_connect,
                          INT *loc_connect, INT map_origin)
{
  INT i, index;

  for (i = 0; i < len_connect; i++) {
    index = bin_search2(domain_connect[i], len_node_map, node_map);
    if (index != -1) { /* found */
      loc_connect[i] = index + map_origin;
    }
    else {
      fprintf(stderr, "error creating local connectivity; i = %" PRId64 "\n", i);
      exit(-1);
    }
  }
}

/*****************************************************************************
 *
 * Searches a monotonic list of values for the value, 'value'.
 * It returns the index (0-based) of the first position found, which
 *   matches 'value'.
 * The list is assumed to be monotonic, and consist of elements
 *   list[0], ..., list[n-1].
 * If no position in list matches value, it returns the value -1.
 *
 *****************************************************************************/

INT bin_search2(INT value, INT num, INT List[])
{
  INT top, bottom = 0, middle, g_mid;

  top = num - 1;
  while (bottom <= top) {
    middle = (bottom + top) >> 1;
    g_mid  = List[middle];
    if (value < g_mid) {
      top = middle - 1;
    }
    else if (value > g_mid) {
      bottom = middle + 1;
    }
    else {
      return middle; /* found */
    }
  }
  return -1;
} /* bin_search2 */

/*****************************************************************************/
void get_file_name(const char *base, const char *ext, int rank, int nprocs, const char *other,
                   char *output)
{
  INT  i1, iTemp1;
  INT  iMaxDigit = 0, iMyDigit = 0;
  char cTemp[128];

  output[0] = '\0';
  ex_copy_string(output, base, MAX_STRING_LEN);
  strcat(output, ".");
  strcat(output, ext);
  if (other != NULL) {
    strcat(output, ".");
    strcat(output, other);
  }
  if (nprocs > 1) {
    /*
     * Find out the number of digits needed to specify the processor ID.
     * This allows numbers like 01-99, i.e., prepending zeros to the
     * name to preserve proper alphabetic sorting of the files.
     */

    iTemp1 = nprocs;
    do {
      iTemp1 /= 10;
      iMaxDigit++;
    } while (iTemp1 >= 1);

    iTemp1 = rank;
    do {
      iTemp1 /= 10;
      iMyDigit++;
    } while (iTemp1 >= 1);

    strcat(output, ".");
    sprintf(cTemp, "%d", nprocs);
    strcat(output, cTemp);
    strcat(output, ".");

    /*
     * Append the proper number of zeros to the filename.
     */
    for (i1 = 0; i1 < iMaxDigit - iMyDigit; i1++) {
      strcat(output, "0");
    }

    sprintf(cTemp, "%d", rank);
    strcat(output, cTemp);
  }
}
