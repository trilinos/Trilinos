/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#undef NDEBUG
#include "exodusII.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  int  exoid, num_dim, num_nodes, num_elem, num_elem_blk;
  int  num_elem_in_block[10], num_nodes_per_elem[10], num_attr[10];
  int  num_node_sets, num_side_sets, error;
  int  i, j, *connect;
  int  ebids[10], ids[10];
  int  CPU_word_size, IO_word_size;
  char title[MAX_LINE_LENGTH + 1], elem_type[MAX_STR_LENGTH + 1];

  float  version;
  float *attrib;
  float  x[100], y[100], z[100];
  char * coord_names[3];

  /* Coordinate Frames */
  int   cf_ids[2]        = {20, 13};
  float pt_coords[9 * 2] = {1, 0, 0, 1, 0, 1, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0};
  char  tags[2]          = {'r', 'c'};

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* Specify compute and i/o word size */

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 4; /* (4 bytes) */

  /* create EXODUS II file */

  exoid = ex_create("test.exo",     /* filename path */
                    EX_CLOBBER,     /* create mode */
                    &CPU_word_size, /* CPU float word size in bytes */
                    &IO_word_size); /* I/O float word size in bytes */
  /* initialize file with parameters */

  num_dim   = 3;
  num_nodes = 19;
  num_elem  = 12;
  ;
  num_elem_blk  = 1;
  num_node_sets = 0;
  num_side_sets = 0;

  error = ex_put_init(exoid, "This is testwt1", num_dim, num_nodes, num_elem, num_elem_blk,
                      num_node_sets, num_side_sets);
  assert(error == 0);

  /* write nodal coordinates values and names to database */

  /* Quad #1 */
  x[0]  = 1.0000000E+00;
  x[1]  = 5.0000000E-01;
  x[2]  = 1.0000000E+00;
  x[3]  = 1.0000000E+00;
  x[4]  = 7.5000000E-01;
  x[5]  = 5.0000000E-01;
  x[6]  = 1.0000000E+00;
  x[7]  = 7.5000000E-01;
  x[8]  = 1.0000000E+00;
  x[9]  = 5.0000000E-01;
  x[10] = 5.0000000E-01;
  x[11] = 5.0000000E-01;
  x[12] = 1.0000000E+00;
  x[13] = 1.0000000E+00;
  x[14] = 7.5000000E-01;
  x[15] = 7.5000000E-01;
  x[16] = 1.0000000E+00;
  x[17] = 7.5000000E-01;
  x[18] = 1.0000000E+00;

  y[0]  = 5.0000000E-01;
  y[1]  = 1.0000000E+00;
  y[2]  = 1.0000000E+00;
  y[3]  = 7.5000000E-01;
  y[4]  = 1.0000000E+00;
  y[5]  = 5.0000000E-01;
  y[6]  = 5.0000000E-01;
  y[7]  = 5.0000000E-01;
  y[8]  = 5.0000000E-01;
  y[9]  = 1.0000000E+00;
  y[10] = 7.5000000E-01;
  y[11] = 1.0000000E+00;
  y[12] = 1.0000000E+00;
  y[13] = 7.5000000E-01;
  y[14] = 1.0000000E+00;
  y[15] = 7.5000000E-01;
  y[16] = 1.0000000E+00;
  y[17] = 1.0000000E+00;
  y[18] = 7.5000000E-01;

  z[0]  = 5.0000000E-01;
  z[1]  = 5.0000000E-01;
  z[2]  = 5.0000000E-01;
  z[3]  = 5.0000000E-01;
  z[4]  = 5.0000000E-01;
  z[5]  = 1.0000000E+00;
  z[6]  = 1.0000000E+00;
  z[7]  = 1.0000000E+00;
  z[8]  = 7.5000000E-01;
  z[9]  = 1.0000000E+00;
  z[10] = 1.0000000E+00;
  z[11] = 7.5000000E-01;
  z[12] = 1.0000000E+00;
  z[13] = 1.0000000E+00;
  z[14] = 1.0000000E+00;
  z[15] = 1.0000000E+00;
  z[16] = 7.5000000E-01;
  z[17] = 7.5000000E-01;
  z[18] = 7.5000000E-01;

  error = ex_put_coord(exoid, x, y, z);
  assert(error == 0);

  coord_names[0] = "xcoor";
  coord_names[1] = "ycoor";
  coord_names[2] = "zcoor";

  error = ex_put_coord_names(exoid, coord_names);
  assert(error == 0);

  /* write element block parameters */
  num_elem_in_block[0]  = 12;
  num_nodes_per_elem[0] = 4;
  ebids[0]              = 10;
  num_attr[0]           = 3;

  error = ex_put_block(exoid, EX_ELEM_BLOCK, ebids[0], "quad", num_elem_in_block[0],
                       num_nodes_per_elem[0], 0, 0, num_attr[0]);
  assert(error == 0);

  /* write element connectivity */
  connect = (int *)calloc(num_elem_in_block[0] * num_nodes_per_elem[0], sizeof(int));

  connect[0]  = 1;
  connect[1]  = 4;
  connect[2]  = 19;
  connect[3]  = 9;
  connect[4]  = 4;
  connect[5]  = 3;
  connect[6]  = 17;
  connect[7]  = 19;
  connect[8]  = 3;
  connect[9]  = 5;
  connect[10] = 18;
  connect[11] = 17;
  connect[12] = 5;
  connect[13] = 2;
  connect[14] = 12;
  connect[15] = 18;
  connect[16] = 9;
  connect[17] = 19;
  connect[18] = 14;
  connect[19] = 7;
  connect[20] = 7;
  connect[21] = 14;
  connect[22] = 16;
  connect[23] = 8;
  connect[24] = 19;
  connect[25] = 17;
  connect[26] = 13;
  connect[27] = 14;
  connect[28] = 17;
  connect[29] = 18;
  connect[30] = 15;
  connect[31] = 13;
  connect[32] = 14;
  connect[33] = 13;
  connect[34] = 15;
  connect[35] = 16;
  connect[36] = 8;
  connect[37] = 16;
  connect[38] = 11;
  connect[39] = 6;
  connect[40] = 18;
  connect[41] = 12;
  connect[42] = 10;
  connect[43] = 15;
  connect[44] = 16;
  connect[45] = 15;
  connect[46] = 10;
  connect[47] = 11;

  error = ex_put_conn(exoid, EX_ELEM_BLOCK, ebids[0], connect, NULL, NULL);
  assert(error == 0);
  free(connect);

  /* write element block attributes  (3 per block) */
  attrib = (float *)calloc(num_elem_in_block[0] * num_attr[0], sizeof(float));

#if 0
  {
    k = 0;
    for (i=0; i < num_elem_in_block[0]; i++) {
      for (j = 0; j < num_attr[0]; j++) {
        attrib[k++] = 10*(i+1) + j+1;
      }
    }
  }

  error = ex_put_attr (exoid, EX_ELEM_BLOCK, ebids[0], &attrib[0]);
  assert(error == 0);
#else
  {
    for (j = 0; j < num_attr[0]; j++) {
      for (i = 0; i < num_elem_in_block[0]; i++) {
        attrib[i] = 10 * (i + 1) + j + 1;
      }
      error = ex_put_one_attr(exoid, EX_ELEM_BLOCK, ebids[0], j + 1, &attrib[0]);
      assert(error == 0);
    }
  }
#endif
  free(attrib);

  /* Add a coordinate frame just to give test coverage... */
  {
    error = ex_put_coordinate_frames(exoid, 2, cf_ids, pt_coords, tags);
    assert(error == 0);
  }

  /* close the EXODUS files
   */
  error = ex_close(exoid);
  assert(error == 0);

  /* Reopen the file and read the attributes to see if they were written correctly */
  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 0; /* use what is stored in file */

  /* open EXODUS II files */

  exoid = ex_open("test.exo",     /* filename path */
                  EX_READ,        /* access mode = READ */
                  &CPU_word_size, /* CPU word size */
                  &IO_word_size,  /* IO word size */
                  &version);      /* ExodusII library version */

  assert(exoid >= 0);
  if (exoid < 0) {
    exit(1);
  }

  error = ex_get_init(exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk, &num_node_sets,
                      &num_side_sets);

  assert(error == 0);

  if (num_elem_blk > 0) {
    error = ex_get_ids(exoid, EX_ELEM_BLOCK, ids);
    assert(error == 0);

    for (i = 0; i < num_elem_blk; i++) {
      error = ex_get_block(exoid, EX_ELEM_BLOCK, ids[i], elem_type, &(num_elem_in_block[i]),
                           &(num_nodes_per_elem[i]), NULL, NULL, &(num_attr[i]));
      assert(error == 0);
    }

    /* read element block attributes */

    attrib = (float *)calloc(num_elem_in_block[0], sizeof(float));
    for (j = 0; j < num_attr[0]; j++) {
      error = ex_get_one_attr(exoid, EX_ELEM_BLOCK, ids[0], j + 1, &attrib[0]);
      assert(error == 0);

      if (error == 0) {
        for (i = 0; i < num_elem_in_block[0]; i++) {
          assert(attrib[i] == 10 * (i + 1) + j + 1);
        }
      }
    }
    free(attrib);
  }

  /* Read the coordinate frame... */
  {
    int   nframes;
    int   rdcf_ids[2];
    float rdpt_coords[9 * 2];
    char  rdtags[2];

    nframes = ex_inquire_int(exoid, EX_INQ_COORD_FRAMES);
    assert(nframes == 2);

    error = ex_get_coordinate_frames(exoid, &nframes, rdcf_ids, rdpt_coords, rdtags);
    assert(error == 0);
    assert(rdtags[0] == tags[0] && rdtags[1] == tags[1]);
    assert(rdcf_ids[0] == cf_ids[0] && rdcf_ids[1] == cf_ids[1]);

    for (i = 0; i < nframes * 9; i++) {
      assert(rdpt_coords[i] == pt_coords[i]);
    }
  }
  error = ex_close(exoid);
  assert(error == 0);
  return 0;
}
