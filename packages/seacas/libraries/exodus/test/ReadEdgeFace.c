/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include "exodusII.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h> /* for malloc(), free() */

#define EX_TEST_FILENAME "edgeFace.exo"

#define EXCHECK(funcall, errmsg)                                                                   \
  if ((funcall) < 0) {                                                                             \
    fprintf(stderr, errmsg);                                                                       \
    return 1;                                                                                      \
  }

#define EXCHKPI(funcall, errmsg, sucmsg, ival)                                                     \
  if ((funcall) < 0) {                                                                             \
    fprintf(stderr, errmsg);                                                                       \
    return 1;                                                                                      \
  }                                                                                                \
  else {                                                                                           \
    fprintf(stdout, sucmsg, ival);                                                                 \
  }

ex_entity_type obj_types[] = {EX_EDGE_BLOCK, EX_FACE_BLOCK, EX_ELEM_BLOCK, EX_NODE_SET,
                              EX_EDGE_SET,   EX_FACE_SET,   EX_SIDE_SET,   EX_ELEM_SET,
                              EX_NODE_MAP,   EX_EDGE_MAP,   EX_FACE_MAP,   EX_ELEM_MAP};

ex_inquiry obj_sizes[] = {
    EX_INQ_EDGE_BLK,  EX_INQ_FACE_BLK,  EX_INQ_ELEM_BLK,  EX_INQ_NODE_SETS,
    EX_INQ_EDGE_SETS, EX_INQ_FACE_SETS, EX_INQ_SIDE_SETS, EX_INQ_ELEM_SETS,
    EX_INQ_NODE_MAP,  EX_INQ_EDGE_MAP,  EX_INQ_FACE_MAP,  EX_INQ_ELEM_MAP,
};

const char *obj_typenames[] = {"   Edge block", "   Face block", "Element block", "    Node set",
                               "    Edge set",  "    Face set",  "    Side set",  " Element set",
                               "    Node map",  "    Edge map",  "    Face map",  " Element map"};

int obj_sizeinq[] = {EX_INQ_EDGE,
                     EX_INQ_FACE,
                     EX_INQ_ELEM,
                     EX_INQ_NS_NODE_LEN,
                     EX_INQ_ES_LEN,
                     EX_INQ_FS_LEN,
                     EX_INQ_SS_ELEM_LEN,
                     EX_INQ_ELS_LEN,
                     -1,
                     -1,
                     -1,
                     -1};

#define OBJECT_IS_BLOCK(i) ((i >= 0) && (i < 3))
#define OBJECT_IS_SET(i) ((i > 2) && (i < 8))

int cReadEdgeFace(int argc, char **argv)
{
  int            exoid;
  int            appWordSize  = 8;
  int            diskWordSize = 8;
  float          exoVersion;
  int            itmp[5];
  int *          ids;
  int            nids;
  int            obj;
  int            i, j;
  int            num_timesteps;
  int            ti;
  char **        obj_names;
  char **        var_names;
  int            have_var_names;
  int            num_vars;    /* number of variables per object */
  int            num_entries; /* number of values per variable per object */
  double *       entry_vals;  /* variable values for each entry of an object */
  ex_init_params modelParams;

  exoid = ex_open(EX_TEST_FILENAME, EX_READ, &appWordSize, &diskWordSize, &exoVersion);
  if (exoid <= 0) {
    fprintf(stderr, "Unable to open \"%s\" for reading.\n", EX_TEST_FILENAME);
    return 1;
  }

  EXCHECK(ex_get_init_ext(exoid, &modelParams), "Unable to read database parameters.\n");

  fprintf(stdout,
          "Title: <%s>\n"
          "Dimension: %" PRId64 "\n"
          "Nodes: %" PRId64 "\n"
          "Edges: %" PRId64 "\n"
          "Faces: %" PRId64 "\n"
          "Elements: %" PRId64 "\n"
          "Edge Blocks: %" PRId64 "\n"
          "Face Blocks: %" PRId64 "\n"
          "Element Blocks: %" PRId64 "\n"
          "Node Sets: %" PRId64 "\n"
          "Edge Sets: %" PRId64 "\n"
          "Face Sets: %" PRId64 "\n"
          "Side Sets: %" PRId64 "\n"
          "Element Sets: %" PRId64 "\n"
          "Node Maps: %" PRId64 "\n"
          "Edge Maps: %" PRId64 "\n"
          "Face Maps: %" PRId64 "\n"
          "Element Maps: %" PRId64 "\n",
          modelParams.title, modelParams.num_dim, modelParams.num_nodes, modelParams.num_edge,
          modelParams.num_face, modelParams.num_elem, modelParams.num_edge_blk,
          modelParams.num_face_blk, modelParams.num_elem_blk, modelParams.num_node_sets,
          modelParams.num_edge_sets, modelParams.num_face_sets, modelParams.num_side_sets,
          modelParams.num_elem_sets, modelParams.num_node_maps, modelParams.num_edge_maps,
          modelParams.num_face_maps, modelParams.num_elem_maps);

  num_timesteps = ex_inquire_int(exoid, EX_INQ_TIME);

  /* *** NEW API *** */
  for (i = 0; i < sizeof(obj_types) / sizeof(obj_types[0]); ++i) {
    int *truth_tab = 0;
    have_var_names = 0;

    EXCHECK(ex_inquire(exoid, obj_sizes[i], &nids, 0, 0),
            "Object ID list size could not be determined.\n");

    if (!nids) {
      fprintf(stdout, "=== %ss: none\n\n", obj_typenames[i]);
      continue;
    }
    else {
      fprintf(stdout, "=== %ss: %d\n", obj_typenames[i], nids);
    }

    ids       = (int *)malloc(nids * sizeof(int));
    obj_names = (char **)malloc(nids * sizeof(char *));
    for (obj = 0; obj < nids; ++obj) {
      obj_names[obj] = (char *)malloc((MAX_STR_LENGTH + 1) * sizeof(char));
    }

    EXCHECK(ex_get_ids(exoid, obj_types[i], ids), "Could not read object ids.\n");
    EXCHECK(ex_get_names(exoid, obj_types[i], obj_names), "Could not read object ids.\n");

    if ((OBJECT_IS_BLOCK(i)) || (OBJECT_IS_SET(i))) {
      int *tp;
      EXCHECK(ex_get_variable_param(exoid, obj_types[i], &num_vars),
              "Could not read number of variables.\n");

      if (num_vars && num_timesteps > 0) {
        truth_tab = (int *)malloc(num_vars * nids * sizeof(int));
        EXCHECK(ex_get_truth_table(exoid, obj_types[i], nids, num_vars, truth_tab),
                "Could not read truth table.\n");
        tp = truth_tab;
        fprintf(stdout, "Truth:");
        for (obj = 0; obj < nids; ++obj) {
          for (j = 0; j < num_vars; ++j, ++tp) {
            fprintf(stdout, " %d", *tp);
          }
          fprintf(stdout, "\n      ");
        }
        fprintf(stdout, "\n");

        var_names = (char **)malloc(num_vars * sizeof(char *));
        for (j = 0; j < num_vars; ++j) {
          var_names[j] = (char *)malloc((MAX_STR_LENGTH + 1) * sizeof(char));
        }

        EXCHECK(ex_get_variable_names(exoid, obj_types[i], num_vars, var_names),
                "Could not read variable names.\n");
        have_var_names = 1;
      }
    }

    if (!have_var_names) {
      var_names = 0;
    }

    for (obj = 0; obj < nids; ++obj) {
      if (obj_names[obj]) {
        fprintf(stdout, "%s %3d (%s): ", obj_typenames[i], ids[obj], obj_names[obj]);
      }
      else {
        fprintf(stdout, "%s %3d: ", obj_typenames[i], ids[obj]);
      }

      if (OBJECT_IS_BLOCK(i)) {
        int *nconn;
        int *econn;
        int *fconn;
        int  ele;
        int  ctr;
        int  num_attrs;
        if (obj_types[i] == EX_ELEM_BLOCK) {
          EXCHECK(ex_get_block(exoid, obj_types[i], ids[obj], 0, itmp, itmp + 1, itmp + 2, itmp + 3,
                               &num_attrs),
                  "Could not read block params.\n");
          fprintf(stdout,
                  "Entries: %3d Nodes/entry: %d Edges/entry: %d Faces/entry: %d Attributes: %d",
                  itmp[0], itmp[1], itmp[2], itmp[3], num_attrs);
        }
        else {
          EXCHECK(ex_get_block(exoid, obj_types[i], ids[obj], 0, itmp, itmp + 1, 0, 0, &num_attrs),
                  "Could not read block params.\n");
          fprintf(stdout, "Entries: %3d Nodes/entry: %d Attributes: %d", itmp[0], itmp[1],
                  num_attrs);
          itmp[2] = itmp[3] = 0;
        }
        fprintf(stdout, "\n   ");
        num_entries = itmp[0];
        nconn       = itmp[1] ? (int *)malloc(itmp[1] * num_entries * sizeof(int)) : 0;
        econn       = itmp[2] ? (int *)malloc(itmp[2] * num_entries * sizeof(int)) : 0;
        fconn       = itmp[3] ? (int *)malloc(itmp[3] * num_entries * sizeof(int)) : 0;
        EXCHECK(ex_get_conn(exoid, obj_types[i], ids[obj], nconn, econn, fconn),
                "Could not read connectivity.\n");
        for (ele = 0; ele < num_entries; ++ele) {
          for (ctr = 0; ctr < itmp[1]; ++ctr) {
            fprintf(stdout, " %2d", nconn[ele * itmp[1] + ctr]);
          }
          if (itmp[2]) {
            fprintf(stdout, "  ++");
            for (ctr = 0; ctr < itmp[2]; ++ctr) {
              fprintf(stdout, " %2d", econn[ele * itmp[2] + ctr]);
            }
          }
          if (itmp[3]) {
            fprintf(stdout, "  ++");
            for (ctr = 0; ctr < itmp[3]; ++ctr) {
              fprintf(stdout, " %2d", fconn[ele * itmp[3] + ctr]);
            }
          }
          fprintf(stdout, "\n   ");
        }
        free(nconn);
        free(econn);
        free(fconn);

        if (num_attrs) {
          char ** attr_names;
          double *attr;
          attr       = (double *)malloc(num_entries * num_attrs * sizeof(double));
          attr_names = (char **)malloc(num_attrs * sizeof(char *));
          for (j = 0; j < num_attrs; ++j) {
            attr_names[j] = (char *)malloc((MAX_STR_LENGTH + 1) * sizeof(char));
          }

          EXCHECK(ex_get_attr_names(exoid, obj_types[i], ids[obj], attr_names),
                  "Could not read attributes names.\n");
          EXCHECK(ex_get_attr(exoid, obj_types[i], ids[obj], attr),
                  "Could not read attribute values.\n");

          fprintf(stdout, "\n      Attributes:\n      ID ");
          for (j = 0; j < num_attrs; ++j) {
            fprintf(stdout, " %s", attr_names[j]);
          }
          fprintf(stdout, "\n");
          for (j = 0; j < num_entries; ++j) {
            int k;
            fprintf(stdout, "      %2d ", j + 1);
            for (k = 0; k < num_attrs; ++k) {
              fprintf(stdout, " %4.1f", attr[j * num_attrs + k]);
            }
            fprintf(stdout, "\n");
          }

          for (j = 0; j < num_attrs; ++j) {
            free(attr_names[j]);
          }
          free(attr_names);
          free(attr);
        }
      }
      else if (OBJECT_IS_SET(i)) {
        int     num_df;
        int *   set_entry;
        int *   set_extra;
        double *set_df;
        EXCHECK(ex_get_set_param(exoid, obj_types[i], ids[obj], &num_entries, &num_df),
                "Could not read set parameters.\n");

        set_entry = (int *)malloc(num_entries * sizeof(int));
        set_extra = (obj_types[i] != EX_NODE_SET && obj_types[i] != EX_ELEM_SET)
                        ? (int *)malloc(num_entries * sizeof(int))
                        : 0;
        EXCHECK(ex_get_set(exoid, obj_types[i], ids[obj], set_entry, set_extra),
                "Could not read set.\n");
        fprintf(stdout, "Entries: %3d Distribution factors: %3d\n", num_entries, num_df);
        if (set_extra) {
          for (j = 0; j < num_entries; ++j) {
            fprintf(stdout, "      %2d %2d\n", set_entry[j], set_extra[j]);
          }
        }
        else {
          for (j = 0; j < num_entries; ++j) {
            fprintf(stdout, "      %2d\n", set_entry[j]);
          }
        }
        free(set_entry);
        free(set_extra);

        set_df = num_df ? (double *)malloc(num_df * sizeof(double)) : 0;
        if (set_df) {
          EXCHECK(ex_get_set_dist_fact(exoid, obj_types[i], ids[obj], set_df),
                  "Could not read set distribution factors.\n");
          fprintf(stdout, "\n    Distribution factors:\n");
          for (j = 0; j < num_df; ++j) {
            fprintf(stdout, "      %4.1f\n", set_df[j]);
          }
          free(set_df);
        }
      }
      else { /* object is map */
        int *map;
        switch (obj_types[i]) {
        case EX_NODE_MAP: num_entries = modelParams.num_nodes; break;
        case EX_EDGE_MAP: num_entries = modelParams.num_edge; break;
        case EX_FACE_MAP: num_entries = modelParams.num_face; break;
        case EX_ELEM_MAP: num_entries = modelParams.num_elem; break;
        default: num_entries = 0;
        }
        if (num_entries) {
          fprintf(stdout, "Entries: %3d\n                :", num_entries);
          map = (int *)malloc(num_entries * sizeof(int));
          EXCHECK(ex_get_num_map(exoid, obj_types[i], ids[obj], map), "Could not read map.\n");
          for (j = 0; j < num_entries; ++j) {
            fprintf(stdout, " %d", map[j]);
          }
          free(map);
        }
        else {
          fprintf(stdout, "Entries: none");
        }
      }
      fprintf(stdout, "\n");

      /* Read results variables */
      if (((OBJECT_IS_BLOCK(i)) || (OBJECT_IS_SET(i))) && num_vars && num_timesteps > 0) {
        /* Print out all the time values to exercise get_var */
        entry_vals = (double *)malloc(num_entries * sizeof(double));
        for (j = 0; j < num_vars; ++j) {
          int k;
          if (!truth_tab[num_vars * obj + j]) {
            continue;
          }

          fprintf(stdout, "      Variable: %s", var_names[j]);
          for (ti = 1; ti <= num_timesteps; ++ti) {
            EXCHECK(ex_get_var(exoid, ti, obj_types[i], 1 + j, ids[obj], num_entries, entry_vals),
                    "Could not read variable values.\n");

            fprintf(stdout, "\n       @t%d ", ti);
            for (k = 0; k < num_entries; ++k) {
              fprintf(stdout, " %4.1f", entry_vals[k]);
            }
          }
          fprintf(stdout, "\n");
        }
        fprintf(stdout, "\n");
        free(entry_vals);
      }
    }

    if (((OBJECT_IS_BLOCK(i)) || (OBJECT_IS_SET(i))) && num_vars && num_timesteps > 0) {
      /* Print out one element's time values to exercise get_var_time */
      entry_vals = (double *)malloc(num_timesteps * sizeof(double));
      EXCHECK(ex_inquire(exoid, obj_sizeinq[i], itmp, 0, 0), "Inquire failed.\n");
      itmp[1] = 11;
      while (itmp[1] > itmp[0]) {
        itmp[1] /= 2;
      }
      for (j = 0; j < num_vars; ++j) {
        /* FIXME: This works for the dataset created by CreateEdgeFace, but not for any dataset in
         * general since
         * NULL truth table entries may mean the referenced elements don't have variable values.
         */
        EXCHECK(ex_get_var_time(exoid, obj_types[i], j + 1, itmp[1], 1, num_timesteps, entry_vals),
                "Could not read variable over time.\n");
        fprintf(stdout, "    Variable over time: %s  Entry: %3d ", var_names[j], itmp[1]);
        for (ti = 1; ti <= num_timesteps; ++ti) {
          fprintf(stdout, " @t%d: %4.1f", ti, entry_vals[ti - 1]);
        }
        fprintf(stdout, "\n");
      }
      free(entry_vals);
    }

    if (var_names) {
      for (j = 0; j < num_vars; ++j) {
        free(var_names[j]);
      }
      free(var_names);
    }
    free(truth_tab);
    free(ids);

    for (obj = 0; obj < nids; ++obj) {
      free(obj_names[obj]);
    }
    free(obj_names);

    fprintf(stdout, "\n");
  }

  EXCHECK(ex_close(exoid), "Unable to close database.\n");

  return 0;
}

#if !defined(USING_CMAKE)
int main(int argc, char *argv[]) { return cReadEdgeFace(argc, argv); }
#endif
