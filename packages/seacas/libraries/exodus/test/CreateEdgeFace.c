/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include "exodusII.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EX_TEST_FILENAME "edgeFace.exo"

/* ======== Coordinates and connectivity ========= */
double coordsX[] = {0., 0., 0., 0., 3., 3., 3., 3., -3., -3., -3., -3.};

double coordsY[] = {-1., 1., 1., -1., -3., 3., 3., -3., -3., 3., 3., -3.};

double coordsZ[] = {-1., -1., 1., 1., -3., -3., 3., 3., -3., -3., 3., 3.};

const char *coordsNames[] = {"X", "Y", "Z"};

int conn1[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4};

int conn2[] = {1, 2, 3, 5};

int econn1[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 17, 18, 19, 20};

int fconn1[] = {4, 5, 7, 6, 3, 2, 8, 9, 11, 10, 1, 3};

int ebconn1[] = {1, 2, 2, 3, 3, 4,  4,  1,  5,  6,  6,  7, 7, 8, 8,  5, 1,  5, 2,  6,
                 4, 8, 3, 7, 9, 10, 10, 11, 11, 12, 12, 9, 9, 1, 10, 2, 12, 4, 11, 3};

int fbconn1[] = {12, 11, 10, 9, 5, 6, 7, 8};

int fbconn2[] = {1, 2, 3, 4};

int fbconn3[] = {1, 5, 6, 2,  3,  7, 8, 4,  2,  6, 7, 3,  4,  8, 5, 1,
                 9, 1, 2, 10, 11, 3, 4, 12, 10, 2, 3, 11, 12, 4, 1, 9};

int nmap1[] = {12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};

int edmap1[] = {1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20};

int famap1[] = {2, 1, 4, 3, 6, 5, 7, 8, 9, 10, 11};

int emap1[] = {1, 2};

const char *eblk_names[]  = {"Eli_WALLACH", "Angelo_NOVI"};
const char *edblk_names[] = {"Aldo_GIUFFRE"};
const char *fablk_names[] = {"Livio_LORENZON", "Claudio_SCARCHILLI", "John_BARTHA"};

const char *nmap_names[]  = {"Luigi_PISTILLI"};
const char *edmap_names[] = {"Antonio_CASALE"};
const char *famap_names[] = {"Sandro_SCARCHILLI"};
const char *emap_names[]  = {"Benito_STEFANELLI"};

/* ======== Sets ========= */
int nset_nodes[] = {5, 6, 9};

int eset_edges[] = {1, 2, 4, 15, 19, 20};

int eset_orient[] = {+1, +1, +1, +1, +1, -1};

double eset_df[] = {2., 2., 0.5, 0.5, 1., 1.};

int fset_faces[] = {3, 9};

int fset_orient[] = {+1, -1};

int sset_elems[] = {1, 1, 1, 2, 2};

int sset_sides[] = {1, 3, 5, 2, 4};

int elset_elems[] = {1, 2};

const char *elset_names[] = {"Clint_EASTWOOD", "Lee_VAN_CLEEF"};

const char *nset_names[] = {"Ennio_MORRICONE"};
const char *eset_names[] = {"Rada_RASSIMOV"};
const char *fset_names[] = {"Enzo_PETITO"};
const char *sset_names[] = {"Luciano_VINCENZONI"};

/* ======== Attributes ========= */
const char *edge_attr_names1[] = {"Sergio_LEONE"};

const char *face_attr_names1[] = {"GOOD"};
const char *face_attr_names2[] = {"BAD"};
const char *face_attr_names3[] = {"UGLY"};

const char *elem_attr_names1[] = {"SPAGHETTI", "WESTERN"};

double edge_attr_values1[] = {1.,  2.,  3.,  5.,  7.,  11., 13., 17., 19., 23.,
                              29., 31., 37., 41., 43., 47., 53., 59., 61., 67.};

double face_attr_values1[] = {71., 73.};

double face_attr_values2[] = {79.};

double face_attr_values3[] = {83., 89., 97., 101., 103., 107., 109., 113.};

double elem_attr_values1[] = {127., 101., 137., 139.};

/* ======== Results variables ========= */
/*           (2 time steps) */

double vals_glo_var[2][2] = {{36., 37.}, {42., 43.}};

double vals_nod_var[2][12] = {{0.1, 0.8, 0.0, 0.4, 0.3, 0.9, 0.8, 0.5, 0.3, 0.7, 0.4, 0.6},
                              {0.7, 0.5, 0.3, 0.5, 0.2, 0.7, 0.9, 0.8, 0.0, 0.2, 0.3, 0.5}};

double vals_edge_var1eb1[2][20] = {
    {20., 19., 18., 17., 16., 15., 14., 13., 12., 11., 10., 9., 8., 7., 6., 5., 4., 3., 2., 1.},
    {21., 20., 19., 18., 17., 16., 15., 14., 13., 12., 11., 10., 9., 8., 7., 6., 5., 4., 3., 2.}};

double vals_edge_var2eb1[2][20] = {
    {1., 1., 0., 0., 1., 1., 2., 0., 2., 0., 1., 1., 1., 1., 0., 0., 2., 2., 2., 2.},
    {1., 1., 0., 0., 1., 1., 2., 0., 2., 0., 1., 1., 1., 1., 0., 0., 2., 2., 2., 2.}};

double vals_face_var1fb1[2][2] = {{0, 1}, {2, 0}};

double vals_face_var1fb3[2][8] = {{1, 0, 2, 0, 3, 0, 4, 0}, {0, 1, 0, 2, 0, 3, 0, 4}};

double vals_elem_var1eb1[2][2] = {{8, 8}, {0, -8}};

double vals_fset_var1fs1[2][2] = {{1., 3.}, {9., 27.}};

#define EXCHECK(funcall, errmsg)                                                                   \
  if ((funcall) < 0) {                                                                             \
    fprintf(stderr, errmsg);                                                                       \
    free(varParams.edge_var_tab);                                                                  \
    free(varParams.face_var_tab);                                                                  \
    free(varParams.elem_var_tab);                                                                  \
    free(varParams.fset_var_tab);                                                                  \
    return 1;                                                                                      \
  }

int ex_have_arg(int argc, char *argv[], const char *aname)
{
  int i;
  for (i = 0; i < argc; ++i) {
    if (!strcmp(argv[i], aname)) {
      return 1;
    }
  }
  return 0;
}

int cCreateEdgeFace(int argc, char *argv[])
{
  int exoid;
  int appWordSize  = 8;
  int diskWordSize = 8;
  /*  int concatBlocks = ex_have_arg( argc, argv, "-pcab" ); */
  int    concatSets   = ex_have_arg(argc, argv, "-pcset");
  int    concatResult = ex_have_arg(argc, argv, "-pvpax");
  double t;

  ex_init_params modelParams = {
      "CreateEdgeFace Test", /* title */
      3,                     /* num_dim */
      12,                    /* num_nodes */
      20,                    /* num_edge */
      1,                     /* num_edge_blk */
      11,                    /* num_face */
      3,                     /* num_face_blk */
      3,                     /* num_elem */
      2,                     /* num_elem_blk */
      1,                     /* num_node_sets */
      1,                     /* num_edge_sets */
      1,                     /* num_face_sets */
      1,                     /* num_side_sets */
      2,                     /* num_elem_sets */
      1,                     /* num_node_map */
      1,                     /* num_edge_map */
      1,                     /* num_face_map */
      1,                     /* num_elem_map */
  };

  ex_block edgeBlocks[1];
  ex_block faceBlocks[3];
  ex_block elemBlocks[2];

  ex_var_params varParams;

  ex_opts(EX_VERBOSE | EX_ABORT);

  exoid = ex_create(EX_TEST_FILENAME, EX_CLOBBER, &appWordSize, &diskWordSize);
  if (exoid <= 0) {
    fprintf(stderr, "Unable to open \"%s\" for writing.\n", EX_TEST_FILENAME);
    return 1;
  }

  edgeBlocks[0].type                = EX_EDGE_BLOCK;
  edgeBlocks[0].id                  = 100;
  edgeBlocks[0].num_entry           = 20;
  edgeBlocks[0].num_nodes_per_entry = 2;
  edgeBlocks[0].num_attribute       = 1;
  ex_copy_string(edgeBlocks[0].topology, "EDGE2", MAX_STR_LENGTH + 1);

  faceBlocks[0].type                = EX_FACE_BLOCK;
  faceBlocks[0].id                  = 500;
  faceBlocks[0].num_entry           = 2;
  faceBlocks[0].num_nodes_per_entry = 4;
  faceBlocks[0].num_attribute       = 1;
  ex_copy_string(faceBlocks[0].topology, "QUAD4", MAX_STR_LENGTH + 1);

  faceBlocks[1].type                = EX_FACE_BLOCK;
  faceBlocks[1].id                  = 600;
  faceBlocks[1].num_entry           = 1;
  faceBlocks[1].num_nodes_per_entry = 4;
  faceBlocks[1].num_attribute       = 1;
  ex_copy_string(faceBlocks[1].topology, "QUAD4", MAX_STR_LENGTH + 1);

  faceBlocks[2].type                = EX_FACE_BLOCK;
  faceBlocks[2].id                  = 700;
  faceBlocks[2].num_entry           = 8;
  faceBlocks[2].num_nodes_per_entry = 4;
  faceBlocks[2].num_attribute       = 1;
  ex_copy_string(faceBlocks[2].topology, "QUAD4", MAX_STR_LENGTH + 1);

  elemBlocks[0].type                = EX_ELEM_BLOCK;
  elemBlocks[0].id                  = 200;
  elemBlocks[0].num_entry           = 2;
  elemBlocks[0].num_nodes_per_entry = 8;
  elemBlocks[0].num_edges_per_entry = 12;
  elemBlocks[0].num_faces_per_entry = 6;
  elemBlocks[0].num_attribute       = 2;
  ex_copy_string(elemBlocks[0].topology, "HEX8", MAX_STR_LENGTH + 1);

  elemBlocks[1].type                = EX_ELEM_BLOCK;
  elemBlocks[1].id                  = 201;
  elemBlocks[1].num_entry           = 1;
  elemBlocks[1].num_nodes_per_entry = 4;
  elemBlocks[1].num_edges_per_entry = 0;
  elemBlocks[1].num_faces_per_entry = 0;
  elemBlocks[1].num_attribute       = 0;
  ex_copy_string(elemBlocks[1].topology, "TET4", MAX_STR_LENGTH + 1);

  varParams.edge_var_tab  = (int *)malloc(2 * sizeof(int));
  varParams.face_var_tab  = (int *)malloc(3 * sizeof(int));
  varParams.elem_var_tab  = (int *)malloc(2 * sizeof(int));
  varParams.nset_var_tab  = (int *)0;
  varParams.eset_var_tab  = (int *)0;
  varParams.fset_var_tab  = (int *)malloc(1 * sizeof(int));
  varParams.sset_var_tab  = (int *)0;
  varParams.elset_var_tab = (int *)0;

  varParams.num_glob        = 2;
  varParams.num_node        = 1;
  varParams.num_edge        = 2;
  varParams.edge_var_tab[0] = 1;
  varParams.edge_var_tab[1] = 1;
  varParams.num_face        = 1;
  varParams.face_var_tab[0] = 1;
  varParams.face_var_tab[1] = 1;
  varParams.face_var_tab[2] = 1;
  varParams.num_elem        = 1;
  varParams.elem_var_tab[0] = 1;
  varParams.elem_var_tab[1] = 0;
  varParams.num_nset        = 0;
  varParams.num_eset        = 0;
  ;
  varParams.num_fset        = 1;
  varParams.fset_var_tab[0] = 1;
  varParams.num_sset        = 0;
  varParams.num_elset       = 0;

  EXCHECK(ex_put_init_ext(exoid, &modelParams), "Unable to initialize database.\n");

  {
    int blk;
    for (blk = 0; blk < modelParams.num_edge_blk; ++blk) {
      EXCHECK(ex_put_block_param(exoid, edgeBlocks[blk]), "Unable to write edge block");
    }
    for (blk = 0; blk < modelParams.num_face_blk; ++blk) {
      EXCHECK(ex_put_block_param(exoid, faceBlocks[blk]), "Unable to write face block");
    }
    for (blk = 0; blk < modelParams.num_elem_blk; ++blk) {
      EXCHECK(ex_put_block_param(exoid, elemBlocks[blk]), "Unable to write elem block");
    }
  }

  EXCHECK(ex_put_coord(exoid, (void *)coordsX, (void *)coordsY, (void *)coordsZ),
          "Unable to write coordinates.\n");

  EXCHECK(ex_put_coord_names(exoid, (char **)coordsNames), "Unable to write coordinate names.\n");

  /*                  =============== Connectivity  ================== */
  /* *** NEW API *** */
  EXCHECK(ex_put_conn(exoid, EX_EDGE_BLOCK, edgeBlocks[0].id, ebconn1, 0, 0),
          "Unable to write edge block connectivity.\n");

  /* *** NEW API *** */
  EXCHECK(ex_put_conn(exoid, EX_FACE_BLOCK, faceBlocks[0].id, fbconn1, 0, 0),
          "Unable to write face block 1 connectivity.\n");
  EXCHECK(ex_put_conn(exoid, EX_FACE_BLOCK, faceBlocks[1].id, fbconn2, 0, 0),
          "Unable to write face block 2 connectivity.\n");
  EXCHECK(ex_put_conn(exoid, EX_FACE_BLOCK, faceBlocks[2].id, fbconn3, 0, 0),
          "Unable to write face block 3 connectivity.\n");

  /* *** NEW API *** */
  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, elemBlocks[0].id, conn1, econn1, fconn1),
          "Unable to write elem block 1 connectivity.\n");

  /* *** NEW API *** */
  EXCHECK(ex_put_conn(exoid, EX_ELEM_BLOCK, elemBlocks[1].id, conn2, 0, 0),
          "Unable to write elem block 2 connectivity.\n");

  /* *** NEW API *** */
  EXCHECK(ex_put_names(exoid, EX_EDGE_BLOCK, (char **)edblk_names),
          "Unable to write edge block names.\n");
  EXCHECK(ex_put_names(exoid, EX_FACE_BLOCK, (char **)fablk_names),
          "Unable to write face block names.\n");
  EXCHECK(ex_put_names(exoid, EX_ELEM_BLOCK, (char **)eblk_names),
          "Unable to write element block names.\n");

  /*                  =============== Number Maps   ================== */
  /* *** NEW API *** */
  EXCHECK(ex_put_num_map(exoid, EX_NODE_MAP, 300, nmap1), "Unable to write node map.\n");
  EXCHECK(ex_put_num_map(exoid, EX_EDGE_MAP, 800, edmap1), "Unable to write edge map.\n");
  EXCHECK(ex_put_num_map(exoid, EX_FACE_MAP, 900, famap1), "Unable to write face map.\n");
  EXCHECK(ex_put_num_map(exoid, EX_ELEM_MAP, 400, emap1), "Unable to write element map.\n");

  /* *** NEW API *** */
  EXCHECK(ex_put_names(exoid, EX_NODE_MAP, (char **)nmap_names),
          "Unable to write node map names.\n");
  EXCHECK(ex_put_names(exoid, EX_EDGE_MAP, (char **)edmap_names),
          "Unable to write edge map names.\n");
  EXCHECK(ex_put_names(exoid, EX_FACE_MAP, (char **)famap_names),
          "Unable to write face map names.\n");
  EXCHECK(ex_put_names(exoid, EX_ELEM_MAP, (char **)emap_names),
          "Unable to write element map names.\n");

  /*                 =============== Attribute names ================ */
  /* *** NEW API *** */
  EXCHECK(ex_put_attr_names(exoid, EX_EDGE_BLOCK, edgeBlocks[0].id, (char **)edge_attr_names1),
          "Unable to write edge block 1 attribute names.\n");

  /* *** NEW API *** */
  EXCHECK(ex_put_attr_names(exoid, EX_FACE_BLOCK, faceBlocks[0].id, (char **)face_attr_names1),
          "Unable to write face block 1 attribute names.\n");
  EXCHECK(ex_put_attr_names(exoid, EX_FACE_BLOCK, faceBlocks[1].id, (char **)face_attr_names2),
          "Unable to write face block 1 attribute names.\n");
  EXCHECK(ex_put_attr_names(exoid, EX_FACE_BLOCK, faceBlocks[2].id, (char **)face_attr_names3),
          "Unable to write face block 1 attribute names.\n");

  /* *** NEW API *** */
  EXCHECK(ex_put_attr_names(exoid, EX_ELEM_BLOCK, elemBlocks[0].id, (char **)elem_attr_names1),
          "Unable to write elem block 1 attribute names.\n");

  /*                  =============== Attribute values =============== */
  /* *** NEW API *** */
  EXCHECK(ex_put_attr(exoid, EX_EDGE_BLOCK, edgeBlocks[0].id, edge_attr_values1),
          "Unable to write edge block 1 attribute values.\n");

  /* *** NEW API *** */
  EXCHECK(ex_put_attr(exoid, EX_FACE_BLOCK, faceBlocks[0].id, face_attr_values1),
          "Unable to write face block 1 attribute values.\n");
  EXCHECK(ex_put_attr(exoid, EX_FACE_BLOCK, faceBlocks[1].id, face_attr_values2),
          "Unable to write face block 1 attribute values.\n");
  EXCHECK(ex_put_attr(exoid, EX_FACE_BLOCK, faceBlocks[2].id, face_attr_values3),
          "Unable to write face block 1 attribute values.\n");

  /* *** NEW API *** */
  EXCHECK(ex_put_attr(exoid, EX_ELEM_BLOCK, elemBlocks[0].id, elem_attr_values1),
          "Unable to write elem block 1 attribute values.\n");

  /*                  =============== Set parameters ================= */
  /* *** NEW API *** */
  EXCHECK(ex_put_names(exoid, EX_NODE_SET, (char **)nset_names),
          "Unable to write node set names.\n");
  EXCHECK(ex_put_names(exoid, EX_EDGE_SET, (char **)eset_names),
          "Unable to write edge set names.\n");
  EXCHECK(ex_put_names(exoid, EX_FACE_SET, (char **)fset_names),
          "Unable to write face set names.\n");
  EXCHECK(ex_put_names(exoid, EX_SIDE_SET, (char **)sset_names),
          "Unable to write side set names.\n");
  EXCHECK(ex_put_names(exoid, EX_ELEM_SET, (char **)elset_names),
          "Unable to write element set names.\n");

  {
    ex_set allSets[1 + 1 + 1 + 1 + 2];

    ex_set *nodeSets = &allSets[0];
    ex_set *edgeSets = &allSets[1];
    ex_set *faceSets = &allSets[2];
    ex_set *sideSets = &allSets[3];
    ex_set *elemSets = &allSets[4];

    nodeSets[0].type                     = EX_NODE_SET;
    nodeSets[0].id                       = 1000;
    nodeSets[0].num_entry                = 3;
    nodeSets[0].num_distribution_factor  = 0;
    nodeSets[0].entry_list               = nset_nodes;
    nodeSets[0].extra_list               = NULL;
    nodeSets[0].distribution_factor_list = NULL;

    edgeSets[0].type                     = EX_EDGE_SET;
    edgeSets[0].id                       = 1200;
    edgeSets[0].num_entry                = 6;
    edgeSets[0].num_distribution_factor  = 6;
    edgeSets[0].entry_list               = eset_edges;
    edgeSets[0].extra_list               = eset_orient;
    edgeSets[0].distribution_factor_list = eset_df;

    faceSets[0].type                     = EX_FACE_SET;
    faceSets[0].id                       = 1400;
    faceSets[0].num_entry                = 2;
    faceSets[0].num_distribution_factor  = 0;
    faceSets[0].entry_list               = fset_faces;
    faceSets[0].extra_list               = fset_orient;
    faceSets[0].distribution_factor_list = NULL;

    sideSets[0].type                     = EX_SIDE_SET;
    sideSets[0].id                       = 1400;
    sideSets[0].num_entry                = 5;
    sideSets[0].num_distribution_factor  = 0;
    sideSets[0].entry_list               = sset_elems;
    sideSets[0].extra_list               = sset_sides;
    sideSets[0].distribution_factor_list = NULL;

    elemSets[0].type                     = EX_ELEM_SET;
    elemSets[0].id                       = 1800;
    elemSets[0].num_entry                = 1;
    elemSets[0].num_distribution_factor  = 0;
    elemSets[0].entry_list               = &elset_elems[0];
    elemSets[0].extra_list               = NULL;
    elemSets[0].distribution_factor_list = NULL;

    elemSets[1].type                     = EX_ELEM_SET;
    elemSets[1].id                       = 1900;
    elemSets[1].num_entry                = 1;
    elemSets[1].num_distribution_factor  = 0;
    elemSets[1].entry_list               = &elset_elems[1];
    elemSets[1].extra_list               = NULL;
    elemSets[1].distribution_factor_list = NULL;

    if (concatSets) {
      EXCHECK(ex_put_sets(exoid, 1 + 2 + 1 + 1 + 1, allSets),
              "Unable to output concatenated sets.\n");
    }
    else {
      EXCHECK(ex_put_sets(exoid, 1, nodeSets), "Unable to write node sets.\n");
      EXCHECK(ex_put_sets(exoid, 1, edgeSets), "Unable to write edge sets.\n");
      EXCHECK(ex_put_sets(exoid, 1, faceSets), "Unable to write face sets.\n");
      EXCHECK(ex_put_sets(exoid, 1, sideSets), "Unable to write side sets.\n");
      EXCHECK(ex_put_sets(exoid, 2, elemSets), "Unable to write element sets.\n");
    }
  }

  /*                  =============== Result variable params ========= */
  /* *** NEW API *** */
  if (concatResult) {
    EXCHECK(ex_put_all_var_param_ext(exoid, &varParams),
            "Unable to write result variable parameter information.\n");
  }
  else {
    EXCHECK(ex_put_variable_param(exoid, EX_GLOBAL, 2),
            "Unable to write global result variable parameters.\n");
    EXCHECK(ex_put_variable_param(exoid, EX_NODAL, 1),
            "Unable to write nodal result variable parameters.\n");
    EXCHECK(ex_put_variable_param(exoid, EX_ELEM_BLOCK, 1),
            "Unable to write element result variable parameters.\n");
    EXCHECK(ex_put_variable_param(exoid, EX_EDGE_BLOCK, 2),
            "Unable to write edge result variable parameters.\n");
    EXCHECK(ex_put_variable_param(exoid, EX_FACE_BLOCK, 1),
            "Unable to write face result variable parameters.\n");
    EXCHECK(ex_put_variable_param(exoid, EX_FACE_SET, 1),
            "Unable to write faceset result variable parameters.\n");
  }

  /*                  =============== Result variable names ========== */
  /* *** NEW API *** */
  EXCHECK(ex_put_variable_name(exoid, EX_GLOBAL, 1, "CALIBER"), "Unable to write variable name.\n");
  EXCHECK(ex_put_variable_name(exoid, EX_GLOBAL, 2, "GUNPOWDER"),
          "Unable to write variable name.\n");
  EXCHECK(ex_put_variable_name(exoid, EX_NODAL, 1, "RHO"), "Unable to write variable name.\n");
  EXCHECK(ex_put_variable_name(exoid, EX_EDGE_BLOCK, 1, "GAMMA1"),
          "Unable to write variable name.\n");
  EXCHECK(ex_put_variable_name(exoid, EX_EDGE_BLOCK, 2, "GAMMA2"),
          "Unable to write variable name.\n");
  EXCHECK(ex_put_variable_name(exoid, EX_FACE_BLOCK, 1, "PHI"), "Unable to write variable name.\n");
  EXCHECK(ex_put_variable_name(exoid, EX_ELEM_BLOCK, 1, "EPSTRN"),
          "Unable to write variable name.\n");
  EXCHECK(ex_put_variable_name(exoid, EX_FACE_SET, 1, "PHI0"), "Unable to write variable name.\n");

  /*                  =============== Result variable values ========= */
  t = 1.;
  /* *** NEW API *** */
  EXCHECK(ex_put_time(exoid, 1, &t), "Unable to write time value.\n");
  EXCHECK(ex_put_var(exoid, 1, EX_GLOBAL, 1, 0 /*N/A*/, 2, vals_glo_var[0]),
          "Unable to write global var 1.\n");
  EXCHECK(ex_put_var(exoid, 1, EX_EDGE_BLOCK, 1, 100, 20, vals_edge_var1eb1[0]),
          "Unable to write edge block 1 var 1.\n");
  EXCHECK(ex_put_var(exoid, 1, EX_EDGE_BLOCK, 2, 100, 20, vals_edge_var2eb1[0]),
          "Unable to write edge block 1 var 2.\n");
  EXCHECK(ex_put_var(exoid, 1, EX_FACE_BLOCK, 1, 500, 2, vals_face_var1fb1[0]),
          "Unable to write face block 1 var 1.\n");
  EXCHECK(ex_put_var(exoid, 1, EX_FACE_BLOCK, 1, 700, 8, vals_face_var1fb3[0]),
          "Unable to write face block 3 var 1.\n");
  EXCHECK(ex_put_var(exoid, 1, EX_ELEM_BLOCK, 1, 200, 2, vals_elem_var1eb1[0]),
          "Unable to write elem block 1 var 1.\n");
  EXCHECK(ex_put_var(exoid, 1, EX_FACE_SET, 1, 1400, 2, vals_fset_var1fs1[0]),
          "Unable to write face set 1 var 1.\n");

  t = 2.;
  EXCHECK(ex_put_time(exoid, 2, &t), "Unable to write time value.\n");
  EXCHECK(ex_put_var(exoid, 2, EX_GLOBAL, 1, 0 /*N/A*/, 2, vals_glo_var[1]),
          "Unable to write global var 1.\n");
  EXCHECK(ex_put_var(exoid, 2, EX_EDGE_BLOCK, 1, 100, 20, vals_edge_var1eb1[1]),
          "Unable to write edge block 1 var 1.\n");
  EXCHECK(ex_put_var(exoid, 2, EX_EDGE_BLOCK, 2, 100, 20, vals_edge_var2eb1[1]),
          "Unable to write edge block 1 var 2.\n");
  EXCHECK(ex_put_var(exoid, 2, EX_FACE_BLOCK, 1, 500, 2, vals_face_var1fb1[1]),
          "Unable to write face block 1 var 1.\n");
  EXCHECK(ex_put_var(exoid, 2, EX_FACE_BLOCK, 1, 700, 8, vals_face_var1fb3[1]),
          "Unable to write face block 3 var 1.\n");
  EXCHECK(ex_put_var(exoid, 2, EX_ELEM_BLOCK, 1, 200, 2, vals_elem_var1eb1[1]),
          "Unable to write elem block 1 var 1.\n");
  EXCHECK(ex_put_var(exoid, 2, EX_FACE_SET, 1, 1400, 2, vals_fset_var1fs1[1]),
          "Unable to write face set 1 var 1.\n");

  EXCHECK(ex_put_var(exoid, 1, EX_NODAL, 1, 1, 12, vals_nod_var[0]),
          "Unable to write nodal var 1.\n");
  EXCHECK(ex_put_var(exoid, 2, EX_NODAL, 1, 1, 12, vals_nod_var[1]),
          "Unable to write nodal var 1.\n");

  EXCHECK(ex_close(exoid), "Unable to close database.\n");

  free(varParams.edge_var_tab);
  free(varParams.face_var_tab);
  free(varParams.elem_var_tab);
  free(varParams.fset_var_tab);
  return 0;
}

#if !defined(USING_CMAKE)
int main(int argc, char *argv[]) { return cCreateEdgeFace(argc, argv); }
#endif
