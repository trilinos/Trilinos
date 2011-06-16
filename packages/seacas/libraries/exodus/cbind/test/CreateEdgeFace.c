#include "exodusII.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define EX_TEST_FILENAME "edgeFace.exo"

/* ======== Coordinates and connectivity ========= */
double coordsX[] = {
   0.,  0.,  0.,  0., 
   3.,  3.,  3.,  3.,
  -3., -3., -3., -3.
};

double coordsY[] = {
  -1.,  1.,  1., -1., 
  -3.,  3.,  3., -3.,
  -3.,  3.,  3., -3.
};

double coordsZ[] = {
  -1., -1.,  1.,  1., 
  -3., -3.,  3.,  3.,
  -3., -3.,  3.,  3.
};

const char* coordsNames[] = { "X", "Y", "Z" };

int conn1[] = {
   1,  2,  3,  4,  5,  6,  7,  8,
   9, 10, 11, 12,  1,  2,  3,  4
};

int conn2[] = { 1, 2, 3, 5 };

int econn1[] = {
   1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13, 14, 15, 16,  1,  2,  3,  4, 17, 18, 19, 20
};

int fconn1[] = {
   4,  5,  7,  6,  3,  2,
   8,  9, 11, 10,  1,  3
};

int ebconn1[] = {
   1,  2,
   2,  3,
   3,  4,
   4,  1,
   5,  6,
   6,  7,
   7,  8,
   8,  5,
   1,  5,
   2,  6,
   4,  8,
   3,  7,
   9, 10,
  10, 11,
  11, 12,
  12,  9,
   9,  1,
  10,  2,
  12,  4,
  11,  3
};

int fbconn1[] = {
  12, 11, 10,  9,
   5,  6,  7,  8
};

int fbconn2[] = {
   1,  2,  3,  4
};

int fbconn3[] = {
   1,  5,  6,  2,
   3,  7,  8,  4,
   2,  6,  7,  3,
   4,  8,  5,  1,
   9,  1,  2, 10,
  11,  3,  4, 12,
  10,  2,  3, 11,
  12,  4,  1,  9
};

int nmap1[] = {
  12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1
};

int edmap1[] = {
  1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20
};

int famap1[] = {
  2, 1, 4, 3, 6, 5, 7, 8, 9, 10, 11
};

int emap1[] = {
  1, 2
};

const char* eblk_names[] = {
  "Eli WALLACH",
  "Angelo NOVI"
};
const char* edblk_names[] = { "Aldo GIUFFRE" };
const char* fablk_names[] = { 
  "Livio LORENZON",
  "Claudio SCARCHILLI",
  "John BARTHA"
};

const char* nmap_names[] = { "Luigi PISTILLI" };
const char* edmap_names[] = { "Antonio CASALE" };
const char* famap_names[] = { "Sandro SCARCHILLI" };
const char* emap_names[] = { "Benito STEFANELLI" };

/* ======== Sets ========= */
int nset_nodes[] = {
  5, 6, 9
};

int eset_edges[] = {
  1, 2, 4, 15, 19, 20
};

int eset_orient[] = {
  +1, +1, +1, +1, +1, -1
};

double eset_df[] = {
  2., 2., 0.5, 0.5, 1., 1.
};

int fset_faces[] = {
  3, 9
};

int fset_orient[] = {
  +1, -1
};

int sset_elems[] = {
  1, 1, 1, 2, 2
};

int sset_sides[] = {
  1, 3, 5, 2, 4
};

int elset_elems[] = {
  1,
  2
};

const char* elset_names[] = {
  "Clint EASTWOOD",
  "Lee VAN CLEEF"
};

const char* nset_names[] = { "Ennio MORRICONE" };
const char* eset_names[] = { "Rada RASSIMOV" };
const char* fset_names[] = { "Enzo PETITO" };
const char* sset_names[] = { "Luciano VINCENZONI" };

/* ======== Attributes ========= */
const char* edge_attr_names1[] = {"Sergio LEONE"};

const char* face_attr_names1[] = {"GOOD"};
const char* face_attr_names2[] = {"BAD"};
const char* face_attr_names3[] = {"UGLY"};

const char* elem_attr_names1[] = {
  "SPAGHETTI",
  "WESTERN"
};

double edge_attr_values1[] = {
   1.,  2.,  3.,  5.,  7., 11., 13., 17., 19., 23.,
  29., 31., 37., 41., 43., 47., 53., 59., 61., 67.
};

double face_attr_values1[] = {
  71., 73.
};

double face_attr_values2[] = {
  79.
};

double face_attr_values3[] = {
  83., 89., 97., 101., 103., 107., 109., 113.
};

double elem_attr_values1[] = {
  127., 101.,
  137., 139.
};

/* ======== Results variables ========= */
/*           (2 time steps) */

double vals_glo_var[2][2] = {
  { 36., 37. },
  { 42., 43. }
};

double vals_nod_var[2][12] = {
  { 0.1, 0.8, 0.0, 0.4, 0.3, 0.9, 0.8, 0.5, 0.3, 0.7, 0.4, 0.6 },
  { 0.7, 0.5, 0.3, 0.5, 0.2, 0.7, 0.9, 0.8, 0.0, 0.2, 0.3, 0.5 }
} ;


double vals_edge_var1eb1[2][20] = {
  { 20., 19., 18., 17., 16., 15., 14., 13., 12., 11., 10., 9., 8., 7., 6., 5., 4., 3., 2., 1. },
  { 21., 20., 19., 18., 17., 16., 15., 14., 13., 12., 11., 10., 9., 8., 7., 6., 5., 4., 3., 2. }
};

double vals_edge_var2eb1[2][20] = {
  { 1., 1., 0., 0., 1., 1., 2., 0., 2., 0., 1., 1., 1., 1., 0., 0., 2., 2., 2., 2. },
  { 1., 1., 0., 0., 1., 1., 2., 0., 2., 0., 1., 1., 1., 1., 0., 0., 2., 2., 2., 2. }
};

double vals_face_var1fb1[2][2] = {
  { 0, 1 },
  { 2, 0 }
};

double vals_face_var1fb3[2][8] = {
  { 1, 0, 2, 0, 3, 0, 4, 0 },
  { 0, 1, 0, 2, 0, 3, 0, 4 }
};

double vals_elem_var1eb1[2][2] = {
  { 8,  8 },
  { 0, -8 }
};

double vals_fset_var1fs1[2][2] = {
  { 1., 3. },
  { 9., 27. }
};

#define EXCHECK(funcall,errmsg)\
  if ( (funcall) < 0 ) \
    { \
      fprintf( stderr, errmsg ); \
      return 1; \
    }

int ex_have_arg( int argc, char* argv[], const char* aname )
{
  int i;
  for ( i = 0; i < argc; ++i )
    if ( ! strcmp( argv[i], aname ) )
      return 1;
  return 0;
}

int cCreateEdgeFace( int argc, char* argv[] )
{
  int exoid;
  int appWordSize = 8;
  int diskWordSize = 8;
  int concatBlocks = ex_have_arg( argc, argv, "-pcab" );
  int concatSets   = ex_have_arg( argc, argv, "-pcset" );
  int concatResult = ex_have_arg( argc, argv, "-pvpax" );
  double t;

  ex_init_params modelParams = {
    "CreateEdgeFace Test", /* title */
    3,  /* num_dim */
    12, /* num_nodes */
    20, /* num_edge */
    1,  /* num_edge_blk */
    11, /* num_face */
    3,  /* num_face_blk */
    2,  /* num_elem */
    2,  /* num_elem_blk */
    1,  /* num_node_sets */
    1,  /* num_edge_sets */
    1,  /* num_face_sets */
    1,  /* num_side_sets */
    2,  /* num_elem_sets */
    1,  /* num_node_map */
    1,  /* num_edge_map */
    1,  /* num_face_map */
    1,  /* num_elem_map */
  };

  ex_block_params blockParams;
  ex_var_params varParams;

  blockParams.edge_blk_id         = (int*)malloc(1 * sizeof(int));
  blockParams.num_edge_this_blk   = (int*)malloc(1 * sizeof(int));
  blockParams.num_nodes_per_edge  = (int*)malloc(1 * sizeof(int));
  blockParams.num_attr_edge       = (int*)malloc(1 * sizeof(int));
  blockParams.face_blk_id         = (int*)malloc(3 * sizeof(int));
  blockParams.num_face_this_blk   = (int*)malloc(3 * sizeof(int));
  blockParams.num_nodes_per_face  = (int*)malloc(3 * sizeof(int));
  blockParams.num_attr_face       = (int*)malloc(3 * sizeof(int));
  blockParams.elem_blk_id         = (int*)malloc(2 * sizeof(int));
  blockParams.num_elem_this_blk   = (int*)malloc(2 * sizeof(int));
  blockParams.num_nodes_per_elem  = (int*)malloc(2 * sizeof(int));
  blockParams.num_edges_per_elem  = (int*)malloc(2 * sizeof(int));
  blockParams.num_faces_per_elem  = (int*)malloc(2 * sizeof(int));
  blockParams.num_attr_elem       = (int*)malloc(2 * sizeof(int));
  
  blockParams.edge_type    = (char**)malloc(1 * sizeof(char*));
  blockParams.edge_type[0] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.face_type    = (char**)malloc(3 * sizeof(char*));
  blockParams.face_type[0] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.face_type[1] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.face_type[2] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.elem_type    = (char**)malloc(2 * sizeof(char*));
  blockParams.elem_type[0] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.elem_type[1] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));

  blockParams.edge_blk_id[0]         = 100;
  blockParams.num_edge_this_blk[0]   = 20;
  blockParams.num_nodes_per_edge[0]  = 2;
  blockParams.num_attr_edge[0]       = 1;

  blockParams.face_blk_id[0]         = 500;
  blockParams.face_blk_id[1]         = 600;
  blockParams.face_blk_id[2]         = 700;
  blockParams.num_face_this_blk[0]   = 2;
  blockParams.num_face_this_blk[1]   = 1;
  blockParams.num_face_this_blk[2]   = 8;
  blockParams.num_nodes_per_face[0]  = 4;
  blockParams.num_nodes_per_face[1]  = 4;
  blockParams.num_nodes_per_face[2]  = 4;
  blockParams.num_attr_face[0]       = 1;
  blockParams.num_attr_face[1]       = 1;
  blockParams.num_attr_face[2]       = 1;

  blockParams.elem_blk_id[0]         = 200;
  blockParams.elem_blk_id[1]         = 201;
  blockParams.num_elem_this_blk[0]   = 2;
  blockParams.num_elem_this_blk[1]   = 1;
  blockParams.num_nodes_per_elem[0]  = 8;
  blockParams.num_nodes_per_elem[1]  = 4;
  blockParams.num_edges_per_elem[0]  = 12;
  blockParams.num_edges_per_elem[1]  = 0;
  blockParams.num_faces_per_elem[0]  = 6;
  blockParams.num_faces_per_elem[1]  = 0;
  blockParams.num_attr_elem[0]       = 2;
  blockParams.num_attr_elem[1]       = 0;

  blockParams.define_maps = 0;

  strcpy(blockParams.edge_type[0], "STRAIGHT2");

  strcpy(blockParams.face_type[0], "QUAD4");
  strcpy(blockParams.face_type[1], "QUAD4");
  strcpy(blockParams.face_type[2], "QUAD4");

  strcpy(blockParams.elem_type[0], "HEX8");
  strcpy(blockParams.elem_type[1], "TET4");

  varParams.edge_var_tab  = (int*)malloc(2 * sizeof(int));
  varParams.face_var_tab  = (int*)malloc(3 * sizeof(int));
  varParams.elem_var_tab  = (int*)malloc(2 * sizeof(int));
  varParams.nset_var_tab  = (int*)0;
  varParams.eset_var_tab  = (int*)0;
  varParams.fset_var_tab  = (int*)malloc(1 * sizeof(int));
  varParams.sset_var_tab  = (int*)0;
  varParams.elset_var_tab = (int*)0;

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
  varParams.num_eset        = 0;;
  varParams.num_fset        = 1;
  varParams.fset_var_tab[0] = 1;
  varParams.num_sset        = 0;
  varParams.num_elset       = 0;

  exoid = ex_create( EX_TEST_FILENAME, EX_CLOBBER, &appWordSize, &diskWordSize );
  if ( exoid <= 0 )
    {
      fprintf( stderr, "Unable to open \"%s\" for writing.\n", EX_TEST_FILENAME );
      return 1;
    }

  /* *** NEW API *** */
  EXCHECK( ex_put_init_ext( exoid, &modelParams ),
	   "Unable to initialize database.\n" );

  /* *** NEW API *** */
  if ( concatBlocks ) {
    EXCHECK( ex_put_concat_all_blocks( exoid, &blockParams ),
	     "Unable to initialize block params.\n" );
  } else {
    int blk;
    for ( blk = 0; blk < modelParams.num_edge_blk; ++blk ) {
      EXCHECK( ex_put_block( exoid, EX_EDGE_BLOCK, blockParams.edge_blk_id[blk], blockParams.edge_type[blk],
			     blockParams.num_edge_this_blk[blk], blockParams.num_nodes_per_edge[blk], 0, 0,
			     blockParams.num_attr_edge[blk] ), "Unable to write edge block" );
    }
    for ( blk = 0; blk < modelParams.num_face_blk; ++blk ) {
      EXCHECK( ex_put_block( exoid, EX_FACE_BLOCK, blockParams.face_blk_id[blk], blockParams.face_type[blk],
			     blockParams.num_face_this_blk[blk], blockParams.num_nodes_per_face[blk], 0, 0,
			     blockParams.num_attr_face[blk] ), "Unable to write face block" );
    }
    for ( blk = 0; blk < modelParams.num_elem_blk; ++blk ) {
      EXCHECK( ex_put_block( exoid, EX_ELEM_BLOCK, blockParams.elem_blk_id[blk], blockParams.elem_type[blk],
			     blockParams.num_elem_this_blk[blk], blockParams.num_nodes_per_elem[blk],
			     blockParams.num_edges_per_elem[blk], blockParams.num_faces_per_elem[blk],
			     blockParams.num_attr_elem[blk] ), "Unable to write elem block" );
    }
  }

  EXCHECK( ex_put_coord( exoid, (void*)coordsX, (void*)coordsY, (void*)coordsZ ),
	   "Unable to write coordinates.\n" );

  EXCHECK( ex_put_coord_names( exoid, (char**)coordsNames ),
	   "Unable to write coordinate names.\n" );

  /*                  =============== Connectivity  ================== */
  /* *** NEW API *** */
  EXCHECK( ex_put_conn( exoid, EX_EDGE_BLOCK, blockParams.edge_blk_id[0], ebconn1, 0, 0 ),
	   "Unable to write edge block connectivity.\n" );

  /* *** NEW API *** */
  EXCHECK( ex_put_conn( exoid, EX_FACE_BLOCK, blockParams.face_blk_id[0], fbconn1, 0, 0 ),
	   "Unable to write face block 1 connectivity.\n" );
  EXCHECK( ex_put_conn( exoid, EX_FACE_BLOCK, blockParams.face_blk_id[1], fbconn2, 0, 0 ),
	   "Unable to write face block 2 connectivity.\n" );
  EXCHECK( ex_put_conn( exoid, EX_FACE_BLOCK, blockParams.face_blk_id[2], fbconn3, 0, 0 ),
	   "Unable to write face block 3 connectivity.\n" );

  /* *** NEW API *** */
  EXCHECK( ex_put_conn( exoid, EX_ELEM_BLOCK, blockParams.elem_blk_id[0], conn1, econn1, fconn1 ),
	   "Unable to write elem block 1 connectivity.\n" );

  /* *** NEW API *** */
  EXCHECK( ex_put_conn( exoid, EX_ELEM_BLOCK, blockParams.elem_blk_id[1], conn2, 0, 0 ),
	   "Unable to write elem block 2 connectivity.\n" );

  /* *** NEW API *** */
  EXCHECK( ex_put_names( exoid, EX_EDGE_BLOCK, (char**)edblk_names ), "Unable to write edge block names.\n" );
  EXCHECK( ex_put_names( exoid, EX_FACE_BLOCK, (char**)fablk_names ), "Unable to write face block names.\n" );
  EXCHECK( ex_put_names( exoid, EX_ELEM_BLOCK, (char**) eblk_names ), "Unable to write element block names.\n" );

  /*                  =============== Number Maps   ================== */
  /* *** NEW API *** */
  EXCHECK( ex_put_num_map( exoid, EX_NODE_MAP, 300, nmap1 ),  "Unable to write node map.\n" );
  EXCHECK( ex_put_num_map( exoid, EX_EDGE_MAP, 800, edmap1 ), "Unable to write edge map.\n" );
  EXCHECK( ex_put_num_map( exoid, EX_FACE_MAP, 900, famap1 ), "Unable to write face map.\n" );
  EXCHECK( ex_put_num_map( exoid, EX_ELEM_MAP, 400, emap1 ),  "Unable to write element map.\n" );

  /* *** NEW API *** */
  EXCHECK( ex_put_names( exoid, EX_NODE_MAP, (char**) nmap_names ), "Unable to write node map names.\n" );
  EXCHECK( ex_put_names( exoid, EX_EDGE_MAP, (char**)edmap_names ), "Unable to write edge map names.\n" );
  EXCHECK( ex_put_names( exoid, EX_FACE_MAP, (char**)famap_names ), "Unable to write face map names.\n" );
  EXCHECK( ex_put_names( exoid, EX_ELEM_MAP, (char**) emap_names ), "Unable to write element map names.\n" );

  /*                 =============== Attribute names ================ */
  /* *** NEW API *** */
  EXCHECK( ex_put_attr_names( exoid, EX_EDGE_BLOCK, blockParams.edge_blk_id[0], (char**)edge_attr_names1 ),
	   "Unable to write edge block 1 attribute names.\n" );

  /* *** NEW API *** */
  EXCHECK( ex_put_attr_names( exoid, EX_FACE_BLOCK, blockParams.face_blk_id[0], (char**)face_attr_names1 ),
	   "Unable to write face block 1 attribute names.\n" );
  EXCHECK( ex_put_attr_names( exoid, EX_FACE_BLOCK, blockParams.face_blk_id[1], (char**)face_attr_names2 ),
	   "Unable to write face block 1 attribute names.\n" );
  EXCHECK( ex_put_attr_names( exoid, EX_FACE_BLOCK, blockParams.face_blk_id[2], (char**)face_attr_names3 ),
	   "Unable to write face block 1 attribute names.\n" );

  /* *** NEW API *** */
  EXCHECK( ex_put_attr_names( exoid, EX_ELEM_BLOCK, blockParams.elem_blk_id[0], (char**)elem_attr_names1 ),
	   "Unable to write elem block 1 attribute names.\n" );

  /*                  =============== Attribute values =============== */
  /* *** NEW API *** */
  EXCHECK( ex_put_attr( exoid, EX_EDGE_BLOCK, blockParams.edge_blk_id[0], edge_attr_values1 ),
	   "Unable to write edge block 1 attribute values.\n" );

  /* *** NEW API *** */
  EXCHECK( ex_put_attr( exoid, EX_FACE_BLOCK, blockParams.face_blk_id[0], face_attr_values1 ),
	   "Unable to write face block 1 attribute values.\n" );
  EXCHECK( ex_put_attr( exoid, EX_FACE_BLOCK, blockParams.face_blk_id[1], face_attr_values2 ),
	   "Unable to write face block 1 attribute values.\n" );
  EXCHECK( ex_put_attr( exoid, EX_FACE_BLOCK, blockParams.face_blk_id[2], face_attr_values3 ),
	   "Unable to write face block 1 attribute values.\n" );

  /* *** NEW API *** */
  EXCHECK( ex_put_attr( exoid, EX_ELEM_BLOCK, blockParams.elem_blk_id[0], elem_attr_values1 ),
	   "Unable to write elem block 1 attribute values.\n" );

  /*                  =============== Set parameters ================= */
  /* *** NEW API *** */
  EXCHECK( ex_put_names( exoid, EX_NODE_SET,  (char**)nset_names ), "Unable to write node set names.\n" );
  EXCHECK( ex_put_names( exoid, EX_EDGE_SET,  (char**)eset_names ), "Unable to write edge set names.\n" );
  EXCHECK( ex_put_names( exoid, EX_FACE_SET,  (char**)fset_names ), "Unable to write face set names.\n" );
  EXCHECK( ex_put_names( exoid, EX_SIDE_SET,  (char**)sset_names ), "Unable to write side set names.\n" );
  EXCHECK( ex_put_names( exoid, EX_ELEM_SET, (char**)elset_names ), "Unable to write element set names.\n" );

  if ( concatSets ) {
    ex_set_specs setParams;

    setParams.sets_ids            = (int*)malloc(2*sizeof(int));
    setParams.num_entries_per_set = (int*)malloc(2*sizeof(int));
    setParams.num_dist_per_set    = (int*)malloc(2*sizeof(int));
    setParams.sets_entry_index    = (int*)malloc(2*sizeof(int));
    setParams.sets_dist_index     = (int*)malloc(2*sizeof(int));
    setParams.sets_entry_list     = (int*)malloc(6*sizeof(int));
    setParams.sets_extra_list     = (int*)malloc(6*sizeof(int));
    setParams.sets_dist_fact      = (double*)malloc(6*sizeof(double));

    setParams.sets_ids[0]            = 1000;
    setParams.num_entries_per_set[0] = 3;
    setParams.num_dist_per_set[0]    = 0;
    setParams.sets_entry_index[0]    = 0;
    setParams.sets_dist_index[0]     = 0;
    setParams.sets_entry_list[0]     = nset_nodes[0];
    setParams.sets_entry_list[1]     = nset_nodes[1];
    setParams.sets_entry_list[2]     = nset_nodes[2];

    EXCHECK( ex_put_concat_sets( exoid, EX_NODE_SET, &setParams ), "Unable to write node sets.\n" );

    setParams.sets_ids[0]            = 1200;
    setParams.num_entries_per_set[0] = 6;
    setParams.num_dist_per_set[0]    = 6;
    setParams.sets_entry_index[0]    = 0;
    setParams.sets_dist_index[0]     = 0;
    setParams.sets_entry_list[0]     = eset_edges[0];
    setParams.sets_entry_list[1]     = eset_edges[1];
    setParams.sets_entry_list[2]     = eset_edges[2];
    setParams.sets_entry_list[3]     = eset_edges[3];
    setParams.sets_entry_list[4]     = eset_edges[4];
    setParams.sets_entry_list[5]     = eset_edges[5];

    setParams.sets_extra_list[0]     = eset_orient[0];
    setParams.sets_extra_list[1]     = eset_orient[1];
    setParams.sets_extra_list[2]     = eset_orient[2];
    setParams.sets_extra_list[3]     = eset_orient[3];
    setParams.sets_extra_list[4]     = eset_orient[4];
    setParams.sets_extra_list[5]     = eset_orient[5];

    memcpy(setParams.sets_dist_fact, eset_df, sizeof(eset_df)/sizeof(eset_df[0]));

    EXCHECK( ex_put_concat_sets( exoid, EX_EDGE_SET, &setParams ), "Unable to write edge sets.\n" );

    setParams.sets_ids[0]            = 1400;
    setParams.num_entries_per_set[0] = 2;
    setParams.num_dist_per_set[0]    = 0;
    setParams.sets_entry_index[0]    = 0;
    setParams.sets_dist_index[0]     = 0;
    setParams.sets_entry_list[0]     = fset_faces[0];
    setParams.sets_entry_list[1]     = fset_faces[1];

    setParams.sets_extra_list[0]     = fset_orient[0];
    setParams.sets_extra_list[1]     = fset_orient[1];

    EXCHECK( ex_put_concat_sets( exoid, EX_FACE_SET, &setParams ), "Unable to write face sets.\n" );

    setParams.sets_ids[0]            = 1400;
    setParams.num_entries_per_set[0] = 5;
    setParams.num_dist_per_set[0]    = 0;
    setParams.sets_entry_index[0]    = 0;
    setParams.sets_dist_index[0]     = 0;
    memcpy(setParams.sets_entry_list, sset_elems, sizeof(sset_elems)/sizeof(sset_elems[0]));
    memcpy(setParams.sets_extra_list, sset_sides, sizeof(sset_sides)/sizeof(sset_sides[0]));

    EXCHECK( ex_put_concat_sets( exoid, EX_SIDE_SET, &setParams ), "Unable to write side sets.\n" );

    setParams.sets_ids[0]            = 1800;
    setParams.sets_ids[1]            = 1900;
    setParams.num_entries_per_set[0] = 1;
    setParams.num_entries_per_set[1] = 1;
    setParams.num_dist_per_set[0]    = 0;
    setParams.num_dist_per_set[1]    = 0;
    setParams.sets_entry_index[0]    = 0;
    setParams.sets_entry_index[1]    = 1;
    setParams.sets_dist_index[0]     = 0;
    setParams.sets_dist_index[1]     = 0;
    memcpy(setParams.sets_entry_list, elset_elems, sizeof(elset_elems)/sizeof(elset_elems[0]));

    EXCHECK( ex_put_concat_sets( exoid, EX_ELEM_SET, &setParams ), "Unable to write element sets.\n" );

  } else {
    /* *** NEW API *** */
    EXCHECK( ex_put_set_param( exoid, EX_NODE_SET, 1000, 3, 0 ), "Unable to write node set params.\n" );
    EXCHECK( ex_put_set( exoid, EX_NODE_SET, 1000, nset_nodes, 0 ), "Unable to write node set.\n" );

    /* *** NEW API *** */
    EXCHECK( ex_put_set_param( exoid, EX_EDGE_SET, 1200, 6, 6 ), "Unable to write edge set params.\n" );
    EXCHECK( ex_put_set( exoid, EX_EDGE_SET, 1200, eset_edges, eset_orient ), "Unable to write edge set.\n" );
    EXCHECK( ex_put_set_dist_fact( exoid, EX_EDGE_SET, 1200, eset_df ), "Unable to write edge set dist factors.\n" );

    /* *** NEW API *** */
    EXCHECK( ex_put_set_param( exoid, EX_FACE_SET, 1400, 2, 0 ), "Unable to write face set params.\n" );
    EXCHECK( ex_put_set( exoid, EX_FACE_SET, 1400, fset_faces, fset_orient ), "Unable to write face set.\n" );

    /* *** NEW API *** */
    EXCHECK( ex_put_set_param( exoid, EX_SIDE_SET, 1600, 5, 0 ), "Unable to write side set params.\n" );
    EXCHECK( ex_put_set( exoid, EX_SIDE_SET, 1600, sset_elems, sset_sides ), "Unable to write side set.\n" );

    /* *** NEW API *** */
    EXCHECK( ex_put_set_param( exoid, EX_ELEM_SET, 1800, 1, 0 ), "Unable to write element set 1 params.\n" );
    EXCHECK( ex_put_set( exoid, EX_ELEM_SET, 1800, elset_elems + 0, 0 ), "Unable to write element set 1.\n" );
    EXCHECK( ex_put_set_param( exoid, EX_ELEM_SET, 1900, 1, 0 ), "Unable to write element set 2 params.\n" );
    EXCHECK( ex_put_set( exoid, EX_ELEM_SET, 1900, elset_elems + 1, 0 ), "Unable to write element set 2.\n" );
  }

  /*                  =============== Result variable params ========= */
  /* *** NEW API *** */
  if ( concatResult ) {
    EXCHECK( ex_put_all_var_param_ext( exoid, &varParams ),
	     "Unable to write result variable parameter information.\n" );
  } else {
    EXCHECK( ex_put_var_param( exoid, "G", 2 ),
	     "Unable to write global result variable parameters.\n" );
    EXCHECK( ex_put_var_param( exoid, "N", 1 ),
	     "Unable to write nodal result variable parameters.\n" );
    EXCHECK( ex_put_var_param( exoid, "E", 1 ),
	     "Unable to write element result variable parameters.\n" );
    EXCHECK( ex_put_var_param( exoid, "L", 2 ),
	     "Unable to write edge result variable parameters.\n" );
    EXCHECK( ex_put_var_param( exoid, "F", 1 ),
	     "Unable to write face result variable parameters.\n" );
    EXCHECK( ex_put_var_param( exoid, "A", 1 ),
	     "Unable to write faceset result variable parameters.\n" );
  }

  /*                  =============== Result variable names ========== */
  /* *** NEW API *** */
  EXCHECK( ex_put_var_name( exoid, "G", 1, "CALIBER" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_var_name( exoid, "g", 2, "GUNPOWDER" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_var_name( exoid, "N", 1, "RHO" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_var_name( exoid, "l", 1, "GAMMA1" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_var_name( exoid, "L", 2, "GAMMA2" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_var_name( exoid, "f", 1, "PHI" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_var_name( exoid, "E", 1, "EPSTRN" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_var_name( exoid, "A", 1, "PHI0" ), "Unable to write variable name.\n" );

  /*                  =============== Result variable values ========= */
  t = 1.;
  /* *** NEW API *** */
  EXCHECK( ex_put_time( exoid, 1, &t ), "Unable to write time value.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_GLOBAL, 1, 0/*N/A*/, 2,      vals_glo_var[0] ), "Unable to write global var 1.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_EDGE_BLOCK, 1, 100, 20, vals_edge_var1eb1[0] ), "Unable to write edge block 1 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_EDGE_BLOCK, 2, 100, 20, vals_edge_var2eb1[0] ), "Unable to write edge block 1 var 2.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_FACE_BLOCK, 1, 500,  2, vals_face_var1fb1[0] ), "Unable to write face block 1 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_FACE_BLOCK, 1, 700,  8, vals_face_var1fb3[0] ), "Unable to write face block 3 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_ELEM_BLOCK, 1, 200,  2, vals_elem_var1eb1[0] ), "Unable to write elem block 1 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_FACE_SET,  1, 1400,  2, vals_fset_var1fs1[0] ), "Unable to write face set 1 var 1.\n" );

  t = 2.;
  EXCHECK( ex_put_time( exoid, 2, &t ), "Unable to write time value.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_GLOBAL, 1, 0/*N/A*/, 2,      vals_glo_var[1] ), "Unable to write global var 1.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_EDGE_BLOCK, 1, 100, 20, vals_edge_var1eb1[1] ), "Unable to write edge block 1 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_EDGE_BLOCK, 2, 100, 20, vals_edge_var2eb1[1] ), "Unable to write edge block 1 var 2.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_FACE_BLOCK, 1, 500,  2, vals_face_var1fb1[1] ), "Unable to write face block 1 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_FACE_BLOCK, 1, 700,  8, vals_face_var1fb3[1] ), "Unable to write face block 3 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_ELEM_BLOCK, 1, 200,  2, vals_elem_var1eb1[1] ), "Unable to write elem block 1 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_FACE_SET,  1, 1400,  2, vals_fset_var1fs1[1] ), "Unable to write face set 1 var 1.\n" );

  EXCHECK( ex_put_nodal_var( exoid, 1, 1, 12, vals_nod_var[0] ), "Unable to write nodal var 1.\n" );
  EXCHECK( ex_put_nodal_var( exoid, 2, 1, 12, vals_nod_var[1] ), "Unable to write nodal var 1.\n" );

  EXCHECK( ex_close( exoid ),
	   "Unable to close database.\n" );

  return 0;
}

#if !defined(USING_CMAKE)
int main( int argc, char* argv[] )
{
  return cCreateEdgeFace(argc, argv);
}
#endif
