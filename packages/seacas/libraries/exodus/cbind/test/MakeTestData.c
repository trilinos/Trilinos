#include "exodusII.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if 0
int mymode = EX_MAPS_INT64_DB|EX_MAPS_INT64_API|EX_BULK_INT64_DB|EX_BULK_INT64_API|EX_IDS_INT64_API|EX_IDS_INT64_DB;
typedef int64_t INT;
#else
int mymode = 0;
typedef int INT;
#endif

#define EX_TEST_FILENAME "ExodusTestData.e"

/* ================ Coordinate Frames ================ */
int cf_ids[2] = {20, 13};
double pt_coords[9*2] = {1,0,0,  1,0,1,  2,0,0,
			0,0,0,  1,0,0,  0,1,0};
char tags[2]={'r', 'c'};

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

INT conn1[] = {
   1,  2,  3,  4,  5,  6,  7,  8,
   9, 10, 11, 12,  1,  2,  3,  4
};

INT conn2[] = { 1, 2, 3, 5 };

INT conn3[] = { 12, 11, 10, 9 };

INT conn4[] = {1, 3, 5, 7, 9, 11}; /* Sphere */

INT conn5[] = {12, 8,   11, 7,    10, 6,    9, 5}; /* Beam */

INT econn1[] = {
   1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
  13, 14, 15, 16,  1,  2,  3,  4, 17, 18, 19, 20
};

INT fconn1[] = {
   4,  5,  7,  6,  3,  2,
   8,  9, 11, 10,  1,  3
};

INT ebconn1[] = {
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

INT fbconn1[] = {
  12, 11, 10,  9,
   5,  6,  7,  8
};

INT fbconn2[] = {
   1,  2,  3,  4
};

INT fbconn3[] = {
   1,  5,  6,  2,
   3,  7,  8,  4,
   2,  6,  7,  3,
   4,  8,  5,  1, 9,  1,  2, 10,
  11,  3,  4, 12,
  10,  2,  3, 11,
  12,  4,  1,  9
};

INT nmap1[] = {
  12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1
};

INT nmap2[] = {
  120, 110, 100, 90, 80, 70, 60, 50, 40, 30, 20, 10
};

INT edmap1[] = {
  1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20
};


INT famap1[] = {
  2, 1, 4, 3, 6, 5, 7, 8, 9, 10, 11
};

INT emap1[] = {
  10, 20, 30, 44, 91,92,93,94,95,96,  88, 87, 86, 85
};

const char* eblk_names[] = {
  "Eli_WALLACH",
  "Angelo_NOVI",
  "A_Shell",
  "Some_Sphere_Elements",
  "Reinforcement"
};

const char* edblk_names[] = { "Aldo_GIUFFRE" };
const char* fablk_names[] = { 
  "Livio_LORENZON",
  "Claudio_SCARCHILLI",
  "John_BARTHA"
};

const char* nmap_names[] = { "Luigi_PISTILLI" };
const char* edmap_names[] = { "Antonio_CASALE" };
const char* famap_names[] = { "Sandro_SCARCHILLI" };
const char* emap_names[] = { "Benito_STEFANELLI" };

/* ======== Sets ========= */
INT nset_nodes[] = {
  5, 6, 9
};

INT eset_edges[] = {
  1, 2, 4, 15, 19, 20
};

INT eset_orient[] = {
  +1, +1, +1, +1, +1, -1
};

double eset_df[] = {
  2., 2., 0.5, 0.5, 1., 1.
};

INT fset_faces[] = {
  3, 9
};

INT fset_orient[] = {
  +1, -1
};

INT sset_elems[] = {
  1, 1, 1, 2, 2
};

INT sset_sides[] = {
  1, 3, 5, 2, 4
};

INT sset1_elems[] = {
  4, 4, 4
};

INT sset1_sides[] = {
  1, 3, 5
};

INT elset_elems[] = {
  1,
  2
};

const char* elset_names[] = {
  "Clint_EASTWOOD",
  "Lee_VAN_CLEEF"
};

const char* nset_names[] = { "Ennio_MORRICONE" };
const char* eset_names[] = { "Rada_RASSIMOV" };
const char* fset_names[] = { "Enzo_PETITO" };
const char* sset_names[] = { "Luciano_VINCENZONI", "A_Very_Long_Name_That_Tests_Capability_for_Long_Names" };

/* ======== Attributes ========= */
const char* node_attr_names[]  = {"Influence_Diameter"};
const char* edge_attr_names1[] = {"Sergio_LEONE"};

const char* face_attr_names1[] = {"GOOD"};
const char* face_attr_names2[] = {"BAD"};
const char* face_attr_names3[] = {"UGLY"};

const char* elem_attr_names1[] = {
  "SPAGHETTI",
  "WESTERN"
};

const char* elem_attr_names3[] = {
  "Thickness"
};

const char* elem_attr_names4[] = {
  "Radius", "Volume"
};

const char* elem_attr_names5[] = {
  "Area", "I_1", "I_2", "J", "V_x", "V_y", "V_z"
};

const char* elem_var_names[] = {"my_stress_xx_1",
				"my_stress_yy_1",
				"my_stress_zz_1",
				"my_stress_xy_1",
				"my_stress_yz_1",
				"my_stress_zx_1",
				"my_stress_xx_2",
				"my_stress_yy_2",
				"my_stress_zz_2",
				"my_stress_xy_2",
				"my_stress_yz_2",
				"my_stress_zx_2",
                                "von_mises_which_is_calculated_the_standard_way_from_stress",
				"tension"
};
  
const char* nset_var_names[] = {"var_name.xx",
				"var_name.yy",
				"var_name.zz",
				"var_name.xy",
				"var_name.yz",
				"var_name.zx"};

const char* sset_var_names[] = {"stressxx",
				"stressyy",
				"stresszz",
				"stressxy",
				"stressyz",
				"stresszx"};

double edge_attr_values1[] = {
   1.,  2.,  3.,  5.,  7., 11., 13., 17., 19., 23.,
  29., 31., 37., 41., 43., 47., 53., 59., 61., 67.
};

double node_attr_values[] = { 12., 11., 10., 9., 8., 7., 6., 5., 4., 3., 2., 1.};

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

double elem_attr_values4[] = {
  .10, 0.0,
  .11, 0.0,
  .12, 0.0,
  .13, 0.0,
  .14, 0.0,
  .15, 0.0
};

double elem_attr_values5[] = { /* 4 elements, 7 attributes/element */
  1.0,   10.0, 11.0, 12.0,   0.0, 0.0, 1.0,
  1.1,   10.1, 11.1, 12.1,   1.0, 0.0, 0.0,
  1.2,   10.2, 11.2, 12.2,   0.0, 1.0, 0.0,
  1.3,   10.3, 11.3, 12.3,   1.0, 1.0, 1.0
};

/* ======== Results variables ========= */
/*           (2 time steps) */

double vals_glo_var[2][3] = {
  { 36., 37., 38.},
  { 42., 43., 44.}
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

double vals_elem_var[2][2*6*2] = {
  { 8.0,  8.0, 8.1, 8.1, 8.2, 8.2, 8.3, 8.3, 8.4, 8.4, 8.5, 8.5, 18.0,  18.0, 18.1, 18.1, 18.2, 18.2, 18.3, 18.3, 18.4, 18.4, 18.5, 18.5 },
  { -7.0,  -7.0, -7.1, -7.1, -7.2, -7.2, -7.3, -7.3, -7.4, -7.4, -7.5, -7.5, -7.0,  -17.0, -17.1, -17.1, -17.2, -17.2, -17.3, -17.3, -17.4, -17.4, -17.5, -17.5 }
};

double vals_elem_var1[2][1] = {
  { 88.8 },
  { 99.9 }
};

double vals_tension[2][4] = {
  { 1000., 2000., 3000., 4000. },
  { 2000., 4000., 6000., 8000. }
};

/* 2 time planes, 6 elements, 3 variables */
double vals_elem_var2[2][6*3] = {
  { 1.01, 1.02, 1.03, 2.01, 2.02, 2.03, 3.01, 3.02, 3.03, 4.01, 4.02, 4.03, 5.01, 5.02, 5.03, 6.01, 6.02, 6.03 },
  {11.01,11.02,11.03,12.01,12.02,12.03,13.01,13.02,13.03,14.01,14.02,14.03,15.01,15.02,15.03,16.01,16.02,16.03 }
};

double vals_nset_var[2][3*6] = {
  { 8.0, 8.1, 8.2,
    7.0, 7.1, 7.2,
    6.0, 6.1, 6.2,
    5.0, 5.1, 5.2,
    4.0, 4.1, 4.2,
    3.0, 3.1, 3.2},
  { -8.0, -8.1, -8.2,
    -7.0, -7.1, -7.2,
    -6.0, -6.1, -6.2,
    -5.0, -5.1, -5.2,
    -4.0, -4.1, -4.2,
    -3.0, -3.1, -3.2},
};

double vals_sset_var[2][5*6] = {
  { 18.0, 18.1, 18.2, 18.3, 18.4,
    17.0, 17.1, 17.2, 17.3, 17.4,
    16.0, 16.1, 16.2, 16.3, 16.4,
    15.0, 15.1, 15.2, 15.3, 15.4,
    14.0, 14.1, 14.2, 14.3, 14.4,
    13.0, 13.1, 13.2, 13.3, 13.4},
  { -18.0, -18.1, -18.2, -18.3, -18.4,
    -17.0, -17.1, -17.2, -17.3, -17.4,
    -16.0, -16.1, -16.2, -16.3, -16.4,
    -15.0, -15.1, -15.2, -15.3, -15.4,
    -14.0, -14.1, -14.2, -14.3, -14.4,
    -13.0, -13.1, -13.2, -13.3, -13.4},
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
  int i;

  ex_opts(EX_VERBOSE);
  
  ex_init_params modelParams = {
    "CreateEdgeFace Test", /* title */
    3,  /* num_dim */
    12, /* num_nodes */
    20, /* num_edge */
    1,  /* num_edge_blk */
    11, /* num_face */
    3,  /* num_face_blk */
    14,  /* num_elem */
    5,  /* num_elem_blk */
    1,  /* num_node_sets */
    1,  /* num_edge_sets */
    1,  /* num_face_sets */
    2,  /* num_side_sets */
    2,  /* num_elem_sets */
    1,  /* num_node_map */
    1,  /* num_edge_map */
    1,  /* num_face_map */
    1,  /* num_elem_map */
  };

  ex_block_params blockParams;
  ex_var_params varParams;

  blockParams.edge_blk_id         = (int*)malloc(1 * sizeof(INT));
  blockParams.num_edge_this_blk   = (int*)malloc(1 * sizeof(int));
  blockParams.num_nodes_per_edge  = (int*)malloc(1 * sizeof(int));
  blockParams.num_attr_edge       = (int*)malloc(1 * sizeof(int));
  blockParams.face_blk_id         = (int*)malloc(3 * sizeof(INT));
  blockParams.num_face_this_blk   = (int*)malloc(3 * sizeof(int));
  blockParams.num_nodes_per_face  = (int*)malloc(3 * sizeof(int));
  blockParams.num_attr_face       = (int*)malloc(3 * sizeof(int));
  blockParams.elem_blk_id         = (int*)malloc(5 * sizeof(INT));
  blockParams.num_elem_this_blk   = (int*)malloc(5 * sizeof(int));
  blockParams.num_nodes_per_elem  = (int*)malloc(5 * sizeof(int));
  blockParams.num_edges_per_elem  = (int*)malloc(5 * sizeof(int));
  blockParams.num_faces_per_elem  = (int*)malloc(5 * sizeof(int));
  blockParams.num_attr_elem       = (int*)malloc(5 * sizeof(int));
  
  blockParams.edge_type    = (char**)malloc(1 * sizeof(char*));
  blockParams.edge_type[0] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.face_type    = (char**)malloc(3 * sizeof(char*));
  blockParams.face_type[0] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.face_type[1] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.face_type[2] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.elem_type    = (char**)malloc(5 * sizeof(char*));
  blockParams.elem_type[0] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.elem_type[1] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.elem_type[2] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.elem_type[3] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));
  blockParams.elem_type[4] = (char*)malloc((MAX_STR_LENGTH+1) * sizeof(char));

  ((INT*)blockParams.edge_blk_id)[0]         = 100;
  blockParams.num_edge_this_blk[0]   = 20;
  blockParams.num_nodes_per_edge[0]  = 2;
  blockParams.num_attr_edge[0]       = 1;

  ((INT*)blockParams.face_blk_id)[0]         = 500;
  ((INT*)blockParams.face_blk_id)[1]         = 600;
  ((INT*)blockParams.face_blk_id)[2]         = 700;
  blockParams.num_face_this_blk[0]   = 2;
  blockParams.num_face_this_blk[1]   = 1;
  blockParams.num_face_this_blk[2]   = 8;
  blockParams.num_nodes_per_face[0]  = 4;
  blockParams.num_nodes_per_face[1]  = 4;
  blockParams.num_nodes_per_face[2]  = 4;
  blockParams.num_attr_face[0]       = 1;
  blockParams.num_attr_face[1]       = 1;
  blockParams.num_attr_face[2]       = 1;

  ((INT*)blockParams.elem_blk_id)[0]         = 200;
  ((INT*)blockParams.elem_blk_id)[1]         = 201;
  ((INT*)blockParams.elem_blk_id)[2]         = 100;
  ((INT*)blockParams.elem_blk_id)[3]         = 500;
  ((INT*)blockParams.elem_blk_id)[4]         = 2147483647;

  blockParams.num_elem_this_blk[0]   = 2;
  blockParams.num_elem_this_blk[1]   = 1;
  blockParams.num_elem_this_blk[2]   = 1;
  blockParams.num_elem_this_blk[3]   = 6;
  blockParams.num_elem_this_blk[4]   = 4;

  blockParams.num_nodes_per_elem[0]  = 8;
  blockParams.num_nodes_per_elem[1]  = 4;
  blockParams.num_nodes_per_elem[2]  = 4;
  blockParams.num_nodes_per_elem[3]  = 1;
  blockParams.num_nodes_per_elem[4]  = 2;

  blockParams.num_edges_per_elem[0]  = 12;
  blockParams.num_edges_per_elem[1]  = 0;
  blockParams.num_edges_per_elem[2]  = 0;
  blockParams.num_edges_per_elem[3]  = 0;
  blockParams.num_edges_per_elem[4]  = 0;

  blockParams.num_faces_per_elem[0]  = 6;
  blockParams.num_faces_per_elem[1]  = 0;
  blockParams.num_faces_per_elem[2]  = 0;
  blockParams.num_faces_per_elem[3]  = 0;
  blockParams.num_faces_per_elem[4]  = 0;

  blockParams.num_attr_elem[0]       = 2;
  blockParams.num_attr_elem[1]       = 0;
  blockParams.num_attr_elem[2]       = 1;
  blockParams.num_attr_elem[3]       = 2;
  blockParams.num_attr_elem[4]       = 7;

  blockParams.define_maps = 0;

  strcpy(blockParams.edge_type[0], "EDGE2");

  strcpy(blockParams.face_type[0], "QUAD4");
  strcpy(blockParams.face_type[1], "QUAD4");
  strcpy(blockParams.face_type[2], "QUAD4");

  strcpy(blockParams.elem_type[0], "HEX8");
  strcpy(blockParams.elem_type[1], "TET4");
  strcpy(blockParams.elem_type[2], "SHELL");
  strcpy(blockParams.elem_type[3], "SPHERE");
  strcpy(blockParams.elem_type[4], "beam");

  varParams.edge_var_tab  = (int*)malloc(2 * sizeof(int));
  varParams.face_var_tab  = (int*)malloc(3 * sizeof(int));
  varParams.elem_var_tab  = (int*)0;
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
  varParams.num_nset        = 0;
  varParams.num_eset        = 0;;
  varParams.num_fset        = 1;
  varParams.fset_var_tab[0] = 1;
  varParams.num_sset        = 0;
  varParams.num_elset       = 0;

  exoid = ex_create( EX_TEST_FILENAME, EX_CLOBBER|mymode, &appWordSize, &diskWordSize );
  if ( exoid <= 0 )
    {
      fprintf( stderr, "Unable to open \"%s\" for writing.\n", EX_TEST_FILENAME );
      return 1;
    }

  ex_set_max_name_length(exoid, 80);
  
  EXCHECK( ex_put_init_ext( exoid, &modelParams ),
	   "Unable to initialize database.\n" );

  /* Add a coordinate frame just to give test coverage... */
  EXCHECK( ex_put_coordinate_frames(exoid, 2, cf_ids, pt_coords, tags),
	   "Unable to output coordinate frame.\n");

  if ( concatBlocks ) {
    EXCHECK( ex_put_concat_all_blocks( exoid, &blockParams ),
	     "Unable to initialize block params.\n" );
  } else {
    int blk;
    for ( blk = 0; blk < modelParams.num_edge_blk; ++blk ) {
      EXCHECK( ex_put_block( exoid, EX_EDGE_BLOCK, ((INT*)blockParams.edge_blk_id)[blk], blockParams.edge_type[blk],
			     blockParams.num_edge_this_blk[blk], blockParams.num_nodes_per_edge[blk], 0, 0,
			     blockParams.num_attr_edge[blk] ), "Unable to write edge block" );
    }
    for ( blk = 0; blk < modelParams.num_face_blk; ++blk ) {
      EXCHECK( ex_put_block( exoid, EX_FACE_BLOCK, ((INT*)blockParams.face_blk_id)[blk], blockParams.face_type[blk],
			     blockParams.num_face_this_blk[blk], blockParams.num_nodes_per_face[blk], 0, 0,
			     blockParams.num_attr_face[blk] ), "Unable to write face block" );
    }
    for ( blk = 0; blk < modelParams.num_elem_blk; ++blk ) {
      EXCHECK( ex_put_block( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[blk], blockParams.elem_type[blk],
			     blockParams.num_elem_this_blk[blk], blockParams.num_nodes_per_elem[blk],
			     blockParams.num_edges_per_elem[blk], blockParams.num_faces_per_elem[blk],
			     blockParams.num_attr_elem[blk] ), "Unable to write elem block" );
    }
  }

  EXCHECK( ex_put_attr_param(exoid, EX_NODAL, 0, 1),
	   "Unable to put nodal attributes.\n" );
  
  EXCHECK( ex_put_coord( exoid, (void*)coordsX, (void*)coordsY, (void*)coordsZ ),
	   "Unable to write coordinates.\n" );

  EXCHECK( ex_put_coord_names( exoid, (char**)coordsNames ),
	   "Unable to write coordinate names.\n" );

  /*                  =============== Connectivity  ================== */
  EXCHECK( ex_put_conn( exoid, EX_EDGE_BLOCK, ((INT*)blockParams.edge_blk_id)[0], ebconn1, 0, 0 ),
	   "Unable to write edge block connectivity.\n" );

  EXCHECK( ex_put_conn( exoid, EX_FACE_BLOCK, ((INT*)blockParams.face_blk_id)[0], fbconn1, 0, 0 ),
	   "Unable to write face block 1 connectivity.\n" );
  EXCHECK( ex_put_conn( exoid, EX_FACE_BLOCK, ((INT*)blockParams.face_blk_id)[1], fbconn2, 0, 0 ),
	   "Unable to write face block 2 connectivity.\n" );
  EXCHECK( ex_put_conn( exoid, EX_FACE_BLOCK, ((INT*)blockParams.face_blk_id)[2], fbconn3, 0, 0 ),
	   "Unable to write face block 3 connectivity.\n" );

  EXCHECK( ex_put_conn( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[0], conn1, econn1, fconn1 ),
	   "Unable to write elem block 1 connectivity.\n" );
  EXCHECK( ex_put_conn( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[1], conn2, 0, 0 ),
	   "Unable to write elem block 2 connectivity.\n" );
  EXCHECK( ex_put_conn( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[2], conn3, 0, 0 ),
	   "Unable to write elem block 3 connectivity.\n" );
  EXCHECK( ex_put_conn( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[3], conn4, 0, 0 ),
	   "Unable to write elem block 4 connectivity.\n" );
  EXCHECK( ex_put_conn( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[4], conn5, 0, 0 ),
	   "Unable to write elem block 5 connectivity.\n" );

  EXCHECK( ex_put_names( exoid, EX_EDGE_BLOCK, (char**)edblk_names ), "Unable to write edge block names.\n" );
  EXCHECK( ex_put_names( exoid, EX_FACE_BLOCK, (char**)fablk_names ), "Unable to write face block names.\n" );
  EXCHECK( ex_put_names( exoid, EX_ELEM_BLOCK, (char**) eblk_names ), "Unable to write element block names.\n" );

  /*                  =============== Number Maps   ================== */
  EXCHECK( ex_put_num_map( exoid, EX_NODE_MAP, 300, nmap1 ),  "Unable to write node map.\n" );
  EXCHECK( ex_put_num_map( exoid, EX_EDGE_MAP, 800, edmap1 ), "Unable to write edge map.\n" );
  EXCHECK( ex_put_num_map( exoid, EX_FACE_MAP, 900, famap1 ), "Unable to write face map.\n" );
  EXCHECK( ex_put_num_map( exoid, EX_ELEM_MAP, 400, emap1 ),  "Unable to write element map.\n" );

  EXCHECK( ex_put_names( exoid, EX_NODE_MAP, (char**) nmap_names ), "Unable to write node map names.\n" );
  EXCHECK( ex_put_names( exoid, EX_EDGE_MAP, (char**)edmap_names ), "Unable to write edge map names.\n" );
  EXCHECK( ex_put_names( exoid, EX_FACE_MAP, (char**)famap_names ), "Unable to write face map names.\n" );
  EXCHECK( ex_put_names( exoid, EX_ELEM_MAP, (char**) emap_names ), "Unable to write element map names.\n" );

  /*                  =============== Id Maps   ================== */
  EXCHECK( ex_put_id_map( exoid, EX_NODE_MAP, nmap2 ),  "Unable to write node id map.\n" );
  EXCHECK( ex_put_id_map( exoid, EX_EDGE_MAP, edmap1 ), "Unable to write edge id map.\n" );
  EXCHECK( ex_put_id_map( exoid, EX_FACE_MAP, famap1 ), "Unable to write face id map.\n" );
  EXCHECK( ex_put_id_map( exoid, EX_ELEM_MAP, emap1 ),  "Unable to write element id map.\n" );

  /*                 =============== Attribute names ================ */
  EXCHECK( ex_put_attr_names( exoid, EX_EDGE_BLOCK, ((INT*)blockParams.edge_blk_id)[0], (char**)edge_attr_names1 ),
	   "Unable to write edge block 1 attribute names.\n" );

  EXCHECK( ex_put_attr_names( exoid, EX_NODAL, 0, (char**)node_attr_names),
	   "Unable to write nodal attribute names.\n" );

  EXCHECK( ex_put_attr_names( exoid, EX_FACE_BLOCK, ((INT*)blockParams.face_blk_id)[0], (char**)face_attr_names1 ),
	   "Unable to write face block 1 attribute names.\n" );
  EXCHECK( ex_put_attr_names( exoid, EX_FACE_BLOCK, ((INT*)blockParams.face_blk_id)[1], (char**)face_attr_names2 ),
	   "Unable to write face block 1 attribute names.\n" );
  EXCHECK( ex_put_attr_names( exoid, EX_FACE_BLOCK, ((INT*)blockParams.face_blk_id)[2], (char**)face_attr_names3 ),
	   "Unable to write face block 1 attribute names.\n" );

  EXCHECK( ex_put_attr_names( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[0], (char**)elem_attr_names1 ),
	   "Unable to write elem block 1 attribute names.\n" );
  EXCHECK( ex_put_attr_names( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[2], (char**)elem_attr_names3 ),
	   "Unable to write elem block 3 attribute names.\n" );
  EXCHECK( ex_put_attr_names( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[3], (char**)elem_attr_names4 ),
	   "Unable to write elem block 4 attribute names.\n" );
  EXCHECK( ex_put_attr_names( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[4], (char**)elem_attr_names5 ),
	   "Unable to write elem block 5 attribute names.\n" );

  /*                  =============== Attribute values =============== */
  EXCHECK( ex_put_attr( exoid, EX_EDGE_BLOCK, ((INT*)blockParams.edge_blk_id)[0], edge_attr_values1 ),
	   "Unable to write edge block 1 attribute values.\n" );
  EXCHECK( ex_put_attr( exoid, EX_NODAL, 0, node_attr_values ),
	   "Unable to write node attribute values.\n" );

  EXCHECK( ex_put_attr( exoid, EX_FACE_BLOCK, ((INT*)blockParams.face_blk_id)[0], face_attr_values1 ),
	   "Unable to write face block 1 attribute values.\n" );
  EXCHECK( ex_put_attr( exoid, EX_FACE_BLOCK, ((INT*)blockParams.face_blk_id)[1], face_attr_values2 ),
	   "Unable to write face block 1 attribute values.\n" );
  EXCHECK( ex_put_attr( exoid, EX_FACE_BLOCK, ((INT*)blockParams.face_blk_id)[2], face_attr_values3 ),
	   "Unable to write face block 1 attribute values.\n" );

  EXCHECK( ex_put_attr( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[0], elem_attr_values1 ),
	   "Unable to write elem block 1 attribute values.\n" );
  EXCHECK( ex_put_attr( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[2], elem_attr_values1 ),
	   "Unable to write elem block 3 attribute values.\n" );

  for (i=0; i < 6; i++) {
    elem_attr_values4[2*i+1] = 4.0 / 3.0 * 3.14 * elem_attr_values4[2*i] * elem_attr_values4[2*i] * elem_attr_values4[2*i];
  }
  EXCHECK( ex_put_attr( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[3], elem_attr_values4 ),
	   "Unable to write elem block 3 attribute values.\n" );
  EXCHECK( ex_put_attr( exoid, EX_ELEM_BLOCK, ((INT*)blockParams.elem_blk_id)[4], elem_attr_values5 ),
	   "Unable to write elem block 4 attribute values.\n" );

  /*                  =============== Set parameters ================= */
  /* *** NEW API *** */
  EXCHECK( ex_put_names( exoid, EX_NODE_SET,  (char**)nset_names ), "Unable to write node set names.\n" );
  EXCHECK( ex_put_names( exoid, EX_EDGE_SET,  (char**)eset_names ), "Unable to write edge set names.\n" );
  EXCHECK( ex_put_names( exoid, EX_FACE_SET,  (char**)fset_names ), "Unable to write face set names.\n" );
  EXCHECK( ex_put_names( exoid, EX_SIDE_SET,  (char**)sset_names ), "Unable to write side set names.\n" );
  EXCHECK( ex_put_names( exoid, EX_ELEM_SET, (char**)elset_names ), "Unable to write element set names.\n" );

  if ( concatSets ) {
    ex_set_specs setParams;

    setParams.sets_ids            = (int*)malloc(2*sizeof(INT));
    setParams.num_entries_per_set = (int*)malloc(2*sizeof(int));
    setParams.num_dist_per_set    = (int*)malloc(2*sizeof(int));
    setParams.sets_entry_index    = (int*)malloc(2*sizeof(int));
    setParams.sets_dist_index     = (int*)malloc(2*sizeof(int));
    setParams.sets_entry_list     = (INT*)malloc(6*sizeof(INT));
    setParams.sets_extra_list     = (INT*)malloc(6*sizeof(INT));
    setParams.sets_dist_fact      = (double*)malloc(6*sizeof(double));

    ((INT*)setParams.sets_ids)[0]            = 1000;
    ((INT*)setParams.num_entries_per_set)[0] = 3;
    ((INT*)setParams.num_dist_per_set)[0]    = 0;
    ((INT*)setParams.sets_entry_index)[0]    = 0;
    ((INT*)setParams.sets_dist_index)[0]     = 0;

    {
      INT* entry_list = setParams.sets_entry_list;
      entry_list[0] = nset_nodes[0];
      entry_list[1] = nset_nodes[1];
      entry_list[2] = nset_nodes[2];
    }

    EXCHECK( ex_put_concat_sets( exoid, EX_NODE_SET, &setParams ), "Unable to write node sets.\n" );

    ((INT*)setParams.sets_ids)[0]            = 1200;
    ((INT*)setParams.num_entries_per_set)[0] = 6;
    ((INT*)setParams.num_dist_per_set)[0]    = 6;
    ((INT*)setParams.sets_entry_index)[0]    = 0;
    ((INT*)setParams.sets_dist_index)[0]     = 0;

    {
      INT* entry_list = setParams.sets_entry_list;
      INT* extra_list = setParams.sets_extra_list;
      
      entry_list[0]     = eset_edges[0];
      entry_list[1]     = eset_edges[1];
      entry_list[2]     = eset_edges[2];
      entry_list[3]     = eset_edges[3];
      entry_list[4]     = eset_edges[4];
      entry_list[5]     = eset_edges[5];

      extra_list[0]     = eset_orient[0];
      extra_list[1]     = eset_orient[1];
      extra_list[2]     = eset_orient[2];
      extra_list[3]     = eset_orient[3];
      extra_list[4]     = eset_orient[4];
      extra_list[5]     = eset_orient[5];
    }

    memcpy(setParams.sets_dist_fact, eset_df, sizeof(eset_df)/sizeof(eset_df[0]));

    EXCHECK( ex_put_concat_sets( exoid, EX_EDGE_SET, &setParams ), "Unable to write edge sets.\n" );

    ((INT*)setParams.sets_ids)[0]            = 1400;
    ((INT*)setParams.num_entries_per_set)[0] = 2;
    ((INT*)setParams.num_dist_per_set)[0]    = 0;
    ((INT*)setParams.sets_entry_index)[0]    = 0;
    ((INT*)setParams.sets_dist_index)[0]     = 0;
    {
      INT *entry_list = setParams.sets_entry_list;
      INT *extra_list = setParams.sets_extra_list;
      
      entry_list[0]     = fset_faces[0];
      entry_list[1]     = fset_faces[1];

      extra_list[0]     = fset_orient[0];
      extra_list[1]     = fset_orient[1];
    }

    EXCHECK( ex_put_concat_sets( exoid, EX_FACE_SET, &setParams ), "Unable to write face sets.\n" );

    ((INT*)setParams.sets_ids)[0]            = 1400;
    ((INT*)setParams.sets_ids)[1]            = 1441;
    ((INT*)setParams.num_entries_per_set)[0] = 5;
    ((INT*)setParams.num_entries_per_set)[1] = 3;
    ((INT*)setParams.num_dist_per_set)[0]    = 0;
    ((INT*)setParams.num_dist_per_set)[1]    = 0;
    ((INT*)setParams.sets_entry_index)[0]    = 0;
    ((INT*)setParams.sets_entry_index)[5]    = 0;
    ((INT*)setParams.sets_dist_index)[0]     = 0;
    memcpy(setParams.sets_entry_list, sset_elems, sizeof(sset_elems)/sizeof(sset_elems[0]));
    memcpy(setParams.sets_extra_list, sset_sides, sizeof(sset_sides)/sizeof(sset_sides[0]));

    EXCHECK( ex_put_concat_sets( exoid, EX_SIDE_SET, &setParams ), "Unable to write side sets.\n" );

    ((INT*)setParams.sets_ids)[0]            = 1800;
    ((INT*)setParams.sets_ids)[1]            = 1900;
    ((INT*)setParams.num_entries_per_set)[0] = 1;
    ((INT*)setParams.num_entries_per_set)[1] = 1;
    ((INT*)setParams.num_dist_per_set)[0]    = 0;
    ((INT*)setParams.num_dist_per_set)[1]    = 0;
    ((INT*)setParams.sets_entry_index)[0]    = 0;
    ((INT*)setParams.sets_entry_index)[1]    = 1;
    ((INT*)setParams.sets_dist_index)[0]     = 0;
    ((INT*)setParams.sets_dist_index)[1]     = 0;
    memcpy(setParams.sets_entry_list, elset_elems, sizeof(elset_elems)/sizeof(elset_elems[0]));

    EXCHECK( ex_put_concat_sets( exoid, EX_ELEM_SET, &setParams ), "Unable to write element sets.\n" );

  } else {
    EXCHECK( ex_put_set_param( exoid, EX_NODE_SET, 1000, 3, 0 ), "Unable to write node set params.\n" );
    EXCHECK( ex_put_set( exoid, EX_NODE_SET, 1000, nset_nodes, 0 ), "Unable to write node set.\n" );

    EXCHECK( ex_put_set_param( exoid, EX_EDGE_SET, 1200, 6, 6 ), "Unable to write edge set params.\n" );
    EXCHECK( ex_put_set( exoid, EX_EDGE_SET, 1200, eset_edges, eset_orient ), "Unable to write edge set.\n" );
    EXCHECK( ex_put_set_dist_fact( exoid, EX_EDGE_SET, 1200, eset_df ), "Unable to write edge set dist factors.\n" );

    EXCHECK( ex_put_set_param( exoid, EX_FACE_SET, 1400, 2, 0 ), "Unable to write face set params.\n" );
    EXCHECK( ex_put_set( exoid, EX_FACE_SET, 1400, fset_faces, fset_orient ), "Unable to write face set.\n" );

    EXCHECK( ex_put_set_param( exoid, EX_SIDE_SET, 1600, 5, 0 ), "Unable to write side set params.\n" );
    EXCHECK( ex_put_set( exoid, EX_SIDE_SET, 1600, sset_elems, sset_sides ), "Unable to write side set.\n" );

    EXCHECK( ex_put_set_param( exoid, EX_SIDE_SET, 1661, 3, 0 ), "Unable to write side set params.\n" );
    EXCHECK( ex_put_set( exoid, EX_SIDE_SET, 1661, sset1_elems, sset1_sides ), "Unable to write side set.\n" );

    EXCHECK( ex_put_set_param( exoid, EX_ELEM_SET, 1800, 1, 0 ), "Unable to write element set 1 params.\n" );
    EXCHECK( ex_put_set( exoid, EX_ELEM_SET, 1800, elset_elems + 0, 0 ), "Unable to write element set 1.\n" );
    EXCHECK( ex_put_set_param( exoid, EX_ELEM_SET, 1900, 1, 0 ), "Unable to write element set 2 params.\n" );
    EXCHECK( ex_put_set( exoid, EX_ELEM_SET, 1900, elset_elems + 1, 0 ), "Unable to write element set 2.\n" );
  }

  /*                  =============== Result variable params ========= */
  if ( concatResult ) {
    EXCHECK( ex_put_all_var_param_ext( exoid, &varParams ),
	     "Unable to write result variable parameter information.\n" );
  } else {
    EXCHECK( ex_put_variable_param( exoid, EX_GLOBAL, 3 ),
	     "Unable to write global result variable parameters.\n" );
    EXCHECK( ex_put_variable_param( exoid, EX_NODAL, 1 ),
	     "Unable to write nodal result variable parameters.\n" );
    EXCHECK( ex_put_variable_param( exoid, EX_ELEM_BLOCK, 14 ),
	     "Unable to write element result variable parameters.\n" );
    EXCHECK( ex_put_variable_param( exoid, EX_EDGE_BLOCK, 2 ),
	     "Unable to write edge result variable parameters.\n" );
    EXCHECK( ex_put_variable_param( exoid, EX_FACE_BLOCK, 1 ),
	     "Unable to write face result variable parameters.\n" );
    EXCHECK( ex_put_variable_param( exoid, EX_FACE_SET, 1 ),
	     "Unable to write faceset result variable parameters.\n" );
    EXCHECK( ex_put_variable_param( exoid, EX_SIDE_SET, 6 ),
	     "Unable to write sideset result variable parameters.\n" );
    EXCHECK( ex_put_variable_param( exoid, EX_NODE_SET, 6 ),
	     "Unable to write nodeset result variable parameters.\n" );
  }

  /*                  =============== Result variable names ========== */
  /* *** NEW API *** */
  EXCHECK( ex_put_variable_name( exoid, EX_GLOBAL,     1, "A_vector_X" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_variable_name( exoid, EX_GLOBAL,     2, "A_vector_Y" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_variable_name( exoid, EX_GLOBAL,     3, "A_vector_Z" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_variable_name( exoid, EX_NODAL,      1, "RHO" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_variable_name( exoid, EX_EDGE_BLOCK, 1, "GAMMA1" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_variable_name( exoid, EX_EDGE_BLOCK, 2, "GAMMA2" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_variable_name( exoid, EX_FACE_BLOCK, 1, "PHI" ), "Unable to write variable name.\n" );
  EXCHECK( ex_put_variable_names( exoid, EX_ELEM_BLOCK, 14, (char**)elem_var_names), "Unable to write variable name.\n" );
  EXCHECK( ex_put_variable_names( exoid, EX_SIDE_SET,   6, (char**)sset_var_names), "Unable to write variable name.\n" );
  EXCHECK( ex_put_variable_names( exoid, EX_NODE_SET,   6, (char**)nset_var_names), "Unable to write variable name.\n" );
  EXCHECK( ex_put_variable_name( exoid, EX_FACE_SET,   1, "PHI0" ), "Unable to write variable name.\n" );

  /*                  =============== Result variable values ========= */
  t = 1.;
  /* *** NEW API *** */
  EXCHECK( ex_put_time( exoid, 1, &t ), "Unable to write time value.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_GLOBAL, 1, 0/*N/A*/, 2,      vals_glo_var[0] ), "Unable to write global var 1.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_EDGE_BLOCK, 1, 100, 20, vals_edge_var1eb1[0] ), "Unable to write edge block 1 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_EDGE_BLOCK, 2, 100, 20, vals_edge_var2eb1[0] ), "Unable to write edge block 1 var 2.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_FACE_BLOCK, 1, 500,  2, vals_face_var1fb1[0] ), "Unable to write face block 1 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_FACE_BLOCK, 1, 700,  8, vals_face_var1fb3[0] ), "Unable to write face block 3 var 1.\n" );

  EXCHECK( ex_put_var( exoid, 1, EX_ELEM_BLOCK, 13,  201,  1, &vals_elem_var1[0] ), "Unable to write elem block 1 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 1, EX_ELEM_BLOCK, 14,  2147483647,  4, &vals_tension[0] ), "Unable to write elem block 1 var 1.\n" );
  for(i=0; i < 6; i++){
    /* There are 2 elements in the block and 12 variables (a composite tensor x 2) */
    EXCHECK( ex_put_var( exoid, 1, EX_ELEM_BLOCK, i+1,  200,  2, &vals_elem_var[0][2*i] ), "Unable to write elem block 1 var 1.\n" );
    EXCHECK( ex_put_var( exoid, 1, EX_ELEM_BLOCK, i+7,  200,  2, &vals_elem_var[0][12+2*i] ), "Unable to write elem block 1 var 1.\n" );
    EXCHECK( ex_put_var( exoid, 1, EX_NODE_SET,   i+1, 1000,  3, &vals_nset_var[0][3*i] ), "Unable to write elem block 1 var 1.\n" );
    EXCHECK( ex_put_var( exoid, 1, EX_SIDE_SET,   i+1, 1600,  5, &vals_sset_var[0][5*i] ), "Unable to write elem block 1 var 1.\n" );
  }
  
  EXCHECK( ex_put_var( exoid, 1, EX_FACE_SET,  1, 1400,  2, vals_fset_var1fs1[0] ), "Unable to write face set 1 var 1.\n" );

  t = 2.;
  EXCHECK( ex_put_time( exoid, 2, &t ), "Unable to write time value.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_GLOBAL, 1, 0/*N/A*/,  2,      vals_glo_var[1] ), "Unable to write global var 1.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_EDGE_BLOCK, 1, 100,  20, vals_edge_var1eb1[1] ), "Unable to write edge block 1 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_EDGE_BLOCK, 2, 100,  20, vals_edge_var2eb1[1] ), "Unable to write edge block 1 var 2.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_FACE_BLOCK, 1, 500,   2, vals_face_var1fb1[1] ), "Unable to write face block 1 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_FACE_BLOCK, 1, 700,   8, vals_face_var1fb3[1] ), "Unable to write face block 3 var 1.\n" );

  EXCHECK( ex_put_var( exoid, 2, EX_ELEM_BLOCK, 13,  201,  1, &vals_elem_var1[1] ), "Unable to write elem block 1 var 1.\n" );
  EXCHECK( ex_put_var( exoid, 2, EX_ELEM_BLOCK, 14,  2147483647,  4, &vals_tension[1] ), "Unable to write elem block 1 var 1.\n" );
  for(i=0; i < 6; i++){
    EXCHECK( ex_put_var( exoid, 2, EX_ELEM_BLOCK, i+1,  200,  2, &vals_elem_var[1][2*i] ), "Unable to write elem block 1 var 1.\n" );
    EXCHECK( ex_put_var( exoid, 2, EX_ELEM_BLOCK, i+7,  200,  2, &vals_elem_var[1][12+2*i] ), "Unable to write elem block 1 var 1.\n" );
    EXCHECK( ex_put_var( exoid, 2, EX_NODE_SET,   i+1, 1000,  3, &vals_nset_var[1][3*i] ), "Unable to write elem block 1 var 1.\n" );
    EXCHECK( ex_put_var( exoid, 2, EX_SIDE_SET,   i+1, 1600,  5, &vals_sset_var[1][5*i] ), "Unable to write elem block 1 var 1.\n" );
  }
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
