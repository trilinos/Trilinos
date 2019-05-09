\page polyhedra Polyhedral Element Support
\section poly Storage of 3D arbitrary polyhedra elements in Exodus.

The 3D polyhedra elements are represented as elements with a variable
number of faces in their connectivity.  The faces can either be
regular faces such as quadrilateral or triangles; or they can be
topologically two-dimensional arbitrary polyhedra themselves.

An arbitrary polyhedra 3D element block will have an element type of
"nfaced" or "NFACED".

The faces that are used in the connectivity of this block should be
defined in one or more face blocks.  If the faces are arbitrary
polyhedra, then they will have a face type of "nsided" or "NSIDED".

An annotated example of defining an arbitrary polyhedral element block
consisting of 3 elements is shown below.

The three elements have the following geometry:

* Element 1: 5 faces.
  * Face 1: triangle with nodes 5, 6, 8
  * Face 2: triangle with nodes 2, 1, 4
  * Face 3: quadrilateral with nodes 6, 2, 4, 8
  * Face 4: quadrilateral with nodes 8, 4, 1, 5
  * Face 5: quadrilateral with nodes 1, 2, 6, 5

* Element 2: 5 faces.
  * Face 6: triangle with nodes 5, 8, 7
  * Face 7: triangle with nodes 1, 3, 4
  * Face 8: quadrilateral with nodes 7, 8, 4, 3
  * Face 9: quadrilateral with nodes 7, 3, 1, 5
  * Face 4: quadrilateral with nodes 8, 4, 1, 5 (shared with element 1)

* Element 3: 7 faces.
  * Face  8: quadrilateral with nodes 7, 8, 4, 3 (shared with element 2)
  * Face 10: pentagonal with nodes 8, 4, 14, 10, 12
  * Face 11: pentagonal with nodes 7, 11, 9, 13, 3
  * Face 12: quadrilateral with nodes 7, 8, 12, 11
  * Face 13: quadrilateral with nodes 11, 12, 10, 9
  * Face 14: quadrilateral with nodes 9, 10, 14, 13
  * Face 15: quadrilateral with nodes 12, 14, 4, 3

The Exodus model is created via the following calls:

* Output the initial information.  Since the model contains faces and
  a face block, the "extended" version of the `ex_put_init_ext()` call must be used:
  ~~~~C
  ex_init_params par;
  strcpy( par.title, "This is the title" );
  par.num_dim = 3;
  par.num_nodes = 14;
  par.num_edge = 0;
  par.num_edge_blk = 0;
  par.num_face = 15;
  par.num_face_blk = 1;
  par.num_elem = 3;
  par.num_elem_blk = 1;
  par.num_node_sets = 0;
  par.num_edge_sets = 0;
  par.num_face_sets = 0;
  par.num_side_sets = 0;
  par.num_elem_sets = 0;
  par.num_node_maps = 0;
  par.num_edge_maps = 0;
  par.num_face_maps = 0;
  par.num_elem_maps = 0;

  ex_put_init_ext (exoid, &par);
  ~~~~

* Coordinate output is normal...

* Define the face block.
  ~~~~C
   block_name = "face_block_1";
   num_face_in_block[0] = 15;
   num_total_nodes_per_blk[0] = 58;
   block_id = 10;

   ex_put_block (exoid, EX_FACE_BLOCK, block_id, "nsided",
                 num_face_in_block[0],
                 num_total_nodes_per_blk[0],
                 0, 0, 0);
   ex_put_name(exoid, EX_FACE_BLOCK, block_id, block_name);
  ~~~~

* Output the face connectivity for "face_block_1".
  The data for the face connectivity is listed above; a portion is shown below...
  ~~~~C
   connect = (int *) calloc(num_total_nodes_per_blk[0], sizeof(int));
   i = 0
   connect[i++] = 5;
   connect[i++] = 6;
   connect[i++] = 8; /* connectivity of face 1 of element 1 */

   connect[i++] = 2;
   connect[i++] = 1;
   connect[i++] = 4; /* face 2 of element 1 */

   connect[i++] = 6;
   connect[i++] = 2;
   connect[i++] = 4;
   connect[i++] = 8; /* face 3 of element 1 */

   connect[i++] = 8;
   connect[i++] = 4;
   connect[i++] = 1;
   connect[i++] = 5; /* face 4 of element 1 */

   connect[i++] = 1;
   connect[i++] = 2;
   connect[i++] = 6;
   connect[i++] = 5; /*  face 5 of element 1 */

   connect[i++] = 5;
   connect[i++] = 8;
   connect[i++] = 7; /* connectivity of face 1 of element 2 */

   ... and so on....
   assert(i == num_total_nodes_per_blk[0]);

   ex_put_conn (exoid, EX_FACE_BLOCK, block_id, connect, NULL, NULL);
  ~~~~

* Output the number of nodes per face count for "face_block_1":
  ~~~~C
   j = 0;
   nnpe[ 1] = 3;   /* Face 1 */
   nnpe[ 2] = 3;
   nnpe[ 3] = 4;
   nnpe[ 4] = 4;
   nnpe[ 5] = 4;
   nnpe[ 6] = 3;
   nnpe[ 7] = 3;
   nnpe[ 8] = 4;
   nnpe[ 9] = 4;
   nnpe[10] = 5;
   nnpe[11] = 5;
   nnpe[12] = 4;
   nnpe[13] = 4;
   nnpe[14] = 4;
   nnpe[15] = 4;

   ex_put_entity_count_per_polyhedra(exoid, EX_FACE_BLOCK, block_id, nnpe);
  ~~~~

* The face block is now fully defined; now define the nfaced element
  block which uses these faces.
  ~~~~C
   block_name = "nfaced_1";

   num_elem_in_block = 3;
   num_total_faces_per_blk = 5 + 5 + 7;
   block_id = 10;

   ex_put_block (exoid, EX_ELEM_BLOCK, block_id, "nfaced",
                 num_elem_in_block,
                 0, /* nodes */
                 0, /* edges  */
                 num_total_faces_per_blk,
                 0); /* attribute count */
   ex_put_name(exoid, EX_ELEM_BLOCK, block_id, block_name);
  ~~~~

   In the `ex_put_block()` function, the element type is "nfaced".  The
   connectivity is defined in terms of the faces, so the node and edge
   arguments are passed zeros.  The nodal connectivity can be defined,
   but it isn't required.  The face connectivity argument for an
   nfaced block is the total number of faces in the connectivity for all
   elements in the nfaced block.

* Write the face connectivity:
  ~~~~C
   /* write element-face connectivity */
   connect = (int *) calloc(num_total_faces_per_blk, sizeof(int));

   i = 0;
   connect[i++] = 1;
   connect[i++] = 2;
   connect[i++] = 3;
   connect[i++] = 4;
   connect[i++] = 5;

   connect[i++] = 4;
   connect[i++] = 6;
   connect[i++] = 7;
   connect[i++] = 8;
   connect[i++] = 9;

   connect[i++] = 8;
   connect[i++] = 10;
   connect[i++] = 11;
   connect[i++] = 12;
   connect[i++] = 13;
   connect[i++] = 14;
   connect[i++] = 15;

   assert(i == num_total_faces_per_blk);
   ex_put_conn (exoid, EX_ELEM_BLOCK, block_id, NULL, NULL, connect);
  ~~~~

* Output the number of faces per element count for "nfaced_1":
  ~~~~C
   nnpe[1] = 5;  /* Number of faces per element 1 */
   nnpe[2] = 5;  /* Number of faces per element 2 */
   nnpe[3] = 7;  /* Number of faces per element 3 */

   ex_put_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK, block_id, nnpe);
  ~~~~

* That's all; the rest of the calls are the same as normal Exodus except:

  * There is a similar `ex_get_entity_count_per_polyhedra()` function for read.

  * The `ex_get_block()` functions return the total number of nodes or
    faces for all faces or element for "nfaced" and "nsided" blocks
    and not the number per element

* An example read/write usage is shown in the
  [testwt-nfaced.c](../test/testwt-nfaced.c) and [testrd-nfaced](../test/testrd-nfaced.c) files.

* These changes are in Exodus version v4.93 and later.
