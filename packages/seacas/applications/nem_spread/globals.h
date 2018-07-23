#ifndef GLOBALS_H
#define GLOBALS_H
/*
 * Copyright (C) 2009-2017 National Technology & Engineering Solutions of
 * Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

#include <cstdlib>
#include <rf_allo.h>

/*---------------------------------------------------------------------------*/
/*      STRUCTURES FOR COMMUNICATION MAPS                                    */
/*---------------------------------------------------------------------------*/
template <typename INT> struct ELEM_COMM_MAP
{
  size_t map_id;
  size_t elem_cnt;
  INT *  elem_ids;
  INT *  side_ids;
  INT *  proc_ids;
};

template <typename INT> struct NODE_COMM_MAP
{
  size_t map_id;
  size_t node_cnt;
  INT *  node_ids;
  INT *  proc_ids;
};

/*---------------------------------------------------------------------------*/
/*	GLOBAL QUANTITITES THAT ARE THE SAME ON ALL PROCESSORS		     */
/*---------------------------------------------------------------------------*/
template <typename T, typename INT> class Globals
{
public:
  int    Num_Dim;    /* Number of physical dimensions in the problem*/
  size_t Num_Node;   /* Total number of nodes in the entire mesh    *
                      * - this is a global quantity                 */
  size_t Num_Elem;   /* Total number of elements in the entire mesh *
                      * - this is a global quantity                 */
  int Num_Elem_Blk;  /* Total number of element blocks              */
  int Num_Node_Set;  /* Total number of node sets defined in the    *
                      * mesh exoII file			     */
  int Num_Side_Set;  /* Total number of side sets defined in the    *
                      * mesh exoII file			     */
  int Max_NP_Elem;   /* Maximum number of nodes in any element      *
                      *  - this is a global quantity                */
  int Num_QA_Recs;   /* Number of QA records in original file       */
  int Num_Info_Recs; /* Number of Info records in original file     */

  int Num_Coordinate_Frames;

  /*---------------------------------------------------------------------------*/
  /*    VARIABLES THAT DEAL WITH SPECIFICATION OF LOAD BALANCE PROPERTIES      */
  /*		THAT ARE THE DIFFERENT ON EACH PROCESSOR		     */
  /*---------------------------------------------------------------------------*/

  INT *Num_Internal_Nodes; /* Number of internal nodes on the current proc*/
  INT *Num_Border_Nodes;   /* Number of border nodes on the current proc  */
  INT *Num_External_Nodes; /* Number of external nodes on the current proc*/
  INT *Num_Internal_Elems; /* Number of Elements on the local processor.  */

  INT *Num_Border_Elems; /* Number of Elements on the local processor *
                          * and shared by other processors. */

  INT *Num_N_Comm_Maps; /* Number of nodal communication maps */

  INT *Num_E_Comm_Maps; /* Number of elemental communication maps */

  ELEM_COMM_MAP<INT> **E_Comm_Map; /* Elemental communication map structure */

  NODE_COMM_MAP<INT> **N_Comm_Map; /* Nodal communication map structure */

  INT **GNodes; /* Data structure which contains the internal, *
                 * border, and external nodes on each processor*
                 * They are structured in that order, and      *
                 * monotonically within each category          *
                 *  Type: int vector of length                 *
                 *       (Num_Internal_Nodes + Num_Border_Nodes*
                 Num_External_Nodes) 		     */

  INT **GElems; /* Data structure which contains the internal  *
                 * elements on each processor.  It is a map    *
                 * from the local element number to the global *
                 * element number.  			     *
                 *  Type: int vector of length                 *
                 *        Num_Internal_Elems                   */

  INT **Proc_Global_Node_Id_Map; /* Data structure which contains the internal  *
                                  * nodes on each processor.  It is a map       *
                                  * from the local node number to the global    *
                                  * node id (as found in node_num_map)	     *
                                  *  Type: int vector of length                 *
                                  *       (Num_Internal_Nodes + Num_Border_Nodes*
                                  *        Num_External_Nodes) 		     */

  INT **Proc_Global_Elem_Id_Map; /* Data structure which contains the internal  *
                                  * elements on each processor.  It is a map    *
                                  * from the local element number to the global *
                                  * element id (as found in elem_num_map)	     *
                                  *  Type: int vector of length                 *
                                  *        Num_Internal_Elems                   */

  INT **Elem_Map; /* Map for Nemesis output */

  INT **GElem_Blks; /* Data structure which contains the mapping   *
                     * from the local element block number to the  *
                     * global element block number                 *
                     *  Type: int vector of length                 *
                     *        Proc_Num_Elem_Blk                    */

  INT **GNode_Sets; /* Data structure which contains the mapping   *
                     * from the local node set number to the       *
                     * global node set number                      *
                     *  Type: int vector of length                 *
                     *        Proc_Num_Node_Sets                   */

  INT **GSide_Sets; /* Data structure which contains the mapping   *
                     * from the local side set number to the       *
                     * global side set number                      *
                     *  Type: int vector of length                 *
                     *        Proc_Num_Side_Sets                   */

  /*---------------------------------------------------------------------------*/
  /*    VARIABLES THAT DEAL WITH SPECIFICATION OF NODAL PROPERTIES   	     */
  /*		THAT ARE THE DIFFERENT ON EACH PROCESSOR		     */
  /*---------------------------------------------------------------------------*/

  T ***Coor; /* 2d dynamically allocated array containing   *
              * the physical space coordinates for each     *
              * node, defined on the local processor.       *
              *  Type: double/float vector of length        *
              *        Num_Dim   by			     *
              *       (Num_Internal_Nodes + Num_Border_Nodes*
              *        Num_External_Nodes) 		     */

  /*---------------------------------------------------------------------------*/
  /*    VARIABLES THAT DEAL WITH SPECIFICATION OF ELEMENT PROPERTIES   	     */
  /*		THAT ARE THE DIFFERENT ON EACH PROCESSOR		     */
  /*---------------------------------------------------------------------------*/

  INT **Proc_Connect_Ptr; /* Vector of pointers to the start of each     *
                           * element in Proc_Elem_Connect                *
                           *  Type: int vector of length                 *
                           *        Num_Internal_Elems                   */

  int **Elem_Type; /* Vector which contains the element code for  *
                    * each defined on the current processor       *
                    *  Type: int vector of length                 *
                    *        Num_Internal_Elems                   */

  INT **Proc_Elem_Connect; /* Connectivity lists for the elements        *
                            * required by the current processor          *
                            *  Type: int vector of variable length       */

  T **Proc_Elem_Attr; /* Attribute list for the elements        *
                       * required by the current processor      *
                       *  Type: float vector of variable length */

  /*---------------------------------------------------------------------------*/
  /*    VARIABLES THAT DEAL WITH SPECIFICATION OF ELEMENT BLOCK PROPERTIES     */
  /*		THAT ARE THE DIFFERENT ON EACH PROCESSOR		     */
  /*---------------------------------------------------------------------------*/

  int *Proc_Num_Elem_Blk; /* Number of element blocks on this processor  */

  INT **Proc_Num_Elem_In_Blk;
  /* Number of elements in the processor's       *
   * element blocks                              *
   *  Type: int vector of length                 *
   *        Proc_Num_Elem_Blk                    */

  INT **Proc_Elem_Blk_Ids;
  /* Element block id's for the processor's      *
   * element blocks                              *
   *  Type: int vector of length                 *
   *        Proc_Num_Elem_Blk                    */

  INT **Proc_Elem_Blk_Types;
  /* Element block types for the processor's     *
   * element blocks                              *
   *  Type: int vector of length                 *
   *        Proc_Num_Elem_Blk                    */

  INT **Proc_Nodes_Per_Elem;
  /* Number of nodes per element for each        *
   * block on the current processor              *
   *  Type: int vector of length                 *
   *        Proc_Num_Elem_Blk                    */

  INT **Proc_Num_Attr; /* Number of attributes for each block on the  *
                        * current processor                           *
                        *  Type: int vector of length                 *
                        *        Proc_Num_Elem_Blk                    */

  INT *Elem_Blk_2_Matls;
  /* Mapping of element block NUMBERS to material*
   * IDs  	                	             *
   *  Type: int vector of length                 *
   *        Proc_Num_Elem_Blk                    */

  /*---------------------------------------------------------------------------*/
  /*	VARIABLES THAT DEAL WITH SPECIFICATION OF NODE SETS		     */
  /*		THAT ARE THE DIFFERENT ON EACH PROCESSOR		     */
  /*---------------------------------------------------------------------------*/

  int *Proc_Num_Node_Sets; /* Number of node sets on the current proc   */

  INT *Proc_NS_List_Length; /* Total length of all the node set lists      */

  INT **Proc_NS_Ids; /* Node sets ids for the node sets in a given  *
                      * processor                                   *
                      *  Type: int vector of length                 *
                      *        Proc_Num_Node_Sets                   */

  INT **Proc_NS_Count; /* Number of nodes in each node set defined on *
                        * the current processor.                      *
                        *  Type: int vector of length                 *
                        *        Proc_Num_Node_Sets                   */

  INT **Proc_NS_DF_Count;
  /* Number of distribution factors associated   *
   * with a node set.                            *
   *  Type: int vector of length                 *
   *        Proc_Num_Node_Sets                   */

  INT **Proc_NS_Pointers; /* Vector of pointers into Proc_NS_List for   *
                           * each node set defined on the current proc  *
                           *  Type: int vector of length                *
                           *        Proc_Num_Node_Sets                  */

  INT **Proc_NS_List; /* Node sets list record for the nodes sets in *
                       * a given processor                           *
                       *  Type: int vector of length                 *
                       *        Proc_NS_List_Length                  */

  T **Proc_NS_Dist_Fact;
  /* Node sets distribution factors for the node *
   * sets on a given processor                   *
   *  Type: float vector of length               *
   *        Proc_NS_List_Length                  */

  INT **Proc_NS_GNMap_List;
  /* Map local node position to global nodeset node list position *
   * in a given processor                        *
   *  Type: int vector of length                 *
   *        Proc_NS_List_Length             */

  /*---------------------------------------------------------------------------*/
  /*	VARIABLES THAT DEAL WITH SPECIFICATION OF SIDE SETS		     */
  /*		THAT ARE THE DIFFERENT ON EACH PROCESSOR		     */
  /*---------------------------------------------------------------------------*/

  int *Proc_Num_Side_Sets; /* Number of side sets on the current proc  */

  INT *Proc_SS_Elem_List_Length;
  /* Total length of all the side set node lists */
  INT *Proc_SS_Side_List_Length;
  /* Total length of all the side set element    *
   * lists                                       */

  INT **Proc_SS_Ids; /* Side sets ids for side sets in a given      *
                      * processor                                   *
                      *  Type: int vector of length                 *
                      *        Proc_Num_Side_Sets                   */
  INT **Proc_SS_Elem_Count;
  /* Side sets count record for elements (see the*
   * EXODUS manual) for the side sets in a given *
   * processor                                   *
   *  Type: int vector of length                 *
   *        Proc_Num_Side_Sets                   */

  INT **Proc_SS_DF_Pointers;
  /* Pointer into the distribution factors for   *
   * the side sets.                              *
   *  Type: int array of dimensions              *
   *    Proc_Num_Side_Sets            */

  T **Proc_SS_Dist_Fact;
  /* Pointer for storage of the distribution     *
   * factors.                                    */

  INT **Proc_SS_DF_Count;
  /* Count of the number of distribution         *
   * factors in this side set.                   *
   *  Type: int vector of length                 *
   *        Proc_Num_Side_Sets                   */

  INT **Proc_SS_Side_Count;
  // Side sets count record for nodes (see the
  // EXODUS manual) for the side sets in a given
  // processor
  //  Type: int vector of length
  //        Proc_Num_Side_Sets

  INT **Proc_SS_Elem_Pointers;
  /* Side sets pointer record for elements (see  *
   * the EXODUS manual) for the side sets in a   *
   * given processor                             *
   *  Type: int vector of length                 *
   *        Proc_Num_Side_Sets                   */

  INT **Proc_SS_Side_Pointers;
  /* Side sets pointer record for nodes (see     *
   * the EXODUS manual) for the side sets in a   *
   * given processor                             *
   *  Type: int vector of length                 *
   *        Proc_Num_Side_Sets                   */

  INT **Proc_SS_Elem_List;
  /* Side sets element list record for the side  *
   * sets in a given processor                   *
   *  Type: int vector of length                 *
   *        Proc_Num_Side_Sets                   */

  INT **Proc_SS_Side_List;
  /* Side sets side list record for the side sets*
   * in a given processor                        *
   *  Type: int vector of length                 *
   *        Proc_SS_Elem_List_Length             */

  INT **Proc_SS_GEMap_List;
  /* Map local element position to global sideset element list position *
   * in a given processor                        *
   *  Type: int vector of length                 *
   *        Proc_SS_Elem_List_Length             */

  /*---------------------------------------------------------------------------*/
  /*		VARIABLES THAT DEAL WITH GENERAL INFORMATION THAT IS         */
  /*			THE SAME ON EVERY PROCESSOR                          */
  /*---------------------------------------------------------------------------*/

  char **QA_Record; /* The QA Records from the original file     */

  char **Info_Record; /* The Information Records from the original *
                       * file                                      */
  INT * Coordinate_Frame_Ids;
  T *   Coordinate_Frame_Coordinates;
  char *Coordinate_Frame_Tags;

  Globals()
      : Num_Dim(0), Num_Node(0), Num_Elem(0), Num_Elem_Blk(0), Num_Node_Set(0), Num_Side_Set(0),
        Max_NP_Elem(0), Num_QA_Recs(0), Num_Info_Recs(0), Num_Coordinate_Frames(0),
        Num_Internal_Nodes(nullptr), Num_Border_Nodes(nullptr), Num_External_Nodes(nullptr),
        Num_Internal_Elems(nullptr), Num_Border_Elems(nullptr), Num_N_Comm_Maps(nullptr),
        Num_E_Comm_Maps(nullptr), E_Comm_Map(nullptr), N_Comm_Map(nullptr), GNodes(nullptr),
        GElems(nullptr), Proc_Global_Node_Id_Map(nullptr), Proc_Global_Elem_Id_Map(nullptr),
        Elem_Map(nullptr), GElem_Blks(nullptr), GNode_Sets(nullptr), GSide_Sets(nullptr),
        Coor(nullptr), Proc_Connect_Ptr(nullptr), Elem_Type(nullptr), Proc_Elem_Connect(nullptr),
        Proc_Elem_Attr(nullptr), Proc_Num_Elem_Blk(nullptr), Proc_Num_Elem_In_Blk(nullptr),
        Proc_Elem_Blk_Ids(nullptr), Proc_Elem_Blk_Types(nullptr), Proc_Nodes_Per_Elem(nullptr),
        Proc_Num_Attr(nullptr), Elem_Blk_2_Matls(nullptr),

        Proc_Num_Node_Sets(nullptr), Proc_NS_List_Length(nullptr), Proc_NS_Ids(nullptr),
        Proc_NS_Count(nullptr), Proc_NS_DF_Count(nullptr), Proc_NS_Pointers(nullptr),
        Proc_NS_List(nullptr), Proc_NS_Dist_Fact(nullptr), Proc_NS_GNMap_List(nullptr),

        Proc_Num_Side_Sets(nullptr), Proc_SS_Elem_List_Length(nullptr),
        Proc_SS_Side_List_Length(nullptr), Proc_SS_Ids(nullptr), Proc_SS_Elem_Count(nullptr),
        Proc_SS_DF_Pointers(nullptr), Proc_SS_Dist_Fact(nullptr), Proc_SS_DF_Count(nullptr),
        Proc_SS_Side_Count(nullptr), Proc_SS_Elem_Pointers(nullptr), Proc_SS_Side_Pointers(nullptr),
        Proc_SS_Elem_List(nullptr), Proc_SS_Side_List(nullptr), Proc_SS_GEMap_List(nullptr),

        QA_Record(nullptr), Info_Record(nullptr), Coordinate_Frame_Ids(nullptr),
        Coordinate_Frame_Coordinates(nullptr), Coordinate_Frame_Tags(nullptr)
  {
  }

  ~Globals()
  {
    safe_free((void **)&Proc_Num_Elem_Blk);
    safe_free((void **)&Proc_Num_Node_Sets);
    safe_free((void **)&Proc_Num_Side_Sets);
    safe_free((void **)&Proc_NS_List_Length);
    safe_free((void **)&Proc_SS_Elem_List_Length);

    safe_free((void **)&Coor);
    safe_free((void **)&Proc_Elem_Attr);
    safe_free((void **)&Num_Internal_Nodes);

    safe_free((void **)&Info_Record);
    safe_free((void **)&GElem_Blks);
    safe_free((void **)&Elem_Type);
    safe_free((void **)&Elem_Type);
  }
};

#endif
