/*
 * Copyright(C) 1999-2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include <cstdlib>
#include <exodusII.h>
#include <rf_allo.h>
#include <vector>

/*---------------------------------------------------------------------------*/
/*      STRUCTURES FOR COMMUNICATION MAPS                                    */
/*---------------------------------------------------------------------------*/
template <typename INT> struct ELEM_COMM_MAP
{
  size_t           map_id{0};
  size_t           elem_cnt{0};
  std::vector<INT> elem_ids{};
  std::vector<INT> side_ids{};
  std::vector<INT> proc_ids{};
};

template <typename INT> struct NODE_COMM_MAP
{
  size_t           map_id{0};
  size_t           node_cnt{0};
  std::vector<INT> node_ids{};
  std::vector<INT> proc_ids{};
};

/*---------------------------------------------------------------------------*/
/*      GLOBAL QUANTITIES THAT ARE THE SAME ON ALL PROCESSORS               */
/*---------------------------------------------------------------------------*/
template <typename T, typename INT> class Globals
{
public:
  int    Num_Dim{0};    /* Number of physical dimensions in the problem*/
  size_t Num_Node{0};   /* Total number of nodes in the entire mesh    *
                         * - this is a global quantity                 */
  size_t Num_Elem{0};   /* Total number of elements in the entire mesh *
                         * - this is a global quantity                 */
  int Num_Elem_Blk{0};  /* Total number of element blocks              */
  int Num_Node_Set{0};  /* Total number of node sets defined in the    *
                         * mesh exoII file                      */
  int Num_Side_Set{0};  /* Total number of side sets defined in the    *
                         * mesh exoII file                      */
  int Max_NP_Elem{0};   /* Maximum number of nodes in any element      *
                         *  - this is a global quantity                */
  int Num_QA_Recs{0};   /* Number of QA records in original file       */
  int Num_Info_Recs{0}; /* Number of Info records in original file     */

  int Num_Coordinate_Frames{0};
  int Num_Assemblies{0};

  /*---------------------------------------------------------------------------*/
  /*    VARIABLES THAT DEAL WITH SPECIFICATION OF LOAD BALANCE PROPERTIES      */
  /*            THAT ARE THE DIFFERENT ON EACH PROCESSOR                     */
  /*---------------------------------------------------------------------------*/

  std::vector<INT> Num_Internal_Nodes{}; /* Number of internal nodes on the current proc*/
  std::vector<INT> Num_Border_Nodes{};   /* Number of border nodes on the current proc  */
  std::vector<INT> Num_External_Nodes{}; /* Number of external nodes on the current proc*/
  std::vector<INT> Num_Internal_Elems{}; /* Number of Elements on the local processor.  */

  std::vector<INT> Num_Border_Elems{}; /* Number of Elements on the local processor *
                                        * and shared by other processors. */

  std::vector<INT> Num_N_Comm_Maps{}; /* Number of nodal communication maps */
  std::vector<INT> Num_E_Comm_Maps{}; /* Number of elemental communication maps */

  std::vector<ELEM_COMM_MAP<INT>> E_Comm_Map{}; /* Elemental communication map structure */
  std::vector<NODE_COMM_MAP<INT>> N_Comm_Map{}; /* Nodal communication map structure */

  std::vector<std::vector<INT>> GNodes{}; /* Data structure which contains the internal, *
                 * border, and external nodes on each processor*
                 * They are structured in that order, and      *
                 * monotonically within each category          *
                 *  Type: int vector of length                 *
                 *       (Num_Internal_Nodes + Num_Border_Nodes*
                 Num_External_Nodes)                 */

  std::vector<std::vector<INT>> GElems{}; /* Data structure which contains the internal  *
                                           * elements on each processor.  It is a map    *
                                           * from the local element number to the global *
                                           * element number.                           *
                                           *  Type: int vector of length                 *
                                           *        Num_Internal_Elems                   */

  std::vector<std::vector<INT>> Elem_Map{}; /* Map for Nemesis output */

  std::vector<std::vector<INT>>
      Proc_Global_Node_Id_Map{}; /* Data structure which contains the internal  *
                                  * nodes on each processor.  It is a map       *
                                  * from the local node number to the global    *
                                  * node id (as found in node_num_map)       *
                                  *  Type: int vector of length                 *
                                  *       (Num_Internal_Nodes + Num_Border_Nodes*
                                  *        Num_External_Nodes)               */

  std::vector<std::vector<INT>>
      Proc_Global_Elem_Id_Map{}; /* Data structure which contains the internal  *
                                  * elements on each processor.  It is a map    *
                                  * from the local element number to the global *
                                  * element id (as found in elem_num_map)            *
                                  *  Type: int vector of length                 *
                                  *        Num_Internal_Elems                   */

  INT **GElem_Blks{nullptr}; /* Data structure which contains the mapping   *
                              * from the local element block number to the  *
                              * global element block number                 *
                              *  Type: int vector of length                 *
                              *        Proc_Num_Elem_Blk                    */

  INT **GNode_Sets{nullptr}; /* Data structure which contains the mapping   *
                              * from the local node set number to the       *
                              * global node set number                      *
                              *  Type: int vector of length                 *
                              *        Proc_Num_Node_Sets                   */

  INT **GSide_Sets{nullptr}; /* Data structure which contains the mapping   *
                              * from the local side set number to the       *
                              * global side set number                      *
                              *  Type: int vector of length                 *
                              *        Proc_Num_Side_Sets                   */

  /*---------------------------------------------------------------------------*/
  /*    VARIABLES THAT DEAL WITH SPECIFICATION OF NODAL PROPERTIES           */
  /*            THAT ARE THE DIFFERENT ON EACH PROCESSOR                     */
  /*---------------------------------------------------------------------------*/

  T ***Coor{nullptr}; /* 2d dynamically allocated array containing   *
                       * the physical space coordinates for each     *
                       * node, defined on the local processor.       *
                       *  Type: double/float vector of length        *
                       *        Num_Dim   by                          *
                       *       (Num_Internal_Nodes + Num_Border_Nodes*
                       *        Num_External_Nodes)                   */

  /*---------------------------------------------------------------------------*/
  /*    VARIABLES THAT DEAL WITH SPECIFICATION OF ELEMENT PROPERTIES         */
  /*            THAT ARE THE DIFFERENT ON EACH PROCESSOR                     */
  /*---------------------------------------------------------------------------*/

  INT **Proc_Connect_Ptr{nullptr}; /* Vector of pointers to the start of each     *
                                    * element in Proc_Elem_Connect                *
                                    *  Type: int vector of length                 *
                                    *        Num_Internal_Elems                   */

  int **Elem_Type{nullptr}; /* Vector which contains the element code for  *
                             * each defined on the current processor       *
                             *  Type: int vector of length                 *
                             *        Num_Internal_Elems                   */

  INT **Proc_Elem_Connect{nullptr}; /* Connectivity lists for the elements        *
                                     * required by the current processor          *
                                     *  Type: int vector of variable length       */

  std::vector<std::vector<T>> Proc_Elem_Attr{}; /* Attribute list for the elements        *
                                                 * required by the current processor      *
                                                 *  Type: float vector of variable length */

  /*---------------------------------------------------------------------------*/
  /*    VARIABLES THAT DEAL WITH SPECIFICATION OF ELEMENT BLOCK PROPERTIES     */
  /*            THAT ARE THE DIFFERENT ON EACH PROCESSOR                     */
  /*---------------------------------------------------------------------------*/

  int *Proc_Num_Elem_Blk{nullptr}; /* Number of element blocks on this processor  */

  INT **Proc_Num_Elem_In_Blk{nullptr};
  /* Number of elements in the processor's       *
   * element blocks                              *
   *  Type: int vector of length                 *
   *        Proc_Num_Elem_Blk                    */

  INT **Proc_Elem_Blk_Ids{nullptr};
  /* Element block id's for the processor's      *
   * element blocks                              *
   *  Type: int vector of length                 *
   *        Proc_Num_Elem_Blk                    */

  INT **Proc_Elem_Blk_Types{nullptr};
  /* Element block types for the processor's     *
   * element blocks                              *
   *  Type: int vector of length                 *
   *        Proc_Num_Elem_Blk                    */

  INT **Proc_Nodes_Per_Elem{nullptr};
  /* Number of nodes per element for each        *
   * block on the current processor              *
   *  Type: int vector of length                 *
   *        Proc_Num_Elem_Blk                    */

  INT **Proc_Num_Attr{nullptr}; /* Number of attributes for each block on the  *
                                 * current processor                           *
                                 *  Type: int vector of length                 *
                                 *        Proc_Num_Elem_Blk                    */

  INT *Elem_Blk_2_Matls{nullptr};
  /* Mapping of element block NUMBERS to material*
   * IDs                                             *
   *  Type: int vector of length                 *
   *        Proc_Num_Elem_Blk                    */

  /*---------------------------------------------------------------------------*/
  /*    VARIABLES THAT DEAL WITH SPECIFICATION OF NODE SETS                  */
  /*            THAT ARE THE DIFFERENT ON EACH PROCESSOR                     */
  /*---------------------------------------------------------------------------*/

  int *Proc_Num_Node_Sets{nullptr}; /* Number of node sets on the current proc   */

  INT *Proc_NS_List_Length{nullptr}; /* Total length of all the node set lists      */

  INT **Proc_NS_Ids{nullptr}; /* Node sets ids for the node sets in a given  *
                               * processor                                   *
                               *  Type: int vector of length                 *
                               *        Proc_Num_Node_Sets                   */

  INT **Proc_NS_Count{nullptr}; /* Number of nodes in each node set defined on *
                                 * the current processor.                      *
                                 *  Type: int vector of length                 *
                                 *        Proc_Num_Node_Sets                   */

  INT **Proc_NS_DF_Count{nullptr};
  /* Number of distribution factors associated   *
   * with a node set.                            *
   *  Type: int vector of length                 *
   *        Proc_Num_Node_Sets                   */

  INT **Proc_NS_Pointers{nullptr}; /* Vector of pointers into Proc_NS_List for   *
                                    * each node set defined on the current proc  *
                                    *  Type: int vector of length                *
                                    *        Proc_Num_Node_Sets                  */

  INT **Proc_NS_List{nullptr}; /* Node sets list record for the nodes sets in *
                                * a given processor                           *
                                *  Type: int vector of length                 *
                                *        Proc_NS_List_Length                  */

  std::vector<std::vector<T>> Proc_NS_Dist_Fact{};
  /* Node sets distribution factors for the node *
   * sets on a given processor                   *
   *  Type: float vector of length               *
   *        Proc_NS_List_Length                  */

  INT **Proc_NS_GNMap_List{nullptr};
  /* Map local node position to global nodeset node list position *
   * in a given processor                        *
   *  Type: int vector of length                 *
   *        Proc_NS_List_Length             */

  /*---------------------------------------------------------------------------*/
  /*    VARIABLES THAT DEAL WITH SPECIFICATION OF SIDE SETS                  */
  /*            THAT ARE THE DIFFERENT ON EACH PROCESSOR                     */
  /*---------------------------------------------------------------------------*/

  int *Proc_Num_Side_Sets{nullptr}; /* Number of side sets on the current proc  */

  INT *Proc_SS_Elem_List_Length{nullptr};
  /* Total length of all the side set node lists */
  INT *Proc_SS_Side_List_Length{nullptr};
  /* Total length of all the side set element    *
   * lists                                       */

  INT **Proc_SS_Ids{nullptr}; /* Side sets ids for side sets in a given      *
                               * processor                                   *
                               *  Type: int vector of length                 *
                               *        Proc_Num_Side_Sets                   */
  INT **Proc_SS_Elem_Count{nullptr};
  /* Side sets count record for elements (see the*
   * EXODUS manual) for the side sets in a given *
   * processor                                   *
   *  Type: int vector of length                 *
   *        Proc_Num_Side_Sets                   */

  INT **Proc_SS_DF_Pointers{nullptr};
  /* Pointer into the distribution factors for   *
   * the side sets.                              *
   *  Type: int array of dimensions              *
   *    Proc_Num_Side_Sets            */

  std::vector<std::vector<T>> Proc_SS_Dist_Fact{};
  /* Pointer for storage of the distribution     *
   * factors.                                    */

  INT **Proc_SS_DF_Count{nullptr};
  /* Count of the number of distribution         *
   * factors in this side set.                   *
   *  Type: int vector of length                 *
   *        Proc_Num_Side_Sets                   */

  INT **Proc_SS_Side_Count{nullptr};
  // Side sets count record for nodes (see the
  // EXODUS manual) for the side sets in a given
  // processor
  //  Type: int vector of length
  //        Proc_Num_Side_Sets

  INT **Proc_SS_Elem_Pointers{nullptr};
  /* Side sets pointer record for elements (see  *
   * the EXODUS manual) for the side sets in a   *
   * given processor                             *
   *  Type: int vector of length                 *
   *        Proc_Num_Side_Sets                   */

  INT **Proc_SS_Side_Pointers{nullptr};
  /* Side sets pointer record for nodes (see     *
   * the EXODUS manual) for the side sets in a   *
   * given processor                             *
   *  Type: int vector of length                 *
   *        Proc_Num_Side_Sets                   */

  INT **Proc_SS_Elem_List{nullptr};
  /* Side sets element list record for the side  *
   * sets in a given processor                   *
   *  Type: int vector of length                 *
   *        Proc_Num_Side_Sets                   */

  INT **Proc_SS_Side_List{nullptr};
  /* Side sets side list record for the side sets*
   * in a given processor                        *
   *  Type: int vector of length                 *
   *        Proc_SS_Elem_List_Length             */

  INT **Proc_SS_GEMap_List{nullptr};
  /* Map local element position to global sideset element list position *
   * in a given processor                        *
   *  Type: int vector of length                 *
   *        Proc_SS_Elem_List_Length             */

  /*---------------------------------------------------------------------------*/
  /*            VARIABLES THAT DEAL WITH GENERAL INFORMATION THAT IS         */
  /*                    THE SAME ON EVERY PROCESSOR                          */
  /*---------------------------------------------------------------------------*/

  char **QA_Record{nullptr}; /* The QA Records from the original file     */

  char **Info_Record{nullptr}; /* The Information Records from the original *
                                * file                                      */
  INT  *Coordinate_Frame_Ids{nullptr};
  T    *Coordinate_Frame_Coordinates{nullptr};
  char *Coordinate_Frame_Tags{nullptr};

  std::vector<ex_assembly> Assemblies{};

  Globals() = default;

  ~Globals()
  {
    safe_free((void **)&Proc_Num_Elem_Blk);
    safe_free((void **)&Proc_Num_Node_Sets);
    safe_free((void **)&Proc_Num_Side_Sets);
    safe_free((void **)&Proc_NS_List_Length);
    safe_free((void **)&Proc_SS_Elem_List_Length);

    safe_free((void **)&Coor);

    safe_free((void **)&Info_Record);
    safe_free((void **)&GElem_Blks);
    safe_free((void **)&Elem_Type);
    safe_free((void **)&Elem_Type);
  }
};
