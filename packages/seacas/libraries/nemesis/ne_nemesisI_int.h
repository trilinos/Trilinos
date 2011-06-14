/*
 * Copyright (c) 1998 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
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
 *     * Neither the name of Sandia Corporation nor the names of its
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
/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: ne_nemesisI_int.h,v $
 *
 * $Author: gdsjaar $
 *
 * $Date: 2007/10/31 21:39:17 $
 *
 * $Revision: 1.21 $
 *
 * $Name:  $
 *====================================================================*/

/**********************************************************************
 * This file contains constant definitions used internally NEMESIS.
 **********************************************************************/

#ifndef _NE_NEMESISI_INT_H
#define _NE_NEMESISI_INT_H

#define NEMESIS_FILE_VERSION		2.6

/* mms: start changing things, no longer need to distinguish **
** between processors for dimensions and variables           */
#define VAR_ELBLK_IDS_GLOBAL            "el_blk_ids_global"
#define VAR_ELBLK_CNT_GLOBAL            "el_blk_cnt_global"
#define VAR_NS_IDS_GLOBAL               "ns_ids_global"
#define VAR_NS_NODE_CNT_GLOBAL          "ns_node_cnt_global"
#define VAR_NS_DF_CNT_GLOBAL            "ns_df_cnt_global"
#define VAR_SS_IDS_GLOBAL               "ss_ids_global"
#define VAR_SS_SIDE_CNT_GLOBAL          "ss_side_cnt_global"
#define VAR_SS_DF_CNT_GLOBAL            "ss_df_cnt_global"
#define VAR_FILE_TYPE                   "nem_ftype"
#define VAR_COMM_MAP                    "comm_map"
#define VAR_NODE_MAP_INT                "node_mapi"
#define VAR_NODE_MAP_INT_IDX            "node_mapi_idx"
#define VAR_NODE_MAP_BOR                "node_mapb"
#define VAR_NODE_MAP_BOR_IDX            "node_mapb_idx"
#define VAR_NODE_MAP_EXT                "node_mape"
#define VAR_NODE_MAP_EXT_IDX            "node_mape_idx"
#define VAR_ELEM_MAP_INT                "elem_mapi"
#define VAR_ELEM_MAP_INT_IDX            "elem_mapi_idx"
#define VAR_ELEM_MAP_BOR                "elem_mapb"
#define VAR_ELEM_MAP_BOR_IDX            "elem_mapb_idx"
#define VAR_INT_N_STAT                  "int_n_stat"
#define VAR_BOR_N_STAT                  "bor_n_stat"
#define VAR_EXT_N_STAT                  "ext_n_stat"
#define VAR_INT_E_STAT                  "int_e_stat"
#define VAR_BOR_E_STAT                  "bor_e_stat"
#define VAR_N_COMM_IDS                  "n_comm_ids"
#define VAR_N_COMM_STAT                 "n_comm_stat"
#define VAR_N_COMM_INFO_IDX             "n_comm_info_idx"
#define VAR_E_COMM_IDS                  "e_comm_ids"
#define VAR_E_COMM_STAT                 "e_comm_stat"
#define VAR_E_COMM_INFO_IDX             "e_comm_info_idx"
#define VAR_N_COMM_NIDS                 "n_comm_nids"
#define VAR_N_COMM_PROC                 "n_comm_proc"
#define VAR_N_COMM_DATA_IDX             "n_comm_data_idx"
#define VAR_E_COMM_EIDS                 "e_comm_eids"
#define VAR_E_COMM_SIDS                 "e_comm_sids"
#define VAR_E_COMM_PROC                 "e_comm_proc"
#define VAR_E_COMM_DATA_IDX             "e_comm_data_idx"

/* variable which allow me to take the tail off some names */
#define END_CMAP_IDX                    "_info_idx"
#define END_CMAP_IDS                    "_ids"

#define DIM_NUM_INT_NODES               "num_int_node"
#define DIM_NUM_BOR_NODES               "num_bor_node"
#define DIM_NUM_EXT_NODES               "num_ext_node"
#define DIM_NUM_INT_ELEMS               "num_int_elem"
#define DIM_NUM_BOR_ELEMS               "num_bor_elem"
#define DIM_NUM_PROCS                   "num_processors"
#define DIM_NUM_PROCS_F                 "num_procs_file"
#define DIM_NUM_NODES_GLOBAL            "num_nodes_global"
#define DIM_NUM_ELEMS_GLOBAL            "num_elems_global"
#define DIM_NUM_NS_GLOBAL               "num_ns_global"
#define DIM_NUM_SS_GLOBAL               "num_ss_global"
#define DIM_NUM_ELBLK_GLOBAL            "num_el_blk_global"
#define DIM_NUM_N_CMAPS                 "num_n_cmaps"
#define DIM_NUM_E_CMAPS                 "num_e_cmaps"
#define DIM_NCNT_CMAP                   "ncnt_cmap"
#define DIM_ECNT_CMAP                   "ecnt_cmap"

#endif /* _NE_NEMESISI_INT_H */
