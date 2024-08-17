// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef CWRAPFMANGLEH_
#define CWRAPFMANGLEH_

#ifndef TRILINOS_NO_CONFIG_H
/*--------------------------------------------------------------------*/
/* Autotools and cmake create a macro to handle the mangling.         */

#include "Zoltan_config.h"

#define Zfw_Initialize                 FC_FUNC_(zfw_initialize,ZFW_INITIALIZE)
#define Zfw_Initialize1                FC_FUNC_(zfw_initialize1,ZFW_INITIALIZE1)
#define Zfw_Create                     FC_FUNC_(zfw_create,ZFW_CREATE)
#define Zfw_Copy                       FC_FUNC_(zfw_copy,ZFW_COPY)
#define Zfw_Copy_To                    FC_FUNC_(zfw_copy_to,ZFW_COPY_TO)
#define Zfw_Destroy                    FC_FUNC_(zfw_destroy,ZFW_DESTROY)
#define Zfw_Align                      FC_FUNC_(zfw_align,ZFW_ALIGN)
#define Zfw_Memory_Stats               FC_FUNC_(zfw_memory_stats,ZFW_MEMORY_STATS)
#define Zfw_Set_Fn0f                   FC_FUNC_(zfw_set_fn0f,ZFW_SET_FN0F)
#define Zfw_Set_Fn1f                   FC_FUNC_(zfw_set_fn1f,ZFW_SET_FN1F)
#define Zfw_Set_Fn2f                   FC_FUNC_(zfw_set_fn2f,ZFW_SET_FN2F)
#define Zfw_Set_Fn3f                   FC_FUNC_(zfw_set_fn3f,ZFW_SET_FN3F)
#define Zfw_Set_Fn4f                   FC_FUNC_(zfw_set_fn4f,ZFW_SET_FN4F)
#define Zfw_Set_Fn5f                   FC_FUNC_(zfw_set_fn5f,ZFW_SET_FN5F)
#define Zfw_Set_Fn6f                   FC_FUNC_(zfw_set_fn6f,ZFW_SET_FN6F)
#define Zfw_Set_Fn7f                   FC_FUNC_(zfw_set_fn7f,ZFW_SET_FN7F)
#define Zfw_Set_Fn8f                   FC_FUNC_(zfw_set_fn8f,ZFW_SET_FN8F)
#define Zfw_Set_Fn9f                   FC_FUNC_(zfw_set_fn9f,ZFW_SET_FN9F)
#define Zfw_Set_FnAf                   FC_FUNC_(zfw_set_fnaf,ZFW_SET_FNAF)
#define Zfw_Set_FnBf                   FC_FUNC_(zfw_set_fnbf,ZFW_SET_FNBF)
#define Zfw_Set_Fn0s                   FC_FUNC_(zfw_set_fn0s,ZFW_SET_FN0S)
#define Zfw_Set_Fn1s                   FC_FUNC_(zfw_set_fn1s,ZFW_SET_FN1S)
#define Zfw_Set_Fn2s                   FC_FUNC_(zfw_set_fn2s,ZFW_SET_FN2S)
#define Zfw_Set_Fn3s                   FC_FUNC_(zfw_set_fn3s,ZFW_SET_FN3S)
#define Zfw_Set_Fn4s                   FC_FUNC_(zfw_set_fn4s,ZFW_SET_FN4S)
#define Zfw_Set_Fn5s                   FC_FUNC_(zfw_set_fn5s,ZFW_SET_FN5S)
#define Zfw_Set_Fn6s                   FC_FUNC_(zfw_set_fn6s,ZFW_SET_FN6S)
#define Zfw_Set_Fn7s                   FC_FUNC_(zfw_set_fn7s,ZFW_SET_FN7S)
#define Zfw_Set_Fn8s                   FC_FUNC_(zfw_set_fn8s,ZFW_SET_FN8S)
#define Zfw_Set_Fn9s                   FC_FUNC_(zfw_set_fn9s,ZFW_SET_FN9S)
#define Zfw_Set_FnAs                   FC_FUNC_(zfw_set_fnas,ZFW_SET_FNAS)
#define Zfw_Set_FnBs                   FC_FUNC_(zfw_set_fnbs,ZFW_SET_FNBS)
#define Zfw_Set_Param                  FC_FUNC_(zfw_set_param,ZFW_SET_PARAM)
#define Zfw_Set_Param_Vec              FC_FUNC_(zfw_set_param_vec,ZFW_SET_PARAM_VEC)
#define Zfw_LB_Partition               FC_FUNC_(zfw_lb_partition,ZFW_LB_PARTITION)
#define Zfw_LB_Eval                    FC_FUNC_(zfw_lb_eval,ZFW_LB_EVAL)
#define Zfw_LB_Set_Part_Sizes          FC_FUNC_(zfw_lb_set_part_sizes,ZFW_LB_SET_PART_SIZES)
#define Zfw_LB_Point_Assign            FC_FUNC_(zfw_lb_point_assign,ZFW_LB_POINT_ASSIGN)
#define Zfw_LB_Point_PP_Assign         FC_FUNC_(zfw_lb_point_pp_assign,ZFW_LB_POINT_PP_ASSIGN)
#define Zfw_LB_Box_Assign              FC_FUNC_(zfw_lb_box_assign,ZFW_LB_BOX_ASSIGN)
#define Zfw_LB_Box_PP_Assign           FC_FUNC_(zfw_lb_box_pp_assign,ZFW_LB_BOX_PP_ASSIGN)
#define Zfw_Invert_Lists               FC_FUNC_(zfw_invert_lists,ZFW_INVERT_LISTS)
#define Zfw_Compute_Destinations       FC_FUNC_(zfw_compute_destinations,ZFW_COMPUTE_DESTINATIONS)
#define Zfw_Migrate                    FC_FUNC_(zfw_migrate,ZFW_MIGRATE)
#define Zfw_Help_Migrate               FC_FUNC_(zfw_help_migrate,ZFW_HELP_MIGRATE)
#define Zfw_Order                      FC_FUNC_(zfw_order,ZFW_ORDER)
#define Zfw_Color                      FC_FUNC_(zfw_color,ZFW_COLOR)
#define Zfw_Color_Test                 FC_FUNC_(zfw_color_test,ZFW_COLOR_TEST)
#define Zfw_Generate_Files             FC_FUNC_(zfw_generate_files,ZFW_GENERATE_FILES)
#define Zfw_RCB_Box                    FC_FUNC_(zfw_rcb_box,ZFW_RCB_BOX)
#define Zfw_Register_Fort_Malloc       FC_FUNC_(zfw_register_fort_malloc,ZFW_REGISTER_FORT_MALLOC)
#define Zfw_Get_Address_int            FC_FUNC_(zfw_get_address_int,ZFW_GET_ADDRESS_INT)
#define Zfw_Get_Address_struct         FC_FUNC_(zfw_get_address_struct,ZFW_GET_ADDRESS_STRUCT)
#define Zfw_Get_Wgt_Dim                FC_FUNC_(zfw_get_wgt_dim,ZFW_GET_WGT_DIM)
#define Zfw_Get_Comm_Dim               FC_FUNC_(zfw_get_comm_dim,ZFW_GET_COMM_DIM)
#define Zfw_Reftree_Get_Child_Order    FC_FUNC_(zfw_reftree_get_child_order,ZFW_REFTREE_GET_CHILD_ORDER)

#else
/*--------------------------------------------------------------------*/
/* procedure name mangling must be done manually                      */

#if defined(FC_FN_LOWER) && defined(FC_FN_NO_UNDER)

#define Zfw_Initialize                 zfw_initialize
#define Zfw_Initialize1                zfw_initialize1
#define Zfw_Create                     zfw_create       
#define Zfw_Copy                       zfw_copy
#define Zfw_Copy_To                    zfw_copy_to
#define Zfw_Destroy                    zfw_destroy       
#define Zfw_Align                      zfw_align       
#define Zfw_Memory_Stats               zfw_memory_stats       
#define Zfw_Set_Fn0f                   zfw_set_fn0f
#define Zfw_Set_Fn1f                   zfw_set_fn1f
#define Zfw_Set_Fn2f                   zfw_set_fn2f
#define Zfw_Set_Fn3f                   zfw_set_fn3f
#define Zfw_Set_Fn4f                   zfw_set_fn4f
#define Zfw_Set_Fn5f                   zfw_set_fn5f
#define Zfw_Set_Fn6f                   zfw_set_fn6f
#define Zfw_Set_Fn7f                   zfw_set_fn7f
#define Zfw_Set_Fn8f                   zfw_set_fn8f
#define Zfw_Set_Fn9f                   zfw_set_fn9f
#define Zfw_Set_FnAf                   zfw_set_fnaf
#define Zfw_Set_FnBf                   zfw_set_fnbf
#define Zfw_Set_Fn0s                   zfw_set_fn0s
#define Zfw_Set_Fn1s                   zfw_set_fn1s
#define Zfw_Set_Fn2s                   zfw_set_fn2s
#define Zfw_Set_Fn3s                   zfw_set_fn3s
#define Zfw_Set_Fn4s                   zfw_set_fn4s
#define Zfw_Set_Fn5s                   zfw_set_fn5s
#define Zfw_Set_Fn6s                   zfw_set_fn6s
#define Zfw_Set_Fn7s                   zfw_set_fn7s
#define Zfw_Set_Fn8s                   zfw_set_fn8s
#define Zfw_Set_Fn9s                   zfw_set_fn9s
#define Zfw_Set_FnAs                   zfw_set_fnas
#define Zfw_Set_FnBs                   zfw_set_fnbs
#define Zfw_Set_Param                  zfw_set_param
#define Zfw_Set_Param_Vec              zfw_set_param_vec
#define Zfw_LB_Partition               zfw_lb_partition
#define Zfw_LB_Eval                    zfw_lb_eval
#define Zfw_LB_Set_Part_Sizes          zfw_lb_set_part_sizes
#define Zfw_LB_Point_Assign            zfw_lb_point_assign
#define Zfw_LB_Point_PP_Assign         zfw_lb_point_pp_assign
#define Zfw_LB_Box_Assign              zfw_lb_box_assign
#define Zfw_LB_Box_PP_Assign           zfw_lb_box_pp_assign
#define Zfw_Invert_Lists               zfw_invert_lists
#define Zfw_Compute_Destinations       zfw_compute_destinations
#define Zfw_Migrate                    zfw_migrate  
#define Zfw_Help_Migrate               zfw_help_migrate  
#define Zfw_Order                      zfw_order  
#define Zfw_Color                      zfw_color  
#define Zfw_Color_Test                 zfw_color_test  
#define Zfw_Generate_Files             zfw_generate_files  
#define Zfw_RCB_Box                    zfw_rcb_box  
#define Zfw_Register_Fort_Malloc       zfw_register_fort_malloc
#define Zfw_Get_Address_int            zfw_get_address_int
#define Zfw_Get_Address_struct         zfw_get_address_struct
#define Zfw_Get_Wgt_Dim                zfw_get_wgt_dim
#define Zfw_Get_Comm_Dim               zfw_get_comm_dim
#define Zfw_Reftree_Get_Child_Order    zfw_reftree_get_child_order

#elif defined(FC_FN_LOWER) && defined(FC_FN_UNDER)

#define Zfw_Initialize                 zfw_initialize_
#define Zfw_Initialize1                zfw_initialize1_
#define Zfw_Create                     zfw_create_
#define Zfw_Copy                       zfw_copy_
#define Zfw_Copy_To                    zfw_copy_to_
#define Zfw_Destroy                    zfw_destroy_
#define Zfw_Align                      zfw_align_
#define Zfw_Memory_Stats               zfw_memory_stats_
#define Zfw_Set_Fn0f                   zfw_set_fn0f_
#define Zfw_Set_Fn1f                   zfw_set_fn1f_
#define Zfw_Set_Fn2f                   zfw_set_fn2f_
#define Zfw_Set_Fn3f                   zfw_set_fn3f_
#define Zfw_Set_Fn4f                   zfw_set_fn4f_
#define Zfw_Set_Fn5f                   zfw_set_fn5f_
#define Zfw_Set_Fn6f                   zfw_set_fn6f_
#define Zfw_Set_Fn7f                   zfw_set_fn7f_
#define Zfw_Set_Fn8f                   zfw_set_fn8f_
#define Zfw_Set_Fn9f                   zfw_set_fn9f_
#define Zfw_Set_FnAf                   zfw_set_fnaf_
#define Zfw_Set_FnBf                   zfw_set_fnbf_
#define Zfw_Set_Fn0s                   zfw_set_fn0s_
#define Zfw_Set_Fn1s                   zfw_set_fn1s_
#define Zfw_Set_Fn2s                   zfw_set_fn2s_
#define Zfw_Set_Fn3s                   zfw_set_fn3s_
#define Zfw_Set_Fn4s                   zfw_set_fn4s_
#define Zfw_Set_Fn5s                   zfw_set_fn5s_
#define Zfw_Set_Fn6s                   zfw_set_fn6s_
#define Zfw_Set_Fn7s                   zfw_set_fn7s_
#define Zfw_Set_Fn8s                   zfw_set_fn8s_
#define Zfw_Set_Fn9s                   zfw_set_fn9s_
#define Zfw_Set_FnAs                   zfw_set_fnas_
#define Zfw_Set_FnBs                   zfw_set_fnbs_
#define Zfw_Set_Param                  zfw_set_param_
#define Zfw_Set_Param_Vec              zfw_set_param_vec_
#define Zfw_LB_Partition               zfw_lb_partition_
#define Zfw_LB_Eval                    zfw_lb_eval_
#define Zfw_LB_Set_Part_Sizes          zfw_lb_set_part_sizes_
#define Zfw_LB_Point_Assign            zfw_lb_point_assign_
#define Zfw_LB_Point_PP_Assign         zfw_lb_point_pp_assign_
#define Zfw_LB_Box_Assign              zfw_lb_box_assign_
#define Zfw_LB_Box_PP_Assign           zfw_lb_box_pp_assign_
#define Zfw_Invert_Lists               zfw_invert_lists_
#define Zfw_Compute_Destinations       zfw_compute_destinations_
#define Zfw_Migrate                    zfw_migrate_
#define Zfw_Help_Migrate               zfw_help_migrate_  
#define Zfw_Order                      zfw_order_  
#define Zfw_Color                      zfw_color_  
#define Zfw_Color_Test                 zfw_color_test_  
#define Zfw_Generate_Files             zfw_generate_files_ 
#define Zfw_RCB_Box                    zfw_rcb_box_
#define Zfw_Register_Fort_Malloc       zfw_register_fort_malloc_
#define Zfw_Get_Address_int            zfw_get_address_int_
#define Zfw_Get_Address_struct         zfw_get_address_struct_
#define Zfw_Get_Wgt_Dim                zfw_get_wgt_dim_
#define Zfw_Get_Comm_Dim               zfw_get_comm_dim_
#define Zfw_Reftree_Get_Child_Order    zfw_reftree_get_child_order_

#elif defined(FC_FN_LOWER) && defined(FC_FN_SECOND_UNDER)

#define Zfw_Initialize                 zfw_initialize__
#define Zfw_Initialize1                zfw_initialize1__
#define Zfw_Create                     zfw_create__
#define Zfw_Copy                       zfw_copy__
#define Zfw_Copy_To                    zfw_copy_to__
#define Zfw_Destroy                    zfw_destroy__
#define Zfw_Align                      zfw_align__
#define Zfw_Memory_Stats               zfw_memory_stats__
#define Zfw_Set_Fn0f                   zfw_set_fn0f__
#define Zfw_Set_Fn1f                   zfw_set_fn1f__
#define Zfw_Set_Fn2f                   zfw_set_fn2f__
#define Zfw_Set_Fn3f                   zfw_set_fn3f__
#define Zfw_Set_Fn4f                   zfw_set_fn4f__
#define Zfw_Set_Fn5f                   zfw_set_fn5f__
#define Zfw_Set_Fn6f                   zfw_set_fn6f__
#define Zfw_Set_Fn7f                   zfw_set_fn7f__
#define Zfw_Set_Fn8f                   zfw_set_fn8f__
#define Zfw_Set_Fn9f                   zfw_set_fn9f__
#define Zfw_Set_FnAf                   zfw_set_fnaf__
#define Zfw_Set_FnBf                   zfw_set_fnbf__
#define Zfw_Set_Fn0s                   zfw_set_fn0s__
#define Zfw_Set_Fn1s                   zfw_set_fn1s__
#define Zfw_Set_Fn2s                   zfw_set_fn2s__
#define Zfw_Set_Fn3s                   zfw_set_fn3s__
#define Zfw_Set_Fn4s                   zfw_set_fn4s__
#define Zfw_Set_Fn5s                   zfw_set_fn5s__
#define Zfw_Set_Fn6s                   zfw_set_fn6s__
#define Zfw_Set_Fn7s                   zfw_set_fn7s__
#define Zfw_Set_Fn8s                   zfw_set_fn8s__
#define Zfw_Set_Fn9s                   zfw_set_fn9s__
#define Zfw_Set_FnAs                   zfw_set_fnas__
#define Zfw_Set_FnBs                   zfw_set_fnbs__
#define Zfw_Set_Param                  zfw_set_param__
#define Zfw_Set_Param_Vec              zfw_set_param_vec__
#define Zfw_LB_Partition               zfw_lb_partition__
#define Zfw_LB_Eval                    zfw_lb_eval__
#define Zfw_LB_Set_Part_Sizes          zfw_lb_set_part_sizes__
#define Zfw_LB_Point_Assign            zfw_lb_point_assign__
#define Zfw_LB_Point_PP_Assign         zfw_lb_point_pp_assign__
#define Zfw_LB_Box_Assign              zfw_lb_box_assign__
#define Zfw_LB_Box_PP_Assign           zfw_lb_box_pp_assign__
#define Zfw_Invert_Lists               zfw_invert_lists__
#define Zfw_Compute_Destinations       zfw_compute_destinations__
#define Zfw_Migrate                    zfw_migrate__
#define Zfw_Help_Migrate               zfw_help_migrate__
#define Zfw_Order                      zfw_order__
#define Zfw_Color                      zfw_color__
#define Zfw_Color_Test                 zfw_color_test__
#define Zfw_Generate_Files             zfw_generate_files__
#define Zfw_RCB_Box                    zfw_rcb_box__
#define Zfw_Register_Fort_Malloc       zfw_register_fort_malloc__
#define Zfw_Get_Address_int            zfw_get_address_int__
#define Zfw_Get_Address_struct         zfw_get_address_struct__
#define Zfw_Get_Wgt_Dim                zfw_get_wgt_dim__
#define Zfw_Get_Comm_Dim               zfw_get_comm_dim__
#define Zfw_Reftree_Get_Child_Order    zfw_reftree_get_child_order__

#elif defined(FC_FN_UPPER) && defined(FC_FN_NO_UNDER)

#define Zfw_Initialize                 ZFW_INITIALIZE
#define Zfw_Initialize1                ZFW_INITIALIZE1
#define Zfw_Create                     ZFW_CREATE       
#define Zfw_Copy                       ZFW_COPY
#define Zfw_Copy_To                    ZFW_COPY_TO
#define Zfw_Destroy                    ZFW_DESTROY       
#define Zfw_Align                      ZFW_ALIGN
#define Zfw_Memory_Stats               ZFW_MEMORY_STATS  
#define Zfw_Set_Fn0f                   ZFW_SET_FN0F
#define Zfw_Set_Fn1f                   ZFW_SET_FN1F
#define Zfw_Set_Fn2f                   ZFW_SET_FN2F
#define Zfw_Set_Fn3f                   ZFW_SET_FN3F
#define Zfw_Set_Fn4f                   ZFW_SET_FN4F
#define Zfw_Set_Fn5f                   ZFW_SET_FN5F
#define Zfw_Set_Fn6f                   ZFW_SET_FN6F
#define Zfw_Set_Fn7f                   ZFW_SET_FN7F
#define Zfw_Set_Fn8f                   ZFW_SET_FN8F
#define Zfw_Set_Fn9f                   ZFW_SET_FN9F
#define Zfw_Set_FnAf                   ZFW_SET_FNAF
#define Zfw_Set_FnBf                   ZFW_SET_FNBF
#define Zfw_Set_Fn0s                   ZFW_SET_FN0S
#define Zfw_Set_Fn1s                   ZFW_SET_FN1S
#define Zfw_Set_Fn2s                   ZFW_SET_FN2S
#define Zfw_Set_Fn3s                   ZFW_SET_FN3S
#define Zfw_Set_Fn4s                   ZFW_SET_FN4S
#define Zfw_Set_Fn5s                   ZFW_SET_FN5S
#define Zfw_Set_Fn6s                   ZFW_SET_FN6S
#define Zfw_Set_Fn7s                   ZFW_SET_FN7S
#define Zfw_Set_Fn8s                   ZFW_SET_FN8S
#define Zfw_Set_Fn9s                   ZFW_SET_FN9S
#define Zfw_Set_FnAs                   ZFW_SET_FNAS
#define Zfw_Set_FnBs                   ZFW_SET_FNBS
#define Zfw_Set_Param                  ZFW_SET_PARAM
#define Zfw_Set_Param_Vec              ZFW_SET_PARAM_VEC
#define Zfw_LB_Partition               ZFW_LB_PARTITION
#define Zfw_LB_Eval                    ZFW_LB_EVAL
#define Zfw_LB_Set_Part_Sizes          ZFW_LB_SET_PART_SIZES
#define Zfw_LB_Point_Assign            ZFW_LB_POINT_ASSIGN
#define Zfw_LB_Point_PP_Assign         ZFW_LB_POINT_PP_ASSIGN
#define Zfw_LB_Box_Assign              ZFW_LB_BOX_ASSIGN
#define Zfw_LB_Box_PP_Assign           ZFW_LB_BOX_PP_ASSIGN
#define Zfw_Invert_Lists               ZFW_INVERT_LISTS
#define Zfw_Compute_Destinations       ZFW_COMPUTE_DESTINATIONS  
#define Zfw_Migrate                    ZFW_MIGRATE  
#define Zfw_Help_Migrate               ZFW_HELP_MIGRATE  
#define Zfw_Order                      ZFW_ORDER  
#define Zfw_Color                      ZFW_COLOR  
#define Zfw_Color_Test                 ZFW_COLOR_TEST  
#define Zfw_Generate_Files             ZFW_GENERATE_FILES  
#define Zfw_RCB_Box                    ZFW_RCB_BOX  
#define Zfw_Register_Fort_Malloc       ZFW_REGISTER_FORT_MALLOC
#define Zfw_Get_Address_int            ZFW_GET_ADDRESS_INT
#define Zfw_Get_Address_struct         ZFW_GET_ADDRESS_STRUCT
#define Zfw_Get_Comm_Dim               ZFW_GET_COMM_DIM
#define Zfw_Reftree_Get_Child_Order    ZFW_REFTREE_GET_CHILD_ORDER

#elif defined(FC_FN_UPPER) && defined(FC_FN_UNDER)

#define Zfw_Initialize                 ZFW_INITIALIZE_
#define Zfw_Initialize1                ZFW_INITIALIZE1_
#define Zfw_Create                     ZFW_CREATE_
#define Zfw_Copy                       ZFW_COPY_
#define Zfw_Copy_To                    ZFW_COPY_TO_
#define Zfw_Destroy                    ZFW_DESTROY_
#define Zfw_Align                      ZFW_ALIGN_
#define Zfw_Memory_Stats               ZFW_MEMORY_STATS_
#define Zfw_Set_Fn0f                   ZFW_SET_FN0F_
#define Zfw_Set_Fn1f                   ZFW_SET_FN1F_
#define Zfw_Set_Fn2f                   ZFW_SET_FN2F_
#define Zfw_Set_Fn3f                   ZFW_SET_FN3F_
#define Zfw_Set_Fn4f                   ZFW_SET_FN4F_
#define Zfw_Set_Fn5f                   ZFW_SET_FN5F_
#define Zfw_Set_Fn6f                   ZFW_SET_FN6F_
#define Zfw_Set_Fn7f                   ZFW_SET_FN7F_
#define Zfw_Set_Fn8f                   ZFW_SET_FN8F_
#define Zfw_Set_Fn9f                   ZFW_SET_FN9F_
#define Zfw_Set_FnAf                   ZFW_SET_FNAF_
#define Zfw_Set_FnBf                   ZFW_SET_FNBF_
#define Zfw_Set_Fn0s                   ZFW_SET_FN0S_
#define Zfw_Set_Fn1s                   ZFW_SET_FN1S_
#define Zfw_Set_Fn2s                   ZFW_SET_FN2S_
#define Zfw_Set_Fn3s                   ZFW_SET_FN3S_
#define Zfw_Set_Fn4s                   ZFW_SET_FN4S_
#define Zfw_Set_Fn5s                   ZFW_SET_FN5S_
#define Zfw_Set_Fn6s                   ZFW_SET_FN6S_
#define Zfw_Set_Fn7s                   ZFW_SET_FN7S_
#define Zfw_Set_Fn8s                   ZFW_SET_FN8S_
#define Zfw_Set_Fn9s                   ZFW_SET_FN9S_
#define Zfw_Set_FnAs                   ZFW_SET_FNAS_
#define Zfw_Set_FnBs                   ZFW_SET_FNBS_
#define Zfw_Set_Param                  ZFW_SET_PARAM_
#define Zfw_Set_Param_Vec              ZFW_SET_PARAM_VEC_
#define Zfw_LB_Partition               ZFW_LB_PARTITION_
#define Zfw_LB_Eval                    ZFW_LB_EVAL_
#define Zfw_LB_Set_Part_Sizes          ZFW_LB_SET_PART_SIZES_
#define Zfw_LB_Point_Assign            ZFW_LB_POINT_ASSIGN_
#define Zfw_LB_Point_PP_Assign         ZFW_LB_POINT_PP_ASSIGN_
#define Zfw_LB_Box_Assign              ZFW_LB_BOX_ASSIGN_
#define Zfw_LB_Box_PP_Assign           ZFW_LB_BOX_PP_ASSIGN_
#define Zfw_Invert_Lists               ZFW_INVERT_LISTS_
#define Zfw_Compute_Destinations       ZFW_COMPUTE_DESTINATIONS_
#define Zfw_Migrate                    ZFW_MIGRATE_
#define Zfw_Help_Migrate               ZFW_HELP_MIGRATE_
#define Zfw_Order                      ZFW_ORDER_
#define Zfw_Color                      ZFW_COLOR_
#define Zfw_Color_Test                 ZFW_COLOR_TEST_
#define Zfw_Generate_Files             ZFW_GENERATE_FILES_
#define Zfw_RCB_Box                    ZFW_RCB_BOX_
#define Zfw_Register_Fort_Malloc       ZFW_REGISTER_FORT_MALLOC_
#define Zfw_Get_Address_int            ZFW_GET_ADDRESS_INT_
#define Zfw_Get_Address_struct         ZFW_GET_ADDRESS_STRUCT_
#define Zfw_Get_Comm_Dim               ZFW_GET_COMM_DIM_
#define Zfw_Reftree_Get_Child_Order    ZFW_REFTREE_GET_CHILD_ORDER_

#elif defined(FC_FN_UPPER) && defined(FC_FN_SECOND_UNDER)

#define Zfw_Initialize                 ZFW_INITIALIZE__
#define Zfw_Initialize1                ZFW_INITIALIZE1__
#define Zfw_Create                     ZFW_CREATE__
#define Zfw_Copy                       ZFW_COPY__
#define Zfw_Copy_To                    ZFW_COPY_TO__
#define Zfw_Destroy                    ZFW_DESTROY__
#define Zfw_Align                      ZFW_ALIGN__
#define Zfw_Memory_Stats               ZFW_MEMORY_STATS__
#define Zfw_Set_Fn0f                   ZFW_SET_FN0F__
#define Zfw_Set_Fn1f                   ZFW_SET_FN1F__
#define Zfw_Set_Fn2f                   ZFW_SET_FN2F__
#define Zfw_Set_Fn3f                   ZFW_SET_FN3F__
#define Zfw_Set_Fn4f                   ZFW_SET_FN4F__
#define Zfw_Set_Fn5f                   ZFW_SET_FN5F__
#define Zfw_Set_Fn6f                   ZFW_SET_FN6F__
#define Zfw_Set_Fn7f                   ZFW_SET_FN7F__
#define Zfw_Set_Fn8f                   ZFW_SET_FN8F__
#define Zfw_Set_Fn9f                   ZFW_SET_FN9F__
#define Zfw_Set_FnAf                   ZFW_SET_FNAF__
#define Zfw_Set_FnBf                   ZFW_SET_FNBF__
#define Zfw_Set_Fn0s                   ZFW_SET_FN0S__
#define Zfw_Set_Fn1s                   ZFW_SET_FN1S__
#define Zfw_Set_Fn2s                   ZFW_SET_FN2S__
#define Zfw_Set_Fn3s                   ZFW_SET_FN3S__
#define Zfw_Set_Fn4s                   ZFW_SET_FN4S__
#define Zfw_Set_Fn5s                   ZFW_SET_FN5S__
#define Zfw_Set_Fn6s                   ZFW_SET_FN6S__
#define Zfw_Set_Fn7s                   ZFW_SET_FN7S__
#define Zfw_Set_Fn8s                   ZFW_SET_FN8S__
#define Zfw_Set_Fn9s                   ZFW_SET_FN9S__
#define Zfw_Set_FnAs                   ZFW_SET_FNAS__
#define Zfw_Set_FnBs                   ZFW_SET_FNBS__
#define Zfw_Set_Param                  ZFW_SET_PARAM__
#define Zfw_Set_Param_Vec              ZFW_SET_PARAM_VEC__
#define Zfw_LB_Partition               ZFW_LB_PARTITION__
#define Zfw_LB_Eval                    ZFW_LB_EVAL__
#define Zfw_LB_Set_Part_Sizes          ZFW_LB_SET_PART_SIZES__
#define Zfw_LB_Point_Assign            ZFW_LB_POINT_ASSIGN__
#define Zfw_LB_Point_PP_Assign         ZFW_LB_POINT_PP_ASSIGN__
#define Zfw_LB_Box_Assign              ZFW_LB_BOX_ASSIGN__
#define Zfw_LB_Box_PP_Assign           ZFW_LB_BOX_PP_ASSIGN__
#define Zfw_Invert_Lists               ZFW_INVERT_LISTS__
#define Zfw_Compute_Destinations       ZFW_COMPUTE_DESTINATIONS__
#define Zfw_Migrate                    ZFW_MIGRATE__
#define Zfw_Help_Migrate               ZFW_HELP_MIGRATE__
#define Zfw_Order                      ZFW_ORDER__
#define Zfw_Color                      ZFW_COLOR__
#define Zfw_Color_Test                 ZFW_COLOR_TEST__
#define Zfw_Generate_Files             ZFW_GENERATE_FILES__
#define Zfw_RCB_Box                    ZFW_RCB_BOX__
#define Zfw_Register_Fort_Malloc       ZFW_REGISTER_FORT_MALLOC__
#define Zfw_Get_Address_int            ZFW_GET_ADDRESS_INT__
#define Zfw_Get_Address_struct         ZFW_GET_ADDRESS_STRUCT__
#define Zfw_Get_Comm_Dim               ZFW_GET_COMM_DIM__
#define Zfw_Reftree_Get_Child_Order    ZFW_REFTREE_GET_CHILD_ORDER__

#else
#error "Unrecognized Fortran Mangling scheme."
#endif

#endif

#endif
