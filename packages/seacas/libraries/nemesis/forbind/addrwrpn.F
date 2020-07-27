C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      subroutine adrdumn
      end

#ifdef Build64
C-----------------------------------------------------------------------

C   Get initial information from nemesis file

      subroutine negii (idne, nproc, nproc_in_f, ftype, ierr)
      implicit none
      integer idne
      integer nproc
      integer nproc_in_f
      character*(*) ftype
      integer ierr

      integer*4 idne4
      integer*4 nproc4
      integer*4 nproc_in_f4
      integer*4 ierr4

      idne4 = idne
      call negii4 (idne4, nproc4, nproc_in_f4, ftype, ierr4)
      nproc = nproc4
      nproc_in_f = nproc_in_f4
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C   Write initial information from nemesis file

      subroutine nepii (idne, nproc, nproc_in_f, ftype, ierr)
      implicit none
      integer idne
      integer nproc
      integer nproc_in_f
      character*(*) ftype
      integer ierr

      integer*4 idne4
      integer*4 nproc4
      integer*4 nproc_in_f4
      integer*4 ierr4

      idne4 = idne
      nproc4 = nproc
      nproc_in_f4 = nproc_in_f
      call nepii4 (idne4, nproc4, nproc_in_f4, ftype, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read initial global information

      subroutine negig (idne, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr)
      implicit none
      integer idne
      integer nnodes_g
      integer nelems_g
      integer nelem_blks_g
      integer nnode_sets_g
      integer nside_sets_g
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negig4 (idne4, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write initial global information

      subroutine nepig (idne, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr)
      implicit none
      integer idne
      integer nnodes_g
      integer nelems_g
      integer nelem_blks_g
      integer nnode_sets_g
      integer nside_sets_g
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepig4 (idne4, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read load balance parameters

      subroutine neglbp (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
      implicit none
      integer idne
      integer nint_nodes
      integer nbor_nodes
      integer next_nodes
      integer nint_elems
      integer nbor_elems
      integer nnode_cmaps
      integer nelem_cmaps
      integer processor
      integer ierr

      integer*4 idne4
      integer*4 processor4
      integer*4 ierr4

      idne4 = idne
      processor4 = processor
      call neglbp4 (idne4, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write load balance parameters

      subroutine neplbp (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
      implicit none
      integer idne
      integer nint_nodes
      integer nbor_nodes
      integer next_nodes
      integer nint_elems
      integer nbor_elems
      integer nnode_cmaps
      integer nelem_cmaps
      integer processor
      integer ierr

      integer*4 idne4
      integer*4 processor4
      integer*4 ierr4

      idne4 = idne
      processor4 = processor
      call neplbp4 (idne4, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write concatenated load balance parameters

      subroutine neplbpc (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps, ierr)
      implicit none
      integer idne
      integer nint_nodes(*)
      integer nbor_nodes(*)
      integer next_nodes(*)
      integer nint_elems(*)
      integer nbor_elems(*)
      integer nnode_cmaps(*)
      integer nelem_cmaps(*)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call neplbpc4 (idne4, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read global node set parameters

      subroutine negnspg (idne, ns_ids_glob, ns_n_cnt_glob,
     $  ns_df_cnt_glob, ierr)
      implicit none
      integer idne
      integer ns_ids_glob(*)
      integer ns_n_cnt_glob(*)
      integer ns_df_cnt_glob(*)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negnspg4 (idne4, ns_ids_glob, ns_n_cnt_glob,
     $  ns_df_cnt_glob, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write global node set parameters

      subroutine nepnspg (idne, global_ids, global_n_cnts,
     $  global_df_cnts, ierr)
      integer idne
      integer global_ids(*)
      integer global_n_cnts(*)
      integer global_df_cnts(*)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepnspg4 (idne4, global_ids, global_n_cnts,
     $  global_df_cnts, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read global side set parameters

      subroutine negsspg (idne, ss_ids_glob, ss_n_cnt_glob,
     $  ss_df_cnt_glob, ierr)
      implicit none
      integer idne
      integer ss_ids_glob(*)
      integer ss_n_cnt_glob(*)
      integer ss_df_cnt_glob(*)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negsspg4 (idne4, ss_ids_glob, ss_n_cnt_glob,
     $  ss_df_cnt_glob, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write global side set parameters

      subroutine nepsspg (idne, global_ids, global_el_cnts,
     $  global_df_cnts, ierr)
      implicit none
      integer idne
      integer global_ids(*)
      integer global_el_cnts(*)
      integer global_df_cnts(*)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepsspg4 (idne4, global_ids, global_el_cnts,
     $  global_df_cnts, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read global element block information

      subroutine negebig (idne, el_blk_ids, el_blk_cnts, ierr)
      implicit none
      integer idne
      integer el_blk_ids(*)
      integer el_blk_cnts(*)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negebig4 (idne4, el_blk_ids, el_blk_cnts, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write global element block information

      subroutine nepebig (idne, el_blk_ids, el_blk_cnts, ierr)
      implicit none
      integer idne
      integer el_blk_ids(*)
      integer el_blk_cnts(*)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepebig4 (idne4, el_blk_ids, el_blk_cnts, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read side set element list and side set side list

      subroutine negnss (idne, ss_id, start_side_num, num_sides,
     $  ss_elem_list, ss_side_list, ierr)
      implicit none
      integer idne
      integer ss_id
      integer start_side_num
      integer num_sides
      integer ss_elem_list(num_sides)
      integer ss_side_list(num_sides)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negnss4 (idne4, ss_id, start_side_num, num_sides,
     $  ss_elem_list, ss_side_list, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write side set element list and side set side list

      subroutine nepnss (idne, ss_id, start_side_num, num_sides,
     $  ss_elem_list, ss_side_list, ierr)
      implicit none
      integer idne
      integer ss_id
      integer start_side_num
      integer num_sides
      integer ss_elem_list(num_sides)
      integer ss_side_list(num_sides)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepnss4 (idne4, ss_id, start_side_num, num_sides,
     $  ss_elem_list, ss_side_list, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read side set distribution factor

      subroutine negnssd (idne, ss_id, start_num, num_df_to_get,
     $  ss_df, ierr)
      implicit none
      integer idne
      integer ss_id
      integer start_num
      integer num_df_to_get
      real ss_df(num_df_to_get)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negnssd4 (idne4, ss_id, start_num, num_df_to_get,
     $  ss_df, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write side set distribution factor

      subroutine nepnssd (idne, ss_id, start_num, num_df_to_get,
     $  ss_df, ierr)
      implicit none
      integer idne
      integer ss_id
      integer start_num
      integer num_df_to_get
      real ss_df(num_df_to_get)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepnssd4 (idne4, ss_id, start_num, num_df_to_get,
     $  ss_df, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read node set list for a single node set

      subroutine negnns (idne, ns_id, start_node_num, num_node,
     $  ns_node_list, ierr)
      implicit none
      integer idne
      integer ns_id
      integer start_node_num
      integer num_node
      integer ns_node_list(num_node)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negnns4 (idne4, ns_id, start_node_num, num_node,
     $  ns_node_list, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write node set list for a single node set

      subroutine nepnns (idne, ns_id, start_node_num, num_node,
     $  ns_node_list, ierr)
      implicit none
      integer idne
      integer ns_id
      integer start_node_num
      integer num_node
      integer ns_node_list(num_node)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepnns4 (idne4, ns_id, start_node_num, num_node,
     $  ns_node_list, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read node set distribution factor

      subroutine negnnsd (idne, ns_id, start_num, num_df_to_get,
     $  ns_df, ierr)
      implicit none
      integer idne
      integer ns_id
      integer start_num
      integer num_df_to_get
      real ns_df(num_df_to_get)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negnnsd4 (idne4, ns_id, start_num, num_df_to_get,
     $  ns_df, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write node set distribution factor

      subroutine nepnnsd (idne, ns_id, start_num, num_df_to_get,
     $  ns_df, ierr)
      implicit none
      integer idne
      integer ns_id
      integer start_num
      integer num_df_to_get
      real ns_df(num_df_to_get)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepnnsd4 (idne4, ns_id, start_num, num_df_to_get,
     $  ns_df, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read coordinates of the nodes

      subroutine negcor (idne, start_node_num, num_nodes, x_coor,
     $  y_coor, z_coor, ierr)
      implicit none
      integer idne
      integer start_node_num
      integer num_nodes
      real x_coor(num_nodes)
      real y_coor(num_nodes)
      real z_coor(num_nodes)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negcor4 (idne4, start_node_num, num_nodes, x_coor,
     $  y_coor, z_coor, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write coordinates of the nodes

      subroutine nepcor (idne, start_node_num, num_nodes, x_coor,
     $  y_coor, z_coor, ierr)
      implicit none
      integer idne
      integer start_node_num
      integer num_nodes
      real x_coor(num_nodes)
      real y_coor(num_nodes)
      real z_coor(num_nodes)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepcor4 (idne4, start_node_num, num_nodes, x_coor,
     $  y_coor, z_coor, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read an element blocks connectivity list

      subroutine negnec (idne, elem_blk_id, start_elem_num,
     $  num_elems, connect, ierr)
      implicit none
      integer idne
      integer elem_blk_id
      integer start_elem_num
      integer num_elems
      integer connect(num_elems)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negnec4 (idne4, elem_blk_id, start_elem_num,
     $  num_elems, connect, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write an element blocks connectivity list

      subroutine nepnec (idne, elem_blk_id, start_elem_num,
     $  num_elems, connect, ierr)
      implicit none
      integer idne
      integer elem_blk_id
      integer start_elem_num
      integer num_elems
      integer connect(num_elems)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepnec4 (idne4, elem_blk_id, start_elem_num,
     $  num_elems, connect, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read an element blocks attributes

      subroutine negneat (idne, elem_blk_id, start_elem_num,
     $  num_elems, attrib, ierr)
      implicit none
      integer idne
      integer elem_blk_id
      integer start_elem_num
      integer num_elems
      real attrib(num_elems)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negneat4 (idne4, elem_blk_id, start_elem_num,
     $  num_elems, attrib, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write an element blocks attributes

      subroutine nepneat (idne, elem_blk_id, start_elem_num,
     $  num_elems, attrib, ierr)
      implicit none
      integer idne
      integer elem_blk_id
      integer start_elem_num
      integer num_elems
      real attrib(num_elems)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepneat4 (idne4, elem_blk_id, start_elem_num,
     $  num_elems, attrib, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------
C  Read the element type for a specific element block

      subroutine negelt (idne, elem_blk_id, elem_type, ierr)
      implicit none
      integer idne
      integer elem_blk_id
      character*(*) elem_type
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negelt4 (idne4, elem_blk_id, elem_type, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read a variable for an element block

      subroutine negnev (idne, time_step, elem_var_index,
     $  elem_blk_id, num_elem_this_blk, start_elem_num,
     $  num_elem, elem_var_vals, ierr)
      implicit none
      integer idne
      integer time_step
      integer elem_var_index
      integer elem_blk_id
      integer num_elem_this_blk
      integer start_elem_num
      integer num_elem
      real elem_var_vals(num_elem)
      integer ierr

      integer*4 idne4
      integer*4 time_step4
      integer*4 elem_var_index4
      integer*4 ierr4

      idne4 = idne
      time_step4 = time_step
      elem_var_index4 = elem_var_index
      call negnev4 (idne4, time_step4, elem_var_index4,
     $  elem_blk_id, num_elem_this_blk, start_elem_num,
     $  num_elem, elem_var_vals, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write a variable slab for an element block

      subroutine nepevs (idne, time_step, elem_var_index,
     $  elem_blk_id, start_pos, num_vals, elem_var_vals,
     $  ierr)
      implicit none
      integer idne
      integer time_step
      integer elem_var_index
      integer elem_blk_id
      integer start_pos
      integer num_vals
      real elem_var_vals(num_vals)
      integer ierr

      integer*4 idne4
      integer*4 time_step4
      integer*4 elem_var_index4
      integer*4 ierr4

      idne4 = idne
      time_step4 = time_step
      elem_var_index4 = elem_var_index
      call nepevs4 (idne4, time_step4, elem_var_index4,
     $  elem_blk_id, start_pos, num_vals, elem_var_vals,
     $  ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the values of a single nodal variable for a single time step

      subroutine negnnv (idne, time_step, nodal_var_index,
     $  start_node_num, num_nodes, nodal_vars, ierr)
      implicit none
      integer idne
      integer time_step
      integer nodal_var_index
      integer start_node_num
      integer num_nodes
      real nodal_vars(num_nodes)
      integer ierr

      integer*4 idne4
      integer*4 time_step4
      integer*4 nodal_var_index4
      integer*4 ierr4

      idne4 = idne
      time_step4 = time_step
      nodal_var_index4 = nodal_var_index
      call negnnv4 (idne4, time_step4, nodal_var_index4,
     $  start_node_num, num_nodes, nodal_vars, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write nodal variable slab

      subroutine nepnvs (idne, time_step, nodal_var_index,
     $  start_pos, num_vals, nodal_var_vals, ierr)
      implicit none
      integer idne
      integer time_step
      integer nodal_var_index
      integer start_pos
      integer num_vals
      real nodal_var_vals(num_vals)
      integer ierr

      integer*4 idne4
      integer*4 time_step4
      integer*4 nodal_var_index4
      integer*4 ierr4

      idne4 = idne
      time_step4 = time_step
      nodal_var_index4 = nodal_var_index
      call nepnvs4 (idne4, time_step4, nodal_var_index4,
     $  start_pos, num_vals, nodal_var_vals, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the element numbering map

      subroutine negnenm (idne, starte, num_ent, elem_map, ierr)
      implicit none
      integer idne
      integer starte
      integer num_ent
      integer elem_map(num_ent)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negnenm4 (idne4, starte, num_ent, elem_map, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the element numbering map

      subroutine nepnenm (idne, starte, num_ent, elem_map, ierr)
      implicit none
      integer idne
      integer starte
      integer num_ent
      integer elem_map(num_ent)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepnenm4 (idne4, starte, num_ent, elem_map, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the node numbering map

      subroutine negnnnm (idne, startn, num_ent, node_map, ierr)
      implicit none
      integer idne
      integer startn
      integer num_ent
      integer node_map(num_ent)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call negnnnm4 (idne4, startn, num_ent, node_map, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the node numbering map

      subroutine nepnnnm (idne, startn, num_ent, node_map, ierr)
      implicit none
      integer idne
      integer startn
      integer num_ent
      integer node_map(num_ent)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepnnnm4 (idne4, startn, num_ent, node_map, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the node map for a processor

      subroutine negnm (idne, node_mapi, node_mapb, node_mape,
     $  processor, ierr)
      implicit none
      integer idne
      integer node_mapi(*)
      integer node_mapb(*)
      integer node_mape(*)
      integer processor
      integer ierr

      integer*4 idne4
      integer*4 processor4
      integer*4 ierr4

      idne4 = idne
      processor4 = processor
      call negnm4 (idne4, node_mapi, node_mapb, node_mape,
     $  processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write a node map for a processor

      subroutine nepnm (idne, node_mapi, node_mapb, node_mape,
     $  processor, ierr)
      implicit none
      integer idne
      integer node_mapi(*)
      integer node_mapb(*)
      integer node_mape(*)
      integer processor
      integer ierr

      integer*4 idne4
      integer*4 processor4
      integer*4 ierr4

      idne4 = idne
      processor4 = processor
      call nepnm4 (idne4, node_mapi, node_mapb, node_mape,
     $  processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the element map for a processor

      subroutine negem (idne, elem_mapi, elem_mapb, processor, ierr)
      implicit none
      integer idne
      integer elem_mapi(*)
      integer elem_mapb(*)
      integer processor
      integer ierr

      integer*4 idne4
      integer*4 processor4
      integer*4 ierr4

      idne4 = idne
      processor4 = processor
      call negem4 (idne4, elem_mapi, elem_mapb, processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the element map for a processor

      subroutine nepem (idne, elem_mapi, elem_mapb, processor, ierr)
      implicit none
      integer idne
      integer elem_mapi(*)
      integer elem_mapb(*)
      integer processor
      integer ierr

      integer*4 idne4
      integer*4 processor4
      integer*4 ierr4

      idne4 = idne
      processor4 = processor
      call nepem4 (idne4, elem_mapi, elem_mapb, processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the communications map parameters for a single processor

      subroutine negcmp (idne, ncmap_ids, ncmap_node_cnts,
     $  ecmap_ids, ecmap_elem_cnts, processor, ierr)
      implicit none
      integer idne
      integer ncmap_ids(*)
      integer ncmap_node_cnts(*)
      integer ecmap_ids(*)
      integer ecmap_elem_cnts(*)
      integer processor
      integer ierr

      integer*4 idne4
      integer*4 processor4
      integer*4 ierr4

      idne4 = idne
      processor4 = processor
      call negcmp4 (idne4, ncmap_ids, ncmap_node_cnts,
     $  ecmap_ids, ecmap_elem_cnts, processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the communications map parameters for a single processor

      subroutine nepcmp (idne, nmap_ids, nmap_node_cnts, emap_ids,
     $  emap_elem_cnts, processor, ierr)
      implicit none
      integer idne
      integer nmap_ids(*)
      integer nmap_node_cnts(*)
      integer emap_ids(*)
      integer emap_elem_cnts(*)
      integer processor
      integer ierr

      integer*4 idne4
      integer*4 processor4
      integer*4 ierr4

      idne4 = idne
      processor4 = processor
      call nepcmp4 (idne4, nmap_ids, nmap_node_cnts, emap_ids,
     $  emap_elem_cnts, processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the communications map parameters for all processors

      subroutine nepcmpc (idne, nmap_ids, nmap_node_cnts,
     $  nproc_ptrs, emap_ids, emap_elem_cnts, eproc_ptrs, ierr)
      implicit none
      integer idne
      integer nmap_ids(*)
      integer nmap_node_cnts(*)
      integer nproc_ptrs(*)
      integer emap_ids(*)
      integer emap_elem_cnts(*)
      integer eproc_ptrs(*)
      integer ierr

      integer*4 idne4
      integer*4 ierr4

      idne4 = idne
      call nepcmpc4 (idne4, nmap_ids, nmap_node_cnts,
     $  nproc_ptrs, emap_ids, emap_elem_cnts, eproc_ptrs, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the nodal communications map for a single processor

      subroutine negncm (idne, map_id, node_ids, proc_ids,
     $  processor, ierr)
      implicit none
      integer idne
      integer map_id
      integer node_ids(*)
      integer proc_ids(*)
      integer processor
      integer ierr

      integer*4 idne4
      integer*4 processor4
      integer*4 ierr4

      idne4 = idne
      processor4 = processor
      call negncm4 (idne4, map_id, node_ids, proc_ids,
     $  processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the nodal communications map for a single processor

      subroutine nepncm (idne, map_id, node_ids, proc_ids,
     $  processor, ierr)
      implicit none
      integer idne
      integer map_id
      integer node_ids(*)
      integer proc_ids(*)
      integer processor
      integer ierr

      integer*4 idne4
      integer*4 processor4
      integer*4 ierr4

      idne4 = idne
      processor4 = processor
      call nepncm4 (idne4, map_id, node_ids, proc_ids,
     $  processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Read the elemental communications map for a single processor

      subroutine negecm (idne, map_id, elem_ids, side_ids,
     $  proc_ids, processor, ierr)
      implicit none
      integer idne
      integer map_id
      integer elem_ids(*)
      integer side_ids(*)
      integer proc_ids(*)
      integer processor
      integer ierr

      integer*4 idne4
      integer*4 processor4
      integer*4 ierr4

      idne4 = idne
      processor4 = processor
      call negecm4 (idne4, map_id, elem_ids, side_ids,
     $  proc_ids, processor4, ierr4)
      ierr = ierr4

      end

C-----------------------------------------------------------------------

C  Write the elemental communications map for a single processor

      subroutine nepecm (idne, map_id, elem_ids, side_ids,
     $  proc_ids, processor, ierr)
      implicit none
      integer idne
      integer map_id
      integer elem_ids(*)
      integer side_ids(*)
      integer proc_ids(*)
      integer processor
      integer ierr

      integer*4 idne4
      integer*4 processor4
      integer*4 ierr4

      idne4 = idne
      processor4 = processor
      call nepecm4 (idne4, map_id, elem_ids, side_ids,
     $  proc_ids, processor4, ierr4)
      ierr = ierr4

         end
C-----------------------------------------------------------------------
#endif
