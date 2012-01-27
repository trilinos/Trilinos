C Copyright (c) 1998 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.  
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

      subroutine adrdumn
      end
      
#ifdef Build64  
C-----------------------------------------------------------------------
C
C   Get initial information from nemesis file
C
      subroutine negii (idne, nproc, nproc_in_f, ftype, ierr)
      implicit none
      integer idne 
      integer nproc 
      integer nproc_in_f 
      character*(*) ftype 
      integer ierr
C
      integer*4 idne4 
      integer*4 nproc4 
      integer*4 nproc_in_f4 
      integer*4 ierr4
C
      idne4 = idne
      call negii4 (idne4, nproc4, nproc_in_f4, ftype, ierr4)
      nproc = nproc4
      nproc_in_f = nproc_in_f4
      ierr = ierr4
C       
      end
      
C-----------------------------------------------------------------------
C
C   Write initial information from nemesis file
C
      subroutine nepii (idne, nproc, nproc_in_f, ftype, ierr)
      implicit none
      integer idne 
      integer nproc 
      integer nproc_in_f 
      character*(*) ftype 
      integer ierr
C
      integer*4 idne4 
      integer*4 nproc4 
      integer*4 nproc_in_f4 
      integer*4 ierr4
C
      idne4 = idne
      nproc4 = nproc
      nproc_in_f4 = nproc_in_f
      call nepii4 (idne4, nproc4, nproc_in_f4, ftype, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read initial global information
C
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
C
      integer*4 idne4 
      integer*4 nnodes_g4 
      integer*4 nelems_g4 
      integer*4 nelem_blks_g4 
      integer*4 nnode_sets_g4 
      integer*4 nside_sets_g4 
      integer*4 ierr4
C
      idne4 = idne
      call negig4 (idne4, nnodes_g4, nelems_g4, nelem_blks_g4,
     $  nnode_sets_g4, nside_sets_g4, ierr4)
      nnodes_g = nnodes_g4
      nelems_g = nelems_g4
      nelem_blks_g = nelem_blks_g4
      nnode_sets_g = nnode_sets_g4
      nside_sets_g = nside_sets_g4
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write initial global information
C
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
C
      integer*4 idne4 
      integer*4 nnodes_g4 
      integer*4 nelems_g4 
      integer*4 nelem_blks_g4 
      integer*4 nnode_sets_g4 
      integer*4 nside_sets_g4 
      integer*4 ierr4
C
      idne4 = idne
      nnodes_g4 = nnodes_g
      nelems_g4 = nelems_g
      nelem_blks_g4 = nelem_blks_g
      nnode_sets_g4 = nnode_sets_g
      nside_sets_g4 = nside_sets_g
      call nepig4 (idne4, nnodes_g4, nelems_g4, nelem_blks_g4,
     $  nnode_sets_g4, nside_sets_g4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read load balance parameters
C
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
C
      integer*4 idne4 
      integer*4 nint_nodes4 
      integer*4 nbor_nodes4 
      integer*4 next_nodes4 
      integer*4 nint_elems4 
      integer*4 nbor_elems4 
      integer*4 nnode_cmaps4 
      integer*4 nelem_cmaps4 
      integer*4 processor4 
      integer*4 ierr4
C
      idne4 = idne
      processor4 = processor
      call neglbp4 (idne4, nint_nodes4, nbor_nodes4, next_nodes4,
     $  nint_elems4, nbor_elems4, nnode_cmaps4, nelem_cmaps4,
     $  processor4, ierr4)
      nint_nodes = nint_nodes4
      nbor_nodes = nbor_nodes4
      next_nodes = next_nodes4
      nint_elems = nint_elems4
      nbor_elems = nbor_elems4
      nnode_cmaps = nnode_cmaps4
      nelem_cmaps = nelem_cmaps4
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write load balance parameters
C
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
C
      integer*4 idne4 
      integer*4 nint_nodes4 
      integer*4 nbor_nodes4 
      integer*4 next_nodes4 
      integer*4 nint_elems4 
      integer*4 nbor_elems4 
      integer*4 nnode_cmaps4 
      integer*4 nelem_cmaps4 
      integer*4 processor4 
      integer*4 ierr4
C
      idne4 = idne
      nint_nodes4 = nint_nodes
      nbor_nodes4 = nbor_nodes
      next_nodes4 = next_nodes
      nint_elems4 = nint_elems
      nbor_elems4 = nbor_elems
      nnode_cmaps4 = nnode_cmaps
      nelem_cmaps4 = nelem_cmaps
      processor4 = processor
      call neplbp4 (idne4, nint_nodes4, nbor_nodes4, next_nodes4,
     $  nint_elems4, nbor_elems4, nnode_cmaps4, nelem_cmaps4,
     $  processor4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write concatenated load balance parameters
C
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
C
      integer nproc 
      integer nproc_in_f 
      character*(80) ftype 
C
      call negii (idne, nproc, nproc_in_f, ftype, ierr)
      if (ierr.ne.0) return
      call neplbpc2 (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps, ierr,
     $  nproc_in_f)
C
      end
C
      subroutine neplbpc2 (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps, ierr,
     $  nproc_in_f)
      implicit none
      integer idne 
      integer nproc_in_f 
      integer nint_nodes(nproc_in_f) 
      integer nbor_nodes(nproc_in_f) 
      integer next_nodes(nproc_in_f) 
      integer nint_elems(nproc_in_f) 
      integer nbor_elems(nproc_in_f) 
      integer nnode_cmaps(nproc_in_f) 
      integer nelem_cmaps(nproc_in_f) 
      integer ierr
C
      integer*4 idne4 
      integer*4 nint_nodes4(nproc_in_f) 
      integer*4 nbor_nodes4(nproc_in_f) 
      integer*4 next_nodes4(nproc_in_f) 
      integer*4 nint_elems4(nproc_in_f) 
      integer*4 nbor_elems4(nproc_in_f) 
      integer*4 nnode_cmaps4(nproc_in_f) 
      integer*4 nelem_cmaps4(nproc_in_f) 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      DO J1X = 1, NPROC_IN_F
        NINT_NODES4(J1X) = NINT_NODES(J1X)
        NBOR_NODES4(J1X) = NBOR_NODES(J1X)
        NEXT_NODES4(J1X) = NEXT_NODES(J1X)
        NINT_ELEMS4(J1X) = NINT_ELEMS(J1X)
        NBOR_ELEMS4(J1X) = NBOR_ELEMS(J1X)
        NNODE_CMAPS4(J1X) = NNODE_CMAPS(J1X)
        NELEM_CMAPS4(J1X) = NELEM_CMAPS(J1X)
      END DO
      call neplbpc4 (idne4, nint_nodes4, nbor_nodes4, next_nodes4,
     $  nint_elems4, nbor_elems4, nnode_cmaps4, nelem_cmaps4,
     $  ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read global node set parameters
C
      subroutine negnspg (idne, ns_ids_glob, ns_n_cnt_glob,
     $  ns_df_cnt_glob, ierr)
      implicit none
      integer idne 
      integer ns_ids_glob(*) 
      integer ns_n_cnt_glob(*) 
      integer ns_df_cnt_glob(*) 
      integer ierr
C
      integer nnodes_g 
      integer nelems_g 
      integer nelem_blks_g 
      integer nnode_sets_g 
      integer nside_sets_g 
      call negig (idne, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr)
      if (ierr.ne.0) return
      call negnspg2 (idne, ns_ids_glob, ns_n_cnt_glob,
     $  ns_df_cnt_glob, ierr, nnode_sets_g)
C
      end
C
      subroutine negnspg2 (idne, ns_ids_glob, ns_n_cnt_glob,
     $  ns_df_cnt_glob, ierr, nnode_sets_g)
      implicit none
      integer idne 
      integer nnode_sets_g 
      integer ns_ids_glob(nnode_sets_g) 
      integer ns_n_cnt_glob(nnode_sets_g) 
      integer ns_df_cnt_glob(nnode_sets_g) 
      integer ierr
C
      integer*4 idne4 
      integer*4 ns_ids_glob4(nnode_sets_g) 
      integer*4 ns_n_cnt_glob4(nnode_sets_g) 
      integer*4 ns_df_cnt_glob4(nnode_sets_g) 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      call negnspg4 (idne4, ns_ids_glob4, ns_n_cnt_glob4,
     $  ns_df_cnt_glob4, ierr4)
      DO J1X = 1, NNODE_SETS_G
        NS_IDS_GLOB(J1X) = NS_IDS_GLOB4(J1X)
        NS_N_CNT_GLOB(J1X) = NS_N_CNT_GLOB4(J1X)
        NS_DF_CNT_GLOB(J1X) = NS_DF_CNT_GLOB4(J1X)
      END DO
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write global node set parameters
C
      subroutine nepnspg (idne, global_ids, global_n_cnts,
     $  global_df_cnts, ierr)
      integer idne 
      integer global_ids(*) 
      integer global_n_cnts(*) 
      integer global_df_cnts(*) 
      integer ierr
C
      integer nnodes_g 
      integer nelems_g 
      integer nelem_blks_g 
      integer nnode_sets_g 
      integer nside_sets_g 
      call negig (idne, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr)
      if (ierr.ne.0) return
      call nepnspg2 (idne, global_ids, global_n_cnts,
     $  global_df_cnts, ierr, nnode_sets_g)
C
      end
C
      subroutine nepnspg2 (idne, global_ids, global_n_cnts,
     $  global_df_cnts, ierr, nnode_sets_g)
      implicit none
      integer idne 
      integer nnode_sets_g
      integer global_ids(nnode_sets_g) 
      integer global_n_cnts(nnode_sets_g) 
      integer global_df_cnts(nnode_sets_g) 
      integer ierr
C
      integer*4 idne4 
      integer*4 global_ids4(nnode_sets_g) 
      integer*4 global_n_cnts4(nnode_sets_g) 
      integer*4 global_df_cnts4(nnode_sets_g) 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      DO J1X = 1, NNODE_SETS_G
        GLOBAL_IDS4(J1X) = GLOBAL_IDS(J1X)
        GLOBAL_N_CNTS4(J1X) = GLOBAL_N_CNTS(J1X)
        GLOBAL_DF_CNTS4(J1X) = GLOBAL_DF_CNTS(J1X)
      END DO
      call nepnspg4 (idne4, global_ids4, global_n_cnts4,
     $  global_df_cnts4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read global side set parameters
C
      subroutine negsspg (idne, ss_ids_glob, ss_n_cnt_glob,
     $  ss_df_cnt_glob, ierr)
      implicit none
      integer idne 
      integer ss_ids_glob(*) 
      integer ss_n_cnt_glob(*) 
      integer ss_df_cnt_glob(*) 
      integer ierr
C
      integer nnodes_g 
      integer nelems_g 
      integer nelem_blks_g 
      integer nnode_sets_g 
      integer nside_sets_g 
C
      call negig (idne, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr)
      if (ierr.ne.0) return
      call negsspg2 (idne, ss_ids_glob, ss_n_cnt_glob,
     $  ss_df_cnt_glob, ierr, nside_sets_g)
C
      end
C
      subroutine negsspg2 (idne, ss_ids_glob, ss_n_cnt_glob,
     $  ss_df_cnt_glob, ierr, nside_sets_g)
      implicit none
      integer idne 
      integer nside_sets_g
      integer ss_ids_glob(nside_sets_g) 
      integer ss_n_cnt_glob(nside_sets_g) 
      integer ss_df_cnt_glob(nside_sets_g) 
      integer ierr
C
      integer*4 idne4 
      integer*4 ss_ids_glob4(nside_sets_g) 
      integer*4 ss_n_cnt_glob4(nside_sets_g) 
      integer*4 ss_df_cnt_glob4(nside_sets_g) 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      call negsspg4 (idne4, ss_ids_glob4, ss_n_cnt_glob4,
     $  ss_df_cnt_glob4, ierr4)
      DO J1X = 1, NSIDE_SETS_G
        SS_IDS_GLOB(J1X) = SS_IDS_GLOB4(J1X)
        SS_N_CNT_GLOB(J1X) = SS_N_CNT_GLOB4(J1X)
        SS_DF_CNT_GLOB(J1X) = SS_DF_CNT_GLOB4(J1X)
      END DO
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write global side set parameters
C
      subroutine nepsspg (idne, global_ids, global_el_cnts,
     $  global_df_cnts, ierr)
      implicit none
      integer idne 
      integer global_ids(*) 
      integer global_el_cnts(*) 
      integer global_df_cnts(*) 
      integer ierr
C
      integer nnodes_g 
      integer nelems_g 
      integer nelem_blks_g 
      integer nnode_sets_g 
      integer nside_sets_g 
C
      call negig (idne, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr)
      if (ierr.ne.0) return
      call nepsspg2 (idne, global_ids, global_el_cnts,
     $  global_df_cnts, ierr, nside_sets_g)
C
      end
C
      subroutine nepsspg2 (idne, global_ids, global_el_cnts,
     $  global_df_cnts, ierr, nside_sets_g)
      implicit none
      integer idne 
      integer nside_sets_g
      integer global_ids(nside_sets_g) 
      integer global_el_cnts(nside_sets_g) 
      integer global_df_cnts(nside_sets_g) 
      integer ierr
C
      integer*4 idne4 
      integer*4 global_ids4(nside_sets_g) 
      integer*4 global_el_cnts4(nside_sets_g) 
      integer*4 global_df_cnts4(nside_sets_g) 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      DO J1X = 1, NSIDE_SETS_G
        GLOBAL_IDS4(J1X) = GLOBAL_IDS(J1X)
        GLOBAL_EL_CNTS4(J1X) = GLOBAL_EL_CNTS(J1X)
        GLOBAL_DF_CNTS4(J1X) = GLOBAL_DF_CNTS(J1X)
      END DO
      call nepsspg4 (idne4, global_ids4, global_el_cnts4,
     $  global_df_cnts4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read global element block information
C
      subroutine negebig (idne, el_blk_ids, el_blk_cnts, ierr)
      implicit none
      integer idne 
      integer el_blk_ids(*) 
      integer el_blk_cnts(*) 
      integer ierr
C
      integer nnodes_g 
      integer nelems_g 
      integer nelem_blks_g 
      integer nnode_sets_g 
      integer nside_sets_g 
C
      call negig (idne, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr)
      if (ierr.ne.0) return
      call negebig2 (idne, el_blk_ids, el_blk_cnts, ierr,
     $  nelem_blks_g)
C
      end
C
      subroutine negebig2 (idne, el_blk_ids, el_blk_cnts, ierr,
     $  nelem_blks_g)
      implicit none
      integer idne 
      integer nelem_blks_g
      integer el_blk_ids(nelem_blks_g) 
      integer el_blk_cnts(nelem_blks_g) 
      integer ierr
C
      integer*4 idne4 
      integer*4 el_blk_ids4(nelem_blks_g) 
      integer*4 el_blk_cnts4(nelem_blks_g) 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      call negebig4 (idne4, el_blk_ids4, el_blk_cnts4, ierr4)
      DO J1X = 1, NELEM_BLKS_G
        EL_BLK_IDS(J1X) = EL_BLK_IDS4(J1X)
        EL_BLK_CNTS(J1X) = EL_BLK_CNTS4(J1X)
      END DO
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write global element block information
C
      subroutine nepebig (idne, el_blk_ids, el_blk_cnts, ierr)
      implicit none
      integer idne 
      integer el_blk_ids(*) 
      integer el_blk_cnts(*) 
      integer ierr
C
      integer nnodes_g 
      integer nelems_g 
      integer nelem_blks_g 
      integer nnode_sets_g 
      integer nside_sets_g 
C
      call negig (idne, nnodes_g, nelems_g, nelem_blks_g,
     $  nnode_sets_g, nside_sets_g, ierr)
      if (ierr.ne.0) return
      call nepebig2 (idne, el_blk_ids, el_blk_cnts, ierr,
     $  nelem_blks_g)
C
      end
C
      subroutine nepebig2 (idne, el_blk_ids, el_blk_cnts, ierr,
     $  nelem_blks_g)
      implicit none
      integer idne 
      integer nelem_blks_g
      integer el_blk_ids(nelem_blks_g) 
      integer el_blk_cnts(nelem_blks_g) 
      integer ierr
C
      integer*4 idne4 
      integer*4 el_blk_ids4(nelem_blks_g) 
      integer*4 el_blk_cnts4(nelem_blks_g) 
      integer*4 ierr4
C
      idne4 = idne
      el_blk_ids4(nelem_blks_g) = el_blk_ids(nelem_blks_g)
      el_blk_cnts4(nelem_blks_g) = el_blk_cnts(nelem_blks_g)
      call nepebig4 (idne4, el_blk_ids4, el_blk_cnts4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read side set element list and side set side list
C
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
C
      integer*4 idne4 
      integer*4 ss_id4 
      integer*4 start_side_num4 
      integer*4 num_sides4 
      integer*4 ss_elem_list4(num_sides) 
      integer*4 ss_side_list4(num_sides) 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      ss_id4 = ss_id
      start_side_num4 = start_side_num
      num_sides4 = num_sides
      call negnss4 (idne4, ss_id4, start_side_num4, num_sides4,
     $  ss_elem_list4, ss_side_list4, ierr4)
      DO J1X = 1, NUM_SIDES
        SS_ELEM_LIST(J1X) = SS_ELEM_LIST4(J1X)
        SS_SIDE_LIST(J1X) = SS_SIDE_LIST4(J1X)
      END DO
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write side set element list and side set side list
C
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
C
      integer*4 idne4 
      integer*4 ss_id4 
      integer*4 start_side_num4 
      integer*4 num_sides4 
      integer*4 ss_elem_list4(num_sides) 
      integer*4 ss_side_list4(num_sides) 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      ss_id4 = ss_id
      start_side_num4 = start_side_num
      num_sides4 = num_sides
      DO J1X = 1, NUM_SIDES
        SS_ELEM_LIST4(J1X) = SS_ELEM_LIST(J1X)
        SS_SIDE_LIST4(J1X) = SS_SIDE_LIST(J1X)
      END DO
      call nepnss4 (idne4, ss_id4, start_side_num4, num_sides4,
     $  ss_elem_list4, ss_side_list4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read side set distribution factor
C
      subroutine negnssd (idne, ss_id, start_num, num_df_to_get,
     $  ss_df, ierr)
      implicit none
      integer idne 
      integer ss_id 
      integer start_num 
      integer num_df_to_get 
      real ss_df(num_df_to_get) 
      integer ierr
C
      integer*4 idne4 
      integer*4 ss_id4 
      integer*4 start_num4 
      integer*4 num_df_to_get4 
      integer*4 ierr4
C
      idne4 = idne
      ss_id4 = ss_id
      start_num4 = start_num
      num_df_to_get4 = num_df_to_get
      call negnssd4 (idne4, ss_id4, start_num4, num_df_to_get4,
     $  ss_df, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write side set distribution factor
C
      subroutine nepnssd (idne, ss_id, start_num, num_df_to_get,
     $  ss_df, ierr)
      implicit none
      integer idne 
      integer ss_id 
      integer start_num 
      integer num_df_to_get 
      real ss_df(num_df_to_get) 
      integer ierr
C
      integer*4 idne4 
      integer*4 ss_id4 
      integer*4 start_num4 
      integer*4 num_df_to_get4 
      integer*4 ierr4
C
      idne4 = idne
      ss_id4 = ss_id
      start_num4 = start_num
      num_df_to_get4 = num_df_to_get
      call nepnssd4 (idne4, ss_id4, start_num4, num_df_to_get4,
     $  ss_df, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read node set list for a single node set
C
      subroutine negnns (idne, ns_id, start_node_num, num_node,
     $  ns_node_list, ierr)
      implicit none
      integer idne 
      integer ns_id 
      integer start_node_num 
      integer num_node 
      integer ns_node_list(num_node) 
      integer ierr
C
      integer*4 idne4 
      integer*4 ns_id4 
      integer*4 start_node_num4 
      integer*4 num_node4 
      integer*4 ns_node_list4(num_node) 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      ns_id4 = ns_id
      start_node_num4 = start_node_num
      num_node4 = num_node
      call negnns4 (idne4, ns_id4, start_node_num4, num_node4,
     $  ns_node_list4, ierr4)
      DO J1X = 1, NUM_NODE
        NS_NODE_LIST(J1X) = NS_NODE_LIST4(J1X)
      END DO
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write node set list for a single node set
C
      subroutine nepnns (idne, ns_id, start_node_num, num_node,
     $  ns_node_list, ierr)
      implicit none
      integer idne 
      integer ns_id 
      integer start_node_num 
      integer num_node 
      integer ns_node_list(num_node) 
      integer ierr
C
      integer*4 idne4 
      integer*4 ns_id4 
      integer*4 start_node_num4 
      integer*4 num_node4 
      integer*4 ns_node_list4(num_node) 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      ns_id4 = ns_id
      start_node_num4 = start_node_num
      num_node4 = num_node
      DO J1X = 1, NUM_NODE
        NS_NODE_LIST4(J1X) = NS_NODE_LIST(J1X)
      END DO
      call nepnns4 (idne4, ns_id4, start_node_num4, num_node4,
     $  ns_node_list4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read node set distribution factor
C
      subroutine negnnsd (idne, ns_id, start_num, num_df_to_get,
     $  ns_df, ierr)
      implicit none
      integer idne 
      integer ns_id 
      integer start_num 
      integer num_df_to_get 
      real ns_df(num_df_to_get) 
      integer ierr
C
      integer*4 idne4 
      integer*4 ns_id4 
      integer*4 start_num4 
      integer*4 num_df_to_get4 
      integer*4 ierr4
C
      idne4 = idne
      ns_id4 = ns_id
      start_num4 = start_num
      num_df_to_get4 = num_df_to_get
      call negnnsd4 (idne4, ns_id4, start_num4, num_df_to_get4,
     $  ns_df, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write node set distribution factor
C
      subroutine nepnnsd (idne, ns_id, start_num, num_df_to_get,
     $  ns_df, ierr)
      implicit none
      integer idne 
      integer ns_id 
      integer start_num 
      integer num_df_to_get 
      real ns_df(num_df_to_get) 
      integer ierr
C
      integer*4 idne4 
      integer*4 ns_id4 
      integer*4 start_num4 
      integer*4 num_df_to_get4 
      integer*4 ierr4
C
      idne4 = idne
      ns_id4 = ns_id
      start_num4 = start_num
      num_df_to_get4 = num_df_to_get
      call nepnnsd4 (idne4, ns_id4, start_num4, num_df_to_get4,
     $  ns_df, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read coordinates of the nodes
C
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
C
      integer*4 idne4 
      integer*4 start_node_num4 
      integer*4 num_nodes4 
      integer*4 ierr4
C
      idne4 = idne
      start_node_num4 = start_node_num
      num_nodes4 = num_nodes
      call negcor4 (idne4, start_node_num4, num_nodes4, x_coor,
     $  y_coor, z_coor, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write coordinates of the nodes
C
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
C
      integer*4 idne4 
      integer*4 start_node_num4 
      integer*4 num_nodes4 
      integer*4 ierr4
C
      idne4 = idne
      start_node_num4 = start_node_num
      num_nodes4 = num_nodes
      call nepcor4 (idne4, start_node_num4, num_nodes4, x_coor,
     $  y_coor, z_coor, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read an element blocks connectivity list
C
      subroutine negnec (idne, elem_blk_id, start_elem_num,
     $  num_elems, connect, ierr)
      implicit none
      integer idne 
      integer elem_blk_id 
      integer start_elem_num 
      integer num_elems 
      integer connect(num_elems) 
      integer ierr
C
      integer*4 idne4 
      integer*4 elem_blk_id4 
      integer*4 start_elem_num4 
      integer*4 num_elems4 
      integer*4 connect4(num_elems) 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      elem_blk_id4 = elem_blk_id
      start_elem_num4 = start_elem_num
      num_elems4 = num_elems
      call negnec4 (idne4, elem_blk_id4, start_elem_num4,
     $  num_elems4, connect4, ierr4)
      DO J1X = 1, NUM_ELEMS
        CONNECT(J1X) = CONNECT4(J1X)
      END DO
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write an element blocks connectivity list
C
      subroutine nepnec (idne, elem_blk_id, start_elem_num,
     $  num_elems, connect, ierr)
      implicit none
      integer idne 
      integer elem_blk_id 
      integer start_elem_num 
      integer num_elems 
      integer connect(num_elems) 
      integer ierr
C
      integer*4 idne4 
      integer*4 elem_blk_id4 
      integer*4 start_elem_num4 
      integer*4 num_elems4 
      integer*4 connect4(num_elems) 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      elem_blk_id4 = elem_blk_id
      start_elem_num4 = start_elem_num
      num_elems4 = num_elems
      DO J1X = 1, NUM_ELEMS
        CONNECT4(J1X) = CONNECT(J1X)
      END DO
      call nepnec4 (idne4, elem_blk_id4, start_elem_num4,
     $  num_elems4, connect4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read an element blocks attributes
C
      subroutine negneat (idne, elem_blk_id, start_elem_num,
     $  num_elems, attrib, ierr)
      implicit none
      integer idne 
      integer elem_blk_id 
      integer start_elem_num 
      integer num_elems 
      real attrib(num_elems) 
      integer ierr
C
      integer*4 idne4 
      integer*4 elem_blk_id4 
      integer*4 start_elem_num4 
      integer*4 num_elems4 
      integer*4 ierr4
C
      idne4 = idne
      elem_blk_id4 = elem_blk_id
      start_elem_num4 = start_elem_num
      num_elems4 = num_elems
      call negneat4 (idne4, elem_blk_id4, start_elem_num4,
     $  num_elems4, attrib, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write an element blocks attributes
C
      subroutine nepneat (idne, elem_blk_id, start_elem_num,
     $  num_elems, attrib, ierr)
      implicit none
      integer idne 
      integer elem_blk_id 
      integer start_elem_num 
      integer num_elems 
      real attrib(num_elems) 
      integer ierr
C
      integer*4 idne4 
      integer*4 elem_blk_id4 
      integer*4 start_elem_num4 
      integer*4 num_elems4 
      integer*4 ierr4
C
      idne4 = idne
      elem_blk_id4 = elem_blk_id
      start_elem_num4 = start_elem_num
      num_elems4 = num_elems
      call nepneat4 (idne4, elem_blk_id4, start_elem_num4,
     $  num_elems4, attrib, ierr4)
      ierr = ierr4
C
      end
C
C-----------------------------------------------------------------------
C  Read the element type for a specific element block
C
      subroutine negelt (idne, elem_blk_id, elem_type, ierr)
      implicit none
      integer idne 
      integer elem_blk_id 
      character*(*) elem_type 
      integer ierr
C
      integer*4 idne4 
      integer*4 elem_blk_id4 
      integer*4 ierr4
C
      idne4 = idne
      elem_blk_id4 = elem_blk_id
      call negelt4 (idne4, elem_blk_id4, elem_type, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read a variable for an element block
C
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
C
      integer*4 idne4 
      integer*4 time_step4 
      integer*4 elem_var_index4 
      integer*4 elem_blk_id4 
      integer*4 num_elem_this_blk4 
      integer*4 start_elem_num4 
      integer*4 num_elem4 
      integer*4 ierr4
C
      idne4 = idne
      time_step4 = time_step
      elem_var_index4 = elem_var_index
      elem_blk_id4 = elem_blk_id
      num_elem_this_blk4 = num_elem_this_blk
      start_elem_num4 = start_elem_num
      num_elem4 = num_elem
      call negnev4 (idne4, time_step4, elem_var_index4,
     $  elem_blk_id4, num_elem_this_blk4, start_elem_num4,
     $  num_elem4, elem_var_vals, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write a variable slab for an element block
C
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
C
      integer*4 idne4 
      integer*4 time_step4 
      integer*4 elem_var_index4 
      integer*4 elem_blk_id4 
      integer*4 start_pos4 
      integer*4 num_vals4 
      integer*4 ierr4
C
      idne4 = idne
      time_step4 = time_step
      elem_var_index4 = elem_var_index
      elem_blk_id4 = elem_blk_id
      start_pos4 = start_pos
      num_vals4 = num_vals
      call nepevs4 (idne4, time_step4, elem_var_index4,
     $  elem_blk_id4, start_pos4, num_vals4, elem_var_vals,
     $  ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read the values of a single nodal variable for a single time step
C
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
C
      integer*4 idne4 
      integer*4 time_step4 
      integer*4 nodal_var_index4 
      integer*4 start_node_num4 
      integer*4 num_nodes4 
      integer*4 ierr4
C
      idne4 = idne
      time_step4 = time_step
      nodal_var_index4 = nodal_var_index
      start_node_num4 = start_node_num
      num_nodes4 = num_nodes
      call negnnv4 (idne4, time_step4, nodal_var_index4,
     $  start_node_num4, num_nodes4, nodal_vars, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write nodal variable slab
C
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
      integer i
C
      integer*4 idne4 
      integer*4 time_step4 
      integer*4 nodal_var_index4 
      integer*4 start_pos4 
      integer*4 num_vals4 
      integer*4 ierr4
C
      idne4 = idne
      time_step4 = time_step
      nodal_var_index4 = nodal_var_index
      start_pos4 = start_pos
      num_vals4 = num_vals
      call nepnvs4 (idne4, time_step4, nodal_var_index4,
     $  start_pos4, num_vals4, nodal_var_vals, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read the element numbering map
C
      subroutine negnenm (idne, starte, num_ent, elem_map, ierr)
      implicit none
      integer idne 
      integer starte 
      integer num_ent 
      integer elem_map(num_ent) 
      integer ierr
      integer i
C
      integer*4 idne4 
      integer*4 starte4 
      integer*4 num_ent4 
      integer*4 elem_map4(num_ent) 
      integer*4 ierr4
C
      idne4 = idne
      starte4 = starte
      num_ent4 = num_ent
      call negnenm4 (idne4, starte4, num_ent4, elem_map4, ierr4)
      do i = 1, num_ent
        elem_map(i) = elem_map4(i)
      enddo
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write the element numbering map
C
      subroutine nepnenm (idne, starte, num_ent, elem_map, ierr)
      implicit none
      integer idne 
      integer starte 
      integer num_ent 
      integer elem_map(num_ent) 
      integer ierr
      integer i
C
      integer*4 idne4 
      integer*4 starte4 
      integer*4 num_ent4 
      integer*4 elem_map4(num_ent) 
      integer*4 ierr4
C
      idne4 = idne
      starte4 = starte
      num_ent4 = num_ent
      do i = 1, num_ent
        elem_map4(i) = elem_map(i)
      enddo
      call nepnenm4 (idne4, starte4, num_ent4, elem_map4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read the node numbering map
C
      subroutine negnnnm (idne, startn, num_ent, node_map, ierr)
      implicit none
      integer idne 
      integer startn 
      integer num_ent 
      integer node_map(num_ent) 
      integer ierr
      integer i
C
      integer*4 idne4 
      integer*4 startn4 
      integer*4 num_ent4 
      integer*4 node_map4(num_ent) 
      integer*4 ierr4
C
      idne4 = idne
      startn4 = startn
      num_ent4 = num_ent
      call negnnnm4 (idne4, startn4, num_ent4, node_map4, ierr4)
      do i = 1, num_ent
        node_map(i) = node_map4(i)
      enddo
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write the node numbering map
C
      subroutine nepnnnm (idne, startn, num_ent, node_map, ierr)
      implicit none
      integer idne 
      integer startn 
      integer num_ent 
      integer node_map(num_ent) 
      integer ierr
      integer i
C
      integer*4 idne4 
      integer*4 startn4 
      integer*4 num_ent4 
      integer*4 node_map4(num_ent) 
      integer*4 ierr4
C
      idne4 = idne
      startn4 = startn
      num_ent4 = num_ent
      do i = 1, num_ent
        node_map4(i) = node_map(i)
      enddo
      call nepnnnm4 (idne4, startn4, num_ent4, node_map4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read the node map for a processor
C
      subroutine negnm (idne, node_mapi, node_mapb, node_mape,
     $  processor, ierr)
      implicit none
      integer idne 
      integer node_mapi(*) 
      integer node_mapb(*) 
      integer node_mape(*) 
      integer processor 
      integer ierr
C
      integer nint_nodes 
      integer nbor_nodes 
      integer next_nodes 
      integer nint_elems 
      integer nbor_elems 
      integer nnode_cmaps 
      integer nelem_cmaps 
      call neglbp (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
      if (ierr.ne.0) return
      call negnm2 (idne, node_mapi, node_mapb, node_mape,
     $  processor, ierr,
     $  nint_nodes, nbor_nodes, next_nodes)
C
      end
C
      subroutine negnm2 (idne, node_mapi, node_mapb, node_mape,
     $  processor, ierr,
     $  nint_nodes, nbor_nodes, next_nodes)
      implicit none
      integer idne 
      integer nint_nodes 
      integer nbor_nodes 
      integer next_nodes 
      integer node_mapi(nint_nodes) 
      integer node_mapb(nbor_nodes) 
      integer node_mape(next_nodes) 
      integer processor 
      integer ierr
      integer i
C
      integer*4 idne4 
      integer*4 node_mapi4(nint_nodes) 
      integer*4 node_mapb4(nbor_nodes) 
      integer*4 node_mape4(next_nodes) 
      integer*4 processor4 
      integer*4 ierr4
C
      idne4 = idne
      processor4 = processor
      call negnm4 (idne4, node_mapi4, node_mapb4, node_mape4,
     $  processor4, ierr4)
      do i = 1, nint_nodes
        node_mapi(i) = node_mapi4(i)
      enddo
      do i = 1, nbor_nodes
        node_mapb(i) = node_mapb4(i)
      enddo
      do i = 1, next_nodes
        node_mape(i) = node_mape4(i)
      enddo
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write a node map for a processor
C
      subroutine nepnm (idne, node_mapi, node_mapb, node_mape,
     $  processor, ierr)
      implicit none
      integer idne 
      integer node_mapi(*) 
      integer node_mapb(*) 
      integer node_mape(*) 
      integer processor 
      integer ierr
C
      integer nint_nodes 
      integer nbor_nodes 
      integer next_nodes 
      integer nint_elems 
      integer nbor_elems 
      integer nnode_cmaps 
      integer nelem_cmaps 
      call neglbp (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
      if (ierr.ne.0) return
      call nepnm2 (idne, node_mapi, node_mapb, node_mape,
     $  processor, ierr,
     $  nint_nodes, nbor_nodes, next_nodes)
C
      end
C
      subroutine nepnm2 (idne, node_mapi, node_mapb, node_mape,
     $  processor, ierr,
     $  nint_nodes, nbor_nodes, next_nodes)
      implicit none
      integer idne 
      integer nint_nodes 
      integer nbor_nodes 
      integer next_nodes 
      integer node_mapi(nint_nodes) 
      integer node_mapb(nbor_nodes) 
      integer node_mape(next_nodes) 
      integer processor 
      integer ierr
C
      integer*4 idne4 
      integer*4 node_mapi4(nint_nodes) 
      integer*4 node_mapb4(nbor_nodes) 
      integer*4 node_mape4(next_nodes) 
      integer*4 processor4 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      DO J1X = 1, NINT_NODES
        NODE_MAPI4(J1X) = NODE_MAPI(J1X)
      END DO
      DO J1X = 1, NBOR_NODES
        NODE_MAPB4(J1X) = NODE_MAPB(J1X)
      END DO
      DO J1X = 1, NEXT_NODES
        NODE_MAPE4(J1X) = NODE_MAPE(J1X)
      END DO
      processor4 = processor
      call nepnm4 (idne4, node_mapi4, node_mapb4, node_mape4,
     $  processor4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read the element map for a processor
C
      subroutine negem (idne, elem_mapi, elem_mapb, processor, ierr)
      implicit none
      integer idne 
      integer elem_mapi(*) 
      integer elem_mapb(*) 
      integer processor 
      integer ierr
C
      integer nint_nodes 
      integer nbor_nodes 
      integer next_nodes 
      integer nint_elems 
      integer nbor_elems 
      integer nnode_cmaps 
      integer nelem_cmaps 
      call neglbp (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
      if (ierr.ne.0) return
      call negem2 (idne, elem_mapi, elem_mapb, processor, ierr,
     $  nint_elems, nbor_elems)
C
      end
C
      subroutine negem2 (idne, elem_mapi, elem_mapb, processor, ierr,
     $  nint_elems, nbor_elems)
      implicit none
      integer idne 
      integer nint_elems 
      integer nbor_elems 
      integer elem_mapi(nint_elems) 
      integer elem_mapb(nbor_elems) 
      integer processor 
      integer ierr
C
      integer*4 idne4 
      integer*4 elem_mapi4(nint_elems) 
      integer*4 elem_mapb4(nbor_elems) 
      integer*4 processor4 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      processor4 = processor
      call negem4 (idne4, elem_mapi4, elem_mapb4, processor4, ierr4)
      DO J1X = 1, NINT_ELEMS
        ELEM_MAPI(J1X) = ELEM_MAPI4(J1X)
      END DO
      DO J1X = 1, NBOR_ELEMS
        ELEM_MAPB(J1X) = ELEM_MAPB4(J1X)
      END DO
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write the element map for a processor
C
      subroutine nepem (idne, elem_mapi, elem_mapb, processor, ierr)
      implicit none
      integer idne 
      integer elem_mapi(*) 
      integer elem_mapb(*) 
      integer processor 
      integer ierr
C
      integer nint_nodes 
      integer nbor_nodes 
      integer next_nodes 
      integer nint_elems 
      integer nbor_elems 
      integer nnode_cmaps 
      integer nelem_cmaps 
      call neglbp (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
      if (ierr.ne.0) return
      call nepem2 (idne, elem_mapi, elem_mapb, processor, ierr,
     $  nint_elems, nbor_elems)
C
      end
C
      subroutine nepem2 (idne, elem_mapi, elem_mapb, processor, ierr,
     $  nint_elems, nbor_elems)
      implicit none
      integer idne 
      integer nint_elems 
      integer nbor_elems 
      integer elem_mapi(nint_elems) 
      integer elem_mapb(nbor_elems) 
      integer processor 
      integer ierr
C
      integer*4 idne4 
      integer*4 elem_mapi4(nint_elems) 
      integer*4 elem_mapb4(nbor_elems) 
      integer*4 processor4 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      DO J1X = 1, NINT_ELEMS
        ELEM_MAPI4(J1X) = ELEM_MAPI(J1X)
      END DO
      DO J1X = 1, NBOR_ELEMS
        ELEM_MAPB4(J1X) = ELEM_MAPB(J1X)
      END DO
      processor4 = processor
      call nepem4 (idne4, elem_mapi4, elem_mapb4, processor4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Read the communications map parameters for a single processor
C
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
C
      integer nint_nodes 
      integer nbor_nodes 
      integer next_nodes 
      integer nint_elems 
      integer nbor_elems 
      integer nnode_cmaps 
      integer nelem_cmaps 
      call neglbp (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
      if (ierr.ne.0) return
      call negcmp2 (idne, ncmap_ids, ncmap_node_cnts,
     $  ecmap_ids, ecmap_elem_cnts, processor, ierr,
     $  nnode_cmaps, nelem_cmaps)
C
      end
C
      subroutine negcmp2 (idne, ncmap_ids, ncmap_node_cnts,
     $  ecmap_ids, ecmap_elem_cnts, processor, ierr,
     $  nnode_cmaps, nelem_cmaps)
      implicit none
      integer idne 
      integer nnode_cmaps 
      integer nelem_cmaps 
      integer ncmap_ids(nnode_cmaps) 
      integer ncmap_node_cnts(nnode_cmaps) 
      integer ecmap_ids(nelem_cmaps) 
      integer ecmap_elem_cnts(nelem_cmaps) 
      integer processor 
      integer ierr
C
      integer*4 idne4 
      integer*4 ncmap_ids4(nnode_cmaps) 
      integer*4 ncmap_node_cnts4(nnode_cmaps) 
      integer*4 ecmap_ids4(nelem_cmaps) 
      integer*4 ecmap_elem_cnts4(nelem_cmaps) 
      integer*4 processor4 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      processor4 = processor
      call negcmp4 (idne4, ncmap_ids4, ncmap_node_cnts4,
     $  ecmap_ids4, ecmap_elem_cnts4, processor4, ierr4)
      DO J1X = 1, NNODE_CMAPS
        NCMAP_IDS(J1X) = NCMAP_IDS4(J1X)
        NCMAP_NODE_CNTS(J1X) = NCMAP_NODE_CNTS4(J1X)
      END DO
      DO J1X = 1, NELEM_CMAPS
        ECMAP_IDS(J1X) = ECMAP_IDS4(J1X)
        ECMAP_ELEM_CNTS(J1X) = ECMAP_ELEM_CNTS4(J1X)
      END DO
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write the communications map parameters for a single processor
C
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
C
      integer nint_nodes 
      integer nbor_nodes 
      integer next_nodes 
      integer nint_elems 
      integer nbor_elems 
      integer nnode_cmaps 
      integer nelem_cmaps 
      call neglbp (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
      if (ierr.ne.0) return
      call nepcmp2 (idne, nmap_ids, nmap_node_cnts, emap_ids,
     $  emap_elem_cnts, processor, ierr,
     $  nnode_cmaps, nelem_cmaps)
C
      end
C
      subroutine nepcmp2 (idne, nmap_ids, nmap_node_cnts, emap_ids,
     $  emap_elem_cnts, processor, ierr,
     $  nnode_cmaps, nelem_cmaps)
      implicit none
      integer idne 
      integer nnode_cmaps 
      integer nelem_cmaps 
      integer nmap_ids(nnode_cmaps) 
      integer nmap_node_cnts(nnode_cmaps) 
      integer emap_ids(nelem_cmaps) 
      integer emap_elem_cnts(nelem_cmaps) 
      integer processor 
      integer ierr
C
      integer*4 idne4 
      integer*4 nmap_ids4(nnode_cmaps) 
      integer*4 nmap_node_cnts4(nnode_cmaps) 
      integer*4 emap_ids4(nelem_cmaps) 
      integer*4 emap_elem_cnts4(nelem_cmaps) 
      integer*4 processor4 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      DO J1X = 1, NNODE_CMAPS
        NMAP_IDS4(J1X) = NMAP_IDS(J1X)
        NMAP_NODE_CNTS4(J1X) = NMAP_NODE_CNTS(J1X)
      END DO
      DO J1X = 1, NELEM_CMAPS
        EMAP_IDS4(J1X) = EMAP_IDS(J1X)
        EMAP_ELEM_CNTS4(J1X) = EMAP_ELEM_CNTS(J1X)
      END DO
      processor4 = processor
      call nepcmp4 (idne4, nmap_ids4, nmap_node_cnts4, emap_ids4,
     $  emap_elem_cnts4, processor4, ierr4)
      ierr = ierr4
C
      end
      
C-----------------------------------------------------------------------
C
C  Write the communications map parameters for all processors
C
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
      
      write (0,*) 'nepcmpc: addrwrap not implemented'
      stop
      end
      
C-----------------------------------------------------------------------
C
C  Read the nodal communications map for a single processor
C
      subroutine negncm (idne, map_id, node_ids, proc_ids,
     $  processor, ierr)
      implicit none
      integer idne 
      integer map_id 
      integer node_ids(*) 
      integer proc_ids(*) 
      integer processor 
      integer ierr
C
      integer nint_nodes 
      integer nbor_nodes 
      integer next_nodes 
      integer nint_elems 
      integer nbor_elems 
      integer nnode_cmaps 
      integer nelem_cmaps 
      call neglbp (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
      if (ierr.ne.0) return
      call negncm2 (idne, map_id, node_ids, proc_ids,
     $  processor, ierr, nnode_cmaps, nelem_cmaps)
C
      end
C
      subroutine negncm2 (idne, map_id, node_ids, proc_ids,
     $  processor, ierr, nnode_cmaps, nelem_cmaps)
      implicit none
      integer idne 
      integer map_id 
      integer node_ids(*) 
      integer proc_ids(*) 
      integer processor 
      integer nnode_cmaps 
      integer nelem_cmaps 
      integer ierr
C
      integer ncmap_ids(nnode_cmaps) 
      integer ncmap_node_cnts(nnode_cmaps) 
      integer ecmap_ids(nelem_cmaps) 
      integer ecmap_elem_cnts(nelem_cmaps) 
      integer node_cnt, i 
C
      call negcmp2 (idne, ncmap_ids, ncmap_node_cnts,
     $  ecmap_ids, ecmap_elem_cnts, processor, ierr,
     $  nnode_cmaps, nelem_cmaps)
      if (ierr.ne.0) return
      node_cnt = -1
      do i = 1, nnode_cmaps
        if (ncmap_ids(i) .eq. map_id) then
          node_cnt = ncmap_node_cnts(i)
        end if
      end do
      if (node_cnt .eq. -1) then
        ierr = -1
        return
      end if
      call negncm3 (idne, map_id, node_ids, proc_ids,
     $  processor, ierr, node_cnt)
C
      end
C
      subroutine negncm3 (idne, map_id, node_ids, proc_ids,
     $  processor, ierr, node_cnt)
      implicit none
      integer idne 
      integer map_id 
      integer node_cnt 
      integer node_ids(node_cnt) 
      integer proc_ids(node_cnt) 
      integer processor 
      integer ierr
C
      integer*4 idne4 
      integer*4 map_id4 
      integer*4 node_ids4(node_cnt) 
      integer*4 proc_ids4(node_cnt) 
      integer*4 processor4 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      map_id4 = map_id
         processor4 = processor
         call negncm4 (idne4, map_id4, node_ids4, proc_ids4,
     $     processor4, ierr4)
         DO J1X = 1, NODE_CNT
           NODE_IDS(J1X) = NODE_IDS4(J1X)
           PROC_IDS(J1X) = PROC_IDS4(J1X)
         END DO
         ierr = ierr4
C
         end
      
C-----------------------------------------------------------------------
C
C  Write the nodal communications map for a single processor
C
      subroutine nepncm (idne, map_id, node_ids, proc_ids,
     $  processor, ierr)
      implicit none
      integer idne 
      integer map_id 
      integer node_ids(*) 
      integer proc_ids(*) 
      integer processor 
      integer ierr
C
      integer nint_nodes 
      integer nbor_nodes 
      integer next_nodes 
      integer nint_elems 
      integer nbor_elems 
      integer nnode_cmaps 
      integer nelem_cmaps 
      call neglbp (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
      if (ierr.ne.0) return
      call nepncm2 (idne, map_id, node_ids, proc_ids,
     $  processor, ierr, nnode_cmaps, nelem_cmaps)
C
      end
C
      subroutine nepncm2 (idne, map_id, node_ids, proc_ids,
     $  processor, ierr, nnode_cmaps, nelem_cmaps)
      implicit none













      integer idne 
      integer map_id 
      integer node_ids(*) 
      integer proc_ids(*) 
      integer processor 
      integer nnode_cmaps 
      integer nelem_cmaps 
      integer ierr
C
      integer ncmap_ids(nnode_cmaps) 
      integer ncmap_node_cnts(nnode_cmaps) 
      integer ecmap_ids(nelem_cmaps) 
      integer ecmap_elem_cnts(nelem_cmaps) 
      integer node_cnt, i 
C
      call negcmp2 (idne, ncmap_ids, ncmap_node_cnts,
     $  ecmap_ids, ecmap_elem_cnts, processor, ierr,
     $  nnode_cmaps, nelem_cmaps)
      if (ierr.ne.0) return
      node_cnt = -1
      do i = 1, nnode_cmaps
        if (ncmap_ids(i) .eq. map_id) then
          node_cnt = ncmap_node_cnts(i)
        end if
      end do
      if (node_cnt .eq. -1) then
        ierr = -1
        return
      end if
      call nepncm3 (idne, map_id, node_ids, proc_ids,
     $  processor, ierr, node_cnt)
C
      end
C
      subroutine nepncm3 (idne, map_id, node_ids, proc_ids,
     $  processor, ierr, node_cnt)
      implicit none
      integer idne 
      integer map_id 
      integer node_cnt
      integer node_ids(node_cnt) 
      integer proc_ids(node_cnt) 
      integer processor 
      integer ierr
C
      integer*4 idne4 
      integer*4 map_id4 
      integer*4 node_ids4(node_cnt) 
      integer*4 proc_ids4(node_cnt) 
      integer*4 processor4 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      map_id4 = map_id
         DO J1X = 1, NODE_CNT
           NODE_IDS4(J1X) = NODE_IDS(J1X)
           PROC_IDS4(J1X) = PROC_IDS(J1X)
         END DO
         processor4 = processor
         call nepncm4 (idne4, map_id4, node_ids4, proc_ids4,
     $     processor4, ierr4)
         ierr = ierr4
C
         end
      
C-----------------------------------------------------------------------
C
C  Read the elemental communications map for a single processor
C
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
C
      integer nint_nodes 
      integer nbor_nodes 
      integer next_nodes 
      integer nint_elems 
      integer nbor_elems 
      integer nnode_cmaps 
      integer nelem_cmaps 
      call neglbp (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
      if (ierr.ne.0) return
      call negecm2 (idne, map_id, elem_ids, side_ids,
     $  proc_ids, processor, ierr, nnode_cmaps, nelem_cmaps)
C
      end
C
      subroutine negecm2 (idne, map_id, elem_ids, side_ids,
     $  proc_ids, processor, ierr, nnode_cmaps, nelem_cmaps)
      implicit none
      integer idne 
      integer map_id 
      integer elem_ids(*) 
      integer side_ids(*) 
      integer proc_ids(*) 
      integer processor 
      integer nnode_cmaps 
      integer nelem_cmaps 
      integer ierr
C
      integer ncmap_ids(nnode_cmaps) 
      integer ncmap_node_cnts(nnode_cmaps) 
      integer ecmap_ids(nelem_cmaps) 
      integer ecmap_elem_cnts(nelem_cmaps) 
      integer elem_cnt,i 
C
      call negcmp2 (idne, ncmap_ids, ncmap_node_cnts,
     $  ecmap_ids, ecmap_elem_cnts, processor, ierr,
     $  nnode_cmaps, nelem_cmaps)
      if (ierr.ne.0) return
      elem_cnt = -1
      do i = 1, nelem_cmaps
        if (ecmap_ids(i) .eq. map_id) then
          elem_cnt = ecmap_elem_cnts(i)
        end if
      end do
      if (elem_cnt .eq. -1) then
        ierr = -1
        return
      end if
      call negecm3 (idne, map_id, elem_ids, side_ids,
     $  proc_ids, processor, ierr, elem_cnt)
C
      end
C
      subroutine negecm3 (idne, map_id, elem_ids, side_ids,
     $  proc_ids, processor, ierr, elem_cnt)
      implicit none
      integer idne 
      integer map_id 
      integer elem_cnt 
      integer elem_ids(elem_cnt) 
      integer side_ids(elem_cnt) 
      integer proc_ids(elem_cnt) 
      integer processor 
      integer ierr
C
      integer*4 idne4 
      integer*4 map_id4 
      integer*4 elem_ids4(elem_cnt) 
      integer*4 side_ids4(elem_cnt) 
      integer*4 proc_ids4(elem_cnt) 
      integer*4 processor4 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      map_id4 = map_id
         processor4 = processor
         call negecm4 (idne4, map_id4, elem_ids4, side_ids4,
     $     proc_ids4, processor4, ierr4)
         DO J1X = 1, ELEM_CNT
           ELEM_IDS(J1X) = ELEM_IDS4(J1X)
           SIDE_IDS(J1X) = SIDE_IDS4(J1X)
           PROC_IDS(J1X) = PROC_IDS4(J1X)
         END DO
         ierr = ierr4
C
         end
      
C-----------------------------------------------------------------------
C
C  Write the elemental communications map for a single processor
C
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
C
      integer nint_nodes 
      integer nbor_nodes 
      integer next_nodes 
      integer nint_elems 
      integer nbor_elems 
      integer nnode_cmaps 
      integer nelem_cmaps 
      call neglbp (idne, nint_nodes, nbor_nodes, next_nodes,
     $  nint_elems, nbor_elems, nnode_cmaps, nelem_cmaps,
     $  processor, ierr)
      if (ierr.ne.0) return
      call nepecm2 (idne, map_id, elem_ids, side_ids,
     $  proc_ids, processor, ierr, nnode_cmaps, nelem_cmaps)
C
      end
C
      subroutine nepecm2 (idne, map_id, elem_ids, side_ids,
     $  proc_ids, processor, ierr, nnode_cmaps, nelem_cmaps)
      implicit none
      integer idne 
      integer map_id 
      integer elem_ids(*) 
      integer side_ids(*) 
      integer proc_ids(*) 
      integer processor 
      integer nnode_cmaps 
      integer nelem_cmaps 
      integer ierr
C
      integer ncmap_ids(nnode_cmaps) 
      integer ncmap_node_cnts(nnode_cmaps) 
      integer ecmap_ids(nelem_cmaps) 
      integer ecmap_elem_cnts(nelem_cmaps) 
      integer elem_cnt,i 
C
      call negcmp2 (idne, ncmap_ids, ncmap_node_cnts,
     $  ecmap_ids, ecmap_elem_cnts, processor, ierr,
     $  nnode_cmaps, nelem_cmaps)
      if (ierr.ne.0) return
      elem_cnt = -1
      do i = 1, nelem_cmaps
        if (ecmap_ids(i) .eq. map_id) then
          elem_cnt = ecmap_elem_cnts(i)
        end if
      end do
      if (elem_cnt .eq. -1) then
        ierr = -1
        return
      end if
      call nepecm3 (idne, map_id, elem_ids, side_ids,
     $  proc_ids, processor, ierr, elem_cnt)
C
      end
C
      subroutine nepecm3 (idne, map_id, elem_ids, side_ids,
     $  proc_ids, processor, ierr, elem_cnt)
      implicit none
      integer idne 
      integer map_id 
      integer elem_cnt 
      integer elem_ids(elem_cnt) 
      integer side_ids(elem_cnt) 
      integer proc_ids(elem_cnt) 
      integer processor 
      integer ierr
C
      integer*4 idne4 
      integer*4 map_id4 
      integer*4 elem_ids4(elem_cnt) 
      integer*4 side_ids4(elem_cnt) 
      integer*4 proc_ids4(elem_cnt) 
      integer*4 processor4 
      integer*4 ierr4
C
      INTEGER J1X
      idne4 = idne
      map_id4 = map_id
         DO J1X = 1, ELEM_CNT
           ELEM_IDS4(J1X) = ELEM_IDS(J1X)
           SIDE_IDS4(J1X) = SIDE_IDS(J1X)
           PROC_IDS4(J1X) = PROC_IDS(J1X)
         END DO
         processor4 = processor
         call nepecm4 (idne4, map_id4, elem_ids4, side_ids4,
     $     proc_ids4, processor4, ierr4)
         ierr = ierr4
C
         end
C-----------------------------------------------------------------------
#endif 
