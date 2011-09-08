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
C*****************************************************************************
C
C NEMESISI_INT - nemesis I internals include file
C
C author - Matthew M. St. John, Sandia National Laboratories
C
C environment - UNIX F77
C
C revision history -
C
C  Based on exodusII_int.inc distributed with ACCESS software
C
C*****************************************************************************
C  --
C  --  This file contains define constants that are used to refer to
C  --  netCDF variables, attributes, or dimensions in which the NEMESIS
C  --  I data are stored.  Using the defined constants will allow the
C  --  names of the netCDF entities to be changed easily in the future
C  --  if needed.
C  --
C  --Constants:
C  --
C  --DIMENSIONS
C  --
C  --   D_NUMINOD  The number of internal nodes
C  --   D_NUMBNOD  The number of border nodes
C  --   D_NUMENOD  The number of external nodes
C  --   D_NUMIELM  The number of internal elements
C  --   D_NUMBELM  The number of border elements
C  --   D_NUMPROC  The total number of processors
C  --   D_NUMPR_F  The number of processors handled in this file
C  --   D_NUMNODG  The number of global nodes
C  --   D_NUMELMG  The number of global elements
C  --   D_NUMNSG   The number of global node sets
C  --   D_NUMSSG   The number of global side sets
C  --   D_NUMEBG   The number of global element blocks
C  --   D_NUMNCM   The number of nodal communications maps per processor
C  --   D_NUMECM   The number of elemental communications maps per processor
C  --   D_NCNTCM   The total number of nodal communications maps
C  --   D_ECNTCM   The total number of elemental communications maps
C  --
C  --VARIABLES
C  --
C  --   V_EBIDSG   The global block ids
C  --   V_NSIDSG   The global node set ids
C  --   V_NSNDCG   The global node set node counts
C  --   V_NSDFCG   The global dist factor for each global node set
C  --   V_SSIDSG   The global side set ids
C  --   V_SSSDCG   The global side counts
C  --   V_SSDFCG   The global dist factor for each global side set
C  --   V_FTYPE    Nemesis file type "p" parallel, "s" scaler
C  --   V_CMAP     Communications map (not currently used
C  --   V_NMAPI    The internal node map
C  --   V_NMAPII   The internal node map index
C  --   V_NMAPB    The border node map
C  --   V_NMAPBI   The border node map index
C  --   V_NMAPE    The external node map
C  --   V_NMAPEI   The external node map index
C  --   V_EMAPI    The internal element map
C  --   V_EMAPII   The internal element map index
C  --   V_EMAPB    The border element map
C  --   V_EMAPBI   The border element map index
C  --   V_INTNST   Status of internal node map
C  --   V_BORNST   Status of border node map
C  --   V_EXTNST   Status of external node map
C  --   V_INTEST   Status of internal element map
C  --   V_BOREST   Status of border element map
C  --   V_NCOMID   Nodal communications map ids
C  --   V_NCOMST   Status of nodal communications map
C  --   V_NCOMII   Nodal communications map information index
C  --   V_ECOMID   Elemental communications map ids
C  --   V_ECOMST   Status of elemental communications map
C  --   V_ECOMII   Elemental communications map information index
C  --   V_NCOMNID  Nodal communications map node ids
C  --   V_NCOMPR   Nodal communications map processors
C  --   V_NCOMDI   Nodal communications map data index
C  --   V_ECOMEID  Elemental communications map element ids
C  --   V_ECOMSID  Elemental communications map side ids
C  --   V_ECOMPR   Elemental communications map processors
C  --   V_ECOMDI   Elemental communications map data index
C  --
C  --MISC
C  --
C  --   END_CMIDX  String found at the end of the index variables
C  --   END_CMIDS  String found at the end of the ids variables
C  --

      character D_NUMINOD*12, D_NUMBNOD*12, D_NUMENOD*12, D_NUMIELM*12,
     1 D_NUMBELM*12, D_NUMPROC*14, D_NUMPR_F*14, D_NUMNODG*16, D_NUMELMG*16,
     1 D_NUMNSG*13, D_NUMSSG*13, D_NUMEBG*17, D_NUMNCM*11, D_NUMECM*11,
     1 D_NCNTCM*9, D_ECNTCM*9, V_EBIDSG*17, V_NSIDSG*13, V_NSNDCG*18,
     1 V_NSDFCG*16, V_SSIDSG*13, V_SSSDCG*18, V_SSDFCG*16, V_FTYPE*9,
     1 V_CMAP*8, V_NMAPI*9, V_NMAPII*13, V_NMAPB*9, V_NMAPBI*13, V_NMAPE*9,
     1 V_NMAPEI*13, V_EMAPI*9, V_EMAPII*13, V_EMAPB*9, V_EMAPBI*13,
     1 V_INTNST*10, V_BORNST*10, V_EXTNST*10, V_INTEST*10, V_BOREST*10,
     1 V_NCOMID*10, V_NCOMST*11, V_NCOMII*15, V_ECOMID*10, V_ECOMST*11,
     1 V_ECOMII*15, V_NCOMNID*11, V_NCOMPR*11, V_NCOMDI*15, V_ECOMEID*11,
     1 V_ECOMSID*11, V_ECOMPR*11, V_ECOMDI*15, END_CMIDX*9, END_CMIDS*4

C
C        Define dimensions
C
      parameter (D_NUMINOD = 'num_int_node', D_NUMBNOD = 'num_bor_node')
      parameter (D_NUMENOD = 'num_ext_node', D_NUMIELM = 'num_int_elem')
      parameter (D_NUMBELM = 'num_bor_elem', D_NUMPROC = 'num_processors')
      parameter (D_NUMPR_F = 'num_procs_file', D_NUMNODG = 'num_nodes_global')
      parameter (D_NUMELMG = 'num_elems_global', D_NUMNSG = 'num_ns_global')
      parameter (D_NUMSSG = 'num_ss_global', D_NUMEBG = 'num_el_blk_global')
      parameter (D_NUMNCM = 'num_n_cmaps', D_NUMECM = 'num_e_cmaps')
      parameter (D_NCNTCM = 'ncnt_cmap', D_ECNTCM = 'ecnt_cmap')

C
C       Define variables
C
      parameter (V_EBIDSG = 'el_blk_ids_global', V_EBCNTG = 'el_blk_cnt_global')
      parameter (V_NSIDSG = 'ns_ids_global', V_NSNDCG = 'ns_node_cnt_global')
      parameter (V_NSDFCG = 'ns_df_cnt_global', V_SSIDSG = 'ss_ids_global')
      parameter (V_SSSDCG = 'ss_side_cnt_global', V_SSDFCG = 'ss_df_cnt_global')
      parameter (V_FTYPE = 'nem_ftype', V_CMAP = 'comm_map')
      parameter (V_NMAPI = 'node_mapi', V_NMAPII = 'node_mapi_idx')
      parameter (V_NMAPB = 'node_mapb', V_NMAPBI = 'node_mapb_idx')
      parameter (V_NMAPE = 'node_mape', V_NMAPEI = 'node_mape_idx')
      parameter (V_EMAPI = 'elem_mapi', V_EMAPII = 'elem_mapi_idx')
      parameter (V_EMAPB = 'elem_mapb', V_EMAPBI = 'elem_mapb_idx')
      parameter (V_INTNST = 'int_n_stat', V_BORNST = 'bor_n_stat')
      parameter (V_EXTNST = 'ext_n_stat', V_INTEST = 'int_e_stat')
      parameter (V_BOREST = 'bor_e_stat', V_NCOMID = 'n_comm_ids')
      parameter (V_NCOMST = 'n_comm_stat', V_NCOMII = 'n_comm_info_idx')
      parameter (V_ECOMID = 'e_comm_ids', V_ECOMST = 'e_comm_stat')
      parameter (V_ECOMII = 'e_comm_info_idx', V_NCOMNID = 'n_comm_nids')
      parameter (V_NCOMPR = 'n_comm_proc', V_NCOMDI = 'n_comm_data_idx')
      parameter (V_ECOMEID = 'e_comm_eids', V_ECOMSID = 'e_comm_sids')
      parameter (V_ECOMPR = 'e_comm_proc', V_ECOMDI = 'e_comm_data_idx')

C
C       Define misc parameters
C
      parameter (END_CMIDX = '_info_idx', END_CMIDS = '_ids')
