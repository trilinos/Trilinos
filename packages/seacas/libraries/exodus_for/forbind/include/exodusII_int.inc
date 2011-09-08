C*****************************************************************************
C Copyright 1994 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
C
C See the COPYRIGHT file for copying and redistribution conditions.
C
C*****************************************************************************
C
C EXODUSII_INT - exodus II internals include file
C
C author - Victor R. Yarberry, Sandia National Laboratories
C
C environment - UNIX F77
C
C revision history - 
C
C  Based on original written by Lynn Clements (RE/SPEC)
C
C  $Id: exodusII_int.inc,v 1.8 2008/07/24 21:40:41 gdsjaar Exp $
C
C*****************************************************************************
C  --
C  --  This file contains define constants that are used to refer to
C  --  netCDF variables, attributes, or dimensions in which the EXODUS
C  --  II data are stored.  Using the defined constants will allow the
C  --  names of the netCDF entities to be changed easily in the future
C  --  if needed.  
C  --
C  --Constants:
C  --
C  --ATTRIBUTES
C  --
C  --   A_TYPE     The EXODUS II file type:
C                    "r" = regular
C                    "h" = history
C  --   A_TITLE    The database title
C  --   A_VERSION  The EXODUS II version number
C  --   A_NAMELB   The element type name for a given element block
C  --
C  --DIMENSIONS
C  --
C  --   D_NUMNOD   The number of nodes
C  --   D_NDIM     The number of coordinates per node
C  --   D_NUMEL    The number of elements
C  --   D_NELBLK   the number of element blocks
C  --   D_NUMELB   The number of elements in a given element block
C  --   D_NELNOD   The number of nodes per element in a given element 
C  --              block
C  --   D_NATRIB   The number of attributes in a given element block
C  --   D_NUMESS   The number of side sets
C  --   D_NEESS    The number of elements in a given side set
C  --   D_NNESS    The number of nodes in a given side set
C  --   D_NUMNPS   The number of node sets
C  --   D_NNNPS    The number of nodes in a given node set
C  --   D_NQAREC   The number of QA records
C  --   D_NINFO    The number of information records
C  --   D_NVARGL   The number of global variables
C  --   D_NVARNP   The number of nodal variables
C  --   D_NVAREL   The number of element variables
C  --   D_NVARHI   The number of history variables
C  --   D_N4       General dimension of length 4
C  --   D_STRLEN   General dimension of length MXSTLN used for lengths 
C  --              of names
C  --   D_LINLEN   General dimension of length MXLNLN used for lengths 
C  --              of long strings
C  --   D_TIMSTP   The number of whole time steps
C  --   D_HTMSTP   The number of history time steps 
C  --
C  --VARIABLES
C  --
C  --   V_COORD    The nodal coordinates
C  --   V_NAMECO   The names of coordinates
C  --   V_IDEBLK   The element block IDs
C  --   V_ICONK    The element connectivity for a given element block
C  --   V_ATRIBK   The list of attributes for a given element block
C  --   V_MAP      The element order map
C  --   V_IDESS    The side set Ids
C  --   V_FACESS   The distribution factor for each node in a given 
C  --              side set
C  --   V_ELSSET   The list of elements in a given side set
C  --   V_NDSSET   The list of nodes in a given side set
C  --   V_IDNPS    The node set IDs
C  --   V_NODSET   The list of nodes in a given node set
C  --   V_FNPS     The list of distribution factors in a given node set
C  --   V_INFO     The information records
C  --   V_TIMEHS   The simulation times for history time steps
C  --   V_TIME     The simulation times for whole time steps
C  --   V_QATITL   The QA records
C  --   V_ISEVOK   The element variable truth table
C  --   V_NAMEGV   The names of global variables
C  --   V_VARGL    The values of global variables
C  --   V_NAMENV   The names of nodal variables
C  --   V_VARNP    The values of nodal variables
C  --   V_NAMEEV   The names of element variables
C  --   V_VAREL    The values of element variables for a given element 
C  --              block and given element within that block
C  --   V_NAMEHV   The names of the history variables
C  --   V_VARHI    The values of the history variables
C  --

      character A_TYPE*4, A_TITLE*5, A_VERSION*7, A_NAMELB*9, 
     1 D_NUMNOD*9, D_NDIM*7, D_NUMEL*8, D_NELBLK*10, D_NUMELB*13, 
     1 D_NELNOD*14, D_NATRIB*14, D_NUMESS*13, D_NEESS*9, D_NNESS*10, 
     1 D_NUMNPS*13, D_NNNPS*10, D_NQAREC*10, D_NINFO*8, D_NVARGL*11, 
     1 D_NVARNP*11, D_NVAREL*12, D_NVARHI*11, D_N4*4, D_STRLEN*10, 
     1 D_LINLEN*8, V_COORD*5, V_NAMECO*10, V_IDEBLK*9, V_ICONK*7, 
     1 V_ATRIBK*6, V_MAP*8, V_IDESS*12, V_FACESS*12, V_ELSSET*7, 
     1 V_NDSSET*7, V_IDNPS*12, V_NODSET*7, V_FNPS*12, V_INFO*12, 
     1 V_TIMEHS*9, V_TIME*10, V_QATITL*10, V_ISEVOK*12, V_NAMEGV*12, 
     1 V_VARGL*12, V_NAMENV*12, V_VARNP*12, V_NAMEEV*13, V_VAREL*13, 
     1 V_NAMEHV*12, V_VARHI*12, D_TIMSTP*9, D_HTMSTP*14 

c
c         Define attributes 
c 
      parameter (A_TYPE = 'type', A_TITLE = 'title') 
      parameter (A_VERSION = 'version') 
c 
c         Define dimensions 
c
      parameter (D_NUMNOD = 'num_nodes', D_NDIM = 'num_dim')
      parameter (D_NUMEL ='num_elem')
      parameter (D_NELBLK = 'num_el_blk', D_NUMELB ='num_el_in_blk')
      parameter (D_NELNOD = 'num_nod_per_el') 
      parameter (D_NATRIB = 'num_att_in_blk')
      parameter (D_NUMESS = 'num_side_sets', D_NEESS = 'num_el_ss')
      parameter (D_NNESS = 'num_nod_ss', D_NUMNPS = 'num_node_sets')
      parameter (D_NNNPS = 'num_nod_ns', D_NQAREC = 'num_qa_rec')
      parameter (D_NINFO = 'num_info', D_NVARGL ='num_glo_var')
      parameter (D_NVARNP = 'num_nod_var', D_NVAREL ='num_elem_var')
      parameter (D_NVARHI = 'num_his_var', D_N4 = 'four')
      parameter (D_STRLEN = 'len_string', D_LINLEN = 'len_line')
      parameter (D_TIMSTP = 'time_step', D_HTMSTP ='hist_time_step')
c
c         Define variables
c
      parameter (V_COORD = 'coord', V_NAMECO = 'coor_names')
      parameter (V_IDEBLK = 'el_blk_id', V_ICONK = 'connect')
      parameter (V_ATRIBK = 'attrib', A_NAMELB = 'elem_type')
      parameter (V_MAP ='elem_map', V_IDESS = 'side_set_ids') 
      parameter (V_FACESS = 'dist_fact_ss')
      parameter (V_ELSSET = 'elem_ss', V_NDSSET = 'node_ss')
      parameter (V_IDNPS = 'node_set_ids', V_NODSET ='node_ns')
      parameter (V_FNPS = 'dist_fact_ns', V_INFO = 'info_records') 
      parameter (V_TIMEHS = 'time_hist')
      parameter (V_TIME = 'time_whole', V_QATITL = 'qa_records')
      parameter (V_ISEVOK = 'elem_var_tab', V_NAMEGV ='name_glo_var')
      parameter (V_VARGL ='vals_glo_var', V_NAMENV = 'name_nod_var')
      parameter (V_VARNP = 'vals_nod_var')
      parameter (V_NAMEEV = 'name_elem_var', V_VAREL = 'vals_elem_var')
      parameter (V_NAMEHV = 'name_his_var', V_VARHI = 'vals_his_var')
c
c        dummy character variables used through out EXODUS
c
      character name*32

c
c        dummy integer variables used throughout EXODUS
c
       integer istart(3), icount(3), id
c
c        declaration for various functions
c
       integer chkfil


