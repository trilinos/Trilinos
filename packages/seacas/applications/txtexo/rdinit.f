C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
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

C $Id: rdinit.f,v 1.3 2007/10/17 18:47:22 gdsjaar Exp $
C=======================================================================
      SUBROUTINE RDINIT (NTXT, VERS, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL, LESSDF, NAMLEN,
     *  *)
C=======================================================================

C   --*** RDINIT *** (TXTEXO) Read database title and initial variables
C   --   Written by Amy Gilkey - revised 12/16/87
C   --
C   --RDINIT reads the title and the initial variables from the text file.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NVERS - OUT - the version number
C   --   TITLE - OUT - the database title
C   --   NDIM - OUT - the number of coordinates per node
C   --   NUMNP - OUT - the number of nodes
C   --   NUMEL - OUT - the number of elements
C   --   NELBLK - OUT - the number of element blocks
C   --   NUMNPS - OUT - the number of node sets
C   --   LNPSNL - OUT - the length of the node sets node list
C   --   NUMESS - OUT - the number of side sets
C   --   LESSEL - OUT - the length of the side sets element list
C   --   LESSNL - OUT - the length of the side sets node list
C   --   * - return statement if end of file or read error
C   --
C   --Database must be rewound before entry; upon exit positioned at end of
C   --initial variables.

      CHARACTER*80 TITLE
      CHARACTER*80 SCRATCH

      READ (NTXT, *, END=110, ERR=110)
      READ (NTXT, '(A)', END=110, ERR=110) TITLE

      READ (NTXT, *, END=120, ERR=120)

C ... Split scratch string into pieces.  The namlen does not exist on older 
C     databases, so we need to see if it is there or not.  The others are always there.
      READ (NTXT, '(A)', END=110, ERR=110) SCRATCH
      READ(SCRATCH( 1:10),*) NDIM
      READ(SCRATCH(11:20),*) VERS
      if (scratch(21:30) .ne. '          ') then
        READ(SCRATCH(21:30),*) NAMLEN
      else
        namlen = 32
      end if

      READ (NTXT, *, END=130, ERR=130) NUMNP, NUMEL, NELBLK
      READ (NTXT, *, END=150, ERR=150) NUMNPS, NUMESS
      READ (NTXT, *, END=160, ERR=160) LNPSNL, LNPSDF
      READ (NTXT, *, END=170, ERR=170) LESSEL, LESSNL, LESSDF

  100 CONTINUE
      RETURN

  110 CONTINUE
      CALL PRTERR ('FATAL', 'Reading TITLE')
      GOTO 180
  120 CONTINUE
      CALL PRTERR ('FATAL', 'Reading NUMBER OF DIMENSIONS or VERSION')
      GOTO 180
  130 CONTINUE
      CALL PRTERR ('FATAL',
     &     'Reading NUMBER OF NODES, ELEMENTS, and ELEMENT BLOCKS')
      GOTO 180
  150 CONTINUE
      CALL PRTERR ('FATAL',
     &   'Reading NUMBER OF NODE SETS and SIDE SETS')
      GOTO 180
  160 CONTINUE
      CALL PRTERR ('FATAL', 'Reading NODE SET NODE LENGTHS')
      GOTO 180
  170 CONTINUE
      CALL PRTERR ('FATAL',
     &     'Reading SIDE SET ELEMENT, NODE, and FACTOR COUNTS')
      GOTO 180
  180 CONTINUE
      RETURN 1
      END
