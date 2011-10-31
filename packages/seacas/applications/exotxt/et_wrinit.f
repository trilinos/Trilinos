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


C=======================================================================
      SUBROUTINE WRINIT (NTXT, VERS, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &                   NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL,
     &                   LESSDF, QAINFO, NAMLEN)
C=======================================================================

C   --*** WRINIT *** (EXOTXT) Write database title and initial variables
C   --   Written by Amy Gilkey - revised 12/04/87
C   --   Modified for ExodusIIv2 database format 10/12/95
C   --
C   --WRINIT writes the title and the initial variables from the database.
C   --
C   --Parameters:
C   --   NTXT   - IN - the text file
C   --   VERS   - IN - the version number
C   --   TITLE  - IN - the database title
C   --   NDIM   - IN - the number of coordinates per node
C   --   NUMNP  - IN - the number of nodes
C   --   NUMEL  - IN - the number of elements
C   --   NELBLK - IN - the number of element blocks
C   --   NUMNPS - IN - the number of nodal point sets
C   --   LNPSNL - IN - the length of the nodal point sets node list
C   --   LNPSDF - IN - the length of the node sets distribution factors list
C   --   NUMESS - IN - the number of side sets
C   --   LESSEL - IN - the length of the concatenated side sets element list
C   --   LESSNL - IN - the length of the concatenated side sets node list
C   --   LESSDF - IN - the length of the side set distribution factors list
C   --   QAINFO - IN - program information array
C   --
C     Header info
      INTEGER NTXT
      REAL    VERS
      CHARACTER*80 TITLE
      INTEGER NDIM, NUMNP, NUMEL, NELBLK
C     Node set info
      INTEGER NUMNPS, LNPSNL, LNPSDF
C     Side set info
      INTEGER NUMESS, LESSEL, LESSNL, LESSDF

      CHARACTER*(*) QAINFO(6)

      WRITE (NTXT, 10030) '! Database Title', (QAINFO(I),I=1,3)
      WRITE (NTXT, '(A)') TITLE

      WRITE (NTXT, '(A)') '! Database initial variables'
      WRITE (NTXT, '(I10, F10.2, I10,5X, A)') NDIM, VERS,NAMLEN,
     &   '! dimensions, version number, name length'
      WRITE (NTXT, 10010) NUMNP, NUMEL, NELBLK,
     &   '! nodes, elements, element blocks'
      WRITE (NTXT, 10020) NUMNPS, NUMESS,
     &   '! #node sets, #side sets'
      WRITE (NTXT, 10020) LNPSNL, LNPSDF,
     &   '! len: node set list, dist fact length'
      WRITE (NTXT, 10010) LESSEL, LESSNL, LESSDF,
     &   '! side sets len: element, node , dist fact'

      RETURN
10010  FORMAT (3I10, 5X, A)
10020  FORMAT (2I10, 15X, A)
10030  FORMAT (A, 5X, 3A32)
      END
