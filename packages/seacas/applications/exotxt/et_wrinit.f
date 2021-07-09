C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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

C ... Some codes are embedding carriage returns in title.
C     Strip them out...
      LTITLE = lenstr(title)
      do i=1, ltitle
        if (ichar(title(i:i)) .eq. 10) title(i:i) = ' '
      end do

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
