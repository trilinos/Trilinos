C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C=======================================================================
      SUBROUTINE DBIEBI (NDB, OPTION, IELB, NUMELB, NUMLNK, NUMATR,
     &                   LINK, ATRIB, NATRDM, NLNKDM, *)
C=======================================================================
C$Log: dbieb1.f,v $
CRevision 1.5  2009/03/25 12:36:43  gdsjaar
CAdd copyright and license notice to all files.
CPermission to assert copyright has been granted; blot is now open source, BSD
C
CRevision 1.4  1998/07/15 14:44:13  gdsjaar
CGeneral cleanup, remove compiler warnings
C
CRevision 1.3  1997/11/11 14:55:54  gdsjaar
CAdded 'external blkdat' to main program to ensure that the block data
Cgets linked into the executable. Wasn't happening on dec alpha
Csystems.
C
CRemoved unreachable lines in several routines
C
CFixed variable name spelling in contor.f
C
CUnsplit strings that were split across lines
C
CRemoved old error variables left over from exodusIIv1
C
CUpped version number
C
CRevision 1.2  1996/06/21 16:07:06  caforsy
CRan ftnchek and removed unused variables.  Reformat output for list
Cvar, list global, and list name.
C
CRevision 1.1  1994/04/07 19:57:50  gdsjaar
CInitial checkin of ACCESS/graphics/blotII2
C
C Revision 1.4  1993/07/28  18:51:58  gdsjaar
C Rename dbieb1 to dbiebi to reduce linker confusion with suplib routines.
C
C Revision 1.3  1993/07/27  20:27:19  gdsjaar
C Added error checking, set to non-verbose exodus II file opening,
C cleaned up dbiv0 routine, removed some unused variables.
C
C Revision 1.2  1993/07/27  14:54:43  gdsjaar
C Cleaned out commented-out code, fixed element block name reading.
C
c Revision 1.1  1993/07/27  13:30:32  gdsjaar
c Initial checkin of ACCESS/graphics/blotII
c
c Revision 1.1.1.1  1990/08/14  16:12:34  gdsjaar
c Testing
c
c Revision 1.1  90/08/14  16:12:33  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:08  gdsjaar
c Initial revision
c

C   --*** DBIEB1 *** (EXOLIB) Read database element block misc.
C   --   Written by Amy Gilkey - revised 10/14/87
C   --   Modified by Greg Sjaardema - 8/8/90
C   --      ---Removed MAX from Dimension statements, Added NATRDM, NLNKDM
C   --
C   --DBIEB1 reads the element block connectivity and attribute information
C   --from the database.  An error message is displayed if the end of file
C   --is read.
C   --
C   --Parameters:
C   --   NDB    - IN  - the database file
C   --   OPTION - IN  - ' ' to not store, '*' to store all, else store options:
C   --                  'C' to store connectivity
C   --                  'A' to store attributes
C   --   IELB   - IN  - the element block number
C   --   NUMELB - IN  - the number of elements in the block
C   --   NUMLNK - IN  - the number of nodes per element;
C   --                  negate to not store connectivity
C   --   NUMATR - IN  - the number of attributes;
C   --                  negate to not store attributes
C   --   LINK   - OUT - the element connectivity for this block
C   --   ATRIB  - OUT - the attributes for this block
C   --	 NATRDM - IN  - dimension of atrib array
C   --   NLNKDM - IN  - dimension of link array
C   --   *      - OUT - return statement if end of file or read error

      include 'exodusII.inc'
      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NUMELB, NUMLNK, NUMATR
      INTEGER LINK(NLNKDM, *)
      REAL ATRIB(NATRDM,*)

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
         call exgelc(ndb, ielb, link, ierr)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0)) THEN
        if (numatr .gt. 0) then
           call exgeat(ndb, ielb, atrib, ierr)
        end if
      END IF

      RETURN

      END
