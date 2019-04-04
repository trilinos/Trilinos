C    Copyright(C) 2008-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C
C=======================================================================
      SUBROUTINE RDSTEP (ISTEP, TIME, NUMELB, IDELB, ISEVOK,
     &                   VISELB, MAXNE, VARVAL, MERR)
C=======================================================================

C   --*** RDSTEP *** (ALGEBRA) Read database time step variables
C   --   Written by Amy Gilkey - revised 11/30/87
C   --   Modified 8/30/95
C   --
C   --RDSTEP reads the input database variables for one time step and
C   --stores the ones used in the equations in array VARVAL.
C   --
C   --Parameters:
C   --   ISTEP  - IN  - the time step number
C   --   TIME   - IN  - the time step time
C   --   NUMELB - IN  - the number of elements per block
C   --   ISEVOK - IN  - the element block variable truth table;
C   --                  variable i of block j exists iff ISEVOK(j,i)
C   --   MAXNE  - IN  - the VARVAL dimension
C   --   VARVAL - OUT - the input data needed
C   --   MERR   - OUT - error flag
C   --
C   --Common Variables:
C   --   Uses ITIME, IGVBEG, INVBEG, IEVBEG, IGVEND, INVEND, IEVEND
C   --      of /DBXVAR/

      include 'ag_namlen.blk'
      include 'ag_var.blk'
      include 'ag_dbase.blk'
      include 'ag_dbnums.blk'
      include 'ag_dbxvar.blk'

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      INTEGER NUMELB(*)
      INTEGER IDELB(*)
      LOGICAL ISEVOK(NELBLK*NVAREL)
      LOGICAL VISELB(NELBLK)
      REAL    VARVAL(MAXNE,*)
      INTEGER MERR
      MERR = 0

C   --Assign time to time/history/globals entry
      VARVAL(IDVAR(ITIME),ISTVAR(ICURTM,ITIME)) = TIME

C      --Read global variables

      IF (IGVBEG .LE. IGVEND) THEN
         CALL STORE (ISTEP, 'G', IGVBEG, IGVEND, NVARGL,
     &        NUMELB, IDELB, ISEVOK, VISELB, MAXNE, VARVAL, MERR)
         IF (MERR .EQ. 1) RETURN
      END IF

C      --Read nodal variables

      IF (INVBEG .LE. INVEND) THEN
         CALL STORE (ISTEP, 'N', INVBEG, INVEND, NUMNP,
     &        NUMELB, IDELB, ISEVOK, VISELB, MAXNE, VARVAL, MERR)
         IF (MERR .EQ. 1) RETURN
      END IF

C      --Read element variables

      IF (IEVBEG .LE. IEVEND) THEN
         CALL STORE (ISTEP, 'E', IEVBEG, IEVEND, NUMEL,
     &        NUMELB, IDELB, ISEVOK, VISELB, MAXNE, VARVAL, MERR)
         IF (MERR .EQ. 1) RETURN
      END IF

      RETURN
      END
