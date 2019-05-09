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
      SUBROUTINE GTMVAR (A, IVAR, IELBLK, INSTEP, LENVAR, VAR)
C=======================================================================

C   --*** GTMVAR *** (BLOT) Read variable
C   --   Written by Amy Gilkey - revised 05/17/88
C   --
C   --GTMVAR returns the values for the requested variable for the
C   --requested time step.  It either reads the values from the sequential
C   --database file and writes them to a direct access scratch file or it
C   --reads the values from the direct access scratch file.
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   NUMELB - the number of elements per element block
C   --   ISEVOK - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   WHOTIM - true iff whole (versus history) time step
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   IVAR - IN - the variable index
C   --      (0 to initialize random file only)
C   --   IELBLK - IN - the element block number, <=0 for all
C   --      (for element blocks only)
C   --   INSTEP - IN - the time step number
C   --      = +n to read time step n
C   --      = -n to transfer time step n to random file only
C   --      =  0 to transfer all time steps to random file
C   --   LENVAR - IN - the length of VAR
C   --   VAR - OUT - the variable values (indeterminate if INSTEP <= 0)
C   --
C   --Common Variables:
C   --   Uses NDB of /DBASE/
C   --   Uses NUMNP, NUMEL, NELBLK, NVARHI, NVARGL, NVARNP, NVAREL,
C   --      NSTEPS, NSTEPW of /DBNUMS/
C   --
C   --Database is rewound upon the first entry of this routine; upon
C   --exit a flag is set to keep track of the database position; the
C   --database should not be moved between calls to this routine.
C   --
C   --A scratch random file is created and read and written in this routine.
C   --It is connected to unit 90.

      DIMENSION A(*)
      REAL VAR(*)

      IF (IVAR .LE. 0) RETURN
      CALL RNDVAR (A, A, A, IVAR, IELBLK, INSTEP, LENVAR, VAR)

      RETURN
      END
