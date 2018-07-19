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

C $Log: msdone.f,v $
C Revision 1.2  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:05:23  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:53:44  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE MSDONE (ANYDEF, ANYUND,
     &   DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS)
C=======================================================================

C   --*** MSDONE *** (MESH) Clean up after mesh plot set
C   --   Written by Amy Gilkey - revised 05/26/88
C   --
C   --MSDONE cleans up the memory after a mesh plot set.  Many variables
C   --may be needed to do "picks" on the mesh plot.
C   --
C   --Parameters:
C   --   ANYDEF - IN - true iff any deformed mesh is to be plotted
C   --   ANYUND - IN - true iff any undeformed mesh is to be plotted
C   --   DOIXF - IN - true iff the IXFAC array is needed
C   --   DON2B - IN - true iff the IN2ELB array is needed
C   --   DOELED - IN - true iff the deformed element quarilateral centers
C   --      are needed
C   --   DOELEU - IN - true iff the undeformed element quarilateral centers
C   --      are needed
C   --   DODEAD - IN - true iff dead nodes are needed
C   --   DONPS - IN - true iff node set information is needed
C   --   DOESS - IN - true iff side set information is needed

      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      LOGICAL ANYDEF, ANYUND
      LOGICAL DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS

C   --Release the memory reserved for plot set

      IF (IS3DIM .AND. DOIXF) THEN
         CALL MDDEL ('IXFAC')
      END IF
      IF (DODEAD) THEN
         CALL MDDEL ('IDN2B')
      END IF
      IF (DONPS) THEN
         CALL MDDEL ('SCRNPS')
      END IF
      IF (DOESS) THEN
         CALL MDDEL ('SCRESS')
      END IF

      RETURN
      END
