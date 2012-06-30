C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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

C $Log: grabrt.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:01:55  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:51:17  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      LOGICAL FUNCTION GRABRT ()
C=======================================================================

C   --*** GRABRT *** (GRPLIB) Check cancel function (PLT)
C   --   Written by Amy Gilkey - revised 04/28/87
C   --
C   --GRABRT checks the cancel flag.  If it is set, it aborts pending
C   --graphics and returns the terminal to alphanumeric mode.
C   --In any case, the value of the cancel flag is returned as the
C   --function value.

C   --Routines Called:
C   --   CPUIFC - (PLTLIB) Check interrupt flag
C   --   PLTBEL - (PLTLIB) Ring bell
C   --   PLTBGN - (PLTLIB) Erase display surface
C   --   PLTFLU - (PLTLIB) Flush graphics buffer and set alphanumeric mode
C   --   GRSNAP - (GRPLIB) Handle device frame snapping

      LOGICAL CPUIFC

C   --Return cancel flag
      GRABRT = CPUIFC (.FALSE.)

      IF (GRABRT) THEN
C      --Abort plot that is being snapped
         CALL GRSNAP ('ABORT', 0)

C      --Flush buffer (do not erase screen), ring bell, set alphanumeric mode
         CALL PLTFLU
         CALL PLTBEL
         CALL PLTFLU

C      --Send abort message
         WRITE (*, '(1X, A)') '*** Plot set aborted ***'
      END IF

      RETURN
      END
