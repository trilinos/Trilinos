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

C $Log: qnpick.f,v $
C Revision 1.2  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:08:46  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:55:44  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE QNPICK (TYPE, ISUP, ISDEF,
     &  A, KXN, KYN, KZN, KHIDEN, KNPSUR)
C=======================================================================

C   --*** QNPICK *** (MESH) Query last mesh display for node information
C   --   Written by Amy Gilkey - revised 04/06/88
C   --
C   --QNPICK returns the address of the last displayed mesh characteristics
C   --for cursor input.  This routine assumes that the characteristics
C   --are in dynamic memory.
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   XN, YN, ZN - the nodal coordinates (ZN for 3D only)
C   --   DHZ, DVT, DPD - the deformed nodal coordinates (DPD for 3D only)
C   --   HIDENP - true iff node i is hidden (3D only)
C   --   NPSURF - the node numbers of the surface nodes (3D only)
C   --
C   --Parameters:
C   --   TYPE - IN - the type of mesh information needed:
C   --      'QUERY'     - return ISUP and ISDEF only
C   --      'ORIGINAL'  - original mesh coordinates
C   --      'DISPLAYED' - last displayed mesh coordinates
C   --   ISUP - OUT - true iff a mesh is displayed
C   --   ISDEF - OUT - true iff the displayed mesh is deformed
C   --   A - IN - the dynamic memory base array
C   --   KXN, KYN, KZN - OUT - the dynamic memory index of the coordinates
C   --      requested (KZN in 3D only)
C   --   KHIDEN - OUT - the dynamic memory index of the node hidden indicator
C   --      (3D only)
C   --   KNPSUR - OUT - the dynamic memory index of the surface node indices
C   --      (3D only)
C   --
C   --Common Variables:
C   --   Uses /PICK/
C   --   Uses IS3DIM of /D3NUMS/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      COMMON /PICK/   INITP, PKDEF,
     &  PKMESH(KTOP), PKRMAT(3,3), PKRCEN(3),
     &  DMESH(KTOP), DVRAT, DXMID, DYMID, DXLAST, DYLAST
      LOGICAL INITP, PKDEF

      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      LOGICAL ISUP, ISDEF
      CHARACTER*(*) TYPE

      ISUP = INITP
      IF (ISUP) THEN
        ISDEF = PKDEF
      ELSE
        ISDEF = .FALSE.
      END IF
      KXN = 1
      KYN = 1
      KZN = 1
      KHIDEN = 1
      KNPSUR = 1

      IF (TYPE .EQ. 'ORIGINAL') THEN
        CALL MDFIND ('XN', KXN, IDUM)
        CALL MDFIND ('YN', KYN, IDUM)
        IF (IS3DIM) CALL MDFIND ('ZN', KZN, IDUM)
      END IF

      IF (ISUP) THEN
        IF (TYPE .EQ. 'DISPLAYED') THEN
          IF (.NOT. PKDEF) THEN
            IF (.NOT. IS3DIM) THEN
              CALL MDFIND ('XN', KXN, IDUM)
              CALL MDFIND ('YN', KYN, IDUM)
              IF (IS3DIM) CALL MDFIND ('ZN', KZN, IDUM)
            ELSE
              CALL MDFIND ('HZ', KXN, IDUM)
              CALL MDFIND ('VT', KYN, IDUM)
              IF (IS3DIM) CALL MDFIND ('PD', KZN, IDUM)
            END IF
          ELSE
            CALL MDFIND ('DHZ', KXN, IDUM)
            CALL MDFIND ('DVT', KYN, IDUM)
            IF (IS3DIM) CALL MDFIND ('DPD', KZN, IDUM)
          END IF
          IF (IS3DIM) CALL MDFIND ('HIDENP', KHIDEN, IDUM)
          IF (IS3DIM) CALL MDFIND ('NPSURF', KNPSUR, IDUM)
        END IF
      END IF

      RETURN
      END
