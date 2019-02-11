C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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
C
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
C

C $Id: pltd2g.f,v 1.1 1993/07/16 16:47:57 gdsjaar Exp $
C $Log: pltd2g.f,v $
C Revision 1.1  1993/07/16 16:47:57  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      LOGICAL FUNCTION PLTD2G(XD,YD,XG,YG)
      DIMENSION UMAP(14)

      CALL PLTGTG(27,UMAP)
      CALL PLTGTG(9,TYPE)
      PLTD2G = .FALSE.
      TOP = XD*UMAP(4) - YD*UMAP(3) - UMAP(5)*UMAP(4) + UMAP(6)*UMAP(3)
      BOTTOM = UMAP(1)*UMAP(4) - UMAP(2)*UMAP(3)
      XG = TOP/BOTTOM
      TOP = XD*UMAP(2) - YD*UMAP(1) - UMAP(5)*UMAP(2) + UMAP(6)*UMAP(1)
      BOTTOM = UMAP(3)*UMAP(2) - UMAP(4)*UMAP(1)
      YG = TOP/BOTTOM
      IF (TYPE.EQ.2.) THEN
         IF (XD.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTD2G',
     *                 'Cannot convert a negative number on log x axis.'
     *                  ,2)
            RETURN

         END IF

         XG = 10.**XG

      ELSE IF (TYPE.EQ.3.) THEN
         IF (YD.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTD2G',
     *                 'Cannot convert a negative number on log y axis.'
     *                  ,2)
            RETURN

         END IF

         YG = 10.**YG

      ELSE IF (TYPE.EQ.4.) THEN
         IF (XD.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTD2G',
     *                 'Cannot convert a negative number on log x axis.'
     *                  ,2)
            RETURN

         END IF

         IF (YD.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTD2G',
     *                 'Cannot convert a negative number on log y axis.'
     *                  ,2)
            RETURN

         END IF

         XG = 10.**XG
         YG = 10.**YG
      END IF

      PLTD2G = .TRUE.
      RETURN

      END
