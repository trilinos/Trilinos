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

C $Log: scolor.f,v $
C Revision 1.3  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:22  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:11:17  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1993/10/04  19:20:04  gdsjaar
c Fixed up message formats
c
c Revision 1.1  1993/09/15  20:44:00  gdsjaar
c Added shading (hidden 5) option to blot.  Simple illumination model.
c
c
C============================================================================
      SUBROUTINE SCOLOR (INIT, INLINE, IFLD, INTYP, RFIELD, IFIELD,
     &   CFIELD, SHDCOL, ISHDCL, IDELB)
C============================================================================

C   --*** BCOLOR ***  (BLOT) Process BLKCOL command
C   --   Written by John Glick - 11/22/88
C   --
C   --Parameters:
C   --   INIT  - IN - .TRUE. iff the purpose of the call is to
C   --            initialize the BLKCOL array and not to modify it.
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   IFLD, INTYP, IFIELD, CFIELD, - IN/OUT - the free-field reader
C   --          index and charcter field.
C   --   BLKCOL - IN/OUT - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
      LOGICAL INIT
      CHARACTER*(*) INLINE(*)
      INTEGER IFLD, INTYP(*), IFIELD(*)
      REAL RFIELD(*)
      CHARACTER*(*) CFIELD(*)
      REAL SHDCOL(7, *)
      INTEGER ISHDCL(3, *)
      INTEGER IDELB(*)

      include 'dbnums.blk'
      COMMON /BCOLR/ BCOLCH
      LOGICAL BCOLCH

      LOGICAL FFMATC, FFEXST
      INTEGER IDCOL, LOCSTR

      CHARACTER*8 COLA, COLF

      include 'shades.blk'

C *****************************************************************

      BCOLCH = .TRUE.

      IF (INIT) THEN

         DO 100 I = 1, NELBLK
            ISHDCL(1,I) = -1
            ISHDCL(2,I) =  0
            ISHDCL(3,I) =  0
  100    CONTINUE

      ELSE

C              If there are no fields on the BLKCOL command line,
C              toggle the ON/OFF flag.

         IF (INTYP(IFLD) .LE. -1) THEN
           DO 110 I = 1, NELBLK
            ISHDCL(1,I) = -ISHDCL(1,I)
 110      CONTINUE

C              Check for ON flag.

         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'OFF', 2)) THEN
           do 120 I = 1, NELBLK
             ISHDCL(1, I) = -1
             ISHDCL(2,I) =  0
             ISHDCL(3,I) =  0
 120       CONTINUE
           CALL FFADDC ('OFF', INLINE(1))

C              Check for RESET flag.

         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'RESET', 2)) THEN
           CALL FFADDC ('RESET', INLINE(1))
           do 130 I = 1, NELBLK
             ISHDCL(1, I) = -1
             ISHDCL(2,I) =  0
             ISHDCL(3,I) =  0
 130       CONTINUE
           
         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'SHOW', 2)) THEN
           CALL FFADDC ('SHOW', INLINE(1))
           CALL SHOCMD ('Valid predefined colors', SHDLST)
           
         ELSE

           ibeg = 0
           iend = 0

           IF (FFMATC (IFLD, INTYP, CFIELD, 'ALL', 2)) THEN
             CALL FFADDC ('ALL', INLINE(1))
             IBEG = 1
             IEND = NELBLK
             
           ELSE IF (FFEXST (IFLD, INTYP)) THEN

             CALL FFINTG (IFLD, INTYP, IFIELD,
     *         'block id', 0, IDBLK, *170)
             CALL FFADDI (IDBLK, INLINE)
             iblk = locint (idblk, nelblk, idelb)
             if (iblk .eq. 0) then
               write (*, 900) idblk
 900           FORMAT (1x, 'Material ID', I5,
     *           ' is not a valid material id.')
               go to 170
             end if
             IBEG = IBLK
             IEND = IBLK
           end if

C                 Get color to assign to the specified blocks.
C                    Check that next field has characters in it.
             
           if (ibeg .ne. 0  .and. iend .ne. 0) then
             IF (INTYP(IFLD) .EQ. 0) THEN
               COLA = CFIELD(IFLD)
               IFLD = IFLD + 1
               CALL ABRSTR (COLF, COLA, SHDLST)
               IF (COLF .EQ. ' ') THEN
                 WRITE (*, 10000) COLA
10000            FORMAT (1X, A, ' not a valid color name.',
     &             ' Rest of command not processed.')
                 GO TO 170
               ELSE
                 IDCOL = LOCSTR (COLF, NCLSHD, SHDLST)
                 RMULT = shades(1,IDCOL)
                 GMULT = shades(2,IDCOL)
                 BMULT = shades(3,IDCOL)
               END IF
               CALL FFADDC (COLF, INLINE)
             ELSE
C ... Color is specified by RGB components.
               CALL FFREAL (IFLD, INTYP, RFIELD,
     *           'Red Multiplier', 1.0, RMULT, *170)
               CALL FFADDR (RMULT, INLINE(1))
               CALL FFREAL (IFLD, INTYP, RFIELD,
     *           'Green Multiplier', 1.0, GMULT, *170)
               CALL FFADDR (GMULT, INLINE(1))
               CALL FFREAL (IFLD, INTYP, RFIELD,
     *           'Blue Multiplier', 1.0, BMULT, *170)
               CALL FFADDR (BMULT, INLINE(1))
             ENDIF
           ELSE
             WRITE (*, 10020)
10020        FORMAT (1X,
     *         'No color specified following block id specifications')
             GO TO 170
           ENDIF

           call ffintg(ifld, intyp, ifield,
     *       'number of colors', 0, NCOL, *170)
           call ffaddi (NCOL, INLINE(1))
C ... Get diffuse and specular values.
           CALL FFREAL (IFLD, INTYP, RFIELD,
     *       'Diffuse Proportion', 0.0, PDIFF, *170)
           CALL FFADDR (PDIFF, INLINE(1))
           CALL FFREAL (IFLD, INTYP, RFIELD,
     *       'Specular Proportion', 0.0, PSPEC, *170)
           CALL FFADDR (PSPEC, INLINE(1))
           CALL FFREAL (IFLD, INTYP, RFIELD,
     *       'Specular Exponent', 0.0, SPEXP, *170)
           CALL FFADDR (SPEXP, INLINE(1))
           
           do 140 iblk = ibeg, iend
             ISHDCL(1, IBLK) = 1
             ISHDCL(2, IBLK) = NCOL
C ... RGB Components of color
             SHDCOL(1, IBLK) = RMULT
             SHDCOL(2, iblk) = GMULT
             SHDCOL(3, iblk) = BMULT
C ... Diffuse and Specular Values, Specular Exponent.
             SHDCOL(4, IBLK) = PDIFF
             SHDCOL(5, iblk) = PSPEC
             SHDCOL(6, iblk) = SPEXP
 140       continue
         ENDIF
         
       ENDIF

 170  CONTINUE
      RETURN
      END
