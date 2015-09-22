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

C---5---10---15---20---25---30---35---40---45---50---55---60---65---70--
C $Log: grcolt.f,v $
C Revision 1.4  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.3  1998/07/07 14:07:30  gdsjaar
C Added variable to SAVE
C
C Modified to pass in a variable instead of a computed value since the
C argument is changed in the routine.
C
C Initialize variable
C
C Revision 1.2  1998/06/12 15:53:12  gdsjaar
C 1. Problem with TIMES array. Blot accesses a dummy timestep even if
C there were no timesteps on the database. Array wasn't allocated, so
C writing off into never-never land.
C
C 2. Inconsistency among some common blocks. Some places weren't using
C the include but had the definition hardwired in. Removed those.
C
C 3. Added 'EXTERNAL BLKDAT' to all routines that used data values set
C in BLKDAT
C
C 4. Cleanup of some A vs. IA argument passing.
C
C Revision 1.1  1994/04/07 20:02:13  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.6  1993/08/24  16:45:35  gdsjaar
c Added entries to common block for additional user-specified colormap
c selection
c
c Revision 1.5  1993/08/16  21:23:05  gdsjaar
c Added inverse color map to default map, added other color map files to
c Imakefile, upped version number.
c
c Revision 1.4  1993/08/16  19:42:41  gdsjaar
c Added new spectrum options: gray, rainbow, hot, terrain, zebra
c
c Revision 1.3  1993/03/01  16:57:21  gdsjaar
c jwswegl - added rainbow spectrum (closer to visible spectrum) specify
c with RAINBOW ncol.  Added CGLOBAL to toggle between global and local
c contour labelling.
c
c Revision 1.2  1990/12/14  08:51:31  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GRCOLT
C=======================================================================

C   --*** GRCOLT *** (GRPLIB) Set alternate color table (PLT)
C   --   Written by Amy Gilkey - revised 04/11/88
C   --
C   --GRCOLT sets the alternate color table (colors 8..) to a specified
C   --color table map.  The following maps are defined:
C   --       spectrum colors blue to red
C   --       rainbow colors violet to red
C   --
C   --Common Variables:
C   --   Uses DEVOK, ICURDV, NUMCOL, MAPALT of /GRPCOM/
C   --   Sets MAPUSE of /GRPCOM/

C   --Routines Called:
C   --   PLTCOL - (PLTLIB) Set color table color
C   --   PLTIQC - (PLTLIB) Get color table color

      EXTERNAL BLKDAT

      PARAMETER (KCOLOR=1)
      PARAMETER (KCOLIN=6, KCOSYM=44)

      include 'grpcom.blk'

C  Flag for rainbow spectrum
      include 'icrnbw.blk'

      REAL SRED(0:4), SGREEN(0:4), SBLUE(0:4)
      SAVE SRED, SGREEN, SBLUE
C      --SRED, SGREEN, SBLUE - color ratios for RGB settings; spectrum is
C      --   divided in 4 parts 

      REAL RRED(32:255), RGREEN(32:255), RBLUE(32:255)
      SAVE RRED, RGREEN, RBLUE
C  RRED, RGREEN, RBLUE - color ratios for RGB settings; rainbow is
C  from TOOPLOT color table; numbering 32-255 has been retained
C  in order to use original coding to generate table.
C  The original routineS used 0-255 for color definition values.
C  These are scaled to the range 0.-1. here.

      INTEGER LSTALT(2), LSTNUM(2)
      SAVE LSTALT, LSTNUM, SATLST, IRNBWL
C      --LSTALT - the color type last set for the device
C      --LSTNUM - the number of colors last set for the device

      DATA SRED   / 0.20, 0.50, 0.80, 0.85, 1.00 /
      DATA SGREEN / 0.60, 0.70, 0.80, 0.60, 0.10 /
      DATA SBLUE  / 1.00, 0.50, 0.25, 0.10, 0.10 /

C  INITIALIZE ALL RAINBOW COLORS TO ZERO - COLOR TABLE HAS NOT BEEN CREATED
      DATA RRED   / 224*0. /
      DATA RGREEN / 224*0. /
      DATA RBLUE  / 224*0. /

      DATA LSTALT / 0, 0 /
      DATA LSTNUM / 0, 0 /
      DATA SATLST / 0.0  /
      DATA IRNBWL / -1 /

      IF (MAPALT(ICURDV) .LE. 0) THEN
         CONTINUE

      ELSE IF ((MAPALT(ICURDV) .NE. LSTALT(ICURDV))
     &   .OR. (NUMCOL(MAPALT(ICURDV),ICURDV) .NE. LSTNUM(ICURDV))
     &   .OR.  IRNBWL.NE.IRAINB  .OR.  SATUR .NE. SATLST) THEN

        IF (MAPALT(ICURDV) .EQ. 1) THEN

C         --Set up a spectrum of colors in colors 8+

          IF(IRAINB.EQ.0 .AND. ISPEC .EQ. 0)THEN
C  Blue-brown-red spectrum (DEFAULT)

            IF (NUMCOL(1,ICURDV) .EQ. 1) THEN
               CALL PLTCOL (8+0, 0.0, 0.0, 1.0)

            ELSE IF (NUMCOL(1,ICURDV) .EQ. 2) THEN
              if (isinv .eq. 0) then
                CALL PLTCOL (8+0, 0.0, 0.0, 1.0)
                CALL PLTCOL (8+1, 1.0, 0.0, 0.0)
              else
                CALL PLTCOL (8+1, 0.0, 0.0, 1.0)
                CALL PLTCOL (8+0, 1.0, 0.0, 0.0)
              end if
            ELSE IF (NUMCOL(1,ICURDV) .EQ. 3) THEN
              if (isinv .eq. 0) then
                CALL PLTCOL (8+0, 0.0, 0.0, 1.0)
                CALL PLTCOL (8+1, 0.0, 1.0, 0.0)
                CALL PLTCOL (8+2, 1.0, 0.0, 0.0)
              else
                CALL PLTCOL (8+2, 0.0, 0.0, 1.0)
                CALL PLTCOL (8+1, 0.0, 1.0, 0.0)
                CALL PLTCOL (8+0, 1.0, 0.0, 0.0)
              end if
            ELSE IF (NUMCOL(1,ICURDV) .EQ. 5) THEN
              if (isinv .eq. 0) then
                CALL PLTCOL (8+0, 0.0, 0.0, 1.0)
                CALL PLTCOL (8+1, 0.0, 1.0, 1.0)
                CALL PLTCOL (8+2, 0.0, 1.0, 0.0)
                CALL PLTCOL (8+3, 1.0, 1.0, 0.0)
                CALL PLTCOL (8+4, 1.0, 0.0, 0.0)
              else
                CALL PLTCOL (8+4, 0.0, 0.0, 1.0)
                CALL PLTCOL (8+3, 0.0, 1.0, 1.0)
                CALL PLTCOL (8+2, 0.0, 1.0, 0.0)
                CALL PLTCOL (8+1, 1.0, 1.0, 0.0)
                CALL PLTCOL (8+0, 1.0, 0.0, 0.0)
              end if
            ELSE
              NEWCOL = NUMCOL(1,ICURDV) - 1
              DO 10 I = 0, NEWCOL
                FRAC = FLOAT(I) / FLOAT(NEWCOL)
                ISEG = INT (FRAC * 4)
                REM = FRAC * 4 - ISEG
                XRED   = SRED(ISEG)
     &            - REM * (SRED(ISEG) - SRED(ISEG+1))
                XGREEN = SGREEN(ISEG)
     &            - REM * (SGREEN(ISEG) - SGREEN(ISEG+1))
                XBLUE  = SBLUE(ISEG)
     &            - REM * (SBLUE(ISEG) - SBLUE(ISEG+1))
                if (isinv .eq. 0) then
                  CALL PLTCOL (8+I, XRED, XGREEN, XBLUE)
                else
                  CALL PLTCOL (8+NEWCOL-I, XRED, XGREEN, XBLUE)
                end if
 10           CONTINUE
            END IF

          ELSE IF(IRAINB.EQ.0 .AND. ISPEC .GE. 1)THEN
            NEWCOL = NUMCOL(1,ICURDV) - 1
            call textur (SATUR, NUMCOL(1,ICURDV), ISPEC, ISINV,
     *        RMULT, GMULT, BMULT)
            SATLST = SATUR
          ELSE 
C  Rainbow spectrum
           IF(RRED(32).EQ.0.)THEN

C  Color table has not been created. Generate it

C     Define colors 32-56 to be a linear variation from red to orange
C     Red    = 255,0,0    (R,G,B)
C     Orange = 255,168,0  (R,G,B)
C
             XIPCOL=168.0/24.0
             DO 50 I=32,56
              RRED(I)=1.
              RGREEN(I)=(I-32)*XIPCOL/255.
              RBLUE(I)=0.
   50        CONTINUE
C
C     Define colors 57-106 to be a linear variation from orange to yello
C     Orange = 255,168,0  (R,G,B)
C     Yellow = 255,255,0  (R,G,B)
C
             XIPCOL=87.0/50.0
             DO 60 I=57,106
              RRED(I)=1.
              RGREEN(I)=(168.0+(I-56)*XIPCOL)/255.
              RBLUE(I)=0.
   60        CONTINUE
C
C     Define colors 107-166 to be a linear variation from yellow to gree
C     Yellow = 255,255,0  (R,G,B)
C     Green  = 0,255,0    (R,G,B)
C
             XIPCOL=255.0/60.0
             DO 70 I=107,166
              RRED(I)=(255.0-(I-106)*XIPCOL)/255.
              RGREEN(I)=1.
              RBLUE(I)=0.
   70        CONTINUE
C
C     Define colors 167-210 to be a linear variation from green to blue
C     Green = 0,255,0  (R,G,B)
C     Blue  = 0,0,255  (R,G,B)
C
             XIPCOL=255.0/44.0
             DO 80 I=167,210
              RRED(I)=0.
              RGREEN(I)=(255.0-(I-166)*XIPCOL)/255.
              RBLUE(I)=(I-166)*XIPCOL/255.
   80        CONTINUE
C
C    Define colors 211-255 to be a linear variation from blue to violet
C    Blue   = 0,0,255    (R,G,B)
C    Purple = 180,0,180  (R,G,B)
C
             XICOLB = 75.0/45.0
             XICOLR = 180.0/45.0
             DO 90 I=211,255
              RRED(I)=(I-210)*XICOLR/255.
              RGREEN(I)=0.
              RBLUE(I)=(255.0-(I-210)*XICOLB)/255.
   90        CONTINUE

            ENDIF


            IF (NUMCOL(1,ICURDV) .EQ. 1) THEN
               CALL PLTCOL (8+0, 0.7059, 0.0, 0.7059)

            ELSE
               NEWCOL = NUMCOL(1,ICURDV) - 1
               DO 100 I = 0, NEWCOL
C  Interpolate on color number, not color values.
C  For a spectrum with less than 224 colors, choose
C  equally-spaced color numbers from the range 255-32
                  FRAC = FLOAT(I) / FLOAT(NEWCOL)
                  NC=NINT(255.-FRAC*223.)
                  CALL PLTCOL (8+I, RRED(NC), RGREEN(NC), RBLUE(NC))
  100          CONTINUE
            END IF

          ENDIF

         END IF

C      --Change old alternate colors to black

         IF (LSTALT(ICURDV) .GT. 0) THEN
            DO 110 I = NEWCOL+1, LSTNUM(ICURDV)-1
               CALL PLTCOL (8+I, 0.0, 0.0, 0.0)
  110       CONTINUE
         END IF
         CALL PLTFLU

         LSTALT(ICURDV) = MAPALT(ICURDV)
         LSTNUM(ICURDV) = NUMCOL(LSTALT(ICURDV),ICURDV)
      END IF

      MAPUSE(ICURDV) = MAPALT(ICURDV)

      IRNBWL=IRAINB

      RETURN
      END
