C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C 

C=======================================================================
*DECK,TRUTBL
      SUBROUTINE TRUTBL(MP,IMP,IDA,IDB,ITRTA,ITRTB)
C
C     ******************************************************************
C
C     SUBROUTINE TO BUILD THE TRUTH TABLE FOR THE RECIPIENT MESH
C     THIS IMPROVES (SPEEDS-UP) THE PROCESS OF WRITING ELEMENT 
C     VARIABLES
C
C     Called by MAPVAR
C
C     ******************************************************************
C
C  MP         INT   Array of element block mapping from donor
C                   mesh to recipient mesh (1:2,1:IMP)
C  IMP        INT   Number of entries in MP
C  IDA        INT   Donor mesh element block I.D. (1:nblksa)
C  IDB        INT   Recipient mesh element block I.D. (1:nblksa)
C  ITRTA      INT   Donor mesh truth table (1:nvarel,1:nblksa)
C  ITRTB      INT   Recipient mesh truth table (1:nvarel,1:nblksb)
C
C     ******************************************************************
C
      include 'exodusII.inc'
C
      include 'aexds1.blk'
      include 'amesh.blk'
      include 'bmesh.blk'
      include 'ex2tp.blk'
      include 'tapes.blk'
      include 'debg.blk'

C
      DIMENSION MP(3,*),IDA(*),IDB(*),ITRTA(NVAREL,*),ITRTB(NVAREL,*)
C
C     ******************************************************************
C
C ... For later sanity check
      do 2 iblkb = 1, nblksb
        do 1 ivar = 1, nvarel
          itrtb(ivar,iblkb) = -999
 1      continue
 2    continue
      
      CALL EXGVTT (NTP2EX,NBLKSA,NVAREL,ITRTA,IERR)
C
C Match element block I.D.'s from donor and recipient mesh
C
      DO 10 IBLKB = 1, NBLKSB
        IDBLKB = IDB(IBLKB)
        DO 20 IBLKA = 1, NBLKSA
          IDBLKA = IDA(IBLKA)
          DO 30 IM = 1, IMP
            IF (IDBLKB .EQ. MP(2,IM) .AND. IDBLKA .EQ. MP(1,IM))THEN
              if (idebug .ge. 3) then
                write (nout,9002) iblkb, idblkb, iblka, idblka, im
                write (ntpout,9002) iblkb, idblkb, iblka, idblka, im
 9002           FORMAT('B:', i5, i5, ' A:', i5, i5, ' IM:', i5)
              end if
              DO 40 IVAR = 1, NVAREL
                ITRTB(IVAR,IBLKB) = ITRTA(IVAR,IBLKA)
   40         CONTINUE
              GO TO 20
            END IF
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C
C Do some error checking - write a warning
C
      DO 50 IBLKB = 1, NBLKSB
        IFND = 0
        IDBLKB = IDB(IBLKB)
        DO 60 IM = 1, IMP
          IF (IDBLKB .EQ. MP(2,IM))THEN
            IFND = 1
          END IF
   60   CONTINUE
        IF (IFND .EQ. 0)THEN
          DO 70 IVAR = 1, NVAREL
            ITRTB(IVAR,IBLKB) = 0
   70     CONTINUE
          WRITE(NOUT,1000)IBLKB,IDBLKB
          WRITE(NTPOUT,1000)IBLKB,IDBLKB
        END IF
   50 CONTINUE
C
C ... Sanity check
      do 90 iblkb = 1, nblksb
        do 80 ivar = 1, nvarel
          if (itrtb(ivar,iblkb) .eq. -999) then
            idblkb = idb(iblkb)
            write(nout, 1010) iblkb, idblkb, ivar
            write(ntpout, 1010) iblkb, idblkb, ivar
          end if
 80     continue
 90   continue

      if (idebug .ge. 3) then
        do 180 iblkb = 1, nblksb
          write (nout,1020) iblkb, idb(iblkb),
     *      (itrtb(ivar,iblkb),ivar=1,nvarel)
 180    continue
        do 190 iblka = 1, nblksa
          write (nout,1020) iblka, ida(iblka),
     *      (itrta(ivar,iblka),ivar=1,nvarel)
 190    continue
      end if
      
      CALL EXPVTT(NTP4EX,NBLKSB,NVAREL,ITRTB,IERR)
C
 1000 FORMAT(5X,'MAPPING BACK INTO DONOR MESH FOR',/,
     1'          RECIPIENT MESH ELEMENT BLOCK NUMBER',I7,/,
     2'          ELEMENT BLOCK I. D.',I7,/,
     3'          WAS NOT FOUND. THIS ELEMENT BLOCK WILL',/,
     4'          NOT BE MAPPED.')
 1010 FORMAT(5X,'Block Number ', i7, /,
     1  'Element Block ID ', i7,/,
     2  'Has unset truth table entry for variable ',i7)
 1020 format('BLK = ', i3, ' Id = ', i4,:,/,(10i2))
      RETURN
      END
