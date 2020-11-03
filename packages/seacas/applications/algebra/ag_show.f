C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE SHOW (TYPE, CMDTBL, NAMECO, BLKTYP, NAMES, TIMES,
     &                 IPTIMS, IDELB, VISELB, SELELB)
C=======================================================================

C   --*** SHOW *** (ALGEBRA) Display ALGEBRA information
C   --   Written by Amy Gilkey - revised 05/18/88
C   --
C   --SHOW prints out the requested information.
C   --
C   --Parameters:
C   --   TYPE      - IN - the type of SHOW requested (none, VARS)
C   --   CMDTBL    - IN - the table of valid SHOW options; if empty, show if
C   --                    not equal to default only
C   --   NAMECO    - IN - the coordinate names
C   --   BLKTYP    - IN - the element block names
C   --   NAMES     - IN - the global, nodal, and element variable names
C   --   TIMES     - IN - the database time steps
C   --   IPTIMS    - IN - the selected times steps
C   --   IDELB     - IN - the element block IDs
C   --   VISELB(i) - IN - true iff element block i is to be written
C   --   SELELB(i) - IN - true iff element block i is selected
C   --
C   --Common Variables:
C   --   Uses IXLHS, NAMVAR, ISTVAR of /VAR../
C   --   Uses NUMALI, NAMALI, NIXALI, IXALI of /ALIAS../
C   --   Uses NDBIN of /DBASE/
C   --   Uses TITLE, TITLEO of /DBTITL/
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK, NVARNP, NVAREL, NVARGL
C   --      of /DBNUMS/
C   --   Uses NPTIMS, TMIN, TMAX, DELT, NINTV, WHONLY of /TIMES/
C   --   Uses ISZOOM, ZMLIM of /ZOOM/

      include 'exodusII.inc'
      include 'ag_namlen.blk'
      include 'ag_numeqn.blk'
      include 'ag_var.blk'
      include 'ag_alias.blk'
      include 'ag_dbase.blk'
      include 'ag_dbtitl.blk'
      include 'ag_dbnums.blk'
      include 'ag_dbnumg.blk'
      include 'ag_dbnumq.blk'
      include 'ag_times.blk'
      include 'ag_zoom.blk'
      include 'ag_filter.blk'
      include 'ag_remove.blk'

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      CHARACTER*(*) TYPE
      CHARACTER*(*) CMDTBL(*)
      CHARACTER*(namlen) NAMECO(*)
      CHARACTER*(MXSTLN) BLKTYP(*)
      CHARACTER*(namlen) NAMES(*)
      REAL TIMES(*)
      INTEGER IPTIMS(*)
      INTEGER IDELB(*)
      LOGICAL VISELB(NELBLK)
      LOGICAL SELELB(NELBLK)

      LOGICAL IFDEF
      CHARACTER*(MXSTLN) SHOTYP
      CHARACTER*(MAXNAM) STRING
      CHARACTER*20 RSTR(6)
      INTEGER INAM(MAXVAR)

      IFDEF = (CMDTBL(1) .NE. ' ')
      IF (CMDTBL(1) .NE. ' ') THEN
         CALL ABRSTR (SHOTYP, TYPE, CMDTBL)
         IF (SHOTYP .EQ. ' ') SHOTYP = TYPE
      ELSE
         SHOTYP = TYPE
      END IF

      IF (SHOTYP .EQ. 'TITLE') THEN
         WRITE (*, 10030) 'Output database title:'
         WRITE (*, 10030) TITLEO(:LENSTR(TITLEO))
      ELSE IF (SHOTYP .EQ. 'SAVE') THEN
         N = 0
         DO 100 IVAR = MAXVAR, IXLHS, -1
            IF (ISTVAR(ICURTM,IVAR) .EQ. -2) THEN
               N = N + 1
               INAM(N) = IVAR
            END IF
  100    CONTINUE
         IF (N .GT. 0) THEN
            WRITE (*, 10030) 'SAVEd variables:'
            WRITE (*, 10000) (NAMVAR(INAM(I))(:namlen), I=1,N)
10000        FORMAT ((12X, 2(A, :, 1X)))
         ELSE IF (IFDEF) THEN
            WRITE (*, 10030) 'There are no SAVEd variables'
         END IF

      ELSE IF (SHOTYP .EQ. 'DELETE') THEN
         N = 0
         DO 110 IVAR = MAXVAR, IXLHS, -1
            IF (ISTVAR(ICURTM,IVAR) .EQ. -1) THEN
               N = N + 1
               INAM(N) = IVAR
            END IF
  110    CONTINUE
         IF (N .GT. 0) THEN
            maxn = 0
            do i=1, n
              k = lenstr(namvar(inam(i)))
              if (k .gt. maxn) maxn = k
            end do
            WRITE (*, 10030) 'Temporary variables:'
            WRITE (*, 10000) (NAMVAR(INAM(I))(:maxn), I=1,N)
         ELSE IF (IFDEF) THEN
            WRITE (*, 10030) 'There are no temporary variables'
         END IF

      ELSE IF (SHOTYP .EQ. 'ALIAS') THEN
         IF (NUMALI .GT. 0) THEN
            WRITE (*, 10030) 'Defined aliases:'
         END IF
         DO 120 IALI = 1, NUMALI
            WRITE (STRING, 10010)
     &         NAMALI(IALI), (NAMES(IXALI(J,IALI)), J=1,NIXALI(IALI))
10010        FORMAT (A, ' = ', 10 (A, ' '))
            CALL SQZSTR (STRING, LSTR)
            WRITE (*, 10030) '   ', STRING(:LSTR)
  120    CONTINUE

      ELSE IF ((SHOTYP .EQ. 'TMIN') .OR. (SHOTYP .EQ. 'TMAX')
     &   .OR. (SHOTYP .EQ. 'DELTIME')
     &   .OR. (SHOTYP .EQ. 'NINTV') .OR. (SHOTYP .EQ. 'ZINTV')
     &   .OR. (SHOTYP .EQ. 'ALLTIMES')
     &   .OR. (SHOTYP .EQ. 'STEPS') .OR. (SHOTYP .EQ. 'TIMES')) THEN

         IF (IFDEF .OR. (NPTIMS .NE. NSTEPS)) THEN
            CALL SHOTSP (TMIN, TMAX, DELT, NINTV, NPTIMS)
         END IF

         IF (((SHOTYP .EQ. 'STEPS') .OR. (SHOTYP .EQ. 'TIMES'))
     &      .AND. (NPTIMS .GT. 0)) THEN
            CALL SHPTIM (NPTIMS, IPTIMS, TIMES)
         END IF
      else if ((shotyp .eq. 'PTIME') .or. (shotyp .eq. 'PTIMES')) then
         call prterr ('CMDREQ', 'Please use SHOW TIMES')

      ELSE IF (SHOTYP .EQ. 'ZOOM') THEN
         IF (ISZOOM) THEN
            CALL NUMSTR (MIN(3,NDIM)*2, 4, ZMLIM, RSTR, LSTR)
            WRITE (STRING, 10020) (RSTR(I)(:LSTR), I=1,MIN(3,NDIM)*2)
10020        FORMAT (6 (A, ' '))
            CALL SQZSTR (STRING, LSTR)
            if (zoomin) then
              rstr(1) = ' Retain Inside'
            else
              rstr(1) = ' Retain Outside'
            end if
            WRITE (*, 10030) 'Zoomed mesh: ', STRING(:LSTR), rstr(1)
         ELSE IF (IFDEF) THEN
            WRITE (*, 10030) 'Zoomed mesh is not defined'
         END IF

      ELSE IF (SHOTYP .EQ. 'FILTER') THEN
         IF (ISFILTER) THEN
           CALL NUMSTR1 (4, VALFLT, RSTR(1), LSTR1)
           CALL NUMSTR1 (4, TIMFLT, RSTR(2), LSTR2)
           if (cmpflt .eq. 1) RSTR(3) = 'lt'
           if (cmpflt .eq. 2) RSTR(3) = 'le'
           if (cmpflt .eq. 3) RSTR(3) = 'eq'
           if (cmpflt .eq. 4) RSTR(3) = 'ne'
           if (cmpflt .eq. 5) RSTR(3) = 'gt'
           if (cmpflt .eq. 6) RSTR(3) = 'ge'

           CALL DBVIX('E', 1, IEV)
           WRITE (*, 10030) 'Remove elements with ',
     *       names(idxflt+iev-1)(:lenstr(names(idxflt+iev-1))),
     *       ' ',RSTR(3)(:2), ' ',RSTR(1)(:LSTR1),' at time ',
     *       RSTR(2)(:LSTR2)
         ELSE IF (IFDEF) THEN
            WRITE (*, 10030) 'Filtering is not defined'
         END IF

      ELSE IF (SHOTYP .EQ. 'REMOVE') THEN
         if (idsglobal) then
           rstr(1) = 'GLOBAL'
         else
           rstr(1) = 'LOCAL'
         end if
         IF (ISREMOVE) THEN
           write (*, 10040) rstr(1)(:lenstr(rstr(1))),
     *       (idsrem(i),i=1,irmcnt)
         ELSE IF (IFDEF) THEN
            WRITE (*, 10030) 'Removal by element id is not defined'
         END IF

      ELSE IF (SHOTYP .EQ. 'VISIBLE') THEN
         NVIS = NUMEQL (.TRUE., NELBLK, VISELB)
         IF (NVIS .LT. NELBLK) THEN
            if (NVIS .GT. 0) then
               WRITE (*, 10030) 'Write elements of element blocks: '
            else
               WRITE(*,*) 'Element blocks are not visible'
            endif
         ELSE IF (IFDEF) THEN
            WRITE (*, 10030) 'Write elements of ALL element blocks: '
         END IF

         if (NVIS .GT. 0) THEN
           N = 0
           STRING = ' '
           DO 130 IELB = 1, NELBLK
             IF (VISELB(IELB)) THEN
               N = N + 1
               WRITE (STRING((N-1)*6+1:N*6), '(I5, 1X)')
     &           IDELB(IELB)
               if ((n+1)*6 .gt. 80) THEN
                 lstr = lenstr(string)
                 write (*,10030) '      ',string(:lstr)
                 N = 0
                 STRING = ' '
               end if
             END IF
 130       CONTINUE
           if (n .gt. 0) then
             lstr = lenstr(string)
             write (*,10030) '      ',string(:lstr)
           end if
         END IF

       ELSE IF ((SHOTYP .EQ. 'BLOCKS')
     &     .OR. (SHOTYP .EQ. 'MATERIAL')) THEN
         NSEL = NUMEQL (.TRUE., NELBLK, SELELB)
         IF (NSEL .LT. NELBLK) THEN
           WRITE (*, 10030) 'Selected element blocks: '
         ELSE IF (IFDEF) THEN
           WRITE (*, 10030) 'ALL element blocks selected: '
         END IF

         if (nsel .gt. 0) then
           N = 0
           STRING = ' '
           DO 140 IELB = 1, NELBLK
             IF (SELELB(IELB)) THEN
               N = N + 1
               WRITE (STRING((N-1)*6+1:N*6), '(I5, 1X)')
     &           IDELB(IELB)
               if ((n+1)*6 .gt. 80) THEN
                 lstr = lenstr(string)
                 write (*,10030) '      ',string(:lstr)
                 N = 0
                 STRING = ' '
               end if
             END IF
 140       CONTINUE
           if (n .gt. 0) then
             lstr = lenstr(string)
             write (*,10030) '      ',string(:lstr)
           end if
         end if

      ELSE IF ((SHOTYP .EQ. 'LOG')
     &   .OR. (SHOTYP .EQ. 'LIST') .OR. (SHOTYP .EQ. 'SHOW')
     &   .OR. (SHOTYP .EQ. 'HELP')
     &   .OR. (SHOTYP .EQ. 'END') .OR. (SHOTYP .EQ. 'EXIT')
     &   .OR. (SHOTYP .EQ. 'QUIT')) THEN
         CALL PRTERR ('CMDERR',
     &      'SHOW is not available for this command')

      ELSE
         CALL SHOCMD ('SHOW Options:', CMDTBL)
      END IF

      RETURN
10030  FORMAT (1X, 10A)
10040  FORMAT (1X, 'Remove elements with ',A,' ids: ',:,/,(10I10))
      END
