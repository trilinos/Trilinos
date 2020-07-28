C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

c sdcgif - FORTRAN shell for handling strings

C  CESC - Escape
      SUBROUTINE CESC (FUNCID, LDR, DATA)
      INTEGER FUNCID
      INTEGER LDR
      CHARACTER*(*) DATA

      CALL CESC1( FUNCID, LDR, DATA, LEN(DATA) )
      RETURN
      END

C  CTX - Text
      SUBROUTINE CTX(X, Y, FLAG, TEXT )
      REAL X,Y
      INTEGER FLAG
      CHARACTER*(*) TEXT

      CALL CTX1( X, Y, FLAG, TEXT, LEN(TEXT) )
      RETURN
      END

C  CGTXX - Get Text Extent
      SUBROUTINE CGTXX( X, Y, STRING, VSTAT, VCONC, XCONC, YCONC,
     1                  X1, Y1, X2, Y2, X3, Y3, X4, Y4)
      REAL X,Y
      CHARACTER*(*) STRING
      INTEGER VSTAT, VCONC
      REAL XCONC, YCONC
      REAL X1, Y1, X2, Y2, X3, Y3, X4, Y4

      CALL CGTXX1( X, Y, STRING, VSTAT, VCONC, XCONC, YCONC,
     1             X1, Y1, X2, Y2, X3, Y3, X4, Y4, LEN(STRING) )
      RETURN
      END

C  CQCHH - Inquire List of Available Character Heights
      SUBROUTINE CQCHH( FONT, TXP, NREQ, FIRST, VSTAT, NTOTAL,
     1                  NLIST, CHHIT )
      CHARACTER*(*)FONT
      INTEGER TXP, NREQ, FIRST, VSTAT, NTOTAL, NLIST
      INTEGER CHHIT(*)

      CALL CQCHH1( FONT, TXP, NREQ, FIRST, VSTAT, NTOTAL,
     1             NLIST, CHHIT, LEN(FONT) )
      RETURN
      END
