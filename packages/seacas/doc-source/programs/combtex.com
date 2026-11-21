$ !!! VERIFY = F$VERIFY("NO")
$!
$! COMBTEX creates a combined .TEX file with all the \include and \input files
$! in the file
$!    P1 is the .TEX file name (no extension)
$!    P8 is the directory
$!
$ IF P8 .NES. "" THEN SET DEFAULT 'P8'
$!
$ BASEFILE = P1
$ BASEFILE = BASEFILE - ".TEX"
$!
$ IF F$SEARCH("''BASEFILE'.TEX") .NES. "" THEN GOTO ENDIF_FILE_OK
$    WRITE SYS$OUTPUT "''BASEFILE'.TEX does not exist"
$    GOTO END
$ ENDIF_FILE_OK:
$!
$ ASSIGN/USER 'BASEFILE'.TEX INFILE
$ ASSIGN/USER 'BASEFILE'.ALLTEX OUTFILE
$ RUN texexe:COMBTEX
$ PURGE 'BASEFILE'.ALLTEX
$!
$ END:
$ !!! IF VERIFY THEN SET VERIFY
