AC_INIT(README)
AM_INIT_AUTOMAKE(siMPI,0.1a9)
AM_CONFIG_HEADER(mpi_config.h)
 
AC_PROG_CC
AC_PROG_RANLIB


AC_HEADER_STDC
AC_CHECK_FUNCS(memcpy)
AC_CHECK_FUNCS(sysconf)

dnl if test "_SC_CLK_TCK" =""; then
dnl  AC_MSG_WARN([_SC_CLK_TCK is being defined to be 60])
dnl  AC_DEFINE(_SC_CLK_TCK, 60, [Clock tick defined to be 60])
dnl fi

AC_MSG_CHECKING(for --with-debug)
AC_ARG_WITH(debug,[  --with-debug          Build with debugging symbols],[
        if test "$withval" = "yes"; then
           withval="-g"
        fi
        CFLAGS="$withval $CFLAGS"
        AC_MSG_RESULT($withval)
],[
        AC_MSG_RESULT(no)
])

AC_MSG_CHECKING(for --with-preallocation)
AC_ARG_WITH(preallocation,[  --with-preallocation      Build with user selected preallocation buffer size],[
        if test "$withval" = "yes"; then
           withval="10"
        fi
	CFLAGS="-DPREALLOCATION_SIZE=$withval $CFLAGS"
        AC_MSG_RESULT($withval)
],[
        AC_MSG_RESULT(no)
])


dnl -----------------------------------------------
dnl FORTRAN specific items
dnl -----------------------------------------------
AC_PROG_F77

if test "$F77" = ""; then
   AC_MSG_WARN([FORTRAN not located])
   F77="/fortran/not/installed"
   MPI_FORTRAN_MANGLING_PATTERNS=""
   MPI_FORTRAN_INTEGER="int"
   MPI_FORTRAN_DOUBLE="double"
else
    dnl -----------------------------------------------
    dnl Figure out which name mangling scheme is used
    dnl -----------------------------------------------
    AC_MSG_CHECKING(Fortran77 mangling)
      cat >conftest.f <<_ACEOF
           PROGRAM STUFF
           CALL CONFTEST()
           CALL CONF_TEST()
           STOP
           END
_ACEOF
      cat >conftest.c <<_ACEOF
#include <stdio.h>
    extern void _CONFTEST();
    void _CONFTEST() {
      fputs("UPPER_%s|",stdout);
    }
    extern void CONFTEST();
    void CONFTEST() {
      fputs("UPPER%s|",stdout);
    }
    extern void CONFTEST_();
    void CONFTEST_() {
      fputs("UPPER%s_|",stdout);
    }
    void _conftest();
    void _conftest() {
      fputs("lower_%s|",stdout);
    }
    void conftest_();
    void conftest_() {
      fputs("lower%s_|",stdout);
    }

    extern void _CONF_TEST();
    void _CONF_TEST() {
      fputs("UPPER_%s",stdout);
    }
    extern void CONF_TEST();
    void CONF_TEST() {
      fputs("UPPER%s",stdout);
    }
    extern void CONF_TEST_();
    void CONF_TEST_() {
      fputs("UPPER%s_",stdout);
    }
    void _conf_test();
    void _conf_test() {
      fputs("lower_%s",stdout);
    }
    void conf_test_();
    void conf_test_() {
      fputs("lower%s_",stdout);
    }
    void conf_test__();
    void conf_test__() {
      fputs("lower%s__",stdout);
    }
_ACEOF
    rm -f conftest.$ac_objext
    (eval $ac_compile) 2>&5
    ac_status=$?
    if test "$ac_status" != "0"; then
      AC_MSG_RESULT(compile error)
      AC_MSG_ERROR([f77 failure])
    fi

    (eval $F77 conftest.$ac_objext conftest.f -o conftest$ac_exeext) 2>&5
    ac_status=$?
    if test "$ac_status" != "0"; then
      AC_MSG_RESULT(linkage error)
      AC_MSG_ERROR([f77 failure])
    fi

    MPI_FORTRAN_MANGLING_PATTERNS=`./conftest$ac_exeext`
    ac_status=$?
    if test "$ac_status" != "0"; then
      AC_MSG_RESULT(runtime error)
      AC_MSG_ERROR([f77 failure])
    fi
    AC_MSG_RESULT($MPI_FORTRAN_MANGLING_PATTERNS)

    dnl -----------------------------------------------
    dnl Detect FORTRAN integer size
    dnl -----------------------------------------------
    AC_MSG_CHECKING(C equivalent of Fortran77 INTEGER )
      cat >conftest.f <<_ACEOF
      PROGRAM FOO
      INTEGER IARRAY(3)
      IARRAY(1) = 11
      IARRAY(2) = 22
      IARRAY(3) = 33
      CALL CINTSIZE(IARRAY)
      STOP
      END
_ACEOF
	
[      cat >conftest.c <<_ACEOF
#include <stdio.h>
#include <stdlib.h>
static void work(void* iarray) {
  short* asShort = (short*)iarray;
  int* asInt = (int*)iarray;
  long* asLong = (long*)iarray;

  if ( asShort[0] == 11 && asShort[1] == 22 && asShort[2] == 33 ) {
    puts("short");
    exit(0);
  }
  
  if ( asInt[0] == 11 && asInt[1] == 22 && asInt[2] == 33 ) {
    puts("int");
    exit(0);
  }
  
  if ( asLong[0] == 11 && asLong[1] == 22 && asLong[2] == 33 ) {
    puts("long");
    exit(0);
  }

  puts("unknown");
  exit(1);
}

void cintsize_(void* iarray) { work(iarray); }
void _cintsize(void* iarray) { work(iarray); }
void CINTSIZE_(void* iarray) { work(iarray); }
void _CINTSIZE(void* iarray) { work(iarray); }
void CINTSIZE(void* iarray)  { work(iarray); }
_ACEOF
]
    rm -f conftest.$ac_objext
    (eval $ac_compile) 2>&5
    ac_status=$?
    if test "$ac_status" != "0"; then
      AC_MSG_RESULT(compile error)
      AC_MSG_ERROR([Size detection failure])
    fi

    (eval $F77 conftest.$ac_objext conftest.f -o conftest$ac_exeext) 2>&5
    ac_status=$?
    if test "$ac_status" != "0"; then
      AC_MSG_RESULT(linkage error)
      AC_MSG_ERROR([Size detection failure])
    fi

    MPI_FORTRAN_INTEGER=`./conftest$ac_exeext`
    ac_status=$?
    if test "$ac_status" != "0"; then
      AC_MSG_RESULT(runtime error)
      AC_MSG_ERROR([Size detection failure])
    fi
    AC_MSG_RESULT($MPI_FORTRAN_INTEGER)

    dnl -----------------------------------------------
    dnl Detect FORTRAN double size
    dnl -----------------------------------------------
    AC_MSG_CHECKING(C equivalent of Fortran77 DOUBLE )
      cat >conftest.f <<_ACEOF
      PROGRAM FOO
      DOUBLE PRECISION DARRAY(3)
      DARRAY(1) = 11.0
      DARRAY(2) = 22.0
      DARRAY(3) = 33.0
      CALL CDOUBLESIZE(DARRAY)
      STOP
      END
_ACEOF
	
[      cat >conftest.c <<_ACEOF
#include <stdio.h>
#include <stdlib.h>
static void work(void* darray) {
  float* asFloat = (float*)darray;
  double* asDouble = (double*)darray;

  if ( asFloat[0] == 11 && asFloat[1] == 22 && asFloat[2] == 33 ) {
    puts("float");
    exit(0);
  }
  
  if ( asDouble[0] == 11 && asDouble[1] == 22 && asDouble[2] == 33 ) {
    puts("double");
    exit(0);
  }

  puts("unknown");
  exit(1);
}

void cdoublesize_(void* darray) { work(darray); }
void _cdoublesize(void* darray) { work(darray); }
void CDOUBLESIZE_(void* darray) { work(darray); }
void _CDOUBLESIZE(void* darray) { work(darray); }
void CDOUBLESIZE(void* darray)  { work(darray); }
_ACEOF
]
    rm -f conftest.$ac_objext
    (eval $ac_compile) 2>&5
    ac_status=$?
    if test "$ac_status" != "0"; then
      AC_MSG_RESULT(compile error)
      AC_MSG_ERROR([Size detection failure])
    fi

    (eval $F77 conftest.$ac_objext conftest.f -o conftest$ac_exeext) 2>&5
    ac_status=$?
    if test "$ac_status" != "0"; then
      AC_MSG_RESULT(linkage error)
      AC_MSG_ERROR([Size detection failure])
    fi

    MPI_FORTRAN_DOUBLE=`./conftest$ac_exeext`
    ac_status=$?
    if test "$ac_status" != "0"; then
      AC_MSG_RESULT(runtime error)
      AC_MSG_ERROR([Size detection failure])
    fi
    AC_MSG_RESULT($MPI_FORTRAN_DOUBLE)

fi
AC_SUBST(MPI_FORTRAN_MANGLING_PATTERNS)
AC_DEFINE_UNQUOTED(MPI_FORTRAN_MANGLING_PATTERNS,"$MPI_FORTRAN_MANGLING_PATTERNS",[Mangling template])
AC_DEFINE_UNQUOTED(MPI_FORTRAN_INTEGER,$MPI_FORTRAN_INTEGER,[C equivalent of FORTRAN INTEGER])
AC_DEFINE_UNQUOTED(MPI_FORTRAN_DOUBLE,$MPI_FORTRAN_DOUBLE,[C equivalent of FORTRAN DOUBLE PRECISION])

#AC_OUTPUT(Makefile)
#AC_OUTPUT(Makefile mpi_path.h)
