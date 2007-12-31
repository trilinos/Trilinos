dnl @synopsis TAC_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/).
dnl On success, it sets the BLAS_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl 	$BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by TAC_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl Many libraries are searched for, from ATLAS to CXML to ESSL.
dnl The user may also use --with-blas=<lib> in order to use some
dnl specific BLAS library <lib>. 
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a BLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_BLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @version $Id$
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
dnl Edited by Jim Willenbring on 5-14-2004 to check for dgemm instead of
dnl sgemm.
dnl Edited by Jim Willenbring on 4-17-2006 to stop looking for BLAS if
dnl a specific BLAS library specified by a user cannot be used.
dnl Edited by Jim Willenbring on 12-5-2007 and renamed from ACX_BLAS
dnl to TAC_BLAS.  Modified the macro to stop using the Fortran 
dnl name-mangling provided by AC_F77_FUNC.  Instead, we now try
dnl lowercase, lowercase_, UPPERCASE, and UPPERCASE_ and then
dnl define a macro F77_BLAS_MANGLE that provides the name-mangling
dnl scheme used by the BLAS library.  This change was made because
dnl the BLAS library on a system does not need to use the same
dnl name-mangling as the compiler and we have an option to build
dnl Trilinos without a Fortran compiler.


AC_DEFUN([TAC_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
tac_blas_ok=no

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) tac_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

tac_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $tac_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"

	AC_MSG_CHECKING([for dgemm_ in $BLAS_LIBS])
	AC_TRY_LINK_FUNC(dgemm_, [tac_blas_ok=yes])
	AC_MSG_RESULT($tac_blas_ok)

        if test $tac_blas_ok = yes; then
		AC_LANG_PUSH(Fortran 77)
		AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name ## _], [Define the Fortran name mangling to be used for the BLAS])
		AC_LANG_POP(Fortran 77)
		ac_cv_use_blas_mangling=lowercaseUnderscore
        else
                AC_MSG_CHECKING([for dgemm in $BLAS_LIBS])
                AC_TRY_LINK_FUNC(dgemm, [tac_blas_ok=yes])
                AC_MSG_RESULT($tac_blas_ok)
        
        	if test $tac_blas_ok = yes; then
                	AC_LANG_PUSH(Fortran 77)
	                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name], [Define the Fortran name mangling to be used for the BLAS])
        	        AC_LANG_POP(Fortran 77)
			ac_cv_use_blas_mangling=lowercaseNoUnderscore
	        else
        	        AC_MSG_CHECKING([for DGEMM in $BLAS_LIBS])
                	AC_TRY_LINK_FUNC(DGEMM, [tac_blas_ok=yes])
	                AC_MSG_RESULT($tac_blas_ok)
        	        if test $tac_blas_ok = yes; then
	        	        AC_LANG_PUSH(Fortran 77)
	        	        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME], [Define the Fortran name mangling to be used for the BLAS])
        	        	AC_LANG_POP(Fortran 77)
				ac_cv_use_blas_mangling=uppercaseNoUnderscore
		        else

		                AC_MSG_CHECKING([for DGEMM_ in $BLAS_LIBS])
		                AC_TRY_LINK_FUNC(DGEMM_, [tac_blas_ok=yes], [user_spec_blas_failed=yes])
		                AC_MSG_RESULT($tac_blas_ok)
			        if test $tac_blas_ok = yes; then
			                AC_LANG_PUSH(Fortran 77)
	                		AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME ## _], [Define the Fortran name mangling to be used for the BLAS])
			                AC_LANG_POP(Fortran 77)
					ac_cv_use_blas_mangling=uppercaseUnderscore
                                fi
			fi
		fi
	fi

        LIBS="$save_LIBS"

fi
fi

# If the user specified a blas library that could not be used we will
# halt the search process rather than risk finding a blas library that
# the user did not specify.

if test "x$user_spec_blas_failed" != xyes; then

# BLAS linked to by default?  (happens on some supercomputers)
if test $tac_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_FUNC(dgemm_, [tac_blas_ok=yes])
	if test $tac_blas_ok = yes; then
                AC_LANG_PUSH(Fortran 77)
                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name ## _], [Define the Fortran name mangling to be used for the BLAS])
                AC_LANG_POP(Fortran 77)
		ac_cv_use_blas_mangling=lowercaseUnderscore
        else
		AC_CHECK_FUNC(dgemm, [tac_blas_ok=yes])
		if test $tac_blas_ok = yes; then
                        AC_LANG_PUSH(Fortran 77)
                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name], [Define the Fortran name mangling to be used for the BLAS])
                        AC_LANG_POP(Fortran 77)
			ac_cv_use_blas_mangling=lowercaseNoUnderscore
                else
			AC_CHECK_FUNC(DGEMM, [tac_blas_ok=yes])
			if test $tac_blas_ok = yes; then
                                AC_LANG_PUSH(Fortran 77)
                                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME], [Define the Fortran name mangling to be used for the BLAS])
                                AC_LANG_POP(Fortran 77)
				ac_cv_use_blas_mangling=uppercaseNoUnderscore
                        else
				AC_CHECK_FUNC(DGEMM_, [tac_blas_ok=yes])
				if test $tac_blas_ok = yes; then
                                        AC_LANG_PUSH(Fortran 77)
                                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME ## _], [Define the Fortran name mangling to be used for the BLAS])
                                        AC_LANG_POP(Fortran 77)
					ac_cv_use_blas_mangling=uppercaseUnderscore
                                fi
                        fi
                fi
        fi

	LIBS="$save_LIBS"
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $tac_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, dgemm_,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[tac_blas_ok=yes
			BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
	if test $tac_blas_ok = yes; then
                AC_LANG_PUSH(Fortran 77)
                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name ## _], [Define the Fortran name mangling to be used for the BLAS])
                AC_LANG_POP(Fortran 77)
		ac_cv_use_blas_mangling=lowercaseUnderscore
        else
		AC_CHECK_LIB(atlas, ATL_xerbla,
	                [AC_CHECK_LIB(f77blas, dgemm,
        	        [AC_CHECK_LIB(cblas, cblas_dgemm,
                	        [tac_blas_ok=yes
				BLAS_LIBS="-lcblas -lf77blas -latlas"],
                       		[], [-lf77blas -latlas])],
                        	[], [-latlas])])
		if test $tac_blas_ok = yes; then
                        AC_LANG_PUSH(Fortran 77)
                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name], [Define the Fortran name mangling to be used for the BLAS])
                        AC_LANG_POP(Fortran 77)
			ac_cv_use_blas_mangling=lowercaseNoUnderscore
                else
			AC_CHECK_LIB(atlas, ATL_xerbla,
	                        [AC_CHECK_LIB(f77blas, DGEMM,
	                        [AC_CHECK_LIB(cblas, cblas_dgemm,
        	                        [tac_blas_ok=yes
		                        BLAS_LIBS="-lcblas -lf77blas -latlas"],
                	                [], [-lf77blas -latlas])],
                        	        [], [-latlas])])
			if test $tac_blas_ok = yes; then
                                AC_LANG_PUSH(Fortran 77)
                                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME], [Define the Fortran name mangling to be used for the BLAS])
                                AC_LANG_POP(Fortran 77)
				ac_cv_use_blas_mangling=uppercaseNoUnderscore
                        else
				AC_CHECK_LIB(atlas, ATL_xerbla,
	                                [AC_CHECK_LIB(f77blas, DGEMM_,
        	                        [AC_CHECK_LIB(cblas, cblas_dgemm,
                	                        [tac_blas_ok=yes
			                        BLAS_LIBS="-lcblas -lf77blas -latlas"],
                        	                [], [-lf77blas -latlas])],
                                	        [], [-latlas])])
				if test $tac_blas_ok = yes; then
                                        AC_LANG_PUSH(Fortran 77)
                                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME ## _], [Define the Fortran name mangling to be used for the BLAS])
                                        AC_LANG_POP(Fortran 77)
					ac_cv_use_blas_mangling=uppercaseUnderscore
                                fi
                        fi
                fi
        fi

fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $tac_blas_ok = no; then
	AC_CHECK_LIB(blas, dgemm_,
		[AC_CHECK_LIB(dgemm, dgemm_,
		[AC_CHECK_LIB(sgemm, sgemm_,
			[tac_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
	if test $tac_blas_ok = yes; then
                AC_LANG_PUSH(Fortran 77)
                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name ## _], [Define the Fortran name mangling to be used for the BLAS])
                AC_LANG_POP(Fortran 77)
		ac_cv_use_blas_mangling=lowercaseUnderscore
        else
		AC_CHECK_LIB(blas, dgemm,
	                [AC_CHECK_LIB(dgemm, dgemm,
        	        [AC_CHECK_LIB(sgemm, sgemm,
                	        [tac_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
                        	[], [-lblas])],
	                        [], [-lblas])])
		if test $tac_blas_ok = yes; then
                        AC_LANG_PUSH(Fortran 77)
                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name], [Define the Fortran name mangling to be used for the BLAS])
                        AC_LANG_POP(Fortran 77)
			ac_cv_use_blas_mangling=lowercaseNoUnderscore
                else
			AC_CHECK_LIB(blas, DGEMM,
	                        [AC_CHECK_LIB(dgemm, DGEMM,
        	                [AC_CHECK_LIB(sgemm, SGEMM,
                	                [tac_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
                        	        [], [-lblas])],
                                	[], [-lblas])])
			if test $tac_blas_ok = yes; then
                                AC_LANG_PUSH(Fortran 77)
                                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME], [Define the Fortran name mangling to be used for the BLAS])
                                AC_LANG_POP(Fortran 77)
				ac_cv_use_blas_mangling=uppercaseNoUnderscore
                        else
				AC_CHECK_LIB(blas, DGEMM_,
	                                [AC_CHECK_LIB(dgemm, DGEMM_,
        	                        [AC_CHECK_LIB(sgemm, SGEMM_,
                	                        [tac_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
                        	                [], [-lblas])],
                                	        [], [-lblas])])
				if test $tac_blas_ok = yes; then
                                        AC_LANG_PUSH(Fortran 77)
                                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME ## _], [Define the Fortran name mangling to be used for the BLAS])
                                        AC_LANG_POP(Fortran 77)
					ac_cv_use_blas_mangling=uppercaseUnderscore
                                fi
                        fi
                fi
        fi

fi

# BLAS in Alpha CXML library?
if test $tac_blas_ok = no; then
	AC_CHECK_LIB(cxml, dgemm_, [tac_blas_ok=yes;BLAS_LIBS="-lcxml"])
	if test $tac_blas_ok = yes; then
                AC_LANG_PUSH(Fortran 77)
                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name ## _], [Define the Fortran name mangling to be used for the BLAS])
                AC_LANG_POP(Fortran 77)
		ac_cv_use_blas_mangling=lowercaseUnderscore
        else
		AC_CHECK_LIB(cxml, dgemm, [tac_blas_ok=yes;BLAS_LIBS="-lcxml"])
		if test $tac_blas_ok = yes; then
                        AC_LANG_PUSH(Fortran 77)
                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name], [Define the Fortran name mangling to be used for the BLAS])
                        AC_LANG_POP(Fortran 77)
			ac_cv_use_blas_mangling=lowercaseNoUnderscore
                else
			AC_CHECK_LIB(cxml, DGEMM, [tac_blas_ok=yes;BLAS_LIBS="-lcxml"])
			if test $tac_blas_ok = yes; then
                                AC_LANG_PUSH(Fortran 77)
                                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME], [Define the Fortran name mangling to be used for the BLAS])
                                AC_LANG_POP(Fortran 77)
				ac_cv_use_blas_mangling=uppercaseNoUnderscore
                        else
				AC_CHECK_LIB(cxml, DGEMM_, [tac_blas_ok=yes;BLAS_LIBS="-lcxml"])
				if test $tac_blas_ok = yes; then
                                        AC_LANG_PUSH(Fortran 77)
                                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME ## _], [Define the Fortran name mangling to be used for the BLAS])
                                        AC_LANG_POP(Fortran 77)
					ac_cv_use_blas_mangling=uppercaseUnderscore
                                fi
                        fi
                fi
        fi

fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $tac_blas_ok = no; then
        AC_CHECK_LIB(dxml, dgemm_, [tac_blas_ok=yes;BLAS_LIBS="-ldxml"])
        if test $tac_blas_ok = yes; then
                AC_LANG_PUSH(Fortran 77)
                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name ## _], [Define the Fortran name mangling to be used for the BLAS])
                AC_LANG_POP(Fortran 77)
		ac_cv_use_blas_mangling=lowercaseUnderscore
        else
                AC_CHECK_LIB(dxml, dgemm, [tac_blas_ok=yes;BLAS_LIBS="-ldxml"])
                if test $tac_blas_ok = yes; then
                        AC_LANG_PUSH(Fortran 77)
                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name], [Define the Fortran name mangling to be used for the BLAS])
                        AC_LANG_POP(Fortran 77)
			ac_cv_use_blas_mangling=lowercaseNoUnderscore
                else
                        AC_CHECK_LIB(dxml, DGEMM, [tac_blas_ok=yes;BLAS_LIBS="-ldxml"])
                        if test $tac_blas_ok = yes; then
                                AC_LANG_PUSH(Fortran 77)
                                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME], [Define the Fortran name mangling to be used for the BLAS])
                                AC_LANG_POP(Fortran 77)
				ac_cv_use_blas_mangling=uppercaseNoUnderscore
                        else
                                AC_CHECK_LIB(dxml, DGEMM_, [tac_blas_ok=yes;BLAS_LIBS="-ldxml"])
                                if test $tac_blas_ok = yes; then
                                        AC_LANG_PUSH(Fortran 77)
                                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME ## _], [Define the Fortran name mangling to be used for the BLAS])
                                        AC_LANG_POP(Fortran 77)
					ac_cv_use_blas_mangling=uppercaseUnderscore
                                fi
                        fi
                fi
        fi

fi

# BLAS in Sun Performance library?
if test $tac_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, dgemm_,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 tac_blas_ok=yes],[],[-lsunmath])])
		if test $tac_blas_ok = yes; then
	                AC_LANG_PUSH(Fortran 77)
        	        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name ## _], [Define the Fortran name mangling to be used for the BLAS])
                	AC_LANG_POP(Fortran 77)
			ac_cv_use_blas_mangling=lowercaseUnderscore
	        else
			AC_CHECK_LIB(sunmath, acosp,
	                        [AC_CHECK_LIB(sunperf, dgemm,
        	                        [BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                	                 tac_blas_ok=yes],[],[-lsunmath])])
			if test $tac_blas_ok = yes; then
	                        AC_LANG_PUSH(Fortran 77)
        	                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name], [Define the Fortran name mangling to be used for the BLAS])
                	        AC_LANG_POP(Fortran 77)
				ac_cv_use_blas_mangling=lowercaseNoUnderscore
                	else
				AC_CHECK_LIB(sunmath, acosp,
	                                [AC_CHECK_LIB(sunperf, DGEMM,
        	                                [BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                	                         tac_blas_ok=yes],[],[-lsunmath])])
				if test $tac_blas_ok = yes; then
	                                AC_LANG_PUSH(Fortran 77)
        	                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME], [Define the Fortran name mangling to be used for the BLAS])
                	                AC_LANG_POP(Fortran 77)
					ac_cv_use_blas_mangling=uppercaseNoUnderscore
                        	else
					AC_CHECK_LIB(sunmath, acosp,
	                                        [AC_CHECK_LIB(sunperf, DGEMM_,
        	                                        [BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                	                                 tac_blas_ok=yes],[],[-lsunmath])])
					if test $tac_blas_ok = yes; then
	                                        AC_LANG_PUSH(Fortran 77)
        	                                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME ## _], [Define the Fortran name mangling to be used for the BLAS])
                	                        AC_LANG_POP(Fortran 77)
						ac_cv_use_blas_mangling=uppercaseUnderscore
                        	        fi
	                        fi
        	        fi
	        fi
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $tac_blas_ok = no; then
	AC_CHECK_LIB(scs, dgemm_, [tac_blas_ok=yes; BLAS_LIBS="-lscs"])
	if test $tac_blas_ok = yes; then
                AC_LANG_PUSH(Fortran 77)
                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name ## _], [Define the Fortran name mangling to be used for the BLAS])
                AC_LANG_POP(Fortran 77)
		ac_cv_use_blas_mangling=lowercaseUnderscore
        else
                AC_CHECK_LIB(scs, dgemm, [tac_blas_ok=yes; BLAS_LIBS="-lscs"])
                if test $tac_blas_ok = yes; then
                        AC_LANG_PUSH(Fortran 77)
                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name], [Define the Fortran name mangling to be used for the BLAS])i
                        AC_LANG_POP(Fortran 77)
			ac_cv_use_blas_mangling=lowercaseNoUnderscore
                else
                        AC_CHECK_LIB(scs, DGEMM, [tac_blas_ok=yes; BLAS_LIBS="-lscs"])
			if test $tac_blas_ok = yes; then
                                AC_LANG_PUSH(Fortran 77)
                                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME], [Define the Fortran name mangling to be used for the BLAS])
                                AC_LANG_POP(Fortran 77)
				ac_cv_use_blas_mangling=uppercaseNoUnderscore
                        else
				AC_CHECK_LIB(scs, DGEMM_, [tac_blas_ok=yes; BLAS_LIBS="-lscs"])
                                if test $tac_blas_ok = yes; then
                                        AC_LANG_PUSH(Fortran 77)
                                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME ## _], [Define the Fortran name mangling to be used for the BLAS])
                                        AC_LANG_POP(Fortran 77)
					ac_cv_use_blas_mangling=uppercaseUnderscore
                                fi
                        fi
                 fi
        fi
fi

# BLAS in SGIMATH library?
if test $tac_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, dgemm_,
		     [tac_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
        if test $tac_blas_ok = yes; then
                AC_LANG_PUSH(Fortran 77)
                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name ## _], [Define the Fortran name mangling to be used for the BLAS])
                AC_LANG_POP(Fortran 77)
		ac_cv_use_blas_mangling=lowercaseUnderscore
        else
		AC_CHECK_LIB(complib.sgimath, dgemm,
                     [tac_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
                if test $tac_blas_ok = yes; then
                        AC_LANG_PUSH(Fortran 77)
                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name], [Define the Fortran name mangling to be used for the BLAS])
                        AC_LANG_POP(Fortran 77)
			ac_cv_use_blas_mangling=lowercaseNoUnderscore
                else
			AC_CHECK_LIB(complib.sgimath, DGEMM,
			     [tac_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
                        if test $tac_blas_ok = yes; then
                                AC_LANG_PUSH(Fortran 77)
                                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME], [Define the Fortran name mangling to be used for the BLAS])
                                AC_LANG_POP(Fortran 77)
				ac_cv_use_blas_mangling=uppercaseNoUnderscore
                        else
				AC_CHECK_LIB(complib.sgimath, DGEMM_,
		                     [tac_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
                                if test $tac_blas_ok = yes; then
                                        AC_LANG_PUSH(Fortran 77)
                                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME ## _], [Define the Fortran name mangling to be used for the BLAS])
                                        AC_LANG_POP(Fortran 77)
					ac_cv_use_blas_mangling=uppercaseUnderscore
                                fi
                        fi
                 fi
        fi
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $tac_blas_ok = no; then
	AC_CHECK_LIB(blas, dgemm_,
		[AC_CHECK_LIB(essl, dgemm_,
			[tac_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
        if test $tac_blas_ok = yes; then
                AC_LANG_PUSH(Fortran 77)
                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name ## _], [Define the Fortran name mangling to be used for the BLAS])
                AC_LANG_POP(Fortran 77)
		ac_cv_use_blas_mangling=lowercaseUnderscore
        else
		AC_CHECK_LIB(blas, dgemm,
	                [AC_CHECK_LIB(essl, dgemm,
        	                [tac_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
                	        [], [-lblas $FLIBS])])
                if test $tac_blas_ok = yes; then
                        AC_LANG_PUSH(Fortran 77)
                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name], [Define the Fortran name mangling to be used for the BLAS])
                        AC_LANG_POP(Fortran 77)
			ac_cv_use_blas_mangling=lowercaseNoUnderscore
                else
			AC_CHECK_LIB(blas, DGEMM,
                	        [AC_CHECK_LIB(essl, DGEMM,
        	                        [tac_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
	                                [], [-lblas $FLIBS])])
			if test $tac_blas_ok = yes; then
                                AC_LANG_PUSH(Fortran 77)
                                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME], [Define the Fortran name mangling to be used for the BLAS])
                                AC_LANG_POP(Fortran 77)
				ac_cv_use_blas_mangling=uppercaseNoUnderscore
                        else
				AC_CHECK_LIB(blas, DGEMM_,
                	                [AC_CHECK_LIB(essl, DGEMM_,
        	                                [tac_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
	                                        [], [-lblas $FLIBS])])
				if test $tac_blas_ok = yes; then
                                        AC_LANG_PUSH(Fortran 77)
                                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME ## _], [Define the Fortran name mangling to be used for the BLAS])
                                        AC_LANG_POP(Fortran 77)
					ac_cv_use_blas_mangling=uppercaseUnderscore
                                fi
                        fi
                 fi
        fi
fi

# Generic BLAS library?
if test $tac_blas_ok = no; then
	AC_CHECK_LIB(blas, dgemm_, [tac_blas_ok=yes; BLAS_LIBS="-lblas"])
        if test $tac_blas_ok = yes; then
                AC_LANG_PUSH(Fortran 77)
                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name ## _], [Define the Fortran name mangling to be used for the BLAS])
                AC_LANG_POP(Fortran 77)
		ac_cv_use_blas_mangling=lowercaseUnderscore
        else
		AC_CHECK_LIB(blas, dgemm, [tac_blas_ok=yes; BLAS_LIBS="-lblas"])
                if test $tac_blas_ok = yes; then
                        AC_LANG_PUSH(Fortran 77)
                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [name], [Define the Fortran name mangling to be used for the BLAS])
                        AC_LANG_POP(Fortran 77)
			ac_cv_use_blas_mangling=lowercaseNoUnderscore
                else
			AC_CHECK_LIB(blas, DGEMM, [tac_blas_ok=yes; BLAS_LIBS="-lblas"])
                        if test $tac_blas_ok = yes; then
                                AC_LANG_PUSH(Fortran 77)
                                AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME], [Define the Fortran name mangling to be used for the BLAS])
                                AC_LANG_POP(Fortran 77)
				ac_cv_use_blas_mangling=uppercaseNoUnderscore
                        else
				AC_CHECK_LIB(blas, DGEMM_, [tac_blas_ok=yes; BLAS_LIBS="-lblas"])
                                if test $tac_blas_ok = yes; then
                                        AC_LANG_PUSH(Fortran 77)
                                        AC_DEFINE([F77_BLAS_MANGLE(name,NAME)], [NAME ## _], [Define the Fortran name mangling to be used for the BLAS])
                                        AC_LANG_POP(Fortran 77)
					ac_cv_use_blas_mangling=uppercaseUnderscore
                                fi
                        fi
                 fi
        fi
fi

AC_SUBST(BLAS_LIBS)

fi # If the user specified library wasn't found, we skipped the remaining
   # checks.

LIBS="$tac_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$tac_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        tac_blas_ok=no
        $2
fi
])dnl TAC_BLAS
