dnl Check case (upper or lower) of F90 module files.  
dnl Also checks module suffix, but we return only ax_cv_f90_modulecase.
AC_DEFUN([AX_F90_MODULE_CASE],[
AC_CACHE_CHECK([fortran 90 module file suffix and case],
ax_cv_f90_modulecase,
[
rm -f conftest*
cat >conftest.f <<EOF
	module conftest
        integer n
        parameter (n=1)
        end module conftest
EOF
# SGI and absoft compilers generates module name in upper case!
testname="conftest"
modcase="lower"
echo "KDDKDD CASE 1" ${FC}
if ${FC} -c conftest.f > conftest.out 2>&1 ; then
    FCMODSUFFIX=`ls conftest* | grep -v conftest.f | grep -v conftest.o`
    echo "KDDKDD CASE 2" ${FCMODSUFFIX}
    FCMODSUFFIX=`echo "${FCMODSUFFIX}" | sed -e 's/conftest\.//g'`
    if test -z "${FCMODSUFFIX}" ; then
	FCMODSUFFIX=`ls CONFTEST* 2>/dev/null \
		| grep -v CONFTEST.f | grep -v CONFTEST.o`
        FCMODSUFFIX=`echo "${FCMODSUFFIX}" | sed -e 's/CONFTEST\.//g'`
	if test -n "${FCMODSUFFIX}" ; then
	    testname="CONFTEST"
	    modcase="upper"
	fi
    fi
    if test -z "${FCMODSUFFIX}" ; then 
        AC_MSG_RESULT(unknown)
	# Use mod if we can't figure it out
	FCMODSUFFIX="mod"   
    else
        AC_MSG_RESULT(${FCMODSUFFIX})
    fi
else
    AC_MSG_RESULT(unknown)
fi
#AC_SUBST(FCMODSUFFIX)
AC_MSG_CHECKING(for case of module names)
if test "${modcase}" = "lower" ; then
    AC_MSG_RESULT(lower)
    ax_cv_f90_modulecase="lower"
else
    AC_MSG_RESULT(upper)
    ax_cv_f90_modulecase="upper"
fi
])])
