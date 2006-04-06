# generated automatically by aclocal 1.7.6 -*- Autoconf -*-

# Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002
# Free Software Foundation, Inc.
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

# Like AC_CONFIG_HEADER, but automatically create stamp file. -*- Autoconf -*-

# Copyright 1996, 1997, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

AC_PREREQ([2.52])

# serial 6

# AM_CONFIG_HEADER is obsolete.  It has been replaced by AC_CONFIG_HEADERS.
AU_DEFUN([AM_CONFIG_HEADER], [AC_CONFIG_HEADERS($@)])

# Add --enable-maintainer-mode option to configure.
# From Jim Meyering

# Copyright 1996, 1998, 2000, 2001, 2002  Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 2

AC_DEFUN([AM_MAINTAINER_MODE],
[AC_MSG_CHECKING([whether to enable maintainer-specific portions of Makefiles])
  dnl maintainer-mode is disabled by default
  AC_ARG_ENABLE(maintainer-mode,
[  --enable-maintainer-mode enable make rules and dependencies not useful
                          (and sometimes confusing) to the casual installer],
      USE_MAINTAINER_MODE=$enableval,
      USE_MAINTAINER_MODE=no)
  AC_MSG_RESULT([$USE_MAINTAINER_MODE])
  AM_CONDITIONAL(MAINTAINER_MODE, [test $USE_MAINTAINER_MODE = yes])
  MAINT=$MAINTAINER_MODE_TRUE
  AC_SUBST(MAINT)dnl
]
)

AU_DEFUN([jm_MAINTAINER_MODE], [AM_MAINTAINER_MODE])

# AM_CONDITIONAL                                              -*- Autoconf -*-

# Copyright 1997, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 5

AC_PREREQ(2.52)

# AM_CONDITIONAL(NAME, SHELL-CONDITION)
# -------------------------------------
# Define a conditional.
AC_DEFUN([AM_CONDITIONAL],
[ifelse([$1], [TRUE],  [AC_FATAL([$0: invalid condition: $1])],
        [$1], [FALSE], [AC_FATAL([$0: invalid condition: $1])])dnl
AC_SUBST([$1_TRUE])
AC_SUBST([$1_FALSE])
if $2; then
  $1_TRUE=
  $1_FALSE='#'
else
  $1_TRUE='#'
  $1_FALSE=
fi
AC_CONFIG_COMMANDS_PRE(
[if test -z "${$1_TRUE}" && test -z "${$1_FALSE}"; then
  AC_MSG_ERROR([conditional "$1" was never defined.
Usually this means the macro was only invoked conditionally.])
fi])])

# Do all the work for Automake.                            -*- Autoconf -*-

# This macro actually does too much some checks are only needed if
# your package does certain things.  But this isn't really a big deal.

# Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
# Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 10

AC_PREREQ([2.54])

# Autoconf 2.50 wants to disallow AM_ names.  We explicitly allow
# the ones we care about.
m4_pattern_allow([^AM_[A-Z]+FLAGS$])dnl

# AM_INIT_AUTOMAKE(PACKAGE, VERSION, [NO-DEFINE])
# AM_INIT_AUTOMAKE([OPTIONS])
# -----------------------------------------------
# The call with PACKAGE and VERSION arguments is the old style
# call (pre autoconf-2.50), which is being phased out.  PACKAGE
# and VERSION should now be passed to AC_INIT and removed from
# the call to AM_INIT_AUTOMAKE.
# We support both call styles for the transition.  After
# the next Automake release, Autoconf can make the AC_INIT
# arguments mandatory, and then we can depend on a new Autoconf
# release and drop the old call support.
AC_DEFUN([AM_INIT_AUTOMAKE],
[AC_REQUIRE([AM_SET_CURRENT_AUTOMAKE_VERSION])dnl
 AC_REQUIRE([AC_PROG_INSTALL])dnl
# test to see if srcdir already configured
if test "`cd $srcdir && pwd`" != "`pwd`" &&
   test -f $srcdir/config.status; then
  AC_MSG_ERROR([source directory already configured; run "make distclean" there first])
fi

# test whether we have cygpath
if test -z "$CYGPATH_W"; then
  if (cygpath --version) >/dev/null 2>/dev/null; then
    CYGPATH_W='cygpath -w'
  else
    CYGPATH_W=echo
  fi
fi
AC_SUBST([CYGPATH_W])

# Define the identity of the package.
dnl Distinguish between old-style and new-style calls.
m4_ifval([$2],
[m4_ifval([$3], [_AM_SET_OPTION([no-define])])dnl
 AC_SUBST([PACKAGE], [$1])dnl
 AC_SUBST([VERSION], [$2])],
[_AM_SET_OPTIONS([$1])dnl
 AC_SUBST([PACKAGE], ['AC_PACKAGE_TARNAME'])dnl
 AC_SUBST([VERSION], ['AC_PACKAGE_VERSION'])])dnl

_AM_IF_OPTION([no-define],,
[AC_DEFINE_UNQUOTED(PACKAGE, "$PACKAGE", [Name of package])
 AC_DEFINE_UNQUOTED(VERSION, "$VERSION", [Version number of package])])dnl

# Some tools Automake needs.
AC_REQUIRE([AM_SANITY_CHECK])dnl
AC_REQUIRE([AC_ARG_PROGRAM])dnl
AM_MISSING_PROG(ACLOCAL, aclocal-${am__api_version})
AM_MISSING_PROG(AUTOCONF, autoconf)
AM_MISSING_PROG(AUTOMAKE, automake-${am__api_version})
AM_MISSING_PROG(AUTOHEADER, autoheader)
AM_MISSING_PROG(MAKEINFO, makeinfo)
AM_MISSING_PROG(AMTAR, tar)
AM_PROG_INSTALL_SH
AM_PROG_INSTALL_STRIP
# We need awk for the "check" target.  The system "awk" is bad on
# some platforms.
AC_REQUIRE([AC_PROG_AWK])dnl
AC_REQUIRE([AC_PROG_MAKE_SET])dnl
AC_REQUIRE([AM_SET_LEADING_DOT])dnl

_AM_IF_OPTION([no-dependencies],,
[AC_PROVIDE_IFELSE([AC_PROG_CC],
                  [_AM_DEPENDENCIES(CC)],
                  [define([AC_PROG_CC],
                          defn([AC_PROG_CC])[_AM_DEPENDENCIES(CC)])])dnl
AC_PROVIDE_IFELSE([AC_PROG_CXX],
                  [_AM_DEPENDENCIES(CXX)],
                  [define([AC_PROG_CXX],
                          defn([AC_PROG_CXX])[_AM_DEPENDENCIES(CXX)])])dnl
])
])


# When config.status generates a header, we must update the stamp-h file.
# This file resides in the same directory as the config header
# that is generated.  The stamp files are numbered to have different names.

# Autoconf calls _AC_AM_CONFIG_HEADER_HOOK (when defined) in the
# loop where config.status creates the headers, so we can generate
# our stamp files there.
AC_DEFUN([_AC_AM_CONFIG_HEADER_HOOK],
[# Compute $1's index in $config_headers.
_am_stamp_count=1
for _am_header in $config_headers :; do
  case $_am_header in
    $1 | $1:* )
      break ;;
    * )
      _am_stamp_count=`expr $_am_stamp_count + 1` ;;
  esac
done
echo "timestamp for $1" >`AS_DIRNAME([$1])`/stamp-h[]$_am_stamp_count])

# Copyright 2002  Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA

# AM_AUTOMAKE_VERSION(VERSION)
# ----------------------------
# Automake X.Y traces this macro to ensure aclocal.m4 has been
# generated from the m4 files accompanying Automake X.Y.
AC_DEFUN([AM_AUTOMAKE_VERSION],[am__api_version="1.7"])

# AM_SET_CURRENT_AUTOMAKE_VERSION
# -------------------------------
# Call AM_AUTOMAKE_VERSION so it can be traced.
# This function is AC_REQUIREd by AC_INIT_AUTOMAKE.
AC_DEFUN([AM_SET_CURRENT_AUTOMAKE_VERSION],
	 [AM_AUTOMAKE_VERSION([1.7.6])])

# Helper functions for option handling.                    -*- Autoconf -*-

# Copyright 2001, 2002  Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 2

# _AM_MANGLE_OPTION(NAME)
# -----------------------
AC_DEFUN([_AM_MANGLE_OPTION],
[[_AM_OPTION_]m4_bpatsubst($1, [[^a-zA-Z0-9_]], [_])])

# _AM_SET_OPTION(NAME)
# ------------------------------
# Set option NAME.  Presently that only means defining a flag for this option.
AC_DEFUN([_AM_SET_OPTION],
[m4_define(_AM_MANGLE_OPTION([$1]), 1)])

# _AM_SET_OPTIONS(OPTIONS)
# ----------------------------------
# OPTIONS is a space-separated list of Automake options.
AC_DEFUN([_AM_SET_OPTIONS],
[AC_FOREACH([_AM_Option], [$1], [_AM_SET_OPTION(_AM_Option)])])

# _AM_IF_OPTION(OPTION, IF-SET, [IF-NOT-SET])
# -------------------------------------------
# Execute IF-SET if OPTION is set, IF-NOT-SET otherwise.
AC_DEFUN([_AM_IF_OPTION],
[m4_ifset(_AM_MANGLE_OPTION([$1]), [$2], [$3])])

#
# Check to make sure that the build environment is sane.
#

# Copyright 1996, 1997, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 3

# AM_SANITY_CHECK
# ---------------
AC_DEFUN([AM_SANITY_CHECK],
[AC_MSG_CHECKING([whether build environment is sane])
# Just in case
sleep 1
echo timestamp > conftest.file
# Do `set' in a subshell so we don't clobber the current shell's
# arguments.  Must try -L first in case configure is actually a
# symlink; some systems play weird games with the mod time of symlinks
# (eg FreeBSD returns the mod time of the symlink's containing
# directory).
if (
   set X `ls -Lt $srcdir/configure conftest.file 2> /dev/null`
   if test "$[*]" = "X"; then
      # -L didn't work.
      set X `ls -t $srcdir/configure conftest.file`
   fi
   rm -f conftest.file
   if test "$[*]" != "X $srcdir/configure conftest.file" \
      && test "$[*]" != "X conftest.file $srcdir/configure"; then

      # If neither matched, then we have a broken ls.  This can happen
      # if, for instance, CONFIG_SHELL is bash and it inherits a
      # broken ls alias from the environment.  This has actually
      # happened.  Such a system could not be considered "sane".
      AC_MSG_ERROR([ls -t appears to fail.  Make sure there is not a broken
alias in your environment])
   fi

   test "$[2]" = conftest.file
   )
then
   # Ok.
   :
else
   AC_MSG_ERROR([newly created file is older than distributed files!
Check your system clock])
fi
AC_MSG_RESULT(yes)])

#  -*- Autoconf -*-


# Copyright 1997, 1999, 2000, 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 3

# AM_MISSING_PROG(NAME, PROGRAM)
# ------------------------------
AC_DEFUN([AM_MISSING_PROG],
[AC_REQUIRE([AM_MISSING_HAS_RUN])
$1=${$1-"${am_missing_run}$2"}
AC_SUBST($1)])


# AM_MISSING_HAS_RUN
# ------------------
# Define MISSING if not defined so far and test if it supports --run.
# If it does, set am_missing_run to use it, otherwise, to nothing.
AC_DEFUN([AM_MISSING_HAS_RUN],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
test x"${MISSING+set}" = xset || MISSING="\${SHELL} $am_aux_dir/missing"
# Use eval to expand $SHELL
if eval "$MISSING --run true"; then
  am_missing_run="$MISSING --run "
else
  am_missing_run=
  AC_MSG_WARN([`missing' script is too old or missing])
fi
])

# AM_AUX_DIR_EXPAND

# Copyright 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# For projects using AC_CONFIG_AUX_DIR([foo]), Autoconf sets
# $ac_aux_dir to `$srcdir/foo'.  In other projects, it is set to
# `$srcdir', `$srcdir/..', or `$srcdir/../..'.
#
# Of course, Automake must honor this variable whenever it calls a
# tool from the auxiliary directory.  The problem is that $srcdir (and
# therefore $ac_aux_dir as well) can be either absolute or relative,
# depending on how configure is run.  This is pretty annoying, since
# it makes $ac_aux_dir quite unusable in subdirectories: in the top
# source directory, any form will work fine, but in subdirectories a
# relative path needs to be adjusted first.
#
# $ac_aux_dir/missing
#    fails when called from a subdirectory if $ac_aux_dir is relative
# $top_srcdir/$ac_aux_dir/missing
#    fails if $ac_aux_dir is absolute,
#    fails when called from a subdirectory in a VPATH build with
#          a relative $ac_aux_dir
#
# The reason of the latter failure is that $top_srcdir and $ac_aux_dir
# are both prefixed by $srcdir.  In an in-source build this is usually
# harmless because $srcdir is `.', but things will broke when you
# start a VPATH build or use an absolute $srcdir.
#
# So we could use something similar to $top_srcdir/$ac_aux_dir/missing,
# iff we strip the leading $srcdir from $ac_aux_dir.  That would be:
#   am_aux_dir='\$(top_srcdir)/'`expr "$ac_aux_dir" : "$srcdir//*\(.*\)"`
# and then we would define $MISSING as
#   MISSING="\${SHELL} $am_aux_dir/missing"
# This will work as long as MISSING is not called from configure, because
# unfortunately $(top_srcdir) has no meaning in configure.
# However there are other variables, like CC, which are often used in
# configure, and could therefore not use this "fixed" $ac_aux_dir.
#
# Another solution, used here, is to always expand $ac_aux_dir to an
# absolute PATH.  The drawback is that using absolute paths prevent a
# configured tree to be moved without reconfiguration.

# Rely on autoconf to set up CDPATH properly.
AC_PREREQ([2.50])

AC_DEFUN([AM_AUX_DIR_EXPAND], [
# expand $ac_aux_dir to an absolute path
am_aux_dir=`cd $ac_aux_dir && pwd`
])

# AM_PROG_INSTALL_SH
# ------------------
# Define $install_sh.

# Copyright 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

AC_DEFUN([AM_PROG_INSTALL_SH],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
install_sh=${install_sh-"$am_aux_dir/install-sh"}
AC_SUBST(install_sh)])

# AM_PROG_INSTALL_STRIP

# Copyright 2001 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# One issue with vendor `install' (even GNU) is that you can't
# specify the program used to strip binaries.  This is especially
# annoying in cross-compiling environments, where the build's strip
# is unlikely to handle the host's binaries.
# Fortunately install-sh will honor a STRIPPROG variable, so we
# always use install-sh in `make install-strip', and initialize
# STRIPPROG with the value of the STRIP variable (set by the user).
AC_DEFUN([AM_PROG_INSTALL_STRIP],
[AC_REQUIRE([AM_PROG_INSTALL_SH])dnl
# Installed binaries are usually stripped using `strip' when the user
# run `make install-strip'.  However `strip' might not be the right
# tool to use in cross-compilation environments, therefore Automake
# will honor the `STRIP' environment variable to overrule this program.
dnl Don't test for $cross_compiling = yes, because it might be `maybe'.
if test "$cross_compiling" != no; then
  AC_CHECK_TOOL([STRIP], [strip], :)
fi
INSTALL_STRIP_PROGRAM="\${SHELL} \$(install_sh) -c -s"
AC_SUBST([INSTALL_STRIP_PROGRAM])])

#                                                          -*- Autoconf -*-
# Copyright (C) 2003  Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 1

# Check whether the underlying file-system supports filenames
# with a leading dot.  For instance MS-DOS doesn't.
AC_DEFUN([AM_SET_LEADING_DOT],
[rm -rf .tst 2>/dev/null
mkdir .tst 2>/dev/null
if test -d .tst; then
  am__leading_dot=.
else
  am__leading_dot=_
fi
rmdir .tst 2>/dev/null
AC_SUBST([am__leading_dot])])

# serial 5						-*- Autoconf -*-

# Copyright (C) 1999, 2000, 2001, 2002, 2003  Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.


# There are a few dirty hacks below to avoid letting `AC_PROG_CC' be
# written in clear, in which case automake, when reading aclocal.m4,
# will think it sees a *use*, and therefore will trigger all it's
# C support machinery.  Also note that it means that autoscan, seeing
# CC etc. in the Makefile, will ask for an AC_PROG_CC use...



# _AM_DEPENDENCIES(NAME)
# ----------------------
# See how the compiler implements dependency checking.
# NAME is "CC", "CXX", "GCJ", or "OBJC".
# We try a few techniques and use that to set a single cache variable.
#
# We don't AC_REQUIRE the corresponding AC_PROG_CC since the latter was
# modified to invoke _AM_DEPENDENCIES(CC); we would have a circular
# dependency, and given that the user is not expected to run this macro,
# just rely on AC_PROG_CC.
AC_DEFUN([_AM_DEPENDENCIES],
[AC_REQUIRE([AM_SET_DEPDIR])dnl
AC_REQUIRE([AM_OUTPUT_DEPENDENCY_COMMANDS])dnl
AC_REQUIRE([AM_MAKE_INCLUDE])dnl
AC_REQUIRE([AM_DEP_TRACK])dnl

ifelse([$1], CC,   [depcc="$CC"   am_compiler_list=],
       [$1], CXX,  [depcc="$CXX"  am_compiler_list=],
       [$1], OBJC, [depcc="$OBJC" am_compiler_list='gcc3 gcc'],
       [$1], GCJ,  [depcc="$GCJ"  am_compiler_list='gcc3 gcc'],
                   [depcc="$$1"   am_compiler_list=])

AC_CACHE_CHECK([dependency style of $depcc],
               [am_cv_$1_dependencies_compiler_type],
[if test -z "$AMDEP_TRUE" && test -f "$am_depcomp"; then
  # We make a subdir and do the tests there.  Otherwise we can end up
  # making bogus files that we don't know about and never remove.  For
  # instance it was reported that on HP-UX the gcc test will end up
  # making a dummy file named `D' -- because `-MD' means `put the output
  # in D'.
  mkdir conftest.dir
  # Copy depcomp to subdir because otherwise we won't find it if we're
  # using a relative directory.
  cp "$am_depcomp" conftest.dir
  cd conftest.dir
  # We will build objects and dependencies in a subdirectory because
  # it helps to detect inapplicable dependency modes.  For instance
  # both Tru64's cc and ICC support -MD to output dependencies as a
  # side effect of compilation, but ICC will put the dependencies in
  # the current directory while Tru64 will put them in the object
  # directory.
  mkdir sub

  am_cv_$1_dependencies_compiler_type=none
  if test "$am_compiler_list" = ""; then
     am_compiler_list=`sed -n ['s/^#*\([a-zA-Z0-9]*\))$/\1/p'] < ./depcomp`
  fi
  for depmode in $am_compiler_list; do
    # Setup a source with many dependencies, because some compilers
    # like to wrap large dependency lists on column 80 (with \), and
    # we should not choose a depcomp mode which is confused by this.
    #
    # We need to recreate these files for each test, as the compiler may
    # overwrite some of them when testing with obscure command lines.
    # This happens at least with the AIX C compiler.
    : > sub/conftest.c
    for i in 1 2 3 4 5 6; do
      echo '#include "conftst'$i'.h"' >> sub/conftest.c
      : > sub/conftst$i.h
    done
    echo "${am__include} ${am__quote}sub/conftest.Po${am__quote}" > confmf

    case $depmode in
    nosideeffect)
      # after this tag, mechanisms are not by side-effect, so they'll
      # only be used when explicitly requested
      if test "x$enable_dependency_tracking" = xyes; then
	continue
      else
	break
      fi
      ;;
    none) break ;;
    esac
    # We check with `-c' and `-o' for the sake of the "dashmstdout"
    # mode.  It turns out that the SunPro C++ compiler does not properly
    # handle `-M -o', and we need to detect this.
    if depmode=$depmode \
       source=sub/conftest.c object=sub/conftest.${OBJEXT-o} \
       depfile=sub/conftest.Po tmpdepfile=sub/conftest.TPo \
       $SHELL ./depcomp $depcc -c -o sub/conftest.${OBJEXT-o} sub/conftest.c \
         >/dev/null 2>conftest.err &&
       grep sub/conftst6.h sub/conftest.Po > /dev/null 2>&1 &&
       grep sub/conftest.${OBJEXT-o} sub/conftest.Po > /dev/null 2>&1 &&
       ${MAKE-make} -s -f confmf > /dev/null 2>&1; then
      # icc doesn't choke on unknown options, it will just issue warnings
      # (even with -Werror).  So we grep stderr for any message
      # that says an option was ignored.
      if grep 'ignoring option' conftest.err >/dev/null 2>&1; then :; else
        am_cv_$1_dependencies_compiler_type=$depmode
        break
      fi
    fi
  done

  cd ..
  rm -rf conftest.dir
else
  am_cv_$1_dependencies_compiler_type=none
fi
])
AC_SUBST([$1DEPMODE], [depmode=$am_cv_$1_dependencies_compiler_type])
AM_CONDITIONAL([am__fastdep$1], [
  test "x$enable_dependency_tracking" != xno \
  && test "$am_cv_$1_dependencies_compiler_type" = gcc3])
])


# AM_SET_DEPDIR
# -------------
# Choose a directory name for dependency files.
# This macro is AC_REQUIREd in _AM_DEPENDENCIES
AC_DEFUN([AM_SET_DEPDIR],
[AC_REQUIRE([AM_SET_LEADING_DOT])dnl
AC_SUBST([DEPDIR], ["${am__leading_dot}deps"])dnl
])


# AM_DEP_TRACK
# ------------
AC_DEFUN([AM_DEP_TRACK],
[AC_ARG_ENABLE(dependency-tracking,
[  --disable-dependency-tracking Speeds up one-time builds
  --enable-dependency-tracking  Do not reject slow dependency extractors])
if test "x$enable_dependency_tracking" != xno; then
  am_depcomp="$ac_aux_dir/depcomp"
  AMDEPBACKSLASH='\'
fi
AM_CONDITIONAL([AMDEP], [test "x$enable_dependency_tracking" != xno])
AC_SUBST([AMDEPBACKSLASH])
])

# Generate code to set up dependency tracking.   -*- Autoconf -*-

# Copyright 1999, 2000, 2001, 2002 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

#serial 2

# _AM_OUTPUT_DEPENDENCY_COMMANDS
# ------------------------------
AC_DEFUN([_AM_OUTPUT_DEPENDENCY_COMMANDS],
[for mf in $CONFIG_FILES; do
  # Strip MF so we end up with the name of the file.
  mf=`echo "$mf" | sed -e 's/:.*$//'`
  # Check whether this is an Automake generated Makefile or not.
  # We used to match only the files named `Makefile.in', but
  # some people rename them; so instead we look at the file content.
  # Grep'ing the first line is not enough: some people post-process
  # each Makefile.in and add a new line on top of each file to say so.
  # So let's grep whole file.
  if grep '^#.*generated by automake' $mf > /dev/null 2>&1; then
    dirpart=`AS_DIRNAME("$mf")`
  else
    continue
  fi
  grep '^DEP_FILES *= *[[^ @%:@]]' < "$mf" > /dev/null || continue
  # Extract the definition of DEP_FILES from the Makefile without
  # running `make'.
  DEPDIR=`sed -n -e '/^DEPDIR = / s///p' < "$mf"`
  test -z "$DEPDIR" && continue
  # When using ansi2knr, U may be empty or an underscore; expand it
  U=`sed -n -e '/^U = / s///p' < "$mf"`
  test -d "$dirpart/$DEPDIR" || mkdir "$dirpart/$DEPDIR"
  # We invoke sed twice because it is the simplest approach to
  # changing $(DEPDIR) to its actual value in the expansion.
  for file in `sed -n -e '
    /^DEP_FILES = .*\\\\$/ {
      s/^DEP_FILES = //
      :loop
	s/\\\\$//
	p
	n
	/\\\\$/ b loop
      p
    }
    /^DEP_FILES = / s/^DEP_FILES = //p' < "$mf" | \
       sed -e 's/\$(DEPDIR)/'"$DEPDIR"'/g' -e 's/\$U/'"$U"'/g'`; do
    # Make sure the directory exists.
    test -f "$dirpart/$file" && continue
    fdir=`AS_DIRNAME(["$file"])`
    AS_MKDIR_P([$dirpart/$fdir])
    # echo "creating $dirpart/$file"
    echo '# dummy' > "$dirpart/$file"
  done
done
])# _AM_OUTPUT_DEPENDENCY_COMMANDS


# AM_OUTPUT_DEPENDENCY_COMMANDS
# -----------------------------
# This macro should only be invoked once -- use via AC_REQUIRE.
#
# This code is only required when automatic dependency tracking
# is enabled.  FIXME.  This creates each `.P' file that we will
# need in order to bootstrap the dependency handling code.
AC_DEFUN([AM_OUTPUT_DEPENDENCY_COMMANDS],
[AC_CONFIG_COMMANDS([depfiles],
     [test x"$AMDEP_TRUE" != x"" || _AM_OUTPUT_DEPENDENCY_COMMANDS],
     [AMDEP_TRUE="$AMDEP_TRUE" ac_aux_dir="$ac_aux_dir"])
])

# Check to see how 'make' treats includes.	-*- Autoconf -*-

# Copyright (C) 2001, 2002, 2003 Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# serial 2

# AM_MAKE_INCLUDE()
# -----------------
# Check to see how make treats includes.
AC_DEFUN([AM_MAKE_INCLUDE],
[am_make=${MAKE-make}
cat > confinc << 'END'
am__doit:
	@echo done
.PHONY: am__doit
END
# If we don't find an include directive, just comment out the code.
AC_MSG_CHECKING([for style of include used by $am_make])
am__include="#"
am__quote=
_am_result=none
# First try GNU make style include.
echo "include confinc" > confmf
# We grep out `Entering directory' and `Leaving directory'
# messages which can occur if `w' ends up in MAKEFLAGS.
# In particular we don't look at `^make:' because GNU make might
# be invoked under some other name (usually "gmake"), in which
# case it prints its new name instead of `make'.
if test "`$am_make -s -f confmf 2> /dev/null | grep -v 'ing directory'`" = "done"; then
   am__include=include
   am__quote=
   _am_result=GNU
fi
# Now try BSD make style include.
if test "$am__include" = "#"; then
   echo '.include "confinc"' > confmf
   if test "`$am_make -s -f confmf 2> /dev/null`" = "done"; then
      am__include=.include
      am__quote="\""
      _am_result=BSD
   fi
fi
AC_SUBST([am__include])
AC_SUBST([am__quote])
AC_MSG_RESULT([$_am_result])
rm -f confinc confmf
])

dnl @synopsis TAC_ARG_CONFIG_MPI
dnl
dnl Test a variety of MPI options:
dnl --enable-mpi       - Turns MPI compiling mode on
dnl --with-mpi         - specify root directory of MPI
dnl --with-mpi-compilers - Turns on MPI compiling mode and sets the MPI C++
dnl                       compiler = mpicxx, mpic++ or mpiCC,
dnl                       the MPI C compiler = mpicc and 
dnl                       the MPI Fortran compiler = mpif77
dnl --with-mpi-incdir - specify include directory for MPI 
dnl --with-mpi-libs    - specify MPI libraries
dnl --with-mpi-libdir  - specify location of MPI libraries
dnl
dnl If any of these options are set, HAVE_MPI will be defined for both
dnl Autoconf and Automake, and HAVE_MPI will be defined in the
dnl generated config.h file
dnl
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_CONFIG_MPI],
[

AC_ARG_ENABLE(mpi,
[AC_HELP_STRING([--enable-mpi],[MPI support])],
[HAVE_PKG_MPI=$enableval],
[HAVE_PKG_MPI=no]
)

AC_ARG_WITH(mpi-compilers,
[AC_HELP_STRING([--with-mpi-compilers=PATH],
[use MPI compilers mpicc, mpif77, and mpicxx, mpic++ or mpiCC in the specified path or in the default path if no path is specified. Enables MPI])],
[
  if test X${withval} != Xno; then
    HAVE_PKG_MPI=yes
    if test X${withval} = Xyes; then
      # Check for mpicxx, if it does not exist, check for mpic++, if it does 
      # not exist, use mpiCC instead.
      AC_CHECK_PROG(MPI_TEMP_CXX, mpicxx, mpicxx, no)
      if test X${MPI_TEMP_CXX} = Xno; then
	AC_CHECK_PROG(MPI_CXX, mpic++, mpic++, mpiCC)
      else 
	MPI_CXX=${MPI_TEMP_CXX}
      fi
      MPI_CC=mpicc
      MPI_F77=mpif77
    else
      if test -f ${withval}/mpicxx; then
        MPI_CXX=${withval}/mpicxx
      elif test -f ${withval}/mpic++; then
	MPI_CXX=${withval}/mpic++
      else
        MPI_CXX=${withval}/mpiCC
      fi
      MPI_CC=${withval}/mpicc
      MPI_F77=${withval}/mpif77
    fi
  fi
]
)

AC_ARG_WITH(mpi,
[AC_HELP_STRING([--with-mpi=MPIROOT],[use MPI root directory (enables MPI)])],
[
  HAVE_PKG_MPI=yes
  MPI_DIR=${withval}
  AC_MSG_CHECKING(MPI directory)
  AC_MSG_RESULT([${MPI_DIR}])
]
)

#AC_ARG_WITH(mpi-include,
#[AC_HELP_STRING([--with-mpi-include],[Obsolete.  Use --with-mpi-incdir=DIR instead.  Do not prefix DIR with '-I'.])],
#[AC_MSG_ERROR([--with-mpi-include is an obsolte option.   Use --with-mpi-incdir=DIR instead.  Do not prefix DIR with '-I'.  For example '--with-mpi-incdir=/usr/lam_path/include'.])]
#)

AC_ARG_WITH(mpi-libs,
[AC_HELP_STRING([--with-mpi-libs="LIBS"],[MPI libraries @<:@"-lmpi"@:>@])],
[
  MPI_LIBS=${withval}
  AC_MSG_CHECKING(user-defined MPI libraries)
  AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_ARG_WITH(mpi-incdir,
[AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@  Do not use -I])],
[
  MPI_INC=${withval}
  AC_MSG_CHECKING(user-defined MPI includes)
  AC_MSG_RESULT([${MPI_INC}])
]
)

AC_ARG_WITH(mpi-libdir,
[AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@  Do not use -L])],
[
  MPI_LIBDIR=${withval}
  AC_MSG_CHECKING(user-defined MPI library directory)
  AC_MSG_RESULT([${MPI_LIBDIR}])
]
)

AC_MSG_CHECKING(whether we are using MPI)
AC_MSG_RESULT([${HAVE_PKG_MPI}])

if test "X${HAVE_PKG_MPI}" = "Xyes"; then
   AC_DEFINE(HAVE_MPI,,[define if we want to use MPI])
fi

dnl Define Automake version of HAVE_MPI if appropriate

AM_CONDITIONAL(HAVE_MPI, [test "X${HAVE_PKG_MPI}" = "Xyes"])


dnl
dnl --------------------------------------------------------------------
dnl Check for MPI compilers (must be done *before* AC_PROG_CXX,
dnl AC_PROG_CC and AC_PROG_F77)
dnl 
dnl --------------------------------------------------------------------

if test -n "${MPI_CXX}"; then
  if test -f ${MPI_CXX}; then
    MPI_CXX_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_CXX_EXISTS, ${MPI_CXX}, yes, no)
  fi

  if test "X${MPI_CXX_EXISTS}" = "Xyes"; then
    CXX=${MPI_CXX}
  else
    echo "-----"
    echo "Cannot find MPI C++ compiler ${MPI_CXX}."
    echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH"
    echo "or specify a C++ compiler using CXX=<compiler>"
    echo "Do not use --with-mpi-compilers if using CXX=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI C++ compiler (${MPI_CXX}) not found.])
  fi
fi

if test -n "${MPI_CC}"; then
  if test -f ${MPI_CC}; then
    MPI_CC_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_CC_EXISTS, ${MPI_CC}, yes, no)
  fi

  if test "X${MPI_CC_EXISTS}" = "Xyes"; then
    CC=${MPI_CC}
  else
    echo "-----"
    echo "Cannot find MPI C compiler ${MPI_CC}."
    echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH"
    echo "or specify a C compiler using CC=<compiler>"
    echo "Do not use --with-mpi-compilers if using CC=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI C compiler (${MPI_CC}) not found.])
  fi
fi

if test -n "${MPI_F77}"; then
  if test -f ${MPI_F77}; then
    MPI_F77_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_F77_EXISTS, ${MPI_F77}, yes, no)
  fi

  if test "X${MPI_F77_EXISTS}" = "Xyes"; then
    F77=${MPI_F77}
  else
    echo "-----"
    echo "Cannot find MPI Fortran compiler ${MPI_F77}."
    echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH"
    echo "or specify a Fortran 77 compiler using F77=<compiler>"
    echo "Do not use --with-mpi-compilers if using F77=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI Fortran 77 compiler (${MPI_F77}) not found.])
  fi
fi
])

dnl @synopsis TAC_ARG_WITH_FLAGS(lcase_name, UCASE_NAME)
dnl
dnl Test for --with-lcase_name="compiler/loader flags".  if defined, prepend 
dnl flags to standard UCASE_NAME definition.
dnl
dnl Use this macro to facilitate additional special flags that should be
dnl passed on to the preprocessor/compilers/loader.
dnl
dnl Example use
dnl 
dnl TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
dnl 
dnl tests for --with-cxxflags and pre-pends to CXXFLAGS
dnl 
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_FLAGS],
[
AC_MSG_CHECKING([whether additional [$2] flags should be added])
AC_ARG_WITH($1,
AC_HELP_STRING([--with-$1], 
[additional [$2] flags to be added: will prepend to [$2]]),
[
$2="${withval} ${$2}"
AC_MSG_RESULT([$2 = ${$2}])
],
AC_MSG_RESULT(no)
)
])


dnl @synopsis TAC_ARG_WITH_LIBS
dnl
dnl Test for --with-libs="name(s)".
dnl 
dnl Prepends the specified name(s) to the list of libraries to link 
dnl with.  
dnl
dnl Example use
dnl
dnl TAC_ARG_WITH_LIBS
dnl 
dnl tests for --with-libs and pre-pends to LIBS
dnl
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_LIBS],
[
AC_MSG_CHECKING([whether additional libraries are needed])
AC_ARG_WITH(libs,
AC_HELP_STRING([--with-libs], 
[List additional libraries here.  For example, --with-libs=-lsuperlu
or --with-libs=/path/libsuperlu.a]),
[
LIBS="${withval} ${LIBS}"
AC_MSG_RESULT([LIBS = ${LIBS}])
],
AC_MSG_RESULT(no)
)
]
)

dnl @synopsis TAC_ARG_WITH_AR
dnl
dnl Test for --with-ar="ar_program ar_flags".
dnl Default is "ar cru"
dnl 
dnl Generates an Automake conditional USE_ALTERNATE_AR that can be tested.  
dnl Generates the user-specified archiver command in @ALTERNATE_AR@.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_AR],
[
AC_ARG_WITH(ar,
AC_HELP_STRING([--with-ar], [override archiver command (default is "ar cru")]),
[
AC_MSG_CHECKING(user-defined archiver)
AC_MSG_RESULT([${withval}])
USE_ALTERNATE_AR=yes
ALTERNATE_AR="${withval}"
]
)

if test -n "${SPECIAL_AR}" && test "X${USE_ALTERNATE_AR}" != "Xyes";
then
  USE_ALTERNATE_AR=yes
  ALTERNATE_AR="${SPECIAL_AR}"
fi

AC_MSG_CHECKING(for special archiver command)
if test "X${USE_ALTERNATE_AR}" = "Xyes"; then
   AC_MSG_RESULT([${ALTERNATE_AR}])
   AM_CONDITIONAL(USE_ALTERNATE_AR, true)
else
   AC_MSG_RESULT([none])
   AM_CONDITIONAL(USE_ALTERNATE_AR, false)
fi
AC_SUBST(ALTERNATE_AR)
])


dnl @synopsis TAC_ARG_CHECK_MPI
dnl
dnl Check to make sure any definitions set in TAC_ARG_CONFIG_MPI
dnl are valid, set the MPI flags.  Test MPI compile using C++ compiler.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_CHECK_MPI],
[

if test "X${HAVE_PKG_MPI}" = "Xyes"; then

  if test -n "${MPI_DIR}" && test -z "${MPI_INC}"; then
    MPI_INC="${MPI_DIR}/include"
  fi

  if test -n "${MPI_INC}"; then
    CPPFLAGS="${CPPFLAGS} -I${MPI_INC}"
  fi

  AC_LANG_CPLUSPLUS 
  AC_MSG_CHECKING(for mpi.h)
  AC_TRY_CPP([#include "mpi.h"],
    [AC_MSG_RESULT(yes)], 
    [
     AC_MSG_RESULT(no)  
     echo "-----"
     echo "Cannot link simple MPI program."
     echo "Try --with-mpi-compilers to specify MPI compilers."
     echo "Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
     echo "to specify all the specific MPI compile options."
     echo "-----"
     AC_MSG_ERROR(MPI cannot link)
    ])

  if test -n "${MPI_DIR}" && test -z "${MPI_LIBDIR}"; then
    MPI_LIBDIR="${MPI_DIR}/lib"
  fi

  if test -n "${MPI_LIBDIR}"; then
    LDFLAGS="${LDFLAGS} -L${MPI_LIBDIR}"
  fi

  if test -z "${MPI_LIBS}" && test -n "${MPI_LIBDIR}"; then
    MPI_LIBS="-lmpi"
  fi

  if test -n "${MPI_LIBS}"; then
    LIBS="${MPI_LIBS} ${LIBS}"
  fi

#   AC_LANG_CPLUSPLUS 
#   AC_MSG_CHECKING(whether MPI will link using C++ compiler)
#   AC_TRY_LINK([#include <mpi.h>],
#   [int c; char** v; MPI_Init(&c,&v);],
#   [AC_MSG_RESULT(yes)], 
#   [AC_MSG_RESULT(no)  
#    echo "-----"
#    echo "Cannot link simple MPI program."
#    echo "Try --with-mpi-cxx to specify MPI C++ compile script."
#    echo "Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
#    echo "to specify all the specific MPI compile options."
#    echo "-----"
#    AC_MSG_ERROR(MPI cannot link)]
#   )

fi
])

dnl Enables export makefile specific code
dnl 
dnl The following AM_CONDITIONALS are set for makefiles to access:
dnl USING_EXPORT_MAKEFILES
dnl USING_PERL via TAC_ARG_WITH_PERL
dnl USING_GNUMAKE
dnl
dnl The following AC_DEFINES are set:
dnl HAVE_EXPORT_MAKEFILES
dnl 
dnl the following variables are set:
dnl PERL_EXE for the perl executable via TAC_ARG_WITH_PERL
dnl 
dnl This file was based on tac_arg_enable_feature.m4 by Mike Heroux
dnl @author Roger Pawlowski <rppawlo@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_EXPORT_MAKEFILES],
[
AC_ARG_ENABLE(export-makefiles,
AC_HELP_STRING([--enable-export-makefiles],[Creates export makefiles in the install (prefix) directory.  This option requires perl to be set in your path or defined with --with-perl=<perl executable>. Note that the export makefiles are always created and used in the build directory, but will not be installable without this option to change the paths. (default is $1)]),
ac_cv_use_export_makefiles=$enableval, 
ac_cv_use_export_makefiles=$1)

AC_MSG_CHECKING(whether to build export makefiles)

if test "X$ac_cv_use_export_makefiles" != "Xno"; then

  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_EXPORT_MAKEFILES],,[Define if you want to build export makefiles.])

else

  AC_MSG_RESULT(no)

fi

AM_CONDITIONAL(USING_EXPORT_MAKEFILES, test X${ac_cv_use_export_makefiles} = Xyes)

# Check for perl to run scripts (Required dependency)
TAC_ARG_WITH_PERL

if test "X$HAVE_PERL" != "Xyes" && 
   test "X$ac_cv_use_export_makefiles" != "Xno"; then
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([Failed to find the perl executable.  The flag --enable-export-makefiles requires perl to be either in your path or explicitly defined by the flag --with-perl=<executable>.  If you do not require the export makefiles to be installed via 'make install', you can disable the export makefiles with --disable-export-makefiles.])
fi

# Check for using gnumake to clean up link lines via 
# gnumake's "shell" command. Optional dependency.
AC_DEFUN([TAC_ARG_WITH_GNUMAKE],
[
AC_ARG_WITH(gnumake,
AC_HELP_STRING([--with-gnumake],[Gnu's make has special functions we can use to eliminate redundant paths in the build and link lines. Enable this if you use gnu-make to build Trilinos. This requires that perl is in your path or that you have specified the perl executable with --with-perl=<perl executable>.  Configure will check for the existence of the perl executable and quit with an error if it is not found. (default is no)]),
ac_cv_use_gnumake=$withval, ac_cv_use_gnumake=no)

AC_MSG_CHECKING(whether gnumake specific code should be enabled)

if test "X$ac_cv_use_gnumake" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([HAVE_GNUMAKE],,[Define if you are using gnumake - this will shorten your link lines.])
else
  AC_MSG_RESULT(no)
fi
AM_CONDITIONAL(USING_GNUMAKE, test "X$ac_cv_use_gnumake" = "Xyes")
])

TAC_ARG_WITH_GNUMAKE

if test "X$HAVE_PERL" != "Xyes" && 
   test "X$ac_cv_use_gnumake" != "Xno"; then
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([The flag --with-gnumake requires perl to be in your path.  The perl executable can alternatively be explicitly defined by the flag --with-perl=<executable>.])
fi

])


dnl @synopsis TAC_ARG_WITH_PERL(DEFAULT_VAL)
dnl
dnl Test for --enable-gnumake and set to DEFAULT_VAL value if feature not specified.
dnl Calls AC_DEFINE to define HAVE_GNUMAKE if value is not equal to "no"
dnl Calls AM_CONDITIONAL to define USING_GNUMAKE to true/false.
dnl 
dnl This file was based on tac_arg_with_ar.m4 by Mike Heroux
dnl @author Roger Pawlowski <rppawlo@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_PERL],
[

AC_ARG_WITH(perl,
AC_HELP_STRING([--with-perl], [supply a perl executable.  For example --with-perl=/usr/bin/perl.]),
[
AC_MSG_CHECKING(for user supplied perl executable)
AC_MSG_RESULT([${withval}])
USER_SPECIFIED_PERL=yes
PERL_EXE="${withval}"
],
[
USER_SPECIFIED_PERL=no
])

if test "X${USER_SPECIFIED_PERL}" = "Xyes"; then
  AC_CHECK_FILE(${PERL_EXE}, [HAVE_PERL=yes], [HAVE_PERL=no])
  AC_SUBST(PERL_EXE, ${PERL_EXE})
else
  AC_CHECK_PROG(HAVE_PERL, perl, yes, no)
  AC_SUBST(PERL_EXE, perl)
fi
AM_CONDITIONAL(USING_PERL, test X${HAVE_PERL} = Xyes)
])


dnl @synopsis AC_CXX_NAMESPACES
dnl
dnl If the compiler can prevent names clashes using namespaces, define
dnl HAVE_NAMESPACES.
dnl
dnl @version $Id$
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_NAMESPACES],
[AC_CACHE_CHECK(whether the compiler implements namespaces,
ac_cv_cxx_namespaces,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([namespace Outer { namespace Inner { int i = 0; }}],
                [using namespace Outer::Inner; return i;],
 ac_cv_cxx_namespaces=yes, ac_cv_cxx_namespaces=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_namespaces" = yes; then
  AC_DEFINE(HAVE_NAMESPACES,,[define if the compiler implements namespaces])
fi
])

dnl @synopsis AC_CXX_HAVE_STL
dnl
dnl If the compiler supports the Standard Template Library, define HAVE_STL.
dnl
dnl @version $Id$
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_HAVE_STL],
[AC_CACHE_CHECK(whether the compiler supports Standard Template Library,
ac_cv_cxx_have_stl,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <list>
#include <deque>
#ifdef HAVE_NAMESPACES
using namespace std;
#endif],[list<int> x; x.push_back(5);
list<int>::iterator iter = x.begin(); if (iter != x.end()) ++iter; return 0;],
 ac_cv_cxx_have_stl=yes, ac_cv_cxx_have_stl=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_have_stl" = yes; then
  AC_DEFINE(HAVE_STL,,[define if the compiler supports Standard Template Library])
fi
])

dnl @synopsis AC_CXX_BOOL
dnl
dnl If the compiler recognizes bool as a separate built-in type,
dnl define HAVE_BOOL. Note that a typedef is not a separate
dnl type since you cannot overload a function such that it accepts either
dnl the basic type or the typedef.
dnl
dnl @version $Id$
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_BOOL],
[AC_CACHE_CHECK(whether the compiler recognizes bool as a built-in type,
ac_cv_cxx_bool,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
int f(int  x){return 1;}
int f(char x){return 1;}
int f(bool x){return 1;}
],[bool b = true; return f(b);],
 ac_cv_cxx_bool=yes, ac_cv_cxx_bool=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_bool" = yes; then
  AC_DEFINE(HAVE_BOOL,,[define if bool is a built-in type])
fi
])

dnl @synopsis AC_CXX_MUTABLE
dnl
dnl If the compiler allows modifying class data members flagged with
dnl the mutable keyword even in const objects (for example in the
dnl body of a const member function), define HAVE_MUTABLE.
dnl
dnl @version $Id$
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_MUTABLE],
[AC_CACHE_CHECK(whether the compiler supports the mutable keyword,
ac_cv_cxx_mutable,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
class A { mutable int i;
          public:
          int f (int n) const { i = n; return i; }
        };
],[A a; return a.f (1);],
 ac_cv_cxx_mutable=yes, ac_cv_cxx_mutable=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_mutable" = yes; then
  AC_DEFINE(HAVE_MUTABLE,,[define if the compiler supports the mutable keyword])
fi
])

dnl @synopsis AC_CXX_NEW_FOR_SCOPING
dnl
dnl If the compiler accepts the new for scoping rules (the scope of a
dnl variable declared inside the parentheses is restricted to the
dnl for-body), define HAVE_NEW_FOR_SCOPING.
dnl
dnl @version $Id$
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_NEW_FOR_SCOPING],
[AC_CACHE_CHECK(whether the compiler accepts the new for scoping rules,
ac_cv_cxx_new_for_scoping,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE(,[
  int z = 0;
  for (int i = 0; i < 10; ++i)
    z = z + i;
  for (int i = 0; i < 10; ++i)
    z = z - i;
  return z;],
 ac_cv_cxx_new_for_scoping=yes, ac_cv_cxx_new_for_scoping=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_new_for_scoping" = yes; then
  AC_DEFINE(HAVE_NEW_FOR_SCOPING,,[define if the compiler accepts the new for scoping rules])
fi
])

dnl @synopsis AC_CXX_STD_SPRINTF
dnl
dnl If the compiler recognizes std::sprintf as a function for IO,
dnl define HAVE_STD_SPRINTF.  If this test fails, use sprintf with no std prefix
dnl Note that we try to compile two versions of this routine, one using cstdio and
dnl another using stdio.h.  This approach is used to eliminate the need to test which
dnl of the two header files is present.  If one or both is usable the test will return true.
dnl

AC_DEFUN([AC_CXX_STD_SPRINTF],
[AC_CACHE_CHECK([[whether the compiler recognizes std::sprintf as supported IO function]],
ac_cv_cxx_std_sprintf,
[ AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
#include <cstdio>
#include <string>
],
[
     int x = 100;
     char s[100];
     std::sprintf(s, "%d", x);
],
 ac_cv_cxx_std_sprintf1=yes, ac_cv_cxx_std_sprintf1=no)

AC_TRY_COMPILE([
#include <stdio.h>
#include <string>
],
[
     int x = 100;
     char s[100];
     std::sprintf(s, "%d", x);
],
 ac_cv_cxx_std_sprintf2=yes, ac_cv_cxx_std_sprintf2=no)

if (test "$ac_cv_cxx_std_sprintf1" = yes || test "$ac_cv_cxx_std_sprintf2" = yes); then
 ac_cv_cxx_std_sprintf=yes
else
 ac_cv_cxx_std_sprintf=no
fi
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_std_sprintf" = yes; then
  AC_DEFINE(HAVE_STD_SPRINTF,,[define if std::sprintf is supported])
fi
])

dnl @synopsis TAC_ARG_WITH_PACKAGE(FEATURE_NAME, FEATURE_DESCRIPTION, HAVE_NAME, DEFAULT_VAL)
dnl
dnl Test for --with-${FEATURE_NAME} and set to DEFAULT_VAL value if feature not specified.
dnl Also calls AC_DEFINE to define HAVE_${HAVE_NAME} if value is not equal to "no"
dnl 
dnl Use this macro to help defining whether or not interfaces for optional 
dnl package should compiled.  For example:
dnl
dnl TAC_ARG_WITH_PACKAGE(zoltan, [Enable Zoltan interface support], ZOLTAN, no)
dnl 
dnl will test for --with-zoltan when configure is run.  If it is defined 
dnl (and not set to "no") then HAVE_ZOLTAN will be defined, 
dnl Otherwise HAVE_ZOLTAN will not be defined.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_PACKAGE],
[
AC_ARG_WITH([$1],
AC_HELP_STRING([--with-$1],[$2 (default is [$4])]),
ac_cv_use_$1=$withval, ac_cv_use_$1=$4)

AC_MSG_CHECKING(whether to use [$1])

if test "X$ac_cv_use_$1" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([HAVE_$3],,[Define if want to build with $1 enabled])
else
  AC_MSG_RESULT(no)
fi
])


dnl @synopsis TAC_ARG_ENABLE_OPTION(FEATURE_NAME, FEATURE_DESCRIPTION, HAVE_NAME, DEFAULT_VAL)
dnl
dnl Test for --enable-${FEATURE_NAME} and set to DEFAULT_VAL value if feature not specified.
dnl Also calls AC_DEFINE to define HAVE_${HAVE_NAME} if value is not equal to "no"
dnl
dnl Use this macro to facilitate definition of options in a package.  For example:
dnl 
dnl TAC_ARG_ENABLE_OPTION(threads, [enable shared memory threads], THREADS, no)
dnl 
dnl will test for --enable-threads when configure is run.  If it is defined (and not set to "no")
dnl then HAVE_THREADS will be defined, Otherwise HAVE_THREADS will not be defined.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_OPTION],
[
AC_ARG_ENABLE([$1],
AC_HELP_STRING([--enable-$1],[$2 (default is [$4])]),
ac_cv_use_$1=$enableval, ac_cv_use_$1=$4)

AC_MSG_CHECKING(whether to use [$1])

if test "X$ac_cv_use_$1" != "Xno"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_$3],1,[Define if want to build with $1 enabled])
else
  AC_MSG_RESULT(no)
fi
])


# @synopsis TAC_ARG_ENABLE_CAN_USE_PACKAGE(PACKAGE_NAME, OPTIONAL_DEPENDENCY_NAME, AC_DEFINE, AM_CONDITIONAL, DEFAULT_BUILD, ALLOW_IMPLICIT_ENABLE, HELP_STRING, IMPLICIT_HELP_STRING)
# 

# Use this macro to facilitate definition of options related to other
# Trilinos packages that a package "can use" but does not depend on.

# This macro supports both explicit and implicit configure options.
# For example to build epetra support into nox (so that nox can call
# epetra) we can explicitly use the flag --enable-nox-epetra.  This
# requires the user to enable both packages and the support between
# them: --enable-nox --enable-epetra --enable-nox-epetra.	

# Some packages, to simplify build requirements, implicitly assume
# that epetra support in nox should be built if both --enable-nox and
# --enable-epetra are supplied.  Users can override this by using the
# explicit command --enable-nox-epetra.

# Usage:
#  TAC_ARG_ENABLE_CAN_USE_PACKAGE(PACKAGE_NAME, 
#                                 OPTIONAL_PACKAGE_DEPENDENCY_NAME, 
#                                 AC_DEFINE, 
#                                 AM_CONDITIONAL, 
#                                 ALLOW_IMPLICIT_ENABLE, 
#                                 HELP_STRING, 
#                                 IMPLICIT_HELP_STRING)
#
# Where:
# PACKAGE_NAME - Name of package currently being enabled.
# OPTIONAL_PACKAGE_DEPENDENCY_NAME - Name of optinal package support that 
#                                    will be built in.
# AC_DEFINE - Name of define that will be put in the 
#             <PackageName>_ConfigDef.h file.  Note that a "HAVE_" 
#             will be automatically prepended to the name.
# AM_CONDITIONAL - Variable that will be set in the makefiles by an 
#                  AM_CONDITIONAL call.
# ALLOW_IMPLICIT_ENABLE - This can be used to turn off implicit enable 
#                         support.  Takes a "yes/no" argument.
# HELP_STRING - Help string for configure option:
#               --enable-<PACKAGE_NAME>-<OPTIONAL_PACKAGE_DEPENDENCY_NAME>
# IMPLICIT_HELP_STRING - Help string for implicit configure option:
#               --enable-<OPTIONAL_PACKAGE_DEPENDENCY_NAME>
# 
# For example, to force explicit configuration of epetra support in nox:
# 
#  TAC_ARG_ENABLE_EXPLICIT_CAN_USE_PACKAGE(nox, epetra, NOX_EPETRA, 
#                                          NOX_USING_EPETRA, no, 
#                                          [Builds epetra support into nox.], 
#                                          [DOES NOTHING!])
# 
# To allow both implicit and explicit configuration of epetra support in nox:
#
# TAC_ARG_ENABLE_EXPLICIT_CAN_USE_PACKAGE(nox, epetra, NOX_EPETRA, 
# NOX_USING_EPETRA, no, yes, 
# [Builds epetra support in nox.], 
# [Builds epetra support in nox.  Can be overridden with --enable-nox-epetra.])
#
# Results of calling this file:
#  1. An AM_CONDITIONAL will be set for the AM_CONDITIONAL argument defined
#     above.
#  2. An AC_DEFINE of the form: HAVE_<AC_DEFINE> where AC_DEFINE is the 
#     argument defined above.
#
# @author Roger Pawlowski <rppawlo@sandia.gov>
# Based on original verison by Jim Willenbring.
#
AC_DEFUN([TAC_ARG_ENABLE_CAN_USE_PACKAGE],
[

dnl Check for implicit enabling of optional package  
AC_ARG_ENABLE([$2],
AC_HELP_STRING([--enable-$2],[$7]),
ac_cv_implicit_use_$2=$enableval, 
ac_cv_implicit_use_$2=no)

dnl If implicit enabling is used, set that as teh default
if test "X$5" != "Xno"; then
  ac_cv_$1_using_$2_default=$ac_cv_implicit_use_$2
else
  ac_cv_$1_using_$2_default=no
fi

AC_ARG_ENABLE([$1-$2],
AC_HELP_STRING([--enable-$1-$2],[$6]),
ac_cv_use_$2=$enableval, 
ac_cv_use_$2=$ac_cv_$1_using_$2_default)

AC_MSG_CHECKING(whether to build optional [$2] dependent code in [$1])

if test "X$ac_cv_use_$2" != "Xno"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_$3],1,[Define if want to build with $1 enabled])
else
  AC_MSG_RESULT(no)
fi

AM_CONDITIONAL($4, test "X$ac_cv_use_$2" = "Xyes")

])

dnl @synopsis TAC_ARG_ENABLE_FEATURE(FEATURE_NAME, FEATURE_DESCRIPTION, HAVE_NAME, DEFAULT_VAL)
dnl
dnl Test for --enable-${FEATURE_NAME} and set to DEFAULT_VAL value if feature not specified.
dnl Also calls AC_DEFINE to define HAVE_${HAVE_NAME} if value is not equal to "no"
dnl 
dnl Use this macro to help defining whether or not optional 
dnl features* should compiled.  For example:
dnl
dnl TAC_ARG_ENABLE_FEATURE(epetra, [Configure and build epetra], EPETRA, yes)
dnl 
dnl will test for --enable-epetra when configure is run.  If it is defined 
dnl and not set to "no" or not defined (default is "yes") then HAVE_EPETRA will
dnl be defined, if --enable-epetra is defined to be "no", HAVE_EPETRA will not
dnl be defined.
dnl
dnl *NOTE: epetra, aztecoo, komplex, ifpack, and other software found in
dnl subdirectories of Trilinos/packages are "packages" in their own right.
dnl However, these packages are also "features" of the larger package
dnl "Trilinos".  Therefore, when configuring from the Trilinos directory,
dnl it is appropriate to refer to these software packages as "features".
dnl
dnl This file was based on tac_arg_with_package.m4 by Mike Heroux
dnl @author James Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_FEATURE],
[
AC_ARG_ENABLE([$1],
AC_HELP_STRING([--enable-$1],[$2 (default is [$4])]),
ac_cv_use_$1=$enableval, ac_cv_use_$1=$4)

AC_MSG_CHECKING(whether to use [$1])

if test "X$ac_cv_use_$1" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([HAVE_$3],,[Define if want to build $1])
else
  AC_MSG_RESULT(no)
fi
])


dnl @synopsis TAC_ARG_ENABLE_FEATURE_SUB_CHECK(FEATURE_NAME, SUB_FEATURE_NAME, FEATURE_DESCRIPTION, HAVE_NAME)
dnl
dnl This hack gets around the fact that TAC_ARG_ENABLE_FEATURE does not support underscores
dnl in its feature names.  TAC_ARG_ENABLE_FEATURE_SUB_CHECK allows exactly one underscore.  Not great,
dnl but arguably better than supporting no underscores.
dnl
dnl TAC_ARG_ENABLE_FEATURE(feature-sub, [Configure and build feature-sub], FEATURE_SUB, yes) 
dnl   fails because tac_arg_enable_feature tests for ac_cv_use_feature-sub which gets 
dnl   rejected because the `-' is not allowed in variables.  (AC_ARG_ENABLE sets ac_cv_use_feature_sub
dnl   to avoid this problem.)  Use:
dnl 
dnl TAC_ARG_ENABLE_FEATURE_SUB_CHECK(feature, sub, [Configure and build feature-sub], FEATURE_SUB) 
dnl   instead.
dnl
dnl This macro will test for --enable-${FEATURE_NAME}-${SUB_FEATURE_NAME} when configure is run.  
dnl If it is defined and not set to "no" or not defined and --disable-${SUB_FEATURE_NAME} is not
dnl specified then HAVE_${HAVE_NAME} will be defined.
dnl
dnl *NOTE: This macro is designed for the use-case when there is an individual Trilinos package 
dnl offering fine-grained control of a Trilinos option.  This way, the individual package 
dnl option is enabled, as long as the Trilinos option is not disabled.  If the Trilinos option is
dnl disabled, then the user must enable each packages option individually.  For instance:
dnl
dnl --disable-tests --enable-teuchos-tests
dnl
dnl *NOTE: epetra, aztecoo, komplex, ifpack, and other software found in
dnl subdirectories of Trilinos/packages are "packages" in their own right.
dnl However, these packages are also "features" of the larger package
dnl "Trilinos".  Therefore, when configuring from the Trilinos directory,
dnl it is appropriate to refer to these software packages as "features".
dnl
dnl This file was based on tac_arg_enable_package.m4 by Jim Willenbring
dnl and tac_arg_enable_package_sub.m4 by Ken Stanley.
dnl
dnl @author Heidi Thornquist <hkthorn@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_FEATURE_SUB_CHECK],
[
AC_ARG_ENABLE([$2],, ac_cv_use_$2=$enableval, ac_cv_use_$2=yes)

AC_ARG_ENABLE([$1-$2],
AC_HELP_STRING([--enable-$1-$2],[$3 (default is yes if --disable-$2 is not specified)]),
ac_cv_use_$1_$2=$enableval, ac_cv_use_$1_$2=${ac_cv_use_$2})

AC_MSG_CHECKING(whether to use [$1-$2])

if test "X$ac_cv_use_$1_$2" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([HAVE_$4],,[Define if want to build $1-$2])
else
  AC_MSG_RESULT(no)
fi
])


dnl @synopsis TAC_ARG_ENABLE_FEATURE_SUB(FEATURE_NAME, SUB_FEATURE_NAME, FEATURE_DESCRIPTION, HAVE_NAME, DEFAULT_VAL)
dnl
dnl This hack gets around the fact that TAC_ARG_ENABLE_FEATURE does not support underscores
dnl in its feature names.  TAC_ARG_ENABLE_FEATURE_SUB allows exactly one underscore.  Not great,
dnl but arguably better than supporting no underscores.
dnl
dnl TAC_ARG_ENABLE_FEATURE(feature-sub, [Configure and build feature-sub], FEATURE_SUB, yes) 
dnl   fails because tac_arg_enable_feature tests for ac_cv_use_feature-sub which gets 
dnl   rejected because the `-' is not allowed in variables.  (AC_ARG_ENABLE sets ac_cv_use_feature_sub
dnl   to avoid this problem.)  Use:
dnl 
dnl TAC_ARG_ENABLE_FEATURE(feature, sub, [Configure and build feature-sub], FEATURE_SUB, yes) 
dnl   instead.
dnl
dnl Test for --enable-${FEATURE_NAME} and set to DEFAULT_VAL value if feature not specified.
dnl Also calls AC_DEFINE to define HAVE_${HAVE_NAME} if value is not equal to "no"
dnl 
dnl Use this macro to help defining whether or not optional 
dnl features* should compiled.  For example:
dnl
dnl TAC_ARG_ENABLE_FEATURE(epetra, [Configure and build epetra], EPETRA, yes)
dnl 
dnl will test for --enable-epetra when configure is run.  If it is defined 
dnl and not set to "no" or not defined (default is "yes") then HAVE_EPETRA will
dnl be defined, if --enable-epetra is defined to be "no", HAVE_EPETRA will not
dnl be defined.
dnl
dnl *NOTE: epetra, aztecoo, komplex, ifpack, and other software found in
dnl subdirectories of Trilinos/packages are "packages" in their own right.
dnl However, these packages are also "features" of the larger package
dnl "Trilinos".  Therefore, when configuring from the Trilinos directory,
dnl it is appropriate to refer to these software packages as "features".
dnl
dnl This file was based on tac_arg_enable_package.m4 by Jim Willenbring
dnl @author Ken Stanley <ksstanl@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_FEATURE_SUB],
[
AC_ARG_ENABLE([$1-$2],
AC_HELP_STRING([--enable-$1-$2],[$3 (default is [$5])]),
ac_cv_use_$1_$2=$enableval, ac_cv_use_$1_$2=$5)

AC_MSG_CHECKING(whether to use [$1-$2])

if test "X$ac_cv_use_$1_$2" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([HAVE_$4],,[Define if want to build $1-$2])
else
  AC_MSG_RESULT(no)
fi
])


AC_DEFUN([TAC_ENABLE_PYTHON],
[
AZ_PYTHON_DEFAULT( )  # Default to no python support
AZ_PYTHON_ENABLE( )   # Check for --enable-python, --disable-python
AZ_PYTHON_WITH( )     # Check for --with-python, --without-python

AC_MSG_CHECKING(whether we should build the python wrappers)
if test -n "$PYTHON"; then
  BUILD_PYTHON=yes
  AC_MSG_RESULT(yes)
  AM_CONDITIONAL(BUILD_PYTHON, true)

  # Ensure that we have python version 2.3 or greater (for distutils)
  # We need 2.3 and not 2.2 because of Boolean types (for parameter lists)
  AZ_PYTHON_VERSION_ENSURE( [2.3] )

  # Python compiler and linker flags
  AZ_PYTHON_CSPEC
  AZ_PYTHON_LSPEC

  # Check that Python.h is available
  save_CPPFLAGS=$CPPFLAGS
  CPPFLAGS="$save_CPPFLAGS $PYTHON_CSPEC"
  AC_LANG([C++])
  AC_CHECK_HEADER(
  [Python.h],
  break,
  AC_MSG_ERROR([You must have Python.h in order to build the Python support!!]))
  CPPFLAGS="$save_CPPFLAGS"

  # Check for Numeric
  AC_PYTHON_MODULE(Numeric,yes)

  # If user specifies prefix, use it for the PYTHON_PREFIX
  if test "$prefix" != "$ac_default_prefix"; then
    PYTHON_PREFIX=$prefix
  fi
  if test "$exec_prefix" != "$ac_default_prefix"; then
    PYTHON_EXECPREFIX=$exec_prefix
  fi

else
  AC_MSG_RESULT(no)
  BUILD_PYTHON=no
  AM_CONDITIONAL(BUILD_PYTHON, false)
fi

# ------------------------------------------------------------------------
# If the python wrappers are to be built, then SWIG (Simple Wrapper
# Interface Generator) is required
# ------------------------------------------------------------------------

if test -n "$PYTHON"; then

  # Check for --with-swig[=path]
  AC_MSG_CHECKING(for --with-swig)
  AC_ARG_WITH(swig,
              [AC_HELP_STRING([--with-swig@<:@=SWIG@:>@],
                              [enable swig and set swig binary])],
              [AC_MSG_RESULT(yes)
               WITH_SWIG=yes
               if test X${withval} != Xyes ; then
                 SWIG=$withval
               fi],
              [AC_MSG_RESULT(no)
               AC_CHECK_PROG(WITH_SWIG,swig,yes,no)])

  # Report error if no swig found
  if test ${WITH_SWIG} = no; then
     AC_MSG_ERROR(
     [Python wrappers require swig (Simple Wrapper Interface Generator).
      See http://www.swig.org])
  fi

  # SWIG configuration
  AC_PROG_SWIG(1.3.23)
  SWIG_ENABLE_CXX
  SWIG_MULTI_MODULE_SUPPORT
  SWIG_PYTHON
fi
AM_CONDITIONAL(HAVE_SWIG,test X${WITH_SWIG} = Xyes)
AC_SUBST(GNU_HAVE_SWIG, ${WITH_SWIG})

])

dnl @synopsis AZ_PYTHON_DEFAULT
dnl @synopsis AZ_PYTHON_ENABLE
dnl @synopsis AZ_PYTHON_WITH
dnl @synopsis AZ_PYTHON_PATH
dnl @synopsis AZ_PYTHON_VERSION_ENSURE( [2.2] )
dnl @synopsis AZ_PYTHON_CSPEC
dnl @synopsis AZ_PYTHON_LSPEC
dnl
dnl @summary New and revised Python support.
dnl
dnl This file provides autoconf support for those applications that
dnl want to embed python. It supports all pythons >= 2.2 which is the
dnl first official release containing distutils. Version 2.2 of python
dnl was released December 21, 2001. Since it actually executes the
dnl python, cross platform configuration will probably not work. Also,
dnl most of the platforms supported are consistent until you look into
dnl MacOSX. The python included with it is installed as a framework
dnl which is a very different environment to set up the normal tools
dnl such as gcc and libtool to deal with. Therefore, once we establish
dnl which python that we are going to use, we use its distutils to
dnl actually compile and link our modules or applications.
dnl
dnl At this time, it does NOT support linking with Python statically.
dnl It does support dynamic linking.
dnl
dnl This set of macros help define $PYTHON, $PYTHON_USE, $PYTHON_CSPEC
dnl and $PYTHON_LSPEC. $PYTHON defines the full executable path for the
dnl Python being linked to and is used within these macros to determine
dnl if that has been specified or found. These macros do execute this
dnl python version so it must be present on the system at configure
dnl time.
dnl
dnl $PYTHON_USE is an automake variable that defines whether Python
dnl support should be included or not in your application.
dnl $PYTHON_CSPEC is a variable that supplies additional CFLAGS for the
dnl compilation of the application/shared library. $PYTHON_LSPEC is a
dnl variable that supplies additional LDFLAGS for linking the
dnl application/shared library.
dnl
dnl The following is an example of how to set up for python usage
dnl within your application in your configure.in:
dnl
dnl   AZ_PYTHON_DEFAULT( )
dnl   AZ_PYTHON_ENABLE( )             # Optional
dnl   AZ_PYTHON_WITH( )               # Optional
dnl   AZ_PYTHON_PATH( )               # or AZ_PYTHON_INSIST( )
dnl   # if $PYTHON is not defined, then the following do nothing.
dnl   AZ_PYTHON_VERSION_ENSURE( [2.2] )
dnl   AZ_PYTHON_CSPEC
dnl   AZ_PYTHON_LSPEC
dnl
dnl The AZ_PYTHON_DEFAULT sets the $PYTHON_USE to false. Thereby,
dnl excluding it if it was optional.
dnl
dnl The AZ_PYTHON_ENABLE looks for the optional configure parameters of
dnl --enable-python/--disable-python and establishes the $PYTHON and
dnl $PYTHON_USE variables accordingly.
dnl
dnl The AZ_PYTHON_WITH looks for the optional configure parameters of
dnl --with-python/--without-python and establishes the $PYTHON and
dnl $PYTHON_USE variables accordingly.
dnl
dnl The AZ_PYTHON_PATH looks for python assuming that none has been
dnl previously found or defined and issues an error if it does not find
dnl it. If it does find it, it establishes the $PYTHON and $PYTHON_USE
dnl variables accordingly. AZ_PYTHON_INSIST could be used here instead
dnl if you want to insist that Python support be included using the
dnl --enable-python or --with-python checks previously done.
dnl
dnl The AZ_PYTHON_VERSION_ENSURE issues an error if the Python
dnl previously found is not of version 2.2 or greater.
dnl
dnl Once that these macros have be run, we can use PYTHON_USE within
dnl the makefile.am file to conditionally add the Python support such
dnl as:
dnl
dnl Makefile.am example showing optional inclusion of directories:
dnl
dnl  if PYTHON_USE
dnl  plugins = plugins
dnl  src = src
dnl  else
dnl  plugins =
dnl  src =
dnl  endif
dnl
dnl  SUBDIRS = . $(plugins) $(src)
dnl
dnl Makefile.am example showing optional shared library build:
dnl
dnl  if PYTHON_USE
dnl  lib_LTLIBRARIES        = libElemList.la
dnl  libElemList_la_SOURCES = libElemList.c
dnl  libElemList_la_CFLAGS  = @PYTHON_CSPEC@
dnl  libElemList_la_LDFLAGS = @PYTHON_LSPEC@
dnl  endif
dnl
dnl Makefile.am example showing optional program build:
dnl
dnl  if PYTHON_USE
dnl  bin_PROGRAMS    = runFunc
dnl  runFunc_SOURCES = runFunc.c
dnl  runFunc_CFLAGS  = @PYTHON_CSPEC@
dnl  runFunc_LDFLAGS = @PYTHON_LSPEC@
dnl  endif
dnl
dnl The above compiles the modules only if PYTHON_USE was specified as
dnl true. Also, the else portion of the if was optional.
dnl
dnl @category InstalledPackages
dnl @author Robert White <kranki@mac.com>
dnl @author Dustin Mitchell <dustin@cs.uchicago.edu>
dnl @version 2005-01-14
dnl @license GPLWithACException

# AZ_PYTHON_DEFAULT( )
# -----------------
# Sets the default to not include Python support.

AC_DEFUN([AZ_PYTHON_DEFAULT],
[
    az_python_use=false
    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
])



# AZ_PYTHON_ENABLE( [path] )
# -----------------------------------------------------------------
# Handles the various --enable-python commands.
# Input:
#   $1 is the optional search path for the python executable if needed
# Ouput:
#   PYTHON_USE (AM_CONDITIONAL) is true if python executable found
#   and --enable-python was requested; otherwise false.
#   $PYTHON contains the full executable path to python if PYTHON_ENABLE_USE
#   is true.
#
# Example:
#   AZ_PYTHON_ENABLE( )
#   or
#   AZ_PYTHON_ENABLE( "/usr/bin" )

AC_DEFUN([AZ_PYTHON_ENABLE],
[
    AC_ARG_VAR([PYTHON],[Python Executable Path])

    # unless PYTHON was supplied to us (as a precious variable),
    # see if --enable-python[=PythonExecutablePath], --enable-python,
    # --disable-python or --enable-python=no was given.
    if test -z "$PYTHON"
    then
        AC_MSG_CHECKING(for --enable-python)
        AC_ARG_ENABLE(
            python,
            AC_HELP_STRING([--enable-python@<:@=PYTHON@:>@],
                [absolute path name of Python executable]
            ),
            [
                if test "$enableval" = "yes"
                then
                    # "yes" was specified, but we don't have a path
                    # for the executable.
                    # So, let's searth the PATH Environment Variable.
                    AC_MSG_RESULT(yes)
                    AC_PATH_PROG(
                        [PYTHON],
                        python,
                        [],
                        $1
                    )
                    if test -z "$PYTHON"
                    then
                        AC_MSG_ERROR(no path to python found)
                    fi
                    az_python_use=true
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                    AZ_PYTHON_PREFIX( )
                elif test "$enableval" = "no"
                then
                    AC_MSG_RESULT(no)
                    az_python_use=false
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                else
                    # $enableval must be the executable path then.
                    AC_SUBST([PYTHON], ["${enableval}"])
                    AC_MSG_RESULT($withval)
                    az_python_use=true
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                    AZ_PYTHON_PREFIX( )
                fi
            ],
            [
                # --with-python was not specified.
                AC_MSG_RESULT(no)
                az_python_use=false
                AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
            ]
        )
    fi

])



# AZ_PYTHON_CSPEC( )
# -----------------
# Set up the c compiler options to compile Python
# embedded programs/libraries in $PYTHON_CSPEC if
# $PYTHON has been defined.

AC_DEFUN([AZ_PYTHON_CSPEC],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -n "$PYTHON"
    then
        az_python_prefix=`${PYTHON} -c "import sys; print sys.prefix"`
        if test -z "$az_python_prefix"
        then
            AC_MSG_ERROR([Python Prefix is not known])
        fi
        az_python_execprefix=`${PYTHON} -c "import sys; print sys.exec_prefix"`
        az_python_version=`$PYTHON -c "import sys; print sys.version[[:3]]"`
        az_python_includespec="-I${az_python_prefix}/include/python${az_python_version}"
        if test x"$python_prefix" != x"$python_execprefix"; then
            az_python_execspec="-I${az_python_execprefix}/include/python${az_python_version}"
            az_python_includespec="${az_python_includespec} $az_python_execspec"
        fi
        az_python_ccshared=`${PYTHON} -c "import distutils.sysconfig; print distutils.sysconfig.get_config_var('CFLAGSFORSHARED')"`
        az_python_cspec="${az_python_ccshared} ${az_python_includespec}"
        AC_SUBST([PYTHON_CSPEC], [${az_python_cspec}])
        AC_MSG_NOTICE([PYTHON_CSPEC=${az_python_cspec}])
    fi
])



# AZ_PYTHON_INSIST( )
# -----------------
# Look for Python and set the output variable 'PYTHON'
# to 'python' if found, empty otherwise.

AC_DEFUN([AZ_PYTHON_PATH],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable not found])
    fi
])



# AZ_PYTHON_LSPEC( )
# -----------------
# Set up the linker options to link Python embedded
# programs/libraries in $PYTHON_LSPEC if $PYTHON
# has been defined.

AC_DEFUN([AZ_PYTHON_LSPEC],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -n "$PYTHON"
    then
        AZ_PYTHON_RUN([
import sys
import distutils.sysconfig
strUseFrameWork = "--enable-framework"
dictConfig = distutils.sysconfig.get_config_vars( )
strConfigArgs = dictConfig.get("CONFIG_ARGS")
strLinkSpec =  dictConfig.get('LDFLAGS')
if -1 ==  strConfigArgs.find(strUseFrameWork):
    strLibPL = dictConfig.get("LIBPL")
    if strLibPL and (strLibPL != ""):
        strLinkSpec += " -L%s" % (strLibPL)
    strSys = dictConfig.get("SYSLIBS")
    if strSys and (strSys != ""):
        strLinkSpec += " %s" % (strSys)
    strSHL = dictConfig.get("SHLIBS")
    if strSHL and (strSHL != ""):
        strLinkSpec += " %s" % (strSHL)
    # Construct the Python Library Name.
    strTmplte = " -lpython%d.%d"
    if (sys.platform == "win32") or (sys.platform == "os2emx"):
        strTmplte = " -lpython%d%d"
    strWrk = strTmplte % ( (sys.hexversion >> 24),
                            ((sys.hexversion >> 16) & 0xff))
    strLinkSpec += strWrk
else:
    # This is not ideal since it changes the search path
    # for Frameworks which could have side-effects on
    # other included Frameworks.  However, it is necessary
    # where someone has installed more than one frameworked
    # Python.  Frameworks are really only used in MacOSX.
    strLibFW = dictConfig.get("PYTHONFRAMEWORKPREFIX")
    if strLibFW and (strLibFW != ""):
        strLinkSpec += " -F%s" % (strLibFW)
strLinkSpec += " %s" % (dictConfig.get('LINKFORSHARED'))
print strLinkSpec
        ])
        AC_SUBST([PYTHON_LSPEC], [${az_python_output}])
        AC_MSG_NOTICE([PYTHON_LSPEC=${az_python_output}])
    fi
])



# AZ_PYTHON_PATH( )
# -----------------
# Look for Python and set the output variable 'PYTHON'
# to 'python' if found, empty otherwise.

AC_DEFUN([AZ_PYTHON_PATH],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    AC_PATH_PROG( PYTHON, python, [], $1 )
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable not found])
    else
        az_python_use=true
    fi
    AM_CONDITIONAL(PYTHON_USE, test "$az_python_use" = "true")
])



# AZ_PYTHON_PREFIX( )
# -------------------
# Use the values of $prefix and $exec_prefix for the corresponding
# values of PYTHON_PREFIX and PYTHON_EXEC_PREFIX.

AC_DEFUN([AZ_PYTHON_PREFIX],
[
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable Path is not known])
    fi
    ax_python_prefix=`${PYTHON} -c "import sys; print sys.prefix"`
    ax_python_execprefix=`${PYTHON} -c "import sys; print sys.exec_prefix"`
    AC_SUBST([PYTHON_PREFIX], ["${ax_python_prefix}"])
    AC_SUBST([PYTHON_EXECPREFIX], ["${ax_python_execprefix}"])
])



# AZ_PYTHON_RUN( PYTHON_PROGRAM )
# -----------------
# Run a Python Test Program saving its output
# in az_python_output and its condition code
# in az_python_cc.

AC_DEFUN([AZ_PYTHON_RUN],
[
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -z "$PYTHON"
    then
        AC_MSG_ERROR([Python Executable not found])
    else
        cat >conftest.py <<_ACEOF
$1
_ACEOF
        az_python_output=`$PYTHON conftest.py`
        az_python_cc=$?
        rm conftest.py
        if test -f "conftest.pyc"
        then
            rm conftest.pyc
        fi
    fi
])



# AZ_PYTHON_VERSION_CHECK( VERSION, [ACTION-IF-TRUE], [ACTION-IF-FALSE] )
# -----------------------------------------------------------------------------
# Run ACTION-IF-TRUE if the Python interpreter has version >= VERSION.
# Run ACTION-IF-FALSE otherwise.
# This test uses sys.hexversion instead of the string equivalant (first
# word of sys.version), in order to cope with versions such as 2.2c1.
# hexversion has been introduced in Python 1.5.2; it's probably not
# worth to support older versions (1.5.1 was released on October 31, 1998).

AC_DEFUN([AZ_PYTHON_VERSION_CHECK],
 [
    AC_ARG_VAR( [PYTHON], [Python Executable Path] )
    if test -n "$PYTHON"
    then
        AC_MSG_CHECKING([whether $PYTHON version >= $1])
        AZ_PYTHON_RUN([
import sys, string
# split strings by '.' and convert to numeric.  Append some zeros
# because we need at least 4 digits for the hex conversion.
minver = map(int, string.split('$1', '.')) + [[0, 0, 0]]
minverhex = 0
for i in xrange(0, 4): minverhex = (minverhex << 8) + minver[[i]]
if sys.hexversion >= minverhex:
    sys.exit( 0 )
else:
    sys.exit( 1 )
        ])
        if test $az_python_cc -eq 0
        then
            $2
        m4_ifvaln(
            [$3],
            [else $3]
        )
        fi
    fi
])



# AZ_PYTHON_VERSION_ENSURE( VERSION )
# -----------------
# Insure that the Python Interpreter Version
# is greater than or equal to the VERSION
# parameter.

AC_DEFUN([AZ_PYTHON_VERSION_ENSURE],
[
    AZ_PYTHON_VERSION_CHECK(
        [$1],
        [AC_MSG_RESULT(yes)],
        [AC_MSG_ERROR(too old)]
    )
])



# AZ_PYTHON_WITH( [path] )
# -----------------------------------------------------------------
# Handles the various --with-python commands.
# Input:
#   $1 is the optional search path for the python executable if needed
# Ouput:
#   PYTHON_USE (AM_CONDITIONAL) is true if python executable found
#   and --with-python was requested; otherwise false.
#   $PYTHON contains the full executable path to python if PYTHON_USE
#   is true.
#
# Example:
#   AZ_PYTHON_WITH( )
#   or
#   AZ_PYTHON_WITH("/usr/bin")

AC_DEFUN([AZ_PYTHON_WITH],
[
    AC_ARG_VAR([PYTHON],[Python Executable Path])

    # unless PYTHON was supplied to us (as a precious variable),
    # see if --with-python[=PythonExecutablePath], --with-python,
    # --without-python or --with-python=no was given.
    if test -z "$PYTHON"
    then
        AC_MSG_CHECKING(for --with-python)
        AC_ARG_WITH(
            python,
            AC_HELP_STRING([--with-python@<:@=PYTHON@:>@],
                [absolute path name of Python executable]
            ),
            [
                if test "$withval" = "yes"
                then
                    # "yes" was specified, but we don't have a path
                    # for the executable.
                    # So, let's searth the PATH Environment Variable.
                    AC_MSG_RESULT(yes)
                    AC_PATH_PROG(
                        [PYTHON],
                        python,
                        [],
                        $1
                    )
                    if test -z "$PYTHON"
                    then
                        AC_MSG_ERROR(no path to python found)
                    fi
                    az_python_use=true
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                    AZ_PYTHON_PREFIX( )
                elif test "$withval" = "no"
                then
                    AC_MSG_RESULT(no)
                    az_python_use=false
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                else
                    # $withval must be the executable path then.
                    AC_SUBST([PYTHON], ["${withval}"])
                    AC_MSG_RESULT($withval)
                    az_python_use=true
                    AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
                    AZ_PYTHON_PREFIX( )
                fi
            ],
            [
                # --with-python was not specified.
                AC_MSG_RESULT(no)
                az_python_use=false
                AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")
            ]
        )
    fi

])

dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_python_module.html
dnl
AC_DEFUN([AC_PYTHON_MODULE],[
	AC_MSG_CHECKING(python module: $1)
	python -c "import $1" 2>/dev/null
	if test $? -eq 0;
	then
		AC_MSG_RESULT(yes)
		eval AS_TR_CPP(HAVE_PYMOD_$1)=yes
	else
		AC_MSG_RESULT(no)
		eval AS_TR_CPP(HAVE_PYMOD_$1)=no
		#
		if test -n "$2"
		then
			AC_MSG_ERROR(failed to find required module $1)
			exit 1
		fi
	fi
])

dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_pkg_swig.html
dnl
AC_DEFUN([AC_PROG_SWIG],[
	if test -z "$SWIG"; then
        	AC_PATH_PROG([SWIG],[swig])
	fi
        if test -z "$SWIG" ; then
                AC_MSG_WARN([cannot find 'swig' program. You should look at http://www.swig.org])
                SWIG='echo "Error: SWIG is not installed. You should look at http://www.swig.org" ; false'
        else
	    AC_MSG_NOTICE([SWIG executable is '$SWIG'])
	    if test -n "$1" ; then
                AC_MSG_CHECKING([for SWIG version])
                [swig_version=`$SWIG -version 2>&1 | grep 'SWIG Version' | sed 's/.*\([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/g'`]
                AC_MSG_RESULT([$swig_version])
                if test -n "$swig_version" ; then
                        # Calculate the required version number components
                        [required=$1]
                        [required_major=`echo $required | sed 's/[^0-9].*//'`]
                        if test -z "$required_major" ; then
                                [required_major=0]
                        fi
                        [required=`echo $required | sed 's/[0-9]*[^0-9]//'`]
                        [required_minor=`echo $required | sed 's/[^0-9].*//'`]
                        if test -z "$required_minor" ; then
                                [required_minor=0]
                        fi
                        [required=`echo $required | sed 's/[0-9]*[^0-9]//'`]
                        [required_patch=`echo $required | sed 's/[^0-9].*//'`]
                        if test -z "$required_patch" ; then
                                [required_patch=0]
                        fi
                        # Calculate the available version number components
                        [available=$swig_version]
                        [available_major=`echo $available | sed 's/[^0-9].*//'`]
                        if test -z "$available_major" ; then
                                [available_major=0]
                        fi
                        [available=`echo $available | sed 's/[0-9]*[^0-9]//'`]
                        [available_minor=`echo $available | sed 's/[^0-9].*//'`]
                        if test -z "$available_minor" ; then
                                [available_minor=0]
                        fi
                        [available=`echo $available | sed 's/[0-9]*[^0-9]//'`]
                        [available_patch=`echo $available | sed 's/[^0-9].*//'`]
                        if test -z "$available_patch" ; then
                                [available_patch=0]
                        fi
                        if test $available_major -ne $required_major \
                                -o $available_minor -ne $required_minor \
                                -o $available_patch -lt $required_patch ; then
                                AC_MSG_WARN([SWIG version >= $1 is required.  You have $swig_version.  You should look at http://www.swig.org])
                                SWIG='echo "Error: SWIG version >= $1 is required.  You have '"$swig_version"'.  You should look at http://www.swig.org" ; false'
                        else
                                SWIG_LIB=`$SWIG -swiglib`
                                AC_MSG_NOTICE([SWIG library directory is '$SWIG_LIB'])
                        fi
                else
                        AC_MSG_WARN([cannot determine SWIG version])
                        SWIG='echo "Error: Cannot determine SWIG version.  You should look at http://www.swig.org" ; false'
                fi
            fi
        fi
        AC_SUBST([SWIG_LIB])
])

# SWIG_ENABLE_CXX()
#
# Enable SWIG C++ support.  This affects all invocations of $(SWIG).
AC_DEFUN([SWIG_ENABLE_CXX],[
        AC_REQUIRE([AC_PROG_SWIG])
        AC_REQUIRE([AC_PROG_CXX])
        SWIG="$SWIG -c++"
])

# SWIG_MULTI_MODULE_SUPPORT()
#
# Enable support for multiple modules.  This effects all invocations
# of $(SWIG).  You have to link all generated modules against the
# appropriate SWIG runtime library.  If you want to build Python
# modules for example, use the SWIG_PYTHON() macro and link the
# modules against $(SWIG_PYTHON_LIBS).
#
AC_DEFUN([SWIG_MULTI_MODULE_SUPPORT],[
        AC_REQUIRE([AC_PROG_SWIG])
#        SWIG="$SWIG -noruntime"
#        SWIG="$SWIG -c"
	SWIG="$SWIG"
])

# SWIG_PYTHON([use-shadow-classes = {no, yes}])
#
# Checks for Python and provides the $(SWIG_PYTHON_CPPFLAGS),
# and $(SWIG_PYTHON_OPT) output variables.
#
# $(SWIG_PYTHON_OPT) contains all necessary SWIG options to generate
# code for Python.  Shadow classes are enabled unless the value of the
# optional first argument is exactly 'no'.  If you need multi module
# support (provided by the SWIG_MULTI_MODULE_SUPPORT() macro) use
# $(SWIG_PYTHON_LIBS) to link against the appropriate library.  It
# contains the SWIG Python runtime library that is needed by the type
# check system for example.
AC_DEFUN([SWIG_PYTHON],[
        AC_REQUIRE([AC_PROG_SWIG])
#        AC_REQUIRE([AC_PYTHON_DEVEL])
#        test "x$1" != "xno" || swig_shadow=" -noproxy"
#        AC_SUBST([SWIG_PYTHON_OPT],[-python$swig_shadow])
        AC_SUBST([SWIG_PYTHON_OPT],[-python])
        AC_SUBST([SWIG_PYTHON_CPPFLAGS],[$PYTHON_CPPFLAGS])
])


dnl @synopsis AC_LIB_WAD
dnl
dnl This macro searches for installed WAD library.
dnl
AC_DEFUN([AC_LIB_WAD],
[
        AC_REQUIRE([AC_PYTHON_DEVEL])
        AC_ARG_ENABLE(wad,
        AC_HELP_STRING([--enable-wad], [enable wad module]),
        [
                case "${enableval}" in
                        no)     ;;
                        *)      if test "x${enableval}" = xyes;
                                then
                                        check_wad="yes"
                                fi ;;
                esac
        ], [])

        if test -n "$check_wad";
        then
                AC_CHECK_LIB(wadpy, _init, [WADPY=-lwadpy], [], $PYTHON_LDFLAGS $PYTHON_EXTRA_LIBS)
                AC_SUBST(WADPY)
        fi
])

dnl @synopsis TAC_ARG_WITH_LIBDIRS
dnl
dnl Test for --with-libdirs="-Llibdir1 -Llibdir2".  if defined, 
dnl prepend "-Llibdir1 -Llibdir2" to LDFLAGS
dnl
dnl Use this macro to facilitate addition of directories to library search path.
dnl 
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_LIBDIRS],
[
AC_MSG_CHECKING([whether additional library search paths defined])
AC_ARG_WITH(libdirs,
AC_HELP_STRING([--with-libdirs], 
[OBSOLETE use --with-ldflags instead.  (ex. --with-ldflags="-L<DIR> -L<DIR2>")]),
[
LDFLAGS="${withval} ${LDFLAGS}"
AC_MSG_RESULT([${withval}])
],
AC_MSG_RESULT(no)
)
])


dnl @synopsis TAC_ARG_WITH_INCDIRS
dnl
dnl Test for --with-incdirs="-Iincdir1 -Iincdir2".  if defined, prepend 
dnl "-Iincdir1 -Iincdir2" to CPPFLAGS
dnl
dnl Use this macro to facilitate addition of directories to include file search path.
dnl 
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_INCDIRS],
[
AC_MSG_CHECKING([whether additional include search paths defined])
AC_ARG_WITH(incdirs,
AC_HELP_STRING([--with-incdirs], 
[additional directories containing include files: will prepend to search here for includes, use -Idir format]),
[
CPPFLAGS="${withval} ${CPPFLAGS}"
AC_MSG_RESULT([${withval}])
],
AC_MSG_RESULT(no)
)
])


dnl @synopsis TAC_ARG_WITH_BLASLIB
dnl
dnl Test for --with-blaslib="name".
dnl 
dnl Prepends the specified name to the list of files to check for BLAS
dnl routines.  
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_BLASLIB],
[
AC_ARG_WITH(blaslib,
AC_HELP_STRING([--with-blaslib], 
[name of library containing BLAS: will search lib directories for
-lname]),
[
USE_BLASLIB=yes
NEWBLASLIB=${withval}
]
)

   BLASLIBS="cxml blas complib.sgimath"

if test "X${USE_BLASLIB}" = "Xyes"; then

   BLASLIBS="${NEWBLASLIB} ${BLASLIBS}"

fi
])


dnl @synopsis TAC_ARG_WITH_LAPACKLIB
dnl
dnl Test for --with-lapacklib="name".
dnl 
dnl Prepends the specified name to the list of files to check for LAPACK
dnl routines.  
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_LAPACKLIB],
[
AC_ARG_WITH(lapacklib,
AC_HELP_STRING([--with-lapacklib], 
[name of library containing LAPACK: will search lib directories for -lname]),
[
USE_LAPACKLIB=yes
NEWLAPACKLIB=${withval}
]
)

   LAPACKLIBS="cxml lapack complib.sgimath"

if test "X${USE_LAPACKLIB}" = "Xyes"; then

   LAPACKLIBS="${NEWLAPACKLIB} ${LAPACKLIBS}"
fi
])


dnl @synopsis ACX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
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
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl Many libraries are searched for, from ATLAS to CXML to ESSL.
dnl The user may also use --with-blas=<lib> in order to use some
dnl specific BLAS library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the BLAS library.
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
AC_DEFUN([ACX_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_blas_ok=no

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_F77_FUNC(sgemm)
AC_F77_FUNC(dgemm)

acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $dgemm in $BLAS_LIBS])
	AC_TRY_LINK_FUNC($dgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($acx_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_FUNC($dgemm, [acx_blas_ok=yes])
	LIBS="$save_LIBS"
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $dgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[acx_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $dgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(cxml, $dgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(dxml, $dgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $dgemm,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 acx_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(scs, $dgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $dgemm,
		     [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $dgemm,
		[AC_CHECK_LIB(essl, $dgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $dgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        acx_blas_ok=no
        $2
fi
])dnl ACX_BLAS

dnl @synopsis ACX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the LAPACK
dnl linear-algebra interface (see http://www.netlib.org/lapack/).
dnl On success, it sets the LAPACK_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with LAPACK, you should link with:
dnl
dnl     $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-lapack=<lib> in order to use some
dnl specific LAPACK library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the LAPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a LAPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_LAPACK.
dnl
dnl @version $Id$
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl edited by Jim Willenbring <jmwille@sandia.gov> to check for sgecon
dnl rather than cheev because by default (as of 8-13-2002) Trilinos
dnl does not build the complex portions of the lapack library.  Edited
dnl again on 5-13-2004 to check for dgecon instead of sgecon.

AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])
acx_lapack_ok=no

AC_ARG_WITH(lapack,
        [AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
        yes | "") ;;
        no) acx_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
AC_F77_FUNC(dgecon)

# We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
        acx_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $dgecon in $LAPACK_LIBS])
        AC_TRY_LINK_FUNC($dgecon, [acx_lapack_ok=yes], [LAPACK_LIBS=""])
        AC_MSG_RESULT($acx_lapack_ok)
        LIBS="$save_LIBS"
        if test acx_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_CHECK_FUNC($dgecon, [acx_lapack_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test $acx_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
                AC_CHECK_LIB($lapack, $dgecon,
                    [acx_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        acx_lapack_ok=no
        $2
fi
])dnl ACX_LAPACK

