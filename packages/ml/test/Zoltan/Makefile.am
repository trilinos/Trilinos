# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

if ML_USING_ZOLTAN
if ML_USING_TEUCHOS
if ML_USING_AZTECOO
if ML_USING_EPETRA
if ML_USING_GALERI

noinst_PROGRAMS = Zoltan

Zoltan_SOURCES      = $(srcdir)/cxx_main.cpp
Zoltan_DEPENDENCIES = $(top_builddir)/src/libml.a
Zoltan_CXXFLAGS     = $(TEST_INCLUDES)
Zoltan_LDADD        = $(TEST_LIBS)

endif
endif
endif
endif
endif
