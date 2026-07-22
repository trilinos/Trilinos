#
# Build a makefile stub that contains Trilinos autoconf-generated
# macros for compiling and building.  We will grab these from
# the tools package Teuchos.
#

TRILINOS_TEUCHOS_BUILD = $(TRILINOS_BUILD_DIR)/packages/teuchos
TRILINOS_MAKE_OPTIONS_FILE = ./trilinos_make_options.mak
TRILINOS_TEUCHOS_MAKEFILE = $(TRILINOS_TEUCHOS_BUILD)/src/Makefile
TRILINOS_TEUCHOS_EXPORT_MAKEFILE = $(TRILINOS_TEUCHOS_BUILD)/Makefile.export.teuchos

# Make rule for this file
$(TRILINOS_MAKE_OPTIONS_FILE) : $(TRILINOS_TEUCHOS_MAKEFILE)
	$(TRILINOS_SRC_DIR)/packages/teuchos/config/generate-makeoptions.pl $(TRILINOS_TEUCHOS_BUILD)/src/Makefile TEUCHOS > $(TRILINOS_MAKE_OPTIONS_FILE)

#
# Read in the Trilinos autoconf-generated macros from the file created above
#

include $(TRILINOS_MAKE_OPTIONS_FILE)
include $(TRILINOS_TEUCHOS_EXPORT_MAKEFILE) # Note, the order is important since TEUCHO_LIBS is redefined here!

#
# File extenstions
#

EXTERNAL_LIB_EXT = a

EXTERNAL_OBJ_EXT = o

# Set makefile buildsystem macros
#

# EXTERNAL_C compiler
EXTERNAL_C = $(TEUCHOS_CC)
EXTERNAL_C_DEP_OPT     = -MM
EXTERNAL_C_COMPILE_OPT = -c
EXTERNAL_C_OUTPUT_OPT = -o \



# EXTERNAL_C++ compiler
EXTERNAL_CXX = $(TEUCHOS_CXX)
EXTERNAL_CXX_DEP_OPT   = -MM
EXTERNAL_CXX_COMPILE_OPT = -c
EXTERNAL_CXX_OUTPUT_OPT = -o \

# Fortran compiler
EXTERNAL_F77 = $(TEUCHOS_F77)
EXTERNAL_F77_COMPILE_OPT = -c
EXTERNAL_F77_OUTPUT_OPT = -o \

# Library creator
EXTERNAL_AR = $(TEUCHOS_libteuchos_a_AR) \

EXTERNAL_RANLIB = $(TEUCHOS_RANLIB) \

# Linker
EXTERNAL_LD = $(TEUCHOS_CXXLD)

# Install directory (taken from Trilinos' --prefix=??? option)
ifneq ($(TEUCHOS_prefix),)
	EXTERNAL_INSTALL_DIR = $(TEUCHOS_prefix)
endif

#
# Program options
#

# Preprocessor macro definitions
EXTERNAL_DEFINES += -D_MIN=min -D_MAX=max $(TEUCHOS_DEFS) $(TEUCHOS_CPPFLAGS)

# Include directories
EXTERNAL_INCLUDES +=  $(EXTERNAL_INCL_DIR)

#EXTERNAL_CPPFLAGS += -v

# Linker Options
EXTERNAL_LDFLAGS = $(TEUCHOS_LDFLAGS) $(TEUCHOS_LIBS)

# EXTERNAL_C, EXTERNAL_C++ and Fortran compiler options

EXTERNAL_CFLAGS   = $(TEUCHOS_CFLAGS)
EXTERNAL_CXXFLAGS = $(TEUCHOS_CXXFLAGS)
EXTERNAL_F77FLAGS = $(TEUCHOS_FFLAGS)

#
# Build Rules
#

# Build object files from EXTERNAL_C source files
%.$(EXTERNAL_OBJ_EXT) : %.c $(TRILINOS_MAKE_OPTIONS_FILE)
	$(EXTERNAL_C) $(EXTERNAL_C_COMPILE_OPT) $(EXTERNAL_CPPFLAGS) $(EXTERNAL_EXTRA_CPPFLAGS) \
	$(EXTERNAL_DEFINES) $(EXTERNAL_INCLUDES) \
	$(EXTERNAL_CFLAGS) $(EXTERNAL_EXTRA_CFLAGS) $(EXTERNAL_C_OUTPUT_OPT)$@ $<

# Build object files from EXTERNAL_C++ source files
%.$(EXTERNAL_OBJ_EXT) : %.cpp $(TRILINOS_MAKE_OPTIONS_FILE)
	$(EXTERNAL_CXX) $(EXTERNAL_CXX_COMPILE_OPT) $(EXTERNAL_CPPFLAGS) $(EXTERNAL_EXTRA_CPPFLAGS) \
	$(EXTERNAL_DEFINES) $(EXTERNAL_INCLUDES) \
	$(EXTERNAL_CXXFLAGS) $(EXTERNAL_EXTRA_CXXFLAGS) $(EXTERNAL_CXX_OUTPUT_OPT)$@ $<

# Build object files from Fotran source files
%.$(EXTERNAL_OBJ_EXT) : %.f $(TRILINOS_MAKE_OPTIONS_FILE)
	$(EXTERNAL_F77) $(EXTERNAL_F77_COMPILE_OPT) $(EXTERNAL_F77FLAGS) $(EXTERNAL_EXTRA_F77FLAGS) $(EXTERNAL_F77_OUTPUT_OPT)$@ $<
#	$(EXTERNAL_F77) $(EXTERNAL_F77_COMPILE_OPT) $(EXTERNAL_CPPFLAGS) $(EXTERNAL_EXTRA_CPPFLAGS) $(EXTERNAL_F77FLAGS) $(EXTERNAL_EXTRA_F77FLAGS) $(EXTERNAL_F77_OUTPUT_OPT)$@ $<

# Build dependency files for EXTERNAL_C source files that include header dependencies
%.d: %.c $(TRILINOS_MAKE_OPTIONS_FILE)
	$(EXTERNAL_C) $(EXTERNAL_C_DEP_OPT) $(EXTERNAL_CPPFLAGS) $(EXTERNAL_EXTRA_CPPFLAGS) \
	$(EXTERNAL_DEFINES) $(EXTERNAL_INCLUDES) \
	$< \
	| sed 's/$(@:.d=\.$(EXTERNAL_OBJ_EXT))/$(@:.d=.$(EXTERNAL_OBJ_EXT)) $@/' | $(EXTERNAL_DEP_POST_PROC) > $@; [ -s $@ ] || rm -f $@
#	| $(EXTERNAL_BASE_DIR)/Moocho/build/dep_post.pl $@ $(EXTERNAL_OBJ_EXT) | $(EXTERNAL_DEP_POST_PROC) > $@; [ -s $@ ] || rm -f $@

# Build dependency files for EXTERNAL_C++ source files that include header dependencies
%.d: %.cpp $(TRILINOS_MAKE_OPTIONS_FILE)
	$(EXTERNAL_CXX) $(EXTERNAL_CXX_DEP_OPT) $(EXTERNAL_CPPFLAGS) $(EXTERNAL_EXTRA_CPPFLAGS) \
	$(EXTERNAL_DEFINES) $(EXTERNAL_INCLUDES) \
	$< \
	| sed 's/$(@:.d=\.$(EXTERNAL_OBJ_EXT))/$(@:.d=.$(EXTERNAL_OBJ_EXT)) $@/' | $(EXTERNAL_DEP_POST_PROC) > $@; [ -s $@ ] || rm -f $@
#	| $(EXTERNAL_BASE_DIR)/Moocho/build/dep_post.pl $@ $(EXTERNAL_OBJ_EXT) | $(EXTERNAL_DEP_POST_PROC) > $@; [ -s $@ ] || rm -f $@

#
# Universal targets
#

clean-obj :
	rm *.o
clean-lib :
	rm *.a
clean-exe :
	rm *.exe
clean : clean-obj clean-lib clean-exe
