#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_user_aggr \
                  ml_filtering \
                  ml_analyze 

# source file definitions, link lines
include ../Makefile-common.include

#make a symbolic link to these files in the build area
.linktosources:
	@if (! [ -h ml_user_aggr.cpp ] ); then ln -s $(ml_user_aggr_SOURCES) .;fi
	@if (! [ -h ml_filtering.cpp ] ); then ln -s $(ml_filtering_SOURCES) .;fi
	@if (! [ -h ml_analyze.cpp ] ); then ln -s $(ml_analyze_SOURCES) .; fi
	touch .linktosources
                                                                                
ml_user_aggr_DEPENDENCIES = \
    .linktosources
ml_filtering_DEPENDENCIES = \
    .linktosources
ml_analyze_DEPENDENCIES = \
    .linktosources
