
# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	echo "Trilinos package trilinoscouplings subdirectory examples built successfully."
	@echo ""


if BUILD_NOX_ML
ML_SUBDIR = ml/NonlinML
else
ML_SUBDIR = 
endif

SUBDIRS = epetraext $(ML_SUBDIR) 

EXTRA_DIST = 
