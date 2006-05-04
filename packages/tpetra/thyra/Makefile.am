
# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package tpetra subdirectory thyra built successfully."
	@echo ""

if BUILD_TESTS
#TEST_SUBDIR=test
endif

if BUILD_EXAMPLES 
EXAMPLE_SUBDIR=example
endif

SUBDIRS = $(EXAMPLE_SUBDIR) $(TEST_SUBDIR)
