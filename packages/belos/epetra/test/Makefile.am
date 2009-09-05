# @HEADER
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:

if BUILD_EPETRA
EPETRA_TEST_SUBDIRS = BlockCG BlockGmres GCRODR PCPG RCG TFQMR MVOPTester
endif

SUBDIRS = $(EPETRA_TEST_SUBDIRS)

EXTRA_DIST = definition
