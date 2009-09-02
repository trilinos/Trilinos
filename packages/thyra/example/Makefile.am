
# The following line helps the test harness recover from build errors.

all-local: generate-sample-output.pl

if USING_EXPORT_MAKEFILES
EXTERNAL_MAKEFILE_EXAMPLE = external
endif

SUBDIRS = operator_vector $(EXTERNAL_MAKEFILE_EXAMPLE)

generate-sample-output.pl : generate-sample-output.pl.stub
	cp generate-sample-output.pl.stub generate-sample-output.pl
	chmod a+x generate-sample-output.pl
