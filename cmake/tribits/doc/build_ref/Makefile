##---------------------------------------------------------------------------##
## MAKEFILE
##---------------------------------------------------------------------------##

EXTRA_ARGS=

.PHONY: clean realclean all

BASE_NAME = TribitsBuildReference
EXTRA_FILES = \
  Makefile \
  TribitsBuildReferenceBody.rst \
  TribitsBuildReferenceTemplate.rst

# NOTE: We *don't* want to rebuild the file if only the version has changed we
# don't add a dependency on TribitsGitVersion.txt!

CREATE_PROJECT_BUILD_REF_CMND = \
	./create-project-build-ref.py \
	--project-name="<Project>" \
	--project-template-file=TribitsBuildReferenceTemplate.rst \
	--file-base=TribitsBuildReference

all: $(EXTRA_FILES)
	$(CREATE_PROJECT_BUILD_REF_CMND) \
	$(EXTRA_ARGS)

$(BASE_NAME).rst: $(EXTRA_FILES)
	$(CREATE_PROJECT_BUILD_REF_CMND) \
	--generate-html= --generate-latex= --generate-pdf= \
	$(EXTRA_ARGS)

clean:
	-rm -f $(BASE_NAME).rst
	-rm -f $(BASE_NAME).html
	-rm -f $(BASE_NAME).pdf

##---------------------------------------------------------------------------##
## end of Makefile
##---------------------------------------------------------------------------##
