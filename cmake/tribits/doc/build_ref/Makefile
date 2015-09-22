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

all: $(BASE_NAME).rst $(BASE_NAME).html $(BASE_NAME).pdf

$(BASE_NAME).pdf: $(BASE_NAME).html 

$(BASE_NAME).html: $(BASE_NAME).rst

$(BASE_NAME).rst: $(EXTRA_FILES)
	./create-project-build-ref.py \
	--project-name="<Project>" \
	--project-template-file=TribitsBuildReferenceTemplate.rst \
	--file-base=TribitsBuildReference \
	$(EXTRA_ARGS)

clean:
	-rm -f $(BASE_NAME).rst
	-rm -f $(BASE_NAME).html
	-rm -f $(BASE_NAME).pdf

##---------------------------------------------------------------------------##
## end of Makefile
##---------------------------------------------------------------------------##
