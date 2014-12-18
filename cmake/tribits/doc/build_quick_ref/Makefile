##---------------------------------------------------------------------------##
## MAKEFILE
##---------------------------------------------------------------------------##

EXTRA_ARGS=

.PHONY: clean realclean all

BASE_NAME = TribitsBuildQuickRef
EXTRA_FILES = \
  Makefile \
  TribitsBuildQuickRefBody.rst \
  TribitsBuildQuickRefTemplate.rst

# NOTE: We *don't* want to rebuild the file if only the version has changed we
# don't add a dependency on TribitsGitVersion.txt!

all: $(BASE_NAME).rst $(BASE_NAME).html $(BASE_NAME).pdf

$(BASE_NAME).pdf: $(BASE_NAME).html 

$(BASE_NAME).html: $(BASE_NAME).rst

$(BASE_NAME).rst: $(EXTRA_FILES)
	./create-project-build-quickref.py \
	--project-name="<Project>" \
	--project-template-file=TribitsBuildQuickRefTemplate.rst \
	--file-base=TribitsBuildQuickRef \
	$(EXTRA_ARGS)

clean:
	-rm -f $(BASE_NAME).rst
	-rm -f $(BASE_NAME).html
	-rm -f $(BASE_NAME).pdf

##---------------------------------------------------------------------------##
## end of Makefile
##---------------------------------------------------------------------------##
