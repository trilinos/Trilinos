##---------------------------------------------------------------------------##
## MAKEFILE
##---------------------------------------------------------------------------##

EXTRA_ARGS=

.PHONY: clean realclean all

TRIBITS_BUILD_REF_DIR=../../cmake/tribits/doc/build_ref

BASE_NAME = TrilinosBuildReference
EXTRA_FILES = \
  create-build-ref.sh \
  Makefile \
  $(TRIBITS_BUILD_REF_DIR)/TribitsBuildReferenceBody.rst \
  $(TRIBITS_BUILD_REF_DIR)/create-project-build-ref.py \
  TrilinosBuildReferenceTemplate.rst

all: $(BASE_NAME).rst $(BASE_NAME).html $(BASE_NAME).pdf

$(BASE_NAME).pdf: $(BASE_NAME).html 

$(BASE_NAME).html: $(BASE_NAME).rst

$(BASE_NAME).rst: $(EXTRA_FILES)
	$(TRIBITS_BUILD_REF_DIR)/create-project-build-ref.py \
        --project-base-dir=../.. \
	--project-template-file=$(PWD)/TrilinosBuildReferenceTemplate.rst \
	--file-base=TrilinosBuildReference \
	$(EXTRA_ARGS)

clean:
	-rm -f $(BASE_NAME).rst
	-rm -f $(BASE_NAME).html
	-rm -f $(BASE_NAME).pdf

##---------------------------------------------------------------------------##
## end of Makefile
##---------------------------------------------------------------------------##
