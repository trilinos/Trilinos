##---------------------------------------------------------------------------##
## MAKEFILE
##---------------------------------------------------------------------------##

EXTRA_ARGS = 

.PHONY: clean all

BASE_NAME = TribitsUsersGuide

include ../Makefile.common_generated_files

EXTRA_FILES = \
  Makefile \
  $(BASE_NAME).rst

all: $(BASE_NAME).html

$(BASE_NAME).html: $(EXTRA_FILES) $(COMMON_DEPENDENT_FILES)
	../../../python_utils/generate-docutils-output.py \
	--file-base=$(BASE_NAME) \
	--generate-latex= --generate-pdf= \
	$(EXTRA_ARGS)

# 	  --generate-latex-options="--stylesheet-path=rst2latex.tex" \

# NOTE: Can't build PDF anymore because it hangs.  Will have to debug later if
# we want PDF again.

clean:
	-rm -f $(BASE_NAME).html
	-rm -f $(BASE_NAME).pdf

##---------------------------------------------------------------------------##
## end of Makefile
##---------------------------------------------------------------------------##
