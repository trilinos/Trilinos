##---------------------------------------------------------------------------##
## MAKEFILE
##---------------------------------------------------------------------------##

EXTRA_ARGS = 

.PHONY: clean all

BASE_NAME = TribitsDevelopersGuide
EXTRA_FILES = \
  generate-dev-guide.sh \
  Makefile \
  TribitsDevelopersGuide.rst \
  TribitsStandardTPLsList.txt \
  TribitsCommonTPLsList.txt \
  checkin-test-help.txt \
  gitdist-help.txt \
  clone_extra_repos-help.txt \
  snapshot-dir-help.txt \
  install_devtools-help.txt \
  TriBITS.README.DIRECTORY_CONTENTS.rst \
  TribitsMacroFunctionDoc.rst \
  UtilsMacroFunctionDoc.rst

# NOTE: We *don't* want to rebuild the file if only the version has changed we
# don't add a dependency on TribitsGitVersion.txt!

all: $(BASE_NAME).rst $(BASE_NAME).html $(BASE_NAME).pdf

$(BASE_NAME).pdf: $(BASE_NAME).html 

$(BASE_NAME).html: $(EXTRA_FILES) 
	../../python_utils/generate-docutils-output.py \
	  --file-base=TribitsDevelopersGuide \
	  --generate-latex-options="--stylesheet-path=rst2latex.tex" \
	  $(EXTRA_ARGS)

clean:
	-rm -f $(BASE_NAME).html
	-rm -f $(BASE_NAME).pdf

##---------------------------------------------------------------------------##
## end of Makefile
##---------------------------------------------------------------------------##
