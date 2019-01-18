##---------------------------------------------------------------------------##
## MAKEFILE
##---------------------------------------------------------------------------##

EXTRA_ARGS = 

.PHONY: clean all

BASE_NAME = TribitsDevelopersGuide
EXTRA_FILES = \
  TribitsGitVersion.txt \
  generate-dev-guide.sh \
  Makefile \
  TribitsDevelopersGuide.rst \
  TribitsStandardTPLsList.txt \
  TribitsCommonTPLsList.txt \
  checkin-test-help.txt \
  gitdist-help.txt \
  gitdist-dist-help-overview.txt \
  gitdist-dist-help-repo-selection-and-setup.txt \
  gitdist-dist-help-dist-repo-status.txt \
  gitdist-dist-help-repo-versions.txt \
  gitdist-dist-help-aliases.txt \
  gitdist-dist-help-usage-tips.txt \
  gitdist-dist-help-script-dependencies.txt \
  gitdist-dist-help-all.txt \
  clone_extra_repos-help.txt \
  snapshot-dir-help.txt \
  is_checkin_tested_commit.txt \
  get-tribits-packages-from-files-list.txt \
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
