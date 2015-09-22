# Define these vars (one at a time) in your input cach to see invalid usage
# checking and the errors and warnings printed!
IF (SPKC_SHOW_TESTONLY_DEPLBIS_ERROR)
  SET(TAL_EXTRALIB_ARGS DEPLIBS b_mixed_lang)
ELSEIF (SPKC_SHOW_TESTONLY_IMPORTEDLIBS_ERROR)
  SET(TAL_EXTRALIB_ARGS IMPORTEDLIBS b_mixed_lang)
ENDIF()
PRINT_NONEMPTY_VAR(TAL_EXTRALIB_ARGS)
