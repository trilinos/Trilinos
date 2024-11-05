# Define these vars (one at a time) in your input cache to see invalid usage
# checking and the errors and warnings printed!
if (SPKC_SHOW_TESTONLY_DEPLBIS_ERROR)
  set(TAL_EXTRALIB_ARGS DEPLIBS b_mixed_lang)
elseif (SPKC_SHOW_TESTONLY_IMPORTEDLIBS_ERROR)
  set(TAL_EXTRALIB_ARGS IMPORTEDLIBS b_mixed_lang)
endif()
print_nonempty_var(TAL_EXTRALIB_ARGS)
