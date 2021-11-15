# Define these vars to see invalid usage checking and errors!
if (SPKB_SHOW_UPSTREAM_DEPLIBS_ERROR)
  set(EXTRA_TAL_ARGS DEPLIBS simplecxx)
elseif (SPKB_SHOW_M_DEPLIBS_ERROR)
  set(EXTRA_TAL_ARGS DEPLIBS m)
elseif (SPKB_SHOW_NONTESTONLY_IMPORTEDLIBS_ERROR)
  set(EXTRA_TAL_ARGS IMPORTEDLIBS simplecxx)
endif()
print_nonempty_var(EXTRA_TAL_ARGS)
