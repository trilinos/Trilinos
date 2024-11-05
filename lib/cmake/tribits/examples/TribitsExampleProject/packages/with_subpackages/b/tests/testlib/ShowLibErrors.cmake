# Define these vars to see invalid usage checking and errors!
if (SPKB_SHOW_PKG_LIB_IMPORTEDLIBS_ERROR)
  set(EXTRA_TAL_ARGS IMPORTEDLIBS pws_b)
elseif (SPKB_SHOW_UPSTREAM_PKG_LIB_IMPORTEDLIBS_ERROR)
  set(EXTRA_TAL_ARGS IMPORTEDLIBS simplecxx)
endif()
print_nonempty_var(EXTRA_TAL_ARGS)
