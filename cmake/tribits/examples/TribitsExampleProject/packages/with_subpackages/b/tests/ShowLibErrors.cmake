# Define these vars (one at a time) in your input cache to see invalid usage
# checking and the errors and warnings printed!
if (SPKB_SHOW_TESTONLY_INSTALLABLE_ERROR)
  set(TAEAT_EXTRALIB_ARGS INSTALLABLE)
elseif (SPKB_SHOW_NON_TESTONLY_LIB_ERROR)
  set(TAEAT_EXTRALIB_ARGS simplecxx)
elseif (SPKB_SHOW_IMPORTED_LIBS_THIS_PKG_ERROR)
  set(TAEAT_EXTRALIB_ARGS IMPORTEDLIBS pws_b)
elseif (SPKB_SHOW_TESTONLY_DEBLIBS_WARNING) # Deprecated
  set(TAEAT_EXTRALIB_ARGS  DEPLIBS b_mixed_lang)
elseif (SPKB_SHOW_NONTESTONLY_DEBLIBS_WARNING) # Deprecated
  set(TAEAT_EXTRALIB_ARGS  DEPLIBS pws_b)
elseif (SPKB_SHOW_EXTERNAL_DEBLIBS_WARNING) # Deprecated
  set(TAEAT_EXTRALIB_ARGS  DEPLIBS m)
endif()
print_nonempty_var(TAEAT_EXTRALIB_ARGS)

