#
# @FUNCTION: TRIBITS_STRIP_QUOTES_FROM_STR()
#
# Remove one set of quotes from the outside of a string if they exist.
#
# Usage::
#
#   TRIBITS_STRIP_QUOTES_FROM_STR(<str_in> <str_var_out>)
#
# If ``<str_in>`` does not contain a quote char ``'"'`` as the first and last
# char, then the original ``<str_in>`` is returned in ``<str_var_out>``.
#
FUNCTION(TRIBITS_STRIP_QUOTES_FROM_STR str_in str_var_out)
  #print_var(str_in)
  string(LENGTH "${str_in}" str_len)
  #print_var(str_len)
  if (str_len LESS 2)
    # Can't have two parath \" chars of does not have 2 chars!
    set(${str_var_out} "${str_in}" PARENT_SCOPE)
    return()
  endif()
  math(EXPR str_len_less_1 "${str_len} - 1")
  math(EXPR str_len_less_2 "${str_len} - 2")
  #print_var(str_len_less_1)
  #print_var(str_len_less_2)
  string(SUBSTRING "${str_in}" 0 1 str_first_char)
  #print_var(str_first_char)
  if (str_len_less_1 GREATER 0)
    string(SUBSTRING "${str_in}" ${str_len_less_1} 1 str_last_char)
    #print_var(str_last_char)
    if (str_first_char STREQUAL "\"" AND str_last_char STREQUAL "\"")
      string(SUBSTRING "${str_in}" 1 ${str_len_less_2} str_parath_removed)
      #print_var(str_parath_removed)
      set(${str_var_out} "${str_parath_removed}" PARENT_SCOPE)
      return()
    endif()
  endif()
  # If we get here, there were not \" chars as the first and last char so
  # return the entire input string!
  set(${str_var_out} "${str_in}" PARENT_SCOPE)
endfunction()