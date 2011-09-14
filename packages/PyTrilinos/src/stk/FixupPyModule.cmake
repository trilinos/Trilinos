
include(BundleUtilities)

function(gp_item_default_embedded_path_override item path)
  set(${path} "@loader_path" PARENT_SCOPE)
endfunction(gp_item_default_embedded_path_override)

# python libs are always a given
function(gp_resolved_file_type_override file type)
  if(file MATCHES Python)
    set(${type} system PARENT_SCOPE)
  endif(file MATCHES Python)
endfunction(gp_resolved_file_type_override)

function(my_set_bundle_key_values keys_var context item exepath dirs copyflag)
  get_filename_component(item_name "${item}" NAME)

  get_item_key("${item}" key)

  list(LENGTH ${keys_var} length_before)
  gp_append_unique(${keys_var} "${key}")
  list(LENGTH ${keys_var} length_after)

  if(NOT length_before EQUAL length_after)
    gp_resolve_item("${context}" "${item}" "${exepath}" "${dirs}" resolved_item)

    gp_item_default_embedded_path("${item}" default_embedded_path)

    if(item MATCHES "[^/]+\\.framework/")
      # For frameworks, construct the name under the embedded path from the
      # opening "${item_name}.framework/" to the closing "/${item_name}":
      #
      string(REGEX REPLACE "^.*(${item_name}.framework/.*/${item_name}).*$" "${default_embedded_path}/\\1" embedded_item "${item}")
    else(item MATCHES "[^/]+\\.framework/")
      # For other items, just use the same name as the original, but in the
      # embedded path:
      #
      set(embedded_item "${default_embedded_path}/${item_name}")
    endif(item MATCHES "[^/]+\\.framework/")

    # Replace @executable_path and resolve ".." references:
    #
    string(REPLACE "@executable_path" "${exepath}" resolved_embedded_item "${embedded_item}")
    # replace loader_path (but should really should use path of using library instead of exepath as a substitute)
    string(REPLACE "@loader_path" "${exepath}" resolved_embedded_item "${resolved_embedded_item}")
    get_filename_component(resolved_embedded_item "${resolved_embedded_item}" ABSOLUTE)

    # *But* -- if we are not copying, then force resolved_embedded_item to be
    # the same as resolved_item. In the case of multiple executables in the
    # original bundle, using the default_embedded_path results in looking for
    # the resolved executable next to the main bundle executable. This is here
    # so that exes in the other sibling directories (like "bin") get fixed up
    # properly...
    #
    if(NOT copyflag)
      set(resolved_embedded_item "${resolved_item}")
    endif(NOT copyflag)

    set(${keys_var} ${${keys_var}} PARENT_SCOPE)
    set(${key}_ITEM "${item}" PARENT_SCOPE)
    set(${key}_RESOLVED_ITEM "${resolved_item}" PARENT_SCOPE)
    set(${key}_DEFAULT_EMBEDDED_PATH "${default_embedded_path}" PARENT_SCOPE)
    set(${key}_EMBEDDED_ITEM "${embedded_item}" PARENT_SCOPE)
    set(${key}_RESOLVED_EMBEDDED_ITEM "${resolved_embedded_item}" PARENT_SCOPE)
    set(${key}_COPYFLAG "${copyflag}" PARENT_SCOPE)
  else(NOT length_before EQUAL length_after)
    #message("warning: item key '${key}' already in the list, subsequent references assumed identical to first")
  endif(NOT length_before EQUAL length_after)
endfunction(my_set_bundle_key_values)

function(fixup_py_module root libs dirs)

  foreach(lib ${libs})
      my_set_bundle_key_values(keys "${lib}" "${lib}" "${root}" "${dirs}" 0)

      set(prereqs "")
      get_prerequisites("${lib}" prereqs 1 1 "${root}" "${dirs}")
      foreach(pr ${prereqs})
        my_set_bundle_key_values(keys "${lib}" "${pr}" "${root}" "${dirs}" 1)
      endforeach(pr)
    endforeach(lib)

    foreach(key ${keys})
      set(${key}_ITEM "${${key}_ITEM}" PARENT_SCOPE)
      set(${key}_RESOLVED_ITEM "${${key}_RESOLVED_ITEM}" PARENT_SCOPE)
      set(${key}_DEFAULT_EMBEDDED_PATH "${${key}_DEFAULT_EMBEDDED_PATH}" PARENT_SCOPE)
      set(${key}_EMBEDDED_ITEM "${${key}_EMBEDDED_ITEM}" PARENT_SCOPE)
      set(${key}_RESOLVED_EMBEDDED_ITEM "${${key}_RESOLVED_EMBEDDED_ITEM}" PARENT_SCOPE)
      set(${key}_COPYFLAG "${${key}_COPYFLAG}" PARENT_SCOPE)
    endforeach(key)

  message(STATUS "fixup_py_module: copying...")
    list(LENGTH keys n)
    math(EXPR n ${n}*2)

    set(i 0)
    foreach(key ${keys})
      math(EXPR i ${i}+1)
      if(${${key}_COPYFLAG})
        message(STATUS "${i}/${n}: copying '${${key}_RESOLVED_ITEM}'")
      else(${${key}_COPYFLAG})
        message(STATUS "${i}/${n}: *NOT* copying '${${key}_RESOLVED_ITEM}'")
      endif(${${key}_COPYFLAG})

      set(show_status 0)
      if(show_status)
        message(STATUS "key='${key}'")
        message(STATUS "item='${${key}_ITEM}'")
        message(STATUS "resolved_item='${${key}_RESOLVED_ITEM}'")
        message(STATUS "default_embedded_path='${${key}_DEFAULT_EMBEDDED_PATH}'")
        message(STATUS "embedded_item='${${key}_EMBEDDED_ITEM}'")
        message(STATUS "resolved_embedded_item='${${key}_RESOLVED_EMBEDDED_ITEM}'")
        message(STATUS "copyflag='${${key}_COPYFLAG}'")
        message(STATUS "")
      endif(show_status)

      if(${${key}_COPYFLAG})
        set(item "${${key}_ITEM}")
        if(item MATCHES "[^/]+\\.framework/")
          copy_resolved_framework_into_bundle("${${key}_RESOLVED_ITEM}"
            "${${key}_RESOLVED_EMBEDDED_ITEM}")
        else()
          message("In else item = " ${item})
          copy_resolved_item_into_bundle("${${key}_RESOLVED_ITEM}"
            "${${key}_RESOLVED_EMBEDDED_ITEM}")
        endif()
      endif(${${key}_COPYFLAG})
    endforeach(key)

    message(STATUS "fixup_py_module: fixing...")
    foreach(key ${keys})
      math(EXPR i ${i}+1)
      if(APPLE)
        message(STATUS "${i}/${n}: fixing up '${${key}_RESOLVED_EMBEDDED_ITEM}'")
        set(prereqs "")
        set(changes "")
        get_prerequisites("${${key}_RESOLVED_EMBEDDED_ITEM}" prereqs 1 0 "${root}" "${dirs}")
        foreach(pr ${prereqs})
          get_item_key("${pr}" rkey)

          if(NOT "${${rkey}_EMBEDDED_ITEM}" STREQUAL "")
            set(changes ${changes} "-change" "${pr}" "${${rkey}_EMBEDDED_ITEM}")
          else(NOT "${${rkey}_EMBEDDED_ITEM}" STREQUAL "")
             message("warning: unexpected reference to '${pr}'")
          endif(NOT "${${rkey}_EMBEDDED_ITEM}" STREQUAL "")
        endforeach(pr)

     if(BU_CHMOD_BUNDLE_ITEMS)
       execute_process(COMMAND chmod u+w "${${key}_RESOLVED_EMBEDDED_ITEM}")
     endif()

     execute_process(COMMAND install_name_tool
       ${changes} -id "${${key}_EMBEDDED_ITEM}" "${${key}_RESOLVED_EMBEDDED_ITEM}"
     )

     #   fixup_bundle_item("${${key}_RESOLVED_EMBEDDED_ITEM}" "${root}" "${dirs}")
      else(APPLE)
        message(STATUS "${i}/${n}: fix-up not required on this platform '${${key}_RESOLVED_EMBEDDED_ITEM}'")
      endif(APPLE)
    endforeach(key)

    message(STATUS "fixup_py_module: cleaning up...")
    clear_bundle_keys(keys)
  

endfunction(fixup_py_module)

