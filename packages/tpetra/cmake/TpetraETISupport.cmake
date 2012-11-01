include(AppendSet)

MACRO(TpetraETI_Explode etivar cs ds lo go n)
  STRING(REGEX REPLACE ".*CS=([^ ]*).*" "\\1" ${cs} ${etivar})
  STRING(REGEX REPLACE ".*DS=([^ ]*).*" "\\1" ${ds} ${etivar})
  STRING(REGEX REPLACE ".*LO=([^ ]*).*" "\\1" ${lo} ${etivar})
  STRING(REGEX REPLACE ".*GO=([^ ]*).*" "\\1" ${go} ${etivar})
  STRING(REGEX REPLACE  ".*N=([^ ]*).*" "\\1" ${n}  ${etivar})
ENDMACRO()

# this macro generates the explicit instantiation macros over the cross product of a number of lists
# it is used for the backwards-compatible method of specifying ETI types
FUNCTION(TpetraExpandTypesetProduct etisetvar datascalarlist lolist golist nodelist)
  foreach(ds ${datascalarlist})
    foreach(lo ${lolist})
      foreach(go ${golist})
        foreach(n ${nodelist})
          SET(etiset "${etiset};CS=${ds} DS=${ds} GO=${go} LO=${lo} N=${n}")
        endforeach()
      endforeach()
    endforeach()
  endforeach()
  SET(${etisetvar} "${${etisetvar}};${etiset}" PARENT_SCOPE)
ENDFUNCTION()

# this macro accepts a master list with elements CS,DS,GO,LO,N 
# it returns the following macros:
# macro_tslgn:  macro to instantiate CS,DS,GO,LO,N
# macro_slgn:  macro to instantiate DS,LO,GO,N  (for CS=DS uses)
# macro_lgn:   macro to instantiate LO,GO,N     (for non-scalar objects)
FUNCTION(GenerateTpetraETIMacros etisetvar tslgn slgn lgn)
  # we make lists of tuples first, because we want to make sure they don't have duplicates
  # this algorithm is O(N^2) in the number of instantiations in etisetvar
  foreach(inst ${etisetvar})
    # strip out the types using regex; we're not assuming that they are ordered
    TpetraETI_Explode(${inst} cs ds lo go n)
    # append tuple to lists
    list(APPEND tuple_tslgn "${cs},${ds},${lo},${go},${n}")
    list(APPEND tuple_slgn "${ds},${lo},${go},${n}")
    list(APPEND tuple_lgn  "${lo},${go},${n}")
  endforeach()
  list(REMOVE_DUPLICATES tuple_tslgn)
  list(REMOVE_DUPLICATES tuple_slgn)
  list(REMOVE_DUPLICATES tuple_lgn)
  # build the tslgn macro string
  set(macro_tslgn "#define TPETRA_INSTANTIATE_TSLGN(INSTMACRO)")
  foreach(f ${tuple_tslgn})
    string(REGEX REPLACE "," ";" f ${f})
    list(GET f 0 cs)
    list(GET f 1 ds)
    list(GET f 2 lo)
    list(GET f 3 go)
    list(GET f 4 n)
    set(macro_tslgn "${macro_tslgn} \\\nINSTMACRO(${cs},${ds},${lo},${go},${n})")
  endforeach()
  set(macro_tslgn "${macro_tslgn}\n")
  # build the slgn macro string
  set(macro_slgn "#define TPETRA_INSTANTIATE_SLGN(INSTMACRO)")
  foreach(f ${tuple_slgn})
    string(REGEX REPLACE "," ";" f ${f})
    list(GET f 0 s)
    list(GET f 1 lo)
    list(GET f 2 go)
    list(GET f 3 n)
    set(macro_slgn "${macro_slgn} \\\nINSTMACRO(${s},${lo},${go},${n})")
  endforeach()
  set(macro_slgn "${macro_slgn}\n")
  # build the lgn macro string
  set(macro_lgn "#define TPETRA_INSTANTIATE_LGN(INSTMACRO)")
  foreach(f ${tuple_lgn})
    string(REGEX REPLACE "," ";" f ${f})
    list(GET f 0 lo)
    list(GET f 1 go)
    list(GET f 2 n)
    set(macro_lgn "${macro_lgn} \\\nINSTMACRO(${lo},${go},${n})")
  endforeach()
  set(macro_lgn "${macro_lgn}\n")
  #
  set(${tslgn} "${macro_tslgn}" PARENT_SCOPE)
  set(${slgn}  "${macro_slgn}"  PARENT_SCOPE)
  set(${lgn}   "${macro_lgn}"   PARENT_SCOPE)
ENDFUNCTION()
