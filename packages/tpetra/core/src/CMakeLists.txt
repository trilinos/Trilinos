INCLUDE(TrilinosCreateClientTemplateHeaders)

# Function that follows the Tpetra convention for mangling C++ types
# so that they can be used as C preprocessor macro arguments.
#
# TYPE_MANGLED_OUT [out] The mangled type name.
#
# TYPE_IN [in] The type to mangle.
FUNCTION(TPETRA_MANGLE_TEMPLATE_PARAMETER TYPE_MANGLED_OUT TYPE_IN)
  STRING(REPLACE "<" "0" TMP0 "${TYPE_IN}")
  STRING(REPLACE ">" "0" TMP1 "${TMP0}")
  STRING(REPLACE "::" "_" TMP2 "${TMP1}")
  # Spaces (as in "long long") get squished out.
  STRING(REPLACE " " "" TMP3 "${TMP2}")
  SET(${TYPE_MANGLED_OUT} ${TMP3} PARENT_SCOPE)
ENDFUNCTION(TPETRA_MANGLE_TEMPLATE_PARAMETER)

# Function that turns a valid Scalar, LocalOrdinal, or GlobalOrdinal
# template parameter into a macro name (all caps, with no white space
# and no punctuation other than underscore).
#
# NAME_OUT [out] The mangled type name.
#
# NAME_IN [in] The type to mangle.
FUNCTION(TPETRA_SLG_MACRO_NAME NAME_OUT NAME_IN)
  STRING(COMPARE EQUAL "${NAME_IN}" "__float128" IS_FLOAT128)
  IF(IS_FLOAT128)
    # __float128 is a special case; we remove the __ from the macro name.
    SET(${NAME_OUT} "FLOAT128" PARENT_SCOPE)
  ELSE()
    STRING(COMPARE EQUAL "${NAME_IN}" "std::complex<float>" IS_COMPLEX_FLOAT)
    IF(IS_COMPLEX_FLOAT)
      SET(${NAME_OUT} "COMPLEX_FLOAT" PARENT_SCOPE)
    ELSE()
      STRING(COMPARE EQUAL "${NAME_IN}" "std::complex<double>" IS_COMPLEX_DOUBLE)
      IF(IS_COMPLEX_DOUBLE)
        SET(${NAME_OUT} "COMPLEX_DOUBLE" PARENT_SCOPE)
      ELSE()
	#long double is a special name; add _ to macro name
        STRING(COMPARE EQUAL "${NAME_IN}" "long double" IS_LONG_DOUBLE)
	IF(IS_LONG_DOUBLE)
	  SET(${NAME_OUT} "LONG_DOUBLE" PARENT_SCOPE)
        ELSE()
          # Convert to upper case, convert double colons to underscores,
          # and hope for the best.
          #
          # It would be nice if CMake were consistent about where output
          # arguments go.  Alas, this is not to be.  TOUPPER puts the
          # output argument last; REPLACE puts it after the search and
          # substitute strings, before the input string.
          STRING(TOUPPER "${NAME_IN}" TMP0)
          STRING(REPLACE "::" "_" TMP1 "${TMP0}")
          STRING(REPLACE " " "_" TMP2 "${TMP1}")
          SET(${NAME_OUT} ${TMP2} PARENT_SCOPE)
        ENDIF()
      ENDIF()
    ENDIF()
  ENDIF()
ENDFUNCTION(TPETRA_SLG_MACRO_NAME)

SET(VALID_GO_TYPES "short;unsigned short;int;unsigned int;long;unsigned long;long long;unsigned long long")

# Whether the input SC (Scalar) type is a valid GO (GlobalOrdinal) type.
FUNCTION(TPETRA_SC_IS_GO IS_GO SC)
  FOREACH(VALID_GO ${VALID_GO_TYPES})
    STRING(COMPARE EQUAL "${VALID_GO}" "${SC}" IS_GO_TMP0)
    IF (IS_GO_TMP0)
      # Now would be a good chance to break from the loop, if I knew
      # how to do that.
      SET(IS_GO_TMP TRUE)
    ENDIF()
  ENDFOREACH()

  SET(${IS_GO} ${IS_GO_TMP} PARENT_SCOPE)
ENDFUNCTION()

# Function that turns a valid Node template parameter into a macro
# name (all caps, with no white space and no punctuation other than
# underscore).
# Upper-cases the Kokkos execution space name, except for Threads,
# which is turned into PTHREAD (cwp: for unclear reasons)
#
# NAME_OUT [out] The mangled type name.
#
# NAME_IN [in] The type to mangle.
FUNCTION(TPETRA_NODE_MACRO_NAME NAME_OUT NAME_IN)
  STRING(REGEX MATCH "Tpetra::KokkosCompat::Kokkos(.*)WrapperNode" TMP0 "${NAME_IN}")
  STRING(COMPARE EQUAL "${TMP0}" "" DOES_NOT_MATCH)

  IF(DOES_NOT_MATCH)
    MESSAGE(FATAL_ERROR "Tpetra: Node $NAME_IN is not a supported Node type.")
  ELSE()
    # Extract the Kokkos execution space (KOKKOS_EXEC_SPACE) from the Node name.
    STRING(REGEX REPLACE "Tpetra::KokkosCompat::Kokkos(.*)WrapperNode" "\\1" KOKKOS_EXEC_SPACE "${NAME_IN}")

    # Special case: Threads.  The macro name unfortunately differs
    # from the execution space name in a way that doesn't fit the
    # pattern of the other execution spaces.
    STRING(COMPARE EQUAL "${KOKKOS_EXEC_SPACE}" "Threads" IS_THREADS)
    IF(IS_THREADS)
      SET(${NAME_OUT} "PTHREAD" PARENT_SCOPE)
    ELSE()
      # The other cases (Cuda, Serial, OpenMP) are easy.
      STRING(TOUPPER "${KOKKOS_EXEC_SPACE}" NAME_OUT_TMP)
      SET(${NAME_OUT} ${NAME_OUT_TMP} PARENT_SCOPE)
    ENDIF()
  ENDIF()
ENDFUNCTION(TPETRA_NODE_MACRO_NAME)

# Function that turns Scalar (SC) and GlobalOrdinal (GO) type names
# into an expression for asking Tpetra whether to build for that
# Scalar type.
#
# SC_MACRO_EXPR [out] Expression for asking Tpetra whether to build
#   for that Scalar type.
#
# SC [in] Original name of the Scalar type.
#
# GO [in] Original name of the GlobalOrdinal type.

# SC_MACRO_NAME [in] Macro-name version of SC.  The
#   TPETRA_SLG_MACRO_NAME function (see above) implements the
#   conversion process from the original name to the macro name.
FUNCTION(TPETRA_SC_MACRO_EXPR SC_MACRO_EXPR SC GO SC_MACRO_NAME)
  # SC = int,char and SC = GO are special cases.  Tpetra doesn't have
  # macros for these cases.  That means the expression is empty.
  STRING(COMPARE EQUAL "${SC}" "int" IS_INT)
  IF(IS_INT)
    SET(SC_MACRO_EXPR_TMP "")
  ELSE()
    STRING(COMPARE EQUAL "${SC}" "char" IS_CHAR)
    IF(IS_CHAR)
      SET(SC_MACRO_EXPR_TMP "")
    ELSE()
      STRING(COMPARE EQUAL "${SC}" "${GO}" IS_GO)
      IF(IS_GO)
        SET(SC_MACRO_EXPR_TMP "")
      ELSE()
        SET(SC_MACRO_EXPR_TMP "&& defined(HAVE_TPETRA_INST_${SC_MACRO_NAME})")
      ENDIF()
    ENDIF()
  ENDIF()

  #MESSAGE(STATUS ">> >> SC = ${SC}, SC_MACRO_EXPR_TMP = ${SC_MACRO_EXPR_TMP}")

  # Set the output argument.
  SET(${SC_MACRO_EXPR} "${SC_MACRO_EXPR_TMP}" PARENT_SCOPE)
ENDFUNCTION(TPETRA_SC_MACRO_EXPR)

# Function to generate one .cpp file for the given (Scalar,
# LocalOrdinal, GlobalOrdinal, Node) template parameter combination,
# for run-time registration of a Tpetra class or function over
# those template parameters.  This is meant to be called by
# TPETRA_PROCESS_ALL_SLGN_TEMPLATES.  This function takes the names
# already mangled, to avoid unnecessary string processing overhead.
#
# OUTPUT_FILE [out] Name of the generated .cpp file.
#
# TEMPLATE_FILE [in] Name of the input .tmpl "template" file.  This
#   function does string substitution in that file, using the input
#   arguments of this function.  For example, @SC_MACRO_EXPR@ (Scalar
#   macro expression) gets substituted for the value of this
#   function's SC_MACRO_EXPR input argument.
#
# CLASS_NAME [in] Name of the Tpetra class (without the Tpetra
#   namespace qualifier; must live in the Tpetra namespace),
#   suitably mangled for use in a file name.
#
# CLASS_MACRO_NAME [in] Name of the Tpetra class, suitably mangled
#   for use in a macro name.
#
# SC_MANGLED_NAME [in] Name of the Scalar (SC) type, mangled for use
#   as a macro argument (e.g., spaces and colons removed).  In the
#   arguments that follow, LO stands for LocalOrdinal, GO for
#   GlobalOrdinal, and NT for Node.
#
# SC_MACRO_EXPR [in] Expression that asks Tpetra
#   whether the given Scalar (SC) type is supported.
#
# LO_MACRO_NAME [in] Name of the LocalOrdinal (LO) type,
#   mangled for use as a macro argument.
#
# GO_MACRO_NAME [in] Name of the GlobalOrdinal (GO) type,
#   mangled for use as a macro argument.
#
# NT_MACRO_NAME [in] Name of the Node (NT) type,
#   mangled for use as a macro argument.
#
FUNCTION(TPETRA_PROCESS_ONE_SLGN_TEMPLATE OUTPUT_FILE TEMPLATE_FILE
  CLASS_NAME CLASS_MACRO_NAME SC_MANGLED_NAME LO_MANGLED_NAME
  GO_MANGLED_NAME NT_MANGLED_NAME SC_MACRO_EXPR LO_MACRO_NAME
  GO_MACRO_NAME NT_MACRO_NAME)

  STRING(REPLACE "ETI_SC_LO_GO_NT.tmpl"
    "${CLASS_NAME}_${SC_MACRO_NAME}_${LO_MACRO_NAME}_${GO_MACRO_NAME}_${NT_MACRO_NAME}.cpp"
    OUT_FILE "${TEMPLATE_FILE}")
  CONFIGURE_FILE("${TEMPLATE_FILE}" "${OUT_FILE}")

  SET(${OUTPUT_FILE} ${OUT_FILE} PARENT_SCOPE)
ENDFUNCTION(TPETRA_PROCESS_ONE_SLGN_TEMPLATE)

# Function to generate one .cpp file for the given (ScalarOut,
# ScalarIn, LocalOrdinal, GlobalOrdinal, Node) template parameter
# combination, for run-time registration of a Tpetra class or function
# over those template parameters. This is meant to be called by
# TPETRA_PROCESS_ALL_CONVERT_TEMPLATES. This function takes the names
# already mangled, to avoid unnecessary string processing overhead.
#
# OUTPUT_FILE [out] Name of the generated .cpp file.
#
# TEMPLATE_FILE [in] Name of the input .tmpl "template" file.  This
#   function does string substitution in that file, using the input
#   arguments of this function.  For example, @SC_MACRO_EXPR@ (Scalar
#   macro expression) gets substituted for the value of this
#   function's SC_MACRO_EXPR input argument.
#
# CLASS_NAME [in] Name of the Tpetra class (without the Tpetra
#   namespace qualifier; must live in the Tpetra namespace),
#   suitably mangled for use in a file name.
#
# CLASS_MACRO_NAME [in] Name of the Tpetra class, suitably mangled
#   for use in a macro name.
#
# SOUT_MANGLED_NAME [in] Name of the output Scalar (SC) type, mangled for use
#   as a macro argument (e.g., spaces and colons removed).  In the
#   arguments that follow, LO stands for LocalOrdinal, GO for
#   GlobalOrdinal, and NT for Node.
#
# SIN_MANGLED_NAME [in] Name of the input Scalar (SC) type, mangled for use
#   as a macro argument (e.g., spaces and colons removed).  In the
#   arguments that follow, LO stands for LocalOrdinal, GO for
#   GlobalOrdinal, and NT for Node.
#
# SC_MACRO_EXPR [in] Expression that asks Tpetra
#   whether the given Scalar (SC) type is supported.
#
# LO_MACRO_NAME [in] Name of the LocalOrdinal (LO) type,
#   mangled for use as a macro argument.
#
# GO_MACRO_NAME [in] Name of the GlobalOrdinal (GO) type,
#   mangled for use as a macro argument.
#
# NT_MACRO_NAME [in] Name of the Node (NT) type,
#   mangled for use as a macro argument.
#
FUNCTION(TPETRA_PROCESS_ONE_CONVERT_TEMPLATE OUTPUT_FILE TEMPLATE_FILE
  CLASS_NAME CLASS_MACRO_NAME SOUT_MANGLED_NAME SIN_MANGLED_NAME LO_MANGLED_NAME
  GO_MANGLED_NAME NT_MANGLED_NAME SOUT_MACRO_EXPR SIN_MACRO_EXPR LO_MACRO_NAME
  GO_MACRO_NAME NT_MACRO_NAME)

  STRING(REPLACE "ETI_SOUT_SIN_LO_GO_NT.tmpl"
    "${CLASS_NAME}_${SOUT_MACRO_NAME}_${SIN_MACRO_NAME}_${LO_MACRO_NAME}_${GO_MACRO_NAME}_${NT_MACRO_NAME}.cpp"
    OUT_FILE "${TEMPLATE_FILE}")
  CONFIGURE_FILE("${TEMPLATE_FILE}" "${OUT_FILE}")

  SET(${OUTPUT_FILE} ${OUT_FILE} PARENT_SCOPE)
ENDFUNCTION(TPETRA_PROCESS_ONE_CONVERT_TEMPLATE)

# Function to generate .cpp files for ETI of a Tpetra class or
# function, over all enabled Scalar, LocalOrdinal, GlobalOrdinal, and
# Node template parameters.  We generate one .cpp file for each
# (Scalar, LocalOrdinal, GlobalOrdinal, Node) type combination over
# which Tpetra does ETI.
#
# OUTPUT_FILES [out] List of the generated .cpp files.
#
# TEMPLATE_FILE [in] Name of the input .tmpl "template" file.  This
#   function does string substitution in that file, using the input
#   arguments of this function.  For example, @SC_MACRO_EXPR@ (Scalar
#   macro expression) gets substituted for the value of this
#   function's SC_MACRO_EXPR input argument.
#
# CLASS_NAME [in] Name of the Tpetra class (without namespace
#   qualifiers; must live in the Tpetra namespace)
#
# CLASS_MACRO_NAME [in] Name of the Tpetra class, suitably mangled for
#   use in a macro name.
#
# SCALAR_TYPES [in] All Scalar types over which to do ETI for the given
#   class.  This may include Scalar = GlobalOrdinal and/or Scalar =
#   int, if appropriate for that class.
#
# LOCALORDINAL_TYPES [in] All LocalOrdinal types over which to do ETI
#   for the given class.
#
# GLOBALORDINAL_TYPES [in] All GlobalOrdinal types over which to do
#   ETI for the given class.
#
# NODE_TYPES [in] All Node types over which to do ETI for the given
#   class.
#
# MUST_HAVE_SCALAR_INT [in] (Boolean) Whether the class must be
#   instantiated with Scalar = int, even if int is not in the set of
#   GlobalOrdinal types.
FUNCTION(TPETRA_PROCESS_ALL_SLGN_TEMPLATES OUTPUT_FILES TEMPLATE_FILE
  CLASS_NAME CLASS_MACRO_NAME SCALAR_TYPES LOCALORDINAL_TYPES
  GLOBALORDINAL_TYPES NODE_TYPES MUST_HAVE_SCALAR_INT)

  SET(OUT_FILES "")
  FOREACH(NT ${NODE_TYPES})
    TPETRA_MANGLE_TEMPLATE_PARAMETER(NT_MANGLED "${NT}")
    TPETRA_NODE_MACRO_NAME(NT_MACRO_NAME "${NT}")
    FOREACH(GO ${GLOBALORDINAL_TYPES})
      TPETRA_MANGLE_TEMPLATE_PARAMETER(GO_MANGLED "${GO}")
      TPETRA_SLG_MACRO_NAME(GO_MACRO_NAME "${GO}")
      FOREACH(LO ${LOCALORDINAL_TYPES})
        TPETRA_MANGLE_TEMPLATE_PARAMETER(LO_MANGLED "${LO}")
        TPETRA_SLG_MACRO_NAME(LO_MACRO_NAME "${LO}")
        FOREACH(SC ${SCALAR_TYPES})
          TPETRA_MANGLE_TEMPLATE_PARAMETER(SC_MANGLED "${SC}")
          TPETRA_SLG_MACRO_NAME(SC_MACRO_NAME "${SC}")
          TPETRA_SC_MACRO_EXPR(SC_MACRO_EXPR "${SC}" "${GO}" "${SC_MACRO_NAME}")

          #MESSAGE(STATUS ">> SC = ${SC}, SC_MACRO_EXPR = ${SC_MACRO_EXPR}")

          # If SC is NOT a GlobalOrdinal type of some kind (not
          # necessarily the current GO), or if it is "int", process
          # it.  Otherwise, then we only have to process it if it
          # equals the current GO.
          TPETRA_SC_IS_GO(IS_GO "${SC}")
          STRING(COMPARE EQUAL "${SC}" "${GO}" IS_CURRENT_GO)
          STRING(COMPARE EQUAL "${SC}" "int" IS_INT)

          IF ((MUST_HAVE_SCALAR_INT AND IS_INT) OR (NOT IS_GO OR IS_CURRENT_GO))
            TPETRA_PROCESS_ONE_SLGN_TEMPLATE(OUT_FILE "${TEMPLATE_FILE}"
              "${CLASS_NAME}" "${CLASS_MACRO_NAME}" "${SC_MANGLED}"
              "${LO_MANGLED}" "${GO_MANGLED}" "${NT_MANGLED}"
              "${SC_MACRO_EXPR}" "${LO_MACRO_NAME}" "${GO_MACRO_NAME}"
              "${NT_MACRO_NAME}")
            LIST(APPEND OUT_FILES ${OUT_FILE})
          ENDIF()
        ENDFOREACH() # SC
      ENDFOREACH() # LO
    ENDFOREACH() # GO
  ENDFOREACH() # NT

  # This is the standard CMake idiom for setting an output variable so
  # that the caller can see the result.
  SET(${OUTPUT_FILES} ${OUT_FILES} PARENT_SCOPE)
ENDFUNCTION(TPETRA_PROCESS_ALL_SLGN_TEMPLATES)

# Function to generate .cpp files for ETI of a Tpetra class or
# function, over all enabled output and input Scalar types,
# LocalOrdinal, GlobalOrdinal, and Node template parameters. We
# generate one .cpp file for each (ScalarOut, ScalarIn, LocalOrdinal,
# GlobalOrdinal, Node) type combination over which Tpetra does ETI.
#
# OUTPUT_FILES [out] List of the generated .cpp files.
#
# TEMPLATE_FILE [in] Name of the input .tmpl "template" file.  This
#   function does string substitution in that file, using the input
#   arguments of this function.  For example, @SC_MACRO_EXPR@ (Scalar
#   macro expression) gets substituted for the value of this
#   function's SC_MACRO_EXPR input argument.
#
# CLASS_NAME [in] Name of the Tpetra class (without namespace
#   qualifiers; must live in the Tpetra namespace)
#
# CLASS_MACRO_NAME [in] Name of the Tpetra class, suitably mangled for
#   use in a macro name.
#
# SOUT_TYPES [in] All output Scalar types over which to do ETI for the given
#   class.  This may include Scalar = GlobalOrdinal and/or Scalar =
#   int, if appropriate for that class.
#
# SIN_TYPES [in] All input Scalar types over which to do ETI for the given
#   class.  This may include Scalar = GlobalOrdinal and/or Scalar =
#   int, if appropriate for that class.
#
# LOCALORDINAL_TYPES [in] All LocalOrdinal types over which to do ETI
#   for the given class.
#
# GLOBALORDINAL_TYPES [in] All GlobalOrdinal types over which to do
#   ETI for the given class.
#
# NODE_TYPES [in] All Node types over which to do ETI for the given
#   class.
#
# MUST_HAVE_SCALAR_INT [in] (Boolean) Whether the class must be
#   instantiated with Scalar = int, even if int is not in the set of
#   GlobalOrdinal types.
FUNCTION(TPETRA_PROCESS_ALL_CONVERT_TEMPLATES OUTPUT_FILES TEMPLATE_FILE
  CLASS_NAME CLASS_MACRO_NAME SOUT_TYPES SIN_TYPES LOCALORDINAL_TYPES
  GLOBALORDINAL_TYPES NODE_TYPES MUST_HAVE_SCALAR_INT)

  SET(OUT_FILES "")
  FOREACH(NT ${NODE_TYPES})
    TPETRA_MANGLE_TEMPLATE_PARAMETER(NT_MANGLED "${NT}")
    TPETRA_NODE_MACRO_NAME(NT_MACRO_NAME "${NT}")
    FOREACH(GO ${GLOBALORDINAL_TYPES})
      TPETRA_MANGLE_TEMPLATE_PARAMETER(GO_MANGLED "${GO}")
      TPETRA_SLG_MACRO_NAME(GO_MACRO_NAME "${GO}")
      FOREACH(LO ${LOCALORDINAL_TYPES})
        TPETRA_MANGLE_TEMPLATE_PARAMETER(LO_MANGLED "${LO}")
        TPETRA_SLG_MACRO_NAME(LO_MACRO_NAME "${LO}")
        FOREACH(SIN ${SIN_TYPES})
          TPETRA_MANGLE_TEMPLATE_PARAMETER(SIN_MANGLED "${SIN}")
          TPETRA_SLG_MACRO_NAME(SIN_MACRO_NAME "${SIN}")
          TPETRA_SC_MACRO_EXPR(SIN_MACRO_EXPR "${SIN}" "${GO}" "${SC_MACRO_NAME}")

          #MESSAGE(STATUS ">> SC = ${SC}, SC_MACRO_EXPR = ${SC_MACRO_EXPR}")

          # If SC is NOT a GlobalOrdinal type of some kind (not
          # necessarily the current GO), or if it is "int", process
          # it.  Otherwise, then we only have to process it if it
          # equals the current GO.
          TPETRA_SC_IS_GO(IS_GO "${SIN}")
          STRING(COMPARE EQUAL "${SIN}" "${GO}" SIN_IS_CURRENT_GO)
          STRING(COMPARE EQUAL "${SIN}" "int" SIN_IS_INT)

          IF ((MUST_HAVE_SCALAR_INT AND SIN_IS_INT) OR (NOT SIN_IS_GO OR SIN_IS_CURRENT_GO))

            FOREACH(SOUT ${SOUT_TYPES})
              TPETRA_MANGLE_TEMPLATE_PARAMETER(SOUT_MANGLED "${SOUT}")
              TPETRA_SLG_MACRO_NAME(SOUT_MACRO_NAME "${SOUT}")
              TPETRA_SC_MACRO_EXPR(SOUT_MACRO_EXPR "${SOUT}" "${GO}" "${SC_MACRO_NAME}")

              #MESSAGE(STATUS ">> SC = ${SC}, SC_MACRO_EXPR = ${SC_MACRO_EXPR}")

              # If SC is NOT a GlobalOrdinal type of some kind (not
              # necessarily the current GO), or if it is "int", process
              # it.  Otherwise, then we only have to process it if it
              # equals the current GO.
              TPETRA_SC_IS_GO(IS_GO "${SOUT}")
              STRING(COMPARE EQUAL "${SOUT}" "${GO}" SOUT_IS_CURRENT_GO)
              STRING(COMPARE EQUAL "${SOUT}" "int" SOUT_IS_INT)

              IF ((MUST_HAVE_SCALAR_INT AND SOUT_IS_INT) OR (NOT SOUT_IS_GO OR SOUT_IS_CURRENT_GO))

                TPETRA_PROCESS_ONE_CONVERT_TEMPLATE(OUT_FILE "${TEMPLATE_FILE}"
                  "${CLASS_NAME}" "${CLASS_MACRO_NAME}" "${SOUT_MANGLED}" "${SIN_MANGLED}"
                  "${LO_MANGLED}" "${GO_MANGLED}" "${NT_MANGLED}"
                  "${SOUT_MACRO_EXPR}" "${SIN_MACRO_EXPR}" "${LO_MACRO_NAME}" "${GO_MACRO_NAME}"
                  "${NT_MACRO_NAME}")
                LIST(APPEND OUT_FILES ${OUT_FILE})
              ENDIF()
            ENDFOREACH() # SOUT
          ENDIF()
        ENDFOREACH() # SIN
      ENDFOREACH() # LO
    ENDFOREACH() # GO
  ENDFOREACH() # NT

  # This is the standard CMake idiom for setting an output variable so
  # that the caller can see the result.
  SET(${OUTPUT_FILES} ${OUT_FILES} PARENT_SCOPE)
ENDFUNCTION(TPETRA_PROCESS_ALL_CONVERT_TEMPLATES)

# Function to generate one .cpp file for the given (LocalOrdinal,
# GlobalOrdinal, Node) template parameter combination, for run-time
# registration of a Tpetra class or function over those template
# parameters.  This is meant to be called by
# TPETRA_PROCESS_ALL_LGN_TEMPLATES.  This function takes the names
# already mangled, to avoid unnecessary string processing overhead.
FUNCTION(TPETRA_PROCESS_ONE_LGN_TEMPLATE OUTPUT_FILE TEMPLATE_FILE
  CLASS_NAME CLASS_MACRO_NAME
  LO_MANGLED_NAME GO_MANGLED_NAME NT_MANGLED_NAME
  LO_MACRO_NAME GO_MACRO_NAME NT_MACRO_NAME)

  STRING(REPLACE "ETI_LO_GO_NT.tmpl"
    "${CLASS_NAME}_${LO_MACRO_NAME}_${GO_MACRO_NAME}_${NT_MACRO_NAME}.cpp"
    OUT_FILE "${TEMPLATE_FILE}")
  CONFIGURE_FILE("${TEMPLATE_FILE}" "${OUT_FILE}")

  SET(${OUTPUT_FILE} ${OUT_FILE} PARENT_SCOPE)
ENDFUNCTION(TPETRA_PROCESS_ONE_LGN_TEMPLATE)

# Function to generate .cpp files for ETI of a Tpetra class or
# function, over all enabled LocalOrdinal, GlobalOrdinal, and Node
# template parameters.  We generate one .cpp file for each
# (LocalOrdinal, GlobalOrdinal, Node) type combination over which
# Tpetra does ETI.
FUNCTION(TPETRA_PROCESS_ALL_LGN_TEMPLATES OUTPUT_FILES TEMPLATE_FILE
  CLASS_NAME CLASS_MACRO_NAME LOCALORDINAL_TYPES GLOBALORDINAL_TYPES
  NODE_TYPES)

  SET(OUT_FILES "")
  FOREACH(NT ${NODE_TYPES})
    TPETRA_MANGLE_TEMPLATE_PARAMETER(NT_MANGLED "${NT}")
    TPETRA_NODE_MACRO_NAME(NT_MACRO_NAME "${NT}")
    FOREACH(GO ${GLOBALORDINAL_TYPES})
      TPETRA_MANGLE_TEMPLATE_PARAMETER(GO_MANGLED "${GO}")
      TPETRA_SLG_MACRO_NAME(GO_MACRO_NAME "${GO}")
      FOREACH(LO ${LOCALORDINAL_TYPES})
        TPETRA_MANGLE_TEMPLATE_PARAMETER(LO_MANGLED "${LO}")
        TPETRA_SLG_MACRO_NAME(LO_MACRO_NAME "${LO}")

        TPETRA_PROCESS_ONE_LGN_TEMPLATE(OUT_FILE "${TEMPLATE_FILE}"
          "${CLASS_NAME}" "${CLASS_MACRO_NAME}"
          "${LO_MANGLED}" "${GO_MANGLED}" "${NT_MANGLED}"
          "${LO_MACRO_NAME}" "${GO_MACRO_NAME}" "${NT_MACRO_NAME}")
        LIST(APPEND OUT_FILES ${OUT_FILE})
      ENDFOREACH() # LO
    ENDFOREACH() # GO
  ENDFOREACH() # NT

  # This is the standard CMake idiom for setting an output variable so
  # that the caller can see the result.
  SET(${OUTPUT_FILES} ${OUT_FILES} PARENT_SCOPE)
ENDFUNCTION(TPETRA_PROCESS_ALL_LGN_TEMPLATES)

# Function to generate one .cpp file for the given (Scalar, Node)
# template parameter combination, for run-time registration of a
# Tpetra class or function over those template parameters.  This is
# meant to be called by TPETRA_PROCESS_ALL_SN_TEMPLATES.  This
# function takes the names already mangled, to avoid unnecessary
# string processing overhead.
#
# OUTPUT_FILE [out] Name of the generated .cpp file.
#
# TEMPLATE_FILE [in] Name of the input .tmpl "template" file.  This
#   function does string substitution in that file, using the input
#   arguments of this function.  For example, @SC_MACRO_EXPR@ (Scalar
#   macro expression) gets substituted for the value of this
#   function's SC_MACRO_EXPR input argument.
#
# CLASS_NAME [in] Name of the Tpetra class (without namespace
#   qualifiers; must live in the Tpetra namespace)
#
# CLASS_MACRO_NAME [in] Name of the Tpetra class, suitably mangled for
#   use in a macro name.
#
# SC_MANGLED_NAME [in] Name of the Scalar (SC) type, mangled for use
#   as a macro argument (e.g., spaces and colons removed).  In the
#   arguments that follow, LO stands for LocalOrdinal, GO for
#   GlobalOrdinal, and NT for Node.
#
# SC_MACRO_EXPR [in] Expression that asks Tpetra whether the given
#   Scalar (SC) type is supported.
#
FUNCTION(TPETRA_PROCESS_ONE_SN_TEMPLATE OUTPUT_FILE TEMPLATE_FILE
  CLASS_NAME CLASS_MACRO_NAME SC_MANGLED_NAME
  GO_MANGLED_NAME NT_MANGLED_NAME SC_MACRO_EXPR
  GO_MACRO_EXPR GO_MACRO_NAME NT_MACRO_NAME)

  STRING(REPLACE "ETI_SC_NT.tmpl"
    "${CLASS_NAME}_${SC_MACRO_NAME}_${NT_MACRO_NAME}.cpp"
    OUT_FILE "${TEMPLATE_FILE}")
  CONFIGURE_FILE("${TEMPLATE_FILE}" "${OUT_FILE}")

  SET(${OUTPUT_FILE} ${OUT_FILE} PARENT_SCOPE)
ENDFUNCTION(TPETRA_PROCESS_ONE_SN_TEMPLATE)

# Function to generate .cpp files for ETI of a Tpetra class, over all
# enabled Scalar and Node template parameters.  We generate one .cpp
# file for each (Scalar, Node) type combination over which Tpetra does
# ETI.
#
# OUTPUT_FILES [out] List of the generated .cpp files.
#
# TEMPLATE_FILE [in] Name of the input .tmpl "template" file.  This
#   function does string substitution in that file, using the input
#   arguments of this function.  For example, @SC_MACRO_EXPR@ (Scalar
#   macro expression) gets substituted for the value of this
#   function's SC_MACRO_EXPR input argument.
#
# CLASS_NAME [in] Name of the Tpetra class (without namespace
#   qualifiers; must live in the Tpetra namespace)
#
# CLASS_MACRO_NAME [in] Name of the Tpetra class, suitably mangled for
#   use in a macro name.
#
# SCALAR_TYPES [in] All Scalar types over which to do ETI for the given
#   class.  This may include Scalar = GlobalOrdinal and/or Scalar =
#   int, if appropriate for that class.
#
# NODE_TYPES [in] All Node types over which to do ETI for the given
#   class.
#
# INSTANTIATE_SCALAR_ORDINAL_TYPES [in] Whether to instantiate for
#   Scalar=GlobalOrdinal types.
#
# MUST_HAVE_SCALAR_INT [in] (Boolean) Whether the class must be
#   instantiated with Scalar = int, even if int is not in the set of
#   GlobalOrdinal types.
FUNCTION(TPETRA_PROCESS_ALL_SN_TEMPLATES OUTPUT_FILES TEMPLATE_FILE
  CLASS_NAME CLASS_MACRO_NAME SCALAR_TYPES GLOBALORDINAL_TYPES
  NODE_TYPES INSTANTIATE_SCALAR_ORDINAL_TYPES MUST_HAVE_SCALAR_INT)

  SET(OUT_FILES "")
  FOREACH(NT ${NODE_TYPES})
    TPETRA_MANGLE_TEMPLATE_PARAMETER(NT_MANGLED "${NT}")
    TPETRA_NODE_MACRO_NAME(NT_MACRO_NAME "${NT}")

    FOREACH(GO ${GLOBALORDINAL_TYPES})
      TPETRA_MANGLE_TEMPLATE_PARAMETER(GO_MANGLED "${GO}")
      TPETRA_SLG_MACRO_NAME(GO_MACRO_NAME "${GO}")

      FOREACH(SC ${SCALAR_TYPES})
        TPETRA_MANGLE_TEMPLATE_PARAMETER(SC_MANGLED "${SC}")
        TPETRA_SLG_MACRO_NAME(SC_MACRO_NAME "${SC}")
        TPETRA_SC_MACRO_EXPR(SC_MACRO_EXPR "${SC}" "${GO}" "${SC_MACRO_NAME}")

        #MESSAGE(STATUS ">> SC = ${SC}, SC_MACRO_EXPR = ${SC_MACRO_EXPR}")

        # If SC is NOT a GlobalOrdinal type of some kind (not
        # necessarily the current GO), or if it is "int", process
        # it.  Otherwise, then we only have to process it if it
        # equals the current GO.
        TPETRA_SC_IS_GO(IS_GO "${SC}")
        STRING(COMPARE EQUAL "${SC}" "${GO}" IS_CURRENT_GO)
        STRING(COMPARE EQUAL "${SC}" "int" IS_INT)

        IF ((MUST_HAVE_SCALAR_INT AND IS_INT) OR
            (NOT IS_GO OR (IS_CURRENT_GO AND INSTANTIATE_SCALAR_ORDINAL_TYPES)))

          IF(IS_CURRENT_GO AND INSTANTIATE_SCALAR_ORDINAL_TYPES)
            # mfh 24 Jul 2019: This is a bit of a hack, because I'm
            # assuming that LO=int.  We've never done ETI for any
            # other LO types, though.
            SET(GO_MACRO_EXPR "&& defined(HAVE_TPETRA_INST_INT_${GO_MACRO_NAME})")
          ELSE()
            SET(GO_MACRO_EXPR "")
          ENDIF()

          #MESSAGE(STATUS ">> SC = ${SC}, GO = ${GO}, GO_MACRO_EXPR = ${GO_MACRO_EXPR}")

          TPETRA_PROCESS_ONE_SN_TEMPLATE(OUT_FILE "${TEMPLATE_FILE}"
            "${CLASS_NAME}" "${CLASS_MACRO_NAME}" "${SC_MANGLED}"
            "${GO_MANGLED}" "${NT_MANGLED}" "${SC_MACRO_EXPR}"
            "${GO_MACRO_EXPR}" "${GO_MACRO_NAME}" "${NT_MACRO_NAME}")
          LIST(APPEND OUT_FILES ${OUT_FILE})
        ENDIF()
      ENDFOREACH() # SC
    ENDFOREACH() # GO
  ENDFOREACH() # NT

  # This is the standard CMake idiom for setting an output variable so
  # that the caller can see the result.
  SET(${OUTPUT_FILES} ${OUT_FILES} PARENT_SCOPE)
ENDFUNCTION(TPETRA_PROCESS_ALL_SN_TEMPLATES)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
IF (NOT Tpetra_ENABLE_DEPRECATED_CODE OR NOT TpetraCore_ENABLE_Epetra)
  LIST(REMOVE_ITEM SOURCES ${DIR}/Epetra_TsqrMessenger.cpp)
  LIST(REMOVE_ITEM SOURCES ${DIR}/Tpetra_EpetraRowMatrix.cpp)
ENDIF()
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Pull in the Kokkos refactor code.
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_refactor)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_config.h )
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_ETIHelperMacros.h)

IF (${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  SET (MultiVector_ETI_SCALARS ${TpetraCore_ETI_SCALARS})
  IF (NOT Tpetra_INST_INT_INT)
    # GO = int is disabled, so the Scalar = GlobalOrdinal case doesn't
    # cover Scalar = int.  We need Scalar = int for communication of
    # process ranks.
    LIST(APPEND MultiVector_ETI_SCALARS "int")
  ENDIF()

  # Generate ETI .cpp files for Tpetra::MultiVector.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(MULTIVECTOR_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "MultiVector" "MULTIVECTOR" "${MultiVector_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${MULTIVECTOR_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Vector.
  # Use the same set of Scalar types as with MultiVector.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(VECTOR_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "Vector" "VECTOR" "${MultiVector_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${VECTOR_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::FEMultiVector.
  # Use the same set of Scalar types as with MultiVector.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(FEMULTIVECTOR_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "FEMultiVector" "FEMULTIVECTOR" "${MultiVector_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${FEMULTIVECTOR_OUTPUT_FILES})

  # Tpetra::CrsMatrix ETI Support
  ASSERT_DEFINED(Tpetra_ENABLE_Reduced_ETI)
  IF(Tpetra_ENABLE_Reduced_ETI)
    # For Reduced_ETI we don't include int
    MESSAGE(STATUS "Tpetra: Using Reduced ETI")
    SET (CrsMatrix_ETI_SCALARS  ${TpetraCore_ETI_SCALARS})
  ELSE()
    # Zoltan2 wants Scalar = int for CrsMatrix (Bug 6298), and thus
    # indirectly for RowMatrix as well.
    MESSAGE(STATUS "Tpetra: Using Full ETI")
    SET (CrsMatrix_ETI_SCALARS ${MultiVector_ETI_SCALARS})
  ENDIF()
    
  # Generate ETI .cpp files for Tpetra::CrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(CRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "CrsMatrix" "CRSMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${CRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::FECrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(FECRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "FECrsMatrix" "FECRSMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${FECRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::RowMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(ROWMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "RowMatrix" "ROWMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${ROWMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::RowMatrixTransposer.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(ROWMATRIXTRANSPOSER_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "RowMatrixTransposer" "ROWMATRIXTRANSPOSER" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${ROWMATRIXTRANSPOSER_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::CrsGraphTransposer.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(CRSGRAPHTRANSPOSER_OUTPUT_FILES "Tpetra_ETI_LO_GO_NT.tmpl" "CrsGraphTransposer" "CRSGRAPHTRANSPOSER" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${CRSGRAPHTRANSPOSER_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::BlockCrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(BLOCKCRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "BlockCrsMatrix" "BLOCKCRSMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${BLOCKCRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::BlockMultiVector.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(BLOCKMULTIVECTOR_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "BlockMultiVector" "BLOCKMULTIVECTOR" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE )
  LIST(APPEND SOURCES ${BLOCKMULTIVECTOR_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::BlockVector.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(BLOCKVECTOR_OUTPUT_FILES
    "Tpetra_ETI_SC_LO_GO_NT.tmpl" "BlockVector" "BLOCKVECTOR"
    "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${BLOCKVECTOR_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::BlockCrsMatrix helpers.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(BLOCKCRSMATRIX_HELPER_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl"
  "BlockCrsMatrix_Helpers" "BLOCKCRSMATRIX_HELPERS" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${BLOCKCRSMATRIX_HELPER_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Details::getDiagCopyWithoutOffsets.
  # Do so for the same set of Scalar types as CrsMatrix (see above),
  # because this function is an implementation detail of CrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(GETDIAGCOPYWITHOUTOFFSETS_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "Details_getDiagCopyWithoutOffsets" "DETAILS_GETDIAGCOPYWITHOUTOFFSETS" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${GETDIAGCOPYWITHOUTOFFSETS_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Details::packCrsMatrix and
  # Tpetra::Details::packCrsMatrixWithOwningPIDs.  Do so for the same
  # set of Scalar types as CrsMatrix (see above), because this
  # function is an implementation detail of CrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(PACKCRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "Details_packCrsMatrix" "DETAILS_PACKCRSMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${PACKCRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for
  # Tpetra::Details::unpackCrsMatrixAndCombine,
  # Tpetra::Details::unpackAndCombineIntoCrsArrays, and
  # Tpetra::Details::unpackAndCombineWithOwningPIDsCount.  Do so for the same
  # set of Scalar types as CrsMatrix (see above), because this
  # function is an implementation detail of CrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(UNPACKCRSMATRIXANDCOMBINE_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "Details_unpackCrsMatrixAndCombine" "DETAILS_UNPACKCRSMATRIXANDCOMBINE" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${UNPACKCRSMATRIXANDCOMBINE_OUTPUT_FILES})

  # DistObject requires the same Scalar type instantiations as
  # MultiVector, plus Scalar = char (for CrsMatrix).
  SET (DistObject_ETI_SCALARS ${MultiVector_ETI_SCALARS})
  LIST(APPEND DistObject_ETI_SCALARS "char")

  # Generate ETI .cpp files for Tpetra::DistObject.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(DISTOBJECT_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "DistObject" "DISTOBJECT" "${DistObject_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${DISTOBJECT_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::leftAndOrRightScaleCrsMatrix.
  # Do this only for non-integer Scalar types, since we really only
  # need this function for linear solvers.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(LEFTANDORRIGHTSCALECRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "leftAndOrRightScaleCrsMatrix" "LEFTANDORRIGHTSCALECRSMATRIX" "${TpetraCore_ETI_SCALARS_NO_ORDS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" FALSE)
  LIST(APPEND SOURCES ${LEFTANDORRIGHTSCALECRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::computeRowAndColumnOneNorms.
  # Do this only for non-integer Scalar types, since we really only
  # need this function for linear solvers.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(COMPUTEROWANDCOLUMNONENORMS_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "computeRowAndColumnOneNorms" "COMPUTEROWANDCOLUMNONENORMS" "${TpetraCore_ETI_SCALARS_NO_ORDS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" FALSE)
  LIST(APPEND SOURCES ${COMPUTEROWANDCOLUMNONENORMS_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::replaceDiagonalCrsMatrix.
  TPETRA_PROCESS_ALL_SLGN_TEMPLATES(REPLACEDIAGONALCRSMATRIX_OUTPUT_FILES "Tpetra_ETI_SC_LO_GO_NT.tmpl" "replaceDiagonalCrsMatrix" "REPLACEDIAGONALCRSMATRIX" "${CrsMatrix_ETI_SCALARS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" FALSE)
  LIST(APPEND SOURCES ${REPLACEDIAGONALCRSMATRIX_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::CrsMatrix::convert.
  TPETRA_PROCESS_ALL_CONVERT_TEMPLATES(CRSMATRIX_CONVERT_OUTPUT_FILES "Tpetra_ETI_SOUT_SIN_LO_GO_NT.tmpl" "CrsMatrix" "CRSMATRIX_CONVERT" "${TpetraCore_ETI_SCALARS_NO_ORDS}" "${TpetraCore_ETI_SCALARS_NO_ORDS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" FALSE)
  LIST(APPEND SOURCES ${CRSMATRIX_CONVERT_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::MultiVector::convert.
  TPETRA_PROCESS_ALL_CONVERT_TEMPLATES(MULTIVECTOR_CONVERT_OUTPUT_FILES "Tpetra_ETI_SOUT_SIN_LO_GO_NT.tmpl" "MultiVector" "MULTIVECTOR_CONVERT" "${TpetraCore_ETI_SCALARS_NO_ORDS}" "${TpetraCore_ETI_SCALARS_NO_ORDS}" "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" FALSE)
  LIST(APPEND SOURCES ${MULTIVECTOR_CONVERT_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::LocalCrsMatrixOperator.
  TPETRA_PROCESS_ALL_SN_TEMPLATES(LOCALCRSMATRIXOPERATOR_OUTPUT_FILES
    "Tpetra_ETI_SC_NT.tmpl" "LocalCrsMatrixOperator"
    "LOCALCRSMATRIXOPERATOR" "${CrsMatrix_ETI_SCALARS}"
    "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" TRUE TRUE)
  LIST(APPEND SOURCES ${LOCALCRSMATRIXOPERATOR_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::LocalCrsMatrixOperator with mixed Scalars.
  ASSERT_DEFINED (Tpetra_INST_DOUBLE)
  ASSERT_DEFINED (Tpetra_INST_FLOAT)
  IF (Tpetra_INST_DOUBLE AND Tpetra_INST_FLOAT)
    TPETRA_PROCESS_ALL_CONVERT_TEMPLATES(LOCALCRSMATRIXOPERATOR_MIXED_OUTPUT_FILES
      "Tpetra_ETI_SOUT_SIN_LO_GO_NT.tmpl" "LocalCrsMatrixOperator"
      "LOCALCRSMATRIXOPERATOR_MIXED"
      "double" "float"
      "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" True True)
    LIST(APPEND SOURCES ${LOCALCRSMATRIXOPERATOR_MIXED_OUTPUT_FILES})
  ENDIF()
  ASSERT_DEFINED (Tpetra_INST_COMPLEX_DOUBLE)
  ASSERT_DEFINED (Tpetra_INST_COMPLEX_FLOAT)
  IF (Tpetra_INST_COMPLEX_DOUBLE AND Tpetra_INST_COMPLEX_FLOAT)
    TPETRA_PROCESS_ALL_CONVERT_TEMPLATES(LOCALCRSMATRIXOPERATOR_MIXED_OUTPUT_FILES
      "Tpetra_ETI_SOUT_SIN_LO_GO_NT.tmpl" "LocalCrsMatrixOperator"
      "LOCALCRSMATRIXOPERATOR_MIXED"
      "std::complex<double>" "std::complex<float>"
      "${TpetraCore_ETI_LORDS}" "${TpetraCore_ETI_GORDS}" "${TpetraCore_ETI_NODES}" True True)
    LIST(APPEND SOURCES ${LOCALCRSMATRIXOPERATOR_MIXED_OUTPUT_FILES})
  ENDIF()

  # Generate ETI .cpp files for Tpetra::CrsGraph.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(CRSGRAPH_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "CrsGraph"
    "CRSGRAPH"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${CRSGRAPH_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::FECrsGraph.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(FECRSGRAPH_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "FECrsGraph"
    "FECRSGRAPH"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${FECRSGRAPH_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Directory.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(DIRECTORY_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "Directory"
    "DIRECTORY"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${DIRECTORY_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Directory's implementation classes.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(DIRECTORYIMPL_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "DirectoryImpl"
    "DIRECTORYIMPL"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${DIRECTORYIMPL_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Export.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(EXPORT_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "Export"
    "EXPORT"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${EXPORT_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Import.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(IMPORT_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "Import"
    "IMPORT"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${IMPORT_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::ImportExportData.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(IMPORTEXPORTDATA_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "ImportExportData"
    "IMPORTEXPORTDATA"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${IMPORTEXPORTDATA_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Details::makeColMap.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(MAKECOLMAP_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "Details_makeColMap"
    "DETAILS_MAKECOLMAP"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${MAKECOLMAP_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Details::packCrsGraph.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(PACKCRSGRAPH_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "Details_packCrsGraph"
    "DETAILS_PACKCRSGRAPH"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${PACKCRSGRAPH_OUTPUT_FILES})

  # Generate ETI .cpp files for Tpetra::Details::unpackCrsGraphAndCombine.
  TPETRA_PROCESS_ALL_LGN_TEMPLATES(UNPACKCRSGRAPHANDCOMBINE_OUTPUT_FILES
    "Tpetra_ETI_LO_GO_NT.tmpl"
    "Details_unpackCrsGraphAndCombine"
    "DETAILS_UNPACKCRSGRAPHANDCOMBINE"
    "${TpetraCore_ETI_LORDS}"
    "${TpetraCore_ETI_GORDS}"
    "${TpetraCore_ETI_NODES}")
  LIST(APPEND SOURCES ${UNPACKCRSGRAPHANDCOMBINE_OUTPUT_FILES})

  #MESSAGE(STATUS "SOURCES = ${SOURCES}")
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  tpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  ADDED_LIB_TARGET_NAME_OUT TPETRA_LIBNAME
  )

# We need to set the linker language explicitly here for CUDA builds.
SET_PROPERTY(
  TARGET ${TPETRA_LIBNAME}
  APPEND PROPERTY LINKER_LANGUAGE CXX
  )

#
# Make a trivial change to this comment if you add / remove a file either to
# / from this directory, or to / from the 'impl' subdirectory.  That ensures
# that running "make" will also rerun CMake in order to regenerate Makefiles.
#
# Here's another change.  Again.

