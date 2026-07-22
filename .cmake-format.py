# https://cmake-format.readthedocs.io/en/latest/

# TRIBITS commands
_exec_kwargs = {"SOURCES": "*",
                "NOEXEPREFIX": "*",
                "NOEXESUFFIX": "*",
                "DIRECTORY": "*",
                "COMM": "*",
                "CATEGORIES": "*",
                "INSTALLABLE": "*",
                "DEPLIBS": "*",
                "TARGET_DEFINES": "*"}

_cmdline = {"pargs": {"nargs": "+", "tags": ["cmdline"]}}
_filelist = {"pargs": {"nargs": "+", "sortable": True}}

_test_kwargs = {"NAME": "*",
                "NAME_POSTFIX": "*",
                "NUM_MPI_PROCS": "*",
                "COMM": "*",
                "NUM_MPI_PROCS": "*",
                "ARGS": _cmdline,
                "POSTFIX_AND_ARGS_0": _cmdline,
                "POSTFIX_AND_ARGS_1": _cmdline,
                "POSTFIX_AND_ARGS_2": _cmdline,
                "POSTFIX_AND_ARGS_3": _cmdline,
                "POSTFIX_AND_ARGS_4": _cmdline,
                "POSTFIX_AND_ARGS_5": _cmdline,
                "POSTFIX_AND_ARGS_6": _cmdline,
                "POSTFIX_AND_ARGS_7": _cmdline,
                "POSTFIX_AND_ARGS_8": _cmdline,
                "POSTFIX_AND_ARGS_9": _cmdline,
                "CATEGORIES": "*",
                "RUN_SERIAL": "*",
                "STANDARD_PASS_OUTPUT": "*",
                "PASS_REGULAR_EXPRESSION": "*",
                "FAIL_REGULAR_EXPRESSION": "*",
                "ENVIRONMENT": "*",
                "WILL_FAIL": "*",
                "EXEC": "*",
                "CMND": "*",
                "OVERALL_NUM_MPI_PROCS": "*",
                "TEST_0": "*",
                "TEST_1": "*",
                "TEST_2": "*",
                "TEST_3": "*"}

_lib_kwargs = {"HEADERS": "*",
               "SOURCES": "*",
               "DEPLIBS": "*",
               "ADDED_LIB_TARGET_NAME_OUT": "*"
               }


# -----------------------------
# Options affecting parsing.
# -----------------------------
with section("parse"):
    additional_commands = {"TRIBITS_ADD_TEST": {"kwargs": _test_kwargs,
                                                "flags": ["NOEXESUFFIX"]},
                           "TRIBITS_ADD_ADVANCED_TEST": {"kwargs": _test_kwargs,
                                                      "flags": ["NOEXEPREFIX", "NOEXESUFFIX"]},
                           "TRIBITS_ADD_EXECUTABLE_AND_TEST": {"kwargs": {**_test_kwargs, **_exec_kwargs},
                                                               "flags": ["NOEXEPREFIX", "NOEXESUFFIX"]},
                           "MUELU_ADD_SERIAL_AND_MPI_TEST": {"kwargs": _test_kwargs,
                                                      "flags": ["NOEXEPREFIX", "NOEXESUFFIX"]},
                           "TRIBITS_ADD_EXECUTABLE": {"kwargs": _exec_kwargs,
                                                      "flags": ["NOEXEPREFIX", "NOEXESUFFIX"]},
                           "TRIBITS_ADD_LIBRARY": {"kwargs": _lib_kwargs},
                           "TRIBITS_COPY_FILES_TO_BINARY_DIR": {"kwargs": {"SOURCE_FILES": _filelist,
                                                                           "SOURCE_DIR": "*",
                                                                           "DEST_FILES": _filelist,
                                                                           "DEST_DIR": "*",
                                                                           "TARGETDEPS": "*",
                                                                           "EXEDEPS": "*",
                                                                           "CATEGORIES": "*",
                                                                           },
                                                                "flags": ["NOEXEPREFIX", "NOEXESUFFIX"]},
                           "TRIBITS_ADD_OPTION_AND_DEFINE": {"pargs": {"nargs": 4}},
                           }

# -----------------------------
# Options affecting formatting.
# -----------------------------
with section("format"):

    # How wide to allow formatted cmake files
    line_width = 200

    # If an argument group contains more than this many sub-groups (parg or kwarg
    # groups) then force it to a vertical layout.
    max_subgroups_hwrap = 2

    # If a positional argument group contains more than this many arguments, then
    # force it to a vertical layout.
    max_pargs_hwrap = 3

    # If a statement is wrapped to more than one line, than dangle the closing
    # parenthesis on its own line.
    dangle_parens = True

    # If the trailing parenthesis must be 'dangled' on its on line, then align it
    # to this reference: `prefix`: the start of the statement,  `prefix-indent`:
    # the start of the statement, plus one indentation  level, `child`: align to
    # the column of the arguments
    dangle_align = 'prefix'

    # Format command names consistently as 'lower' or 'upper' case
    command_case = 'upper'

    # If true, separate function names from parentheses with a space
    separate_fn_name_with_space = False

    # If true, separate flow control names from their parentheses with a space
    separate_ctrl_name_with_space = False

    # If true, the argument lists which are known to be sortable will be sorted lexicographically.
    enable_sort = True

    # If true, the parsers may infer whether or not an argument list is sortable (without annotation).
    autosort = True


# ------------------------------------------------
# Options affecting comment reflow and formatting.
# ------------------------------------------------
with section("markup"):
    # enable comment markup parsing and reflow
    enable_markup = False
