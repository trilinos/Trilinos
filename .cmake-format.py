# TRIBITS commands

_exec_kwargs = {"SOURCES": ".*",
                "NOEXEPREFIX": ".*",
                "NOEXESUFFIX": ".*",
                "DIRECTORY": ".*",
                "COMM": ".*",
                "CATEGORIES": ".*",
                "INSTALLABLE": ".*",
                "DEPLIBS": ".*"}

_test_kwargs = {"NAME": ".*",
                "NAME_POSTFIX": ".*",
                "NUM_MPI_PROCS": ".*",
                "COMM": ".*",
                "NUM_MPI_PROCS": ".*",
                "ARGS": ".*",
                "POSTFIX_AND_ARGS_0": ".*",
                "POSTFIX_AND_ARGS_1": ".*",
                "POSTFIX_AND_ARGS_2": ".*",
                "POSTFIX_AND_ARGS_3": ".*",
                "POSTFIX_AND_ARGS_4": ".*",
                "POSTFIX_AND_ARGS_5": ".*",
                "POSTFIX_AND_ARGS_6": ".*",
                "POSTFIX_AND_ARGS_7": ".*",
                "POSTFIX_AND_ARGS_8": ".*",
                "CATEGORIES": ".*",
                "RUN_SERIAL": ".*",
                "STANDARD_PASS_OUTPUT": ".*",
                "PASS_REGULAR_EXPRESSION": ".*",
                "FAIL_REGULAR_EXPRESSION": ".*",
                "ENVIRONMENT": ".*",
                "WILL_FAIL": ".*",
                "EXEC": ".*",
                "CMND": ".*",
                "OVERALL_NUM_MPI_PROCS": ".*",
                "TEST_0": ".*",
                "TEST_1": ".*",
                "TEST_2": ".*",
                "TEST_3": ".*"}

_lib_kwargs = {"HEADERS": ".*",
               "SOURCES": ".*"}


# -----------------------------
# Options affecting parsing.
# -----------------------------
with section("parse"):
    additional_commands = {"TRIBITS_ADD_TEST": {"kwargs": _test_kwargs},
                           "TRIBITS_ADD_ADVANCED_TEST": {"kwargs": _test_kwargs},
                           "TRIBITS_ADD_EXECUTABLE_AND_TEST": {"kwargs": {**_test_kwargs, **_exec_kwargs}},
                           "MUELU_ADD_SERIAL_AND_MPI_TEST": {"kwargs": _test_kwargs},
                           "TRIBITS_ADD_EXECUTABLE": {"kwargs": _exec_kwargs},
                           "TRIBITS_ADD_LIBRARY": {"kwargs": _lib_kwargs},
                           "TRIBITS_COPY_FILES_TO_BINARY_DIR": {"kwargs": {"SOURCE_FILES": ".*",
                                                                           "SOURCE_DIR": ".*",
                                                                           "DEST_DIR": ".*"}}
                           }

# -----------------------------
# Options affecting formatting.
# -----------------------------
with section("format"):

    # How wide to allow formatted cmake files
    line_width = 200

    # If an argument group contains more than this many sub-groups (parg or kwarg
    # groups) then force it to a vertical layout.
    max_subgroups_hwrap = 3

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

    # If an argument group contains more than this many sub-groups (parg or kwarg
    # groups) then force it to a vertical layout.
    max_subgroups_hwrap = 1

    # If a positional argument group contains more than this many arguments, then
    # force it to a vertical layout.
    max_pargs_hwrap = 20

    # If an argument group contains more than this many sub-groups (parg or kwarg
    # groups) then force it to a vertical layout.
    max_subgroups_hwrap = 2


# ------------------------------------------------
# Options affecting comment reflow and formatting.
# ------------------------------------------------
with section("markup"):
    # enable comment markup parsing and reflow
    enable_markup = False
