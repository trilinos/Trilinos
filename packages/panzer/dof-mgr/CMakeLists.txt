TRIBITS_SUBPACKAGE(DofMgr)

# No CMake check that Phalanx' and Tpetra's spaces agree.  CMake cannot
# evaluate C++ types, so any check here can only approximate: comparing the
# option values gets the booleans wrong (ON and TRUE are unequal strings, as
# are OFF and FALSE), and gating on a GPU backend being enabled still rejects a
# host execution space on a GPU build, where Tpetra's host node types ignore
# Tpetra_ALLOCATE_IN_SHARED_SPACE entirely and there is nothing to conflict.
# The static_assert in Panzer_NodeType.hpp compares the actual node types and
# is exact.  See doc/design_notes/PanzerTpetraNodeTypePlan.txt, decision N6.

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_SUBPACKAGE_POSTPROCESS()

# TOUCH (to trigger reconfigure)
