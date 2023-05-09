# FIXME (mfh 11 Feb 2015) I am removing KokkosClassic::ThrustGPUNode,
# so the example in this subdirectory won't build.  We need to test
# this example, though!  We should figure out how to test it in a
# sensible way without needing ThrustGPUNode.
#
#IF(${PACKAGE_NAME}_ENABLE_Ifpack2 AND
#   TpetraClassic_ENABLE_ThrustGPUNode AND TpetraClassic_ENABLE_TPINode)
#  ADD_SUBDIRECTORY(clone)
#ENDIF()

ADD_SUBDIRECTORY(multiplesolve)
ADD_SUBDIRECTORY(memory)
ADD_SUBDIRECTORY(separatesetups)
ADD_SUBDIRECTORY(levelwrap)
ADD_SUBDIRECTORY(blockcrs)
