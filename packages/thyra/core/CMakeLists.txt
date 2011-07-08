INCLUDE(SubPackageMacros)

SUBPACKAGE(Core)

ADD_SUBDIRECTORY(src)

# Add the targets for building the Thyra test library so that other packages
# can use the same test library without requiring Thyra tests be enabled.
ADD_SUBDIRECTORY(test/nonlinear/models)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_TEST_DIRECTORIES(example)

SUBPACKAGE_POSTPROCESS()
