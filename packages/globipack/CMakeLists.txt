

TRIBITS_PACKAGE( GlobiPack ENABLE_SHADOWING_WARNINGS CLEANED )

# GlobiPack package is deprecated April 2019; to be removed May 2019
MESSAGE(WARNING "***** WARNING:  Package GlobiPack is deprecated in April 2019; it will be removed in May 2019 *****")

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_PACKAGE_POSTPROCESS()
