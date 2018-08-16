# Sumbit to Trilinos CDash server on casl-dev until can
# address VRI Kanban story #2404
#SET_DEFAULT(CTEST_DROP_SITE "casl-dev.ornl.gov")
SET_DEFAULT(CTEST_DROP_SITE "testing-vm.sandia.gov")
SET_DEFAULT(CTEST_DROP_LOCATION "/cdash/submit.php?project=Trilinos")

# Must overridde the Trilinos defaults to send to /extended/cdash
SET_DEFAULT(CTEST_DROP_SITE_COVERAGE ${CTEST_DROP_SITE})
SET_DEFAULT(CTEST_DROP_LOCATION_COVERAGE ${CTEST_DROP_LOCATION})
