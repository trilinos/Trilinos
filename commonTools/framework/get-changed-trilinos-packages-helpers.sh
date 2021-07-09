#
# Determine paths
#

if [ "$TRILINOS_DIR" == "" ] ; then
  # Grab from the symlink (only works on Linux)
  _ABS_FILE_PATH=`readlink -f $0` || \
   echo "Could not follow symlink to set TRILINOS_DIR!"
  if [ "$_ABS_FILE_PATH" != "" ] ; then
    _SCRIPT_DIR=`dirname $_ABS_FILE_PATH`
    TRILINOS_DIR=$_SCRIPT_DIR/../..
  fi
fi

if [ "$TRILINOS_DIR" == "" ] ; then
  echo "ERROR: Cannot determine TRILINOS_DIR!  Please set env var TRILINOS_DIR!"
  exit 4
fi

echo "TRILINOS_DIR=$TRILINOS_DIR"

# Allow a different source tree for the Trilinos scripts
if [ "$TRILINOS_SCRIPTS_DIR" == "" ] ; then
  TRILINOS_SCRIPTS_DIR=${TRILINOS_DIR}
fi

echo "TRILINOS_SCRIPTS_DIR=$TRILINOS_SCRIPTS_DIR"

ORIG_CWD=$PWD

# Allow override of TriBITS for testing purposes
if [ "${GCTP_TRIBITS_DIR_OVERRIDE}" != "" ] ; then
  TRIBITS_DIR=$GCTP_TRIBITS_DIR_OVERRIDE
else
  TRIBITS_DIR=$TRILINOS_SCRIPTS_DIR/cmake/tribits
fi
echo "TRIBITS_DIR=$TRIBITS_DIR"


#
# Functions
#


function comma_list_to_list() {
  echo "$1" | sed "s|,| |g"
}


function list_to_comma_list() {
  echo "$@" | sed "s| |,|g"
}


# list_contains_ele <ele> <a0> <a1> ...
# 
# Returns if <ele> is contained in list (0 is success)
function list_contains_ele() {
  ele_to_find="$1" ; shift
  list="$@"
  #echo "ele_to_find='${ele_to_find}'"
  for ELE in ${list} ; do
    #echo "ELE='${ELE}'"
    if [[ "${ELE}" == "${ele_to_find}" ]] ; then
      #echo "Contains ${ele_to_find}!"
      return 0
    fi
  done
  return 1
}


# comma_list_contains_ele <ele> <a0>,<a1>,...
# 
# Returns if <ele> is contained in list (0 is success)
function comma_list_contains_ele() {
  ele="$1"
  comma_list="$2"
  list=$(comma_list_to_list "${comma_list}")
  if list_contains_ele "${ele}" "${list}"; then
    return 0
  else
    return 1
  fi
}


# Generates TrilinosPackageDependencies.xml
function generate_trilinos_package_dependencies_xml_file() {
  cmake \
    -D Trilinos_DEPS_XML_OUTPUT_FILE=TrilinosPackageDependencies.xml \
    -P $TRIBITS_DIR/ci_support/TribitsDumpDepsXmlScript.cmake \
    &> TribitsDumpDepsXmlScript.log
  echo "Wrote the file 'TrilinosPackageDependencies.xml'"
}


# Take in and return a filtered comma-seprated list of packages
#
# Non PT and ST packages as well as packages listed in
# TRILINOS_EXCLUDE_PACKAGES_FROM_PR_TESTING are filtered out.
function trilinos_filter_packages_to_test() {
  input_packages_comma_list="$1"
  fullFilteredPackagesCommaList=$(
    ${TRIBITS_DIR}/ci_support/filter-packages-list.py \
    --deps-xml-file=TrilinosPackageDependencies.xml \
    --input-packages-list="${input_packages_comma_list}" \
    --keep-test-test-categories=PT,ST)
  fullFilteredPackagesList=$(comma_list_to_list "$fullFilteredPackagesCommaList")
  filteredPackagesList=()
  for pkg in ${fullFilteredPackagesList} ; do
    if ! list_contains_ele "${pkg}" "${TRILINOS_EXCLUDE_PACKAGES_FROM_PR_TESTING[@]}";then
      filteredPackagesList+=($pkg)
    fi
  done
  list_to_comma_list "${filteredPackagesList[@]}"
}


function trilinos_get_all_toplevel_packages() {
  $TRIBITS_DIR/ci_support/get-tribits-packages.py \
    --deps-xml-file=TrilinosPackageDependencies.xml
}
