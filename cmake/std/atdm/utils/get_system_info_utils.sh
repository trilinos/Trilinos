# Utilities for system matching

function get_knownSystemNameInBuildName() {
  knownSystemNameInBuildName=
  for known_build_name in ${ATDM_KNOWN_SYSTEM_NAMES_LIST[@]} ; do
    if [[ $ATDM_CONFIG_BUILD_NAME == *"${known_build_name}"* ]] ; then
      knownSystemNameInBuildName=${known_build_name}
      break
    fi
  done
  echo ${knownSystemNameInBuildName}
}


function list_contains() {
  list_name=$1  # Name of bash array of items is being searched
  item=$2   # Item we are looking for in that array
  list=$(eval echo \${${list_name}[@]})  # Get the array elements
  for list_item in ${list}; do
    if [[ "${list_item}" == "${item}" ]] ; then
      #echo "List contains '${item}'"
      return 0
    fi
  done
  return 1
}


function assert_selected_system_matches_known_system_type_matches() {

  if [[ "${systemNameTypeMatchedList}" == "" ]] ; then
    return 0
  fi

  local selectedSystemNameMatchesSystemTypeMatches=false
  list_contains systemNameTypeMatchedList ${ATDM_SYSTEM_NAME} \
    && selectedSystemNameMatchesSystemTypeMatches=true

  if [[ "${selectedSystemNameMatchesSystemTypeMatches}" == "false" ]]; then
    echo
    echo "***"
    echo "*** Error, the system name '${knownSystemNameInBuildName}' given in the build name:"
    echo "***"
    echo "***   ${ATDM_CONFIG_BUILD_NAME}"
    echo "***"
    echo "*** does not match the current host:"
    echo "***"
    echo "***   $realHostname"
    echo "***"
    echo "**** and does not match one of the supported system types on this machine which includes:"
    echo "***"
    echo "***   (${systemNameTypeMatchedList[@]})"
    echo "***"
    echo "*** To address this, either remove '${knownSystemNameInBuildName}' from the build name"
    echo "*** or change it to one of the above suppported system types."
    echo "***"
    echo
    return 1
  fi
  return 0

}
