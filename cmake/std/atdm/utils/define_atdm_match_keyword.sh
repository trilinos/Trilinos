function atdm_match_keyword() {
  input_string=$1
  keyword=$2

  lower_input_string=$(echo "$input_string" | tr '[:upper:]' '[:lower:]')
  lower_keyword=$(echo "$keyword" | tr '[:upper:]' '[:lower:]')

  #echo "lower_input_string='${lower_input_string}'"
  #echo "lower_keyword='${lower_keyword}'"

  if   [[ "${lower_input_string}" == "${lower_keyword}" ]] ; then
    return 0
  elif [[ "${lower_input_string}" == *"-${lower_keyword}" ]] ; then
    return 0
  elif [[ "${lower_input_string}" == *"_${lower_keyword}" ]] ; then
    return 0
  elif [[ "${lower_input_string}" == "${lower_keyword}-"* ]] ; then
    return 0
  elif [[ "${lower_input_string}" == "${lower_keyword}_"* ]] ; then
    return 0
  elif [[ "${lower_input_string}" == *"-${lower_keyword}-"* ]] ; then
    return 0
  elif [[ "${lower_input_string}" == *"-${lower_keyword}_"* ]] ; then
    return 0
  elif [[ "${lower_input_string}" == *"_${lower_keyword}-"* ]] ; then
    return 0
  elif [[ "${lower_input_string}" == *"_${lower_keyword}_"* ]] ; then
    return 0
  fi

  return 1

}


function atdm_match_any_keyword() {
  input_string=$1 ; shift
  keyword_array=$@

  #echo "input_string=$input_string"
  #echo "keyword_array=$keyword_array"

  for keyword in ${keyword_array[@]} ; do
    if atdm_match_keyword "${input_string}" ${keyword}; then
      return 0
    fi
  done

  return 1
}


function atdm_match_buildname_keyword() {
  keyword=$1
  atdm_match_keyword "$ATDM_CONFIG_BUILD_NAME" ${keyword}
}


function atdm_match_any_buildname_keyword() {
  atdm_match_any_keyword "$ATDM_CONFIG_BUILD_NAME" $@
}
