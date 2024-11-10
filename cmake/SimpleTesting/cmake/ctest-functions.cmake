include_guard()
message("+--------------------------------------+")
message("| ctest-functions.cmake START          |")
message("+--------------------------------------+")



macro(banner arg_banner_message)
    message("+----------------------------------------------------------+")
    message("+ ${arg_banner_message}")
    message("+----------------------------------------------------------+")
endmacro()



macro(submit_by_parts arg_parts_value)
    banner("submit_by_parts() START")
    message(">>> arg_parts_value: ${arg_parts_value}")
    message(">>> CTEST_DROP_METHOD : ${CTEST_DROP_METHOD}")
    message(">>> CTEST_DROP_LOCATION : ${CTEST_DROP_LOCATION}")
    message(">>> CDash URL1 = ${build_url1}")
    message(">>> CDash URL2 = ${build_url2}")
    message(">>> CDash URL3 = ${build_url3}")
    if(NOT skip_by_parts_submit)
        message(">>> ctest_submit(PARTS        ${arg_parts_value}")
        message("                 RETRY_COUNT  ${ctest_submit_retry_count}")
        message("                 RETRY_DELAY  ${ctest_submit_retry_delay}")
        message("                 RETURN_VALUE ctest_submit_error)")

        # https://cmake.org/cmake/help/latest/command/ctest_submit.html
        ctest_submit(PARTS ${arg_parts_value}
                     RETRY_COUNT ${ctest_submit_retry_count}
                     RETRY_DELAY ${ctest_submit_retry_delay}
                     BUILD_ID    cdash_build_id
                     RETURN_VALUE ctest_submit_error)

        if(ctest_submit_error EQUAL 0)
            message(">>> ${arg_parts_value} submit: OK")
            message(">>> CDash Build-ID : ${cdash_build_id}")
        else()
            message(">>> ${arg_parts_value} submit: FAILED")
            message(">>> - The ERROR code is ${ctest_submit_error}")
        endif()
    else()
        message(">>> SKIPPED")
        message(">>> skip_by_parts_submit : ${skip_by_parts_submit}")
    endif()
    banner("submit_by_parts() FINISH")
endmacro()



macro(submit_single_submit)
    banner("submit_single_submit() START")
    if(NOT skip_single_submit)
        message(">>> ctest_submit(RETRY_COUNT  ${ctest_submit_retry_count}")
        message("                 RETRY_DELAY  ${ctest_submit_retry_delay}")
        message("                 RETURN_VALUE error_code)")

        ctest_submit(RETRY_COUNT ${ctest_submit_retry_count}
                     RETRY_DELAY ${ctest_submit_retry_delay}
                     RETURN_VALUE error_code)

        if(error_code EQUAL 0)
            message(">>> Single submit: OK")
        else()
            message(">>> Single submit: FAILED")
            message(">>> - The ERROR code is ${error_code}")
        endif()
    else()
        message(">>> SKIPPED")
        message(">>> skip_single_submit : ${skip_single_submit}")
    endif()
    banner("submit_single_submit() FINISH")
endmacro()



macro(submit_upload_config_files)
    banner("submit_upload_config_files() START")
    if( NOT skip_upload_config_files )
        if( NOT (skip_single_submit AND skip_by_parts_submit) )
            message(">>> ctest_upload(FILES ${configure_command_file}")
            message("                 ${configure_file}")
            message("                 ${package_enables_file}")
            message("                 ${genconfig_build_name_file})")

            ctest_upload(FILES ${configure_command_file}
                               ${configure_file}
                               ${package_enables_file}
                               ${genconfig_build_name_file})

            message(">>> ctest_submit(PARTS upload")
            message("                 RETRY_COUNT  ${ctest_submit_retry_count}")
            message("                 RETRY_DELAY  ${ctest_submit_retry_delay}")
            message("                 RETURN_VALUE file_upload_erorr)")

            ctest_submit(PARTS Upload
                         RETRY_COUNT ${ctest_submit_retry_count}
                         RETRY_DELAY ${ctest_submit_retry_delay}
                         RETURN_VALUE file_upload_error)

            if(file_upload_error EQUAL 0)
                message(">>> Config Files Upload: OK")
            else()
                message(">>> Config Files Upload: FAILED")
                message(">>> - The ERROR code is ${file_upload_error}")
            endif()
        else()
            message(">>> SKIPPED")
            message(">>> skip_single_submit   : ${skip_single_submit}")
            message(">>> skip_by_parts_submit : ${skip_by_parts_submit}")
        endif()
    else()
        message(">>> SKIPPED")
        message(">>> skip_upload_config_files : ${skip_upload_config_files}")
    endif()
    banner("submit_upload_config_files() FINISH")
endmacro()



macro(print_options_list)
    banner("OPTIONS")
    message(">>> CTEST_BUILD_NAME         = ${CTEST_BUILD_NAME}")
    message(">>> PARALLEL_LEVEL           = ${PARALLEL_LEVEL}")
    message(">>> TEST_PARALLEL_LEVEL      = ${TEST_PARALLEL_LEVEL}")
    message(">>> skip_by_parts_submit     = ${skip_by_parts_submit}")
    message(">>> skip_single_submit       = ${skip_single_submit}")
    message(">>> skip_update_step         = ${skip_update_step}")
    message(">>> skip_upload_config_files = ${skip_upload_config_files}")
    message(">>> skip_clean_build_dir     = ${skip_clean_build_dir}")
    message(">>> subproject_count         = ${subproject_count}")
    message(">>> dashboard_model          = ${dashboard_model}")
    message(">>> dashboard_track          = ${dashboard_track}")
    message(">>> genconfig_build_name_file= ${genconfig_build_name_file}")
    message(">>> configure_command_file   = ${configure_command_file}")
    message(">>> configure_file           = ${configure_file}")
    message(">>> build_root               = ${build_root}")
    message(">>> build_dir                = ${build_dir}")
    message(">>> SOURCE_DIR               = ${SOURCE_DIR}")
    message(">>> CMAKE_CURRENT_LIST_DIR   = ${CMAKE_CURRENT_LIST_DIR}")
    message(">>> CTEST_SOURCE_DIRECTORY   = ${CTEST_SOURCE_DIRECTORY}")
    message(">>> CTEST_BINARY_DIRECTORY   = ${CTEST_BINARY_DIRECTORY}")
endmacro()



function(generate_build_url1 url_output cdash_site cdash_location project_name build_name build_stamp machine_name)
    banner("generate_build_url1() START")
    message(">>> cdash_site    : ${cdash_site}")
    message(">>> cdash_location: ${cdash_location}")
    message(">>> project_name  : ${project_name}")
    message(">>> build_name    : ${build_name}")
    message(">>> build_stamp   : ${build_stamp}")
    message(">>> machine_name  : ${machine_name}")
    string(REPLACE " " "%20" url_output_tmp
    	"https://${cdash_site}${cdash_location}index.php?project=${project_name}&display=project&filtercount=3&showfilters=1&filtercombine=and&field1=site&compare1=61&value1=${machine_name}&field2=buildname&compare2=61&value2=${build_name}&field3=buildstamp&compare3=61&value3=${build_stamp}"
    )
    set(${url_output} ${url_output_tmp} PARENT_SCOPE)
    banner("generate_build_url1() FINISH")
endfunction()



function(generate_build_url2 url_output cdash_site cdash_location project_name build_name build_stamp)
    banner("generate_build_url2() START")
    message(">>> cdash_site    : ${cdash_site}")
    message(">>> cdash_location: ${cdash_location}")
    message(">>> project_name  : ${project_name}")
    message(">>> build_name    : ${build_name}")
    message(">>> build_stamp   : ${build_stamp}")
    string(REPLACE " " "%20" url_output_tmp
        "https://${cdash_site}${cdash_location}index.php?project=${project_name}&display=project&filtercount=2&showfilters=0&filtercombine=and&field1=buildname&compare1=61&value1=${build_name}&field2=buildstamp&compare2=61&value2=${build_stamp}"
    )
    set(${url_output} ${url_output_tmp} PARENT_SCOPE)
    banner("generate_build_url2() FINISH")
endfunction()



function(generate_build_url3 url_output cdash_site cdash_location project_name build_name build_stamp)
    banner("generate_build_url2() START")
    message(">>> cdash_site    : ${cdash_site}")
    message(">>> cdash_location: ${cdash_location}")
    message(">>> project_name  : ${project_name}")
    message(">>> build_name    : ${build_name}")
    message(">>> build_stamp   : ${build_stamp}")
    string(REPLACE " " "%20" url_output_tmp
        "https://${cdash_site}${cdash_location}index.php?project=${project_name}&filtercount=2&showfilters=0&filtercombine=and&field1=buildname&compare1=61&value1=${build_name}&field2=buildstamp&compare2=61&value2=${build_stamp}"
    )
    set(${url_output} ${url_output_tmp} PARENT_SCOPE)
    banner("generate_build_url3() FINISH")
endfunction()


# generate_build_url4 generates a link to view all builds for a particular PR
function(generate_build_url4 url_output cdash_site cdash_location project_name pr_num)
    banner("generate_build_url4() START")
    message(">>> cdash_site    : ${cdash_site}")
    message(">>> cdash_location: ${cdash_location}")
    message(">>> project_name  : ${project_name}")
    message(">>> pr_num        : ${pr_num}")
    string(REPLACE " " "%20" url_output_tmp
        "https://${cdash_site}${cdash_location}index.php?project=${project_name}&display=project&begin=2024-01-01&end=now&filtercount=1&showfilters=1&field1=buildname&compare1=65&value1=PR-${pr_num}"
    )
    set(${url_output} ${url_output_tmp} PARENT_SCOPE)
    banner("generate_build_url4() FINISH")
endfunction()


message("+--------------------------------------+")
message("| ctest-functions.cmake FINISH         |")
message("+--------------------------------------+")
