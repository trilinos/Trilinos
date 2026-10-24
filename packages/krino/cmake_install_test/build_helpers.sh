#!/usr/bin/env bash

function execute()
{
    if [[ $# -eq 0 ]]; then
        return 0
    fi

    printf '%%'
    printf ' %q' "$@"
    printf '\n'

    set +e
    "$@"
    local status=$?
    set -e

    if [[ $status -ne 0 ]]; then
        echo "ERROR: command failed with status ${status}: $*"
        exit "$status"
    fi
}

function cd_to_clean_dir()
{
    local dir="$1"

    if [[ -z "$dir" || "$dir" == "/" ]]; then
        echo "ERROR: refusing to clean invalid directory: '${dir}'"
        exit 1
    fi

    if [[ -e "$dir" && ! -d "$dir" ]]; then
        echo "ERROR: path exists but is not a directory: '${dir}'"
        exit 1
    fi

    if [[ -d "$dir" ]]; then
        execute find "$dir" -mindepth 1 -maxdepth 1 -exec rm -rf -- {} +
    else
        execute mkdir -p -- "$dir"
    fi

    execute cd "$dir"
}

function cd_to_possibly_new_dir()
{
    local dir="$1"

    if [[ -z "$dir" || "$dir" == "/" ]]; then
        echo "ERROR: invalid directory: '${dir}'"
        exit 1
    fi

    if [[ -e "$dir" && ! -d "$dir" ]]; then
        echo "ERROR: path exists but is not a directory: '${dir}'"
        exit 1
    fi

    execute mkdir -p -- "$dir"
    execute cd "$dir"
}

function exit_with_message()
{
    echo "$1"
    exit 1
}

function make_and_install()
{
    make -j 16 || exit_with_message "Failed building $1"
    make install || exit_with_message "Failed installing $1"
}

function build_yaml()
{
    local product_name=yaml
    local build_type=$1

    execute cd "${output_dir}/${product_name}"
    cd_to_clean_dir ${product_name}_build
    execute cmake \
        -DCMAKE_BUILD_TYPE="${build_type}" \
        -DYAML_CPP_BUILD_TESTS=false \
        -DCMAKE_INSTALL_PREFIX="../${product_name}_install" \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        ../yaml-cpp
    make_and_install "$product_name"
}

function link_current_krino_and_stk()
{
    execute rm -rf -- Trilinos/packages/krino
    execute rm -rf -- Trilinos/packages/stk
    execute ln -s "${sierra_proj}/krino" Trilinos/packages/krino
    execute ln -s "${sierra_proj}/stk" Trilinos/packages/stk
}

function setup_build_with_krino()
{
    cd_to_possibly_new_dir "${output_dir}"
    cd_to_possibly_new_dir yaml
    if [ ! -d yaml-cpp ] ; then
      execute tar -xzf ${sierra_proj}/TPLs_src/spack/spack_tpls/yaml-cpp/*.tar.gz
    fi
    
    cd_to_possibly_new_dir "../trilinos"
    if [ ! -d Trilinos ] ; then
      execute git clone -b develop https://github.com/trilinos/Trilinos.git Trilinos
    else
      execute cd Trilinos
      execute git checkout develop
      execute git clean -fd
      execute git reset --hard
      execute git pull
      execute cd ..
    fi
    
    link_current_krino_and_stk
}

function get_krino_variant_name()
{
    if [[ "$enable_pykrino" == "ON" ]]; then
        echo "pykrino"
    else
        echo "krino"
    fi
}

function get_trilinos_build_name()
{
    local build_type="$1"
    local krino_name
    krino_name="$(get_krino_variant_name)"

    echo "trilinos_${krino_name}_${build_type}"
}

function build_trilinos_with_krino()
{
    local product_name="trilinos"
    local build_type="$1"

    if [[ -z "$build_type" ]]; then
        echo "ERROR: build_trilinos_with_krino requires a build type argument"
        return 1
    fi

    local product_build_name
    product_build_name="$(get_trilinos_build_name "$build_type")"

    local build_dir="${product_build_name}_build"
    local install_dir="../${product_build_name}_install"

    execute cd "${output_dir}/${product_name}"
    cd_to_clean_dir "${build_dir}"

    execute rm -rf -- "${install_dir}"

    execute env \
        TRILINOS_INSTALL_DIR="${install_dir}" \
        CMAKE_BUILD_TYPE="${build_type}" \
        ENABLE_PYKRINO="${enable_pykrino}" \
        ENABLE_YAML="${enable_yaml}" \
        "${sierra_proj}/krino/cmake_install_test/run_cmake_krino"

    make_and_install "${product_name}"
}

function get_modules_file()
{
    # KRINO_MODULES_FILE overrides the choice by OS release.
    if [[ -n "${KRINO_MODULES_FILE}" ]]; then
        echo "${KRINO_MODULES_FILE}"
        return 0
    fi

    local modules_dir="${sierra_proj}/krino/cmake_install_test"
    local rhel_major
    rhel_major="$(. /etc/os-release && echo "${VERSION_ID%%.*}")"

    case "${rhel_major}" in
        8) echo "${modules_dir}/load_gcc_modules_rhel8" ;;
        9) echo "${modules_dir}/load_gcc_modules_rhel9" ;;
        *) echo "ERROR: no modules file for OS version '${rhel_major}'; set KRINO_MODULES_FILE" >&2
           return 1 ;;
    esac
}

function setup_environment()
{
    execute source /etc/profile.d/modules.sh
    execute module purge

    local modules_file
    modules_file="$(get_modules_file)" || exit 1
    execute source "${modules_file}"

    execute module list
    execute env
}

function run_tests()
{
    local build_type="$1"

    local product_name="trilinos"
    local product_build_name="$(get_trilinos_build_name "$build_type")"

    execute cd "${output_dir}/${product_name}/${product_build_name}_build"
    execute ctest -j 16 --output-on-failure
}
