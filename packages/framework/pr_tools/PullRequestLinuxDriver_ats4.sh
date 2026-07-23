#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

# Environment variables for configuration (can be set externally):
#   FLUX_ALLOC_ARGS: Additional arguments for flux alloc
#                    (default: --nodes=2 --cores-per-task=32 --gpus=0)
#   FLUX_OUTPUT_LOG: Path for output log
#                    (default: ${WORKSPACE}/flux-output-${JOB_NAME}.log)
#   FLUX_ERROR_LOG:  Path for error log
#                    (default: ${WORKSPACE}/flux-error-${JOB_NAME}.log)

# Set default output and error log paths if not provided
if [ "${FLUX_OUTPUT_LOG:-}" == "" ]; then
  FLUX_OUTPUT_LOG="${WORKSPACE}/flux-output-${JOB_NAME}.log"
fi

if [ "${FLUX_ERROR_LOG:-}" == "" ]; then
  FLUX_ERROR_LOG="${WORKSPACE}/flux-error-${JOB_NAME}.log"
fi

export TRILINOS_MAX_CORES=96
export TRILINOS_NUM_CONCURRENT_TESTS=96


flux alloc --job-name=trilinos_nightly --time-limit=12h \
  --nodes 1 ${FLUX_ALLOC_ARGS} \
  --output=${FLUX_OUTPUT_LOG} \
  --error=${FLUX_ERROR_LOG} \
  ${WORKSPACE}/Trilinos/packages/framework/pr_tools/PullRequestLinuxDriver.sh $@
