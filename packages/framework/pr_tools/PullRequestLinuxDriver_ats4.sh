#!/bin/bash -l

# Use flux to allocate resources and run the job
# flux alloc waits for completion (like bsub -Is)
flux alloc --name=trilinos_nightly --time=12:00 \
  --nodes=1 \
  --output=${WORKSPACE}/flux-output-${JOB_NAME}.log \
  --error=${WORKSPACE}/flux-error-${JOB_NAME}.log \
  ${WORKSPACE}/Trilinos/packages/framework/pr_tools/PullRequestLinuxDriver.sh
