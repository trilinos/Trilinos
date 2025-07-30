# Trilinos Pull Request Testing Tools

## Entry Point Scripts

### PullRequestLinuxDriverTest.py
This script is intended to drive the "core" configure/build/test/report part of the process.  As such, this is the tool that is called as part of most of the GitHub Actions (seen in `.github/workflows/AT2.yml`).

### PullRequestLinuxDriverMerge.py
This script does the merge-related parts of the CI testing process.

#### Behavior
1. Merge a passed source ref (branch or commit) into a passed target branch
