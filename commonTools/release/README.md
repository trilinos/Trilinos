## Trilinos Release Script

Requirements:
- Python 3.9+
- Local copy of Trilinos repository that has a remote `origin` to the SSH authenticated Trilinos
- Valid SSH key that is registered with your GitHub account.
- `GITHUB_TOKEN` environment variable set to a token for the target repository.
  - GitHub fine-grain token with:
    - Admin (read/write), needed for modifying branch protection rules.
    - Pull Request (read/write)

Install dependencies
```
python3 -m pip install -r requirements.txt
```

## Usage

The release process is split into subcommands so each phase can be run independently.

```
release.py pre-release <rel_version> <dev_version> [--dir DIR] [--debug]
release.py release     <rel_branch>                [--dir DIR] [--debug]
release.py patch       <rel_branch>                [--dir DIR] [--debug]
```

- `rel_version` will be desired release version in semver form (i.e "17.1.0")
- `dev_version` will be the version associated with the next upcoming release (i.e "17.2.0-dev")
- `rel_branch` refers to the name of an existing release branch that exists in the project's repository

NOTE: The standard Trilinos development versioning adds a "-dev" which is semver compliant.

### `pre-release` Mode

Pre-release mode will prepare everything for a Trilinos release up until the point
where a release is tagged and announced (see `release` mode below).

The set of pre-release operations are:
- creates and pushes the release branch to GitHub
- applies branch protection rules to the pushed release branch on GitHub
- opens a GitHub pull request updating the release branch's `Version.cmake`
- opens a GitHub pull request updating the  `develop` branch's `Version.cmake`
- creates the GitHub issue label for the next release's release notes

```
release.py pre-release 17.1.0 17.2.0-dev
```

### `release` Mode

Given the name of a release branch, fetch that branch and create a release tag at the
tip of the release branch. Release mode will validate that the `pre-release` mode's `Version.cmake`
has merged by checking if dev-mode is turned off in branch's `Version.cmake`.

```
release.py release trilinos-release-17-1-branch
```

### `patch` Mode

Given the name of a release branch, fetch that branch and create a new branch from it that
bumps the patch number in `Version.cmake`. Once the version number has been bumped, push version
update branch and create a pull request into the target release branch.

The patch mode is intended to behave similarly to the `pre-release` mode in that it does not perform
any of the operations that formalizes the release, such as creating a new tag or announcing the release.
The user should use the `release` mode after using the `patch` mode when the patch is ready for release.

```
release.py patch trilinos-release-17-1-branch
# Version.cmake PR gets merged...
release.py release trilinos-release-17-1-branch
```
