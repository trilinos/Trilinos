#!/usr/bin/env python3
import argparse
import sys
import os
import shutil
from pathlib import Path
from utils import *
import logging
logger = logging.getLogger(__name__)

script_path = Path(__file__).resolve()
git_cmd = shutil.which("git")


def add_common_args(parser):
    """Add options shared by every subcommand."""
    parser.add_argument(
        "--dir",
        help="Path to another local Trilinos directory."
    )
    parser.add_argument(
        "--debug",
        action='store_true',
        help="Debug outputs."
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Trilinos Release Automation Script")
    subparsers = parser.add_subparsers(dest="command", required=True)

    pre = subparsers.add_parser(
        "pre-release",
        help="Create, push, and protect the release branch. Update Version.cmakes "
             "for release branch and dev branch and open PRs for those changes. Create GitHub issue label for the next release notes."
    )
    pre.add_argument(
        "rel_version",
        help="Version string in format MAJOR.MINOR.PATCH (e.g., 17.1.0)"
    )
    pre.add_argument(
        "dev_version",
        help="Development version string in format MAJOR.MINOR.PATCH-dev (e.g., 17.1.0-dev)"
    )
    add_common_args(pre)
    pre.set_defaults(func=pre_release)

    rel = subparsers.add_parser(
        "release",
        help="Finalize an existing release branch for release. Create a git tag for the current version of the release branch and pushes the tag to GitHub."
    )
    rel.add_argument(
        "rel_branch",
        help="Name of an existing release branch (e.g., trilinos-release-17-1-branch)"
    )
    add_common_args(rel)
    rel.set_defaults(func=release)

    pat = subparsers.add_parser(
        "patch",
        help="Bump the patch number of an existing release "
             "branch's Version.cmake via a pull request. User is expected to use 'release' subcommand "
             "to finalize the patch release."
    )
    pat.add_argument(
        "rel_branch",
        help="Name of an existing release branch (e.g., trilinos-release-17-1-branch)"
    )
    add_common_args(pat)
    pat.set_defaults(func=patch)

    args = parser.parse_args()
    return args


def pre_checks(args):
    git_root = get_git_root(args.dir if args.dir else script_path)
    logger.debug(f"git = {git_cmd}")

    if git_cmd is None:
        raise FileNotFoundError(f"Cannot find git executable")

    if (not does_remote_exists("fork", git_root) or not does_remote_exists("origin", git_root)):
        raise RuntimeError(f"Missing remotes 'fork' or 'origin'")

    if not is_git_workspace_clean(git_root):
        raise RuntimeError(f"Git workspace not clean: {git_root}")

    token = os.getenv('GITHUB_TOKEN')
    if not token:
        raise RuntimeError(f"GITHUB_TOKEN environment variable not set.")


def pre_release(args):
    """Pre-release phase.

    Idea is to perform every step that is necessary before tagging and announcing the release.
    This should allow developers to iterating on the pushed release-branch before it's ready.
    """
    git_root = get_git_root(args.dir if args.dir else script_path)

    rel = parse_semver(args.rel_version)
    rel_branch = f"trilinos-release-{rel['major']}-{rel['minor']}-branch"
    dev = parse_semver(args.dev_version)
    dev_update_branch = f"update-dev-version-{dev['major']}-{dev['minor']}"

    main_branch = "master"
    dev_branch = "develop"

    fork_owner = get_remote_owner("fork", git_root)

    #####################################
    # Fetch latest master and develop

    fetch_branch(dev_branch, git_root, merge=True)
    fetch_branch(main_branch, git_root, merge=True)

    #####################################
    # Checkout and push new release branch

    checkout_branch(rel_branch, git_root)
    push(rel_branch, git_root)
    set_release_branch_protection(rel_branch)
    print(f"> Pushed created {rel_branch} to origin repository")

    #####################################
    # Checkout and push new branch to
    # update release branch's Version.cmake

    rel_update_branch = f"update-rel-version-{rel['major']}-{rel['minor']}"
    checkout_branch(rel_update_branch, git_root)
    dev_mode = False
    update_version_cmake(args.rel_version, rel_branch, dev_mode, git_root)
    release_commit_msg = "Update release Version.cmake"
    commit_tracked(release_commit_msg, git_root)
    push(rel_update_branch, git_root, remote="fork")

    #####################################
    # Create Pull Request for release branch version update
    title = f"Framework: Update {args.rel_version} release Version.cmake"
    body = "@trilinos/framework"
    pr = create_pull_request(rel_branch, f"{fork_owner}:{rel_update_branch}", title, body)
    print(f"> Created release branch update PR: {pr.html_url}")

    #####################################
    # Checkout and push new branch from
    # develop to update develop's Version.cmake

    checkout_branch(dev_branch, git_root, remote=True)
    checkout_branch(dev_update_branch, git_root)
    dev_mode = True
    update_version_cmake(args.dev_version, dev_update_branch, dev_mode, git_root)
    dev_commit_msg = "Update develop Version.cmake"
    commit_tracked(dev_commit_msg, git_root)
    push(dev_update_branch, git_root, remote="fork")

    #####################################
    # Create PR to update Version.cmake in develop branch
    title = "Framework: Update develop release Version.cmake"
    pr = create_pull_request(dev_branch, f"{fork_owner}:{dev_update_branch}", title, body)
    print(f"> Created dev branch update PR: {pr.html_url}")

    #####################################
    # Create a GitHub issue label for the next set of release notes

    label_name = f"{dev['major']}.{dev['minor']} release note"
    create_release_label(label_name)
    print(f"> Created GitHub label: {label_name}")

    return 0


def release(args):
    """Release phase

    Complete the release process for a given release branch that exists in the remote repository.
    This can be done for a release branch created from the 'pre-release' phase or after a 'patch' bump
    has been applied to an existing release branch.
    """
    git_root = get_git_root(args.dir if args.dir else script_path)
    rel_branch = args.rel_branch

    #####################################
    # Ensure the release branch exists on origin and sync to the remote tip.
    verify_remote_branch_exists(rel_branch, git_root)
    fetch_branch(rel_branch, git_root)
    checkout_branch(rel_branch, git_root, remote=True)

    #####################################
    # Validate that the branch is in release mode before tagging. If the
    # Version.cmake update PR has not been merged yet, refuse to tag.

    version_info = read_version_cmake(git_root)
    if not is_release_mode(version_info):
        raise RuntimeError(f"Release branch '{rel_branch}' is not in release mode.")

    #####################################
    # Create and push the release tag using the version from Version.cmake

    version = version_info["version"]
    rel = parse_semver(version)
    tag_name = f"trilinos-release-{rel['major']}-{rel['minor']}-{rel['patch']}"
    tag_message = f"Trilinos release {version}"
    create_tag(tag_name, tag_message, rel_branch, git_root)
    push_tag(tag_name, git_root)
    print(f"> Created and pushed tag: {tag_name}")

    return 0


def patch(args):
    """Patch mode

    Bump the patch version number from the given release branch.
    """
    git_root = get_git_root(args.dir if args.dir else script_path)
    rel_branch = args.rel_branch

    fork_owner = get_remote_owner("fork", git_root)

    #####################################
    # Ensure the release branch exists on origin and sync to the remote tip.

    verify_remote_branch_exists(rel_branch, git_root)
    fetch_branch(rel_branch, git_root)
    checkout_branch(rel_branch, git_root, remote=True)

    #####################################
    # Determine the new patch version by bumping the current patch number.

    version_info = read_version_cmake(git_root)
    current = parse_semver(version_info["version"])
    major = current['major']
    minor = current['minor']
    new_patch = int(current['patch']) + 1
    new_version = f"{major}.{minor}.{new_patch}"
    print(f"> Bumping {rel_branch} from {version_info['version']} to {new_version}")

    #####################################
    # Checkout and push a new branch to update the release branch's Version.cmake

    patch_update_branch = f"update-patch-version-{major}-{minor}-{new_patch}"
    checkout_branch(patch_update_branch, git_root)
    dev_mode = False
    update_version_cmake(new_version, rel_branch, dev_mode, git_root)
    patch_commit_msg = "Update patch release Version.cmake"
    commit_tracked(patch_commit_msg, git_root)
    push(patch_update_branch, git_root, remote="fork")

    #####################################
    # Create Pull Request for the patch version update

    title = f"Framework: Bump {new_version} patch version"
    body = "@trilinos/framework"
    pr = create_pull_request(rel_branch, f"{fork_owner}:{patch_update_branch}", title, body)
    print(f"> Created patch release update PR: {pr.html_url}")

    return 0


def main():
    args = parse_args()
    if args.debug:
        logging.basicConfig(level=logging.INFO)
        logger.parent.setLevel(logging.DEBUG)
    logger.debug(f"args: {args}")

    pre_checks(args)

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
