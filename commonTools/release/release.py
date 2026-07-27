#!/usr/bin/env python3
import argparse
import sys
import shutil
from pathlib import Path
from utils import *
import logging
logger = logging.getLogger(__name__)

script_path = Path(__file__).resolve()
git_cmd = shutil.which("git")


def parse_args():
    parser = argparse.ArgumentParser(description="Trilinos Release Automation Script")
    parser.add_argument(
        "rel_version",
        help="Version string in format MAJOR.MINOR.PATCH (e.g., 17.1.0)"
    )
    parser.add_argument(
        "dev_version",
        help="Development version string in format MAJOR.MINOR.PATCH-dev (e.g., 17.1.0-dev)"
    )
    parser.add_argument(
        "--dir",
        help="Path to another local Trilinos directory."
    )
    parser.add_argument(
        "--debug",
        action='store_true',
        help="Debug outputs."
    )
    args = parser.parse_args()
    return args


def pre_checks(args):
    logger.debug(f"git = {git_cmd}")
    if git_cmd is None:
        raise FileNotFoundError(f"Cannot find git executable")
        exit(1)


release_commit_msg = "Update release Version.cmake"
dev_commit_msg = "Update dev Version.cmake"


def main():
    args = parse_args()
    if args.debug:
        logging.basicConfig(level=logging.INFO)
        logger.parent.setLevel(logging.DEBUG)
    logger.debug(f"args: {args}")

    pre_checks(args)
    git_root = get_git_root(args.dir if args.dir else script_path)


    rel = parse_semver(args.rel_version)
    rel_branch= f"trilinos-release-{rel['major']}-{rel['minor']}-branch"
    dev = parse_semver(args.dev_version)
    dev_update_branch= f"update-dev-version-{dev['major']}-{dev['minor']}"

    main_branch = "master"
    dev_branch = "develop"

    # Fetch latest master and develop

    fetch_branch(dev_branch, git_root, merge=True)
    fetch_branch(main_branch, git_root, merge=True)

    # Checkout new release branch

    checkout_branch(rel_branch, git_root)

    # Push new release branch to origin

    push(rel_branch, git_root)

    # Update Version.cmake in release branch
    rel_update_branch = f"update-rel-version-{dev['major']}-{dev['minor']}"
    checkout_branch(rel_update_branch, git_root)
    dev_mode = False
    update_version_cmake(args.rel_version, rel_branch, dev_mode, git_root)
    commit_tracked(release_commit_msg, git_root)
    push(rel_update_branch, git_root)

    # Create branch protection rules for pushed release branch

    # TODO

    # Make PR from Version.cmake update branch to release branch
    # Probably want to abstract this into another fucntion with templated title/body
    create_pull_request(rel_branch, rel_update_branch, "Update release Version.cmake", "")

   # Checkout remote "develop" and update
    checkout_branch(dev_branch, git_root, remote=True)

    # Checkout new branch for dev version update
    checkout_branch(dev_update_branch, git_root)
    dev_mode = True
    update_version_cmake(args.dev_version, dev_update_branch, dev_mode, git_root)
    commit_tracked(dev_commit_msg, git_root)
    push(dev_update_branch, git_root)

    # Create PR to update Version.cmake in develop branch
    create_pull_request(dev_branch, dev_update_branch, "Update develop Version.cmake", "")

    return 0


if __name__ == "__main__":
    sys.exit(main())
