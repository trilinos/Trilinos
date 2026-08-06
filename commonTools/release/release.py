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

    #####################################
    # Fetch latest master and develop

    fetch_branch(dev_branch, git_root, merge=True)
    fetch_branch(main_branch, git_root, merge=True)

    #####################################
    # Checkout and push new release branch

    checkout_branch(rel_branch, git_root)
    push(rel_branch, git_root)
    set_release_branch_protection(rel_branch)
    print(f"Pushed created {rel_branch} to origin repository")

    #####################################
    # Checkout and push new branch to
    # update release branch's Version.cmake

    rel_update_branch = f"update-rel-version-{rel['major']}-{rel['minor']}"
    checkout_branch(rel_update_branch, git_root)
    dev_mode = False
    update_version_cmake(args.rel_version, rel_branch, dev_mode, git_root)
    release_commit_msg = "Update release Version.cmake"
    commit_tracked(release_commit_msg, git_root)
    push(rel_update_branch, git_root)

    #####################################
    # Create Pull Request for release branch version update
    title = f"Framework: Update {args.rel_version} release Version.cmake"
    body = "@trilinos/framework"
    pr = create_pull_request(rel_branch, rel_update_branch, title, body)
    print(f"Created release branch update PR: {pr.html_url}")

    #####################################
    # Checkout and push new branch from 
    # develop to update develop's Version.cmake

    checkout_branch(dev_branch, git_root, remote=True)
    checkout_branch(dev_update_branch, git_root)
    dev_mode = True
    update_version_cmake(args.dev_version, dev_update_branch, dev_mode, git_root)
    dev_commit_msg = "Update develop Version.cmake"
    commit_tracked(dev_commit_msg, git_root)
    push(dev_update_branch, git_root)

    #####################################
    # Create PR to update Version.cmake in develop branch
    title = "Framework: Update develop release Version.cmake"
    pr = create_pull_request(dev_branch, dev_update_branch, title, body)
    print(f"Created dev branch update PR: {pr.html_url}")


    #####################################
    # Create tags for the release branch and push tag to GitHub

    tag_name = f"trilinos-release-{rel['major']}-{rel['minor']}-{rel['patch']}"
    tag_message = f"Trilinos release {args.rel_version}"
    create_tag(tag_name, tag_message, rel_branch, git_root)
    push_tag(tag_name, git_root)
    print(f"Created and pushed tag: {tag_name}")

    #####################################
    # Create a GitHub issue label for the next set of release notes

    label_name = f"{dev['major']}.{dev['minor']} release note"
    create_release_label(label_name)
    print(f"Created GitHub label: {label_name}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

