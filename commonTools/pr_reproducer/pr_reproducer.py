"""
Helper script for reproducing PR builds

This script uses the Github API to retrieve information about PRs and parses
.github/workflows/AT2.yml for build information.
"""
from pathlib import Path
import questionary
import argparse
from shutil import which
import os
from pprint import pformat
import subprocess
from github import Github
from tools import (is_valid_pr_number,
                   convert_to_pr_number,
                   is_valid_source_dir,
                   convert_to_valid_source_dir,
                   get_github_upstream_of_local_branch,
                   get_pr_number,
                   parse_AT2_workflows,
                   generate_package_enables,
                   image_available,
                   pull_image,
                   launch_container)
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

remote = "origin"
repo = "trilinos/Trilinos"

##################################################

parser = argparse.ArgumentParser(description="""This script helps with reproducining PR builds for the Trilinos project.
Run the script without any arguments in interactive mode or specify required information using the command line flags.""")
parser.add_argument("--pr", type=convert_to_pr_number, help="Number of the PR that should be reproduced. Allows to skip the interactive prompt.")
parser.add_argument("--source", type=convert_to_valid_source_dir, help="Trilinos source directory. Allows to skip the interactive prompt.")
parser.add_argument("--build", help="PR build that should be reproduced. Allows to skip the interactive prompt.")
parser.add_argument("--debug", action="store_true", help="debug mode")
args = parser.parse_args()

if args.debug:
    logger.parent.setLevel(logging.DEBUG)

logo = """

         %%%%%%%%%%
     %%%      %%% %%%%
   %%%      %%% %%%%  %%%
  %%      %%% %%%%      %%
 %%  %%%%%%%%%%%%%%%%%%%%%%%       %%%%%%%%%%%%  %%%%%%%%      %%%  %%%%         %%%%  %%       %%%      %%%%%%         %%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%       %%%%%%%%%%%%  %%%%%%%%%%%%%%  %%%  %%%%         %%%% %%%%%%    %%%  %%%%%%%%%%%%%%  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  %%         %%%%                %%%%  %%%  %%%%         %%%% %%%%%%%%  %%%  %%%        %%%% %%%%%%%%%%%
%%%%%  %%% %%%%%          %%         %%%%      %%%%%%%%%%%%%   %%%  %%%%         %%%% %%%% %%%%%%%% %%%%        %%%%   %%%%%%%%%%%%
%%%  %%%   %%%%%          %          %%%%      %%%     %%%%%   %%%  %%%%%%%%%%%  %%%% %%%%   %%%%%%  %%%%%%%%%%%%%%   %%%%%%%%%%%%%
 %%%%%     %%%%%         %%          %%%%      %%%      %%%%%  %%%  %%%%%%%%%%   %%%% %%%%      %%      %%%%%%%%%    %%%%%%%%%%%%
  %%       %%%%%        %%
   %%%     %%%%%      %%
     %%%%  %%%%%   %%%
         %%%%%%%%%%

"""

questionary.print(logo)

##################################################
# Check that podman is available
podman_cmd = which("podman")
logger.debug(f"podman_cmd = {podman_cmd}")
assert podman_cmd is not None, "No command \"podman\" in path"

##################################################
# Check that cmake is available
cmake_cmd = which("cmake")
logger.debug(f"cmake_cmd = {cmake_cmd}")
assert cmake_cmd is not None, "No command \"cmake\" in path"

##################################################
# Trilinos source directory

if args.source is not None:
    trilinos_source = args.source
else:
    # relies on this file being in commonTools/or_reproducer/
    trilinos_default_source = Path(os.path.abspath(__file__)).parent.parent.parent
    logger.debug(f"trilinos_default_source = {trilinos_default_source}")
    trilinos_source = questionary.path("Where is the Trilinos source code checked out?",
                                       default=str(trilinos_default_source) if trilinos_default_source.exists() else "",
                                       validate=is_valid_source_dir).ask()
    trilinos_source = Path(trilinos_source)
logger.debug(f"trilinos_source = {trilinos_source}")

local_git_sha = subprocess.run(["git", "rev-parse", "HEAD"],
                               capture_output=True, cwd=trilinos_source, universal_newlines=True).stdout[:-1]
logger.debug(f"local_git_sha = {local_git_sha}")

##################################################
# Available builds

pr_builds = parse_AT2_workflows(trilinos_source)
logger.debug("pr_builds = \n" + pformat(pr_builds))
assert pr_builds is not None and len(pr_builds) > 0

##################################################
# Get PR number

if args.pr is not None:
    pr_number = args.pr
else:
    upstream = get_github_upstream_of_local_branch(trilinos_source)
    logger.debug(f"upstream = {upstream}")
    if upstream is not None:
        pr_number_default = get_pr_number(repo, upstream)
    else:
        pr_number_default = None
    logger.debug(f"pr_number_default = {pr_number_default}")

    pr_number = questionary.text("What is the number of the PR that should be reproduced?",
                                 default=str(pr_number_default) if pr_number_default is not None else "",
                                 validate=is_valid_pr_number).ask()
    pr_number = convert_to_pr_number(pr_number)
logger.debug(f"pr_number = {pr_number}")

##################################################
# Choose build

if args.build is not None:
    pr_build = args.build
    assert pr_build in pr_builds, f"Specified build \"{pr_build}\" is not a known build. Known builds: {list(pr_builds.keys())}"
else:
    pr_build = questionary.select("Which PR build do you want to reproduce?", choices=pr_builds.keys()).ask()
logger.debug(f"pr_build = {pr_build}")

image = pr_builds[pr_build]["image"]
tag = pr_builds[pr_build]["tag"]
genconfig_build_id = pr_builds[pr_build]["genconfig_build_id"]
logger.debug(f"image = {image}")
logger.debug(f"tag = {tag}")
logger.debug(f"genconfig_build_id = {genconfig_build_id}")

##################################################
# Retrieve PR information from Github

g = Github()
repo = g.get_repo(repo)
pr = repo.get_pull(pr_number)
pr_head_ref = pr.head.ref
logger.debug(f"pr.head.sha = {pr.head.sha}")
logger.debug(f"pr.merge_commit_sha = {pr.merge_commit_sha}")
pr_base_ref = pr.base.ref
logger.debug(f"pr_base_ref = {pr_base_ref}")

if pr.head.sha != local_git_sha:
    do_continue = questionary.confirm(f"The current commit SHA {local_git_sha} of the local repo does not match the SHA {pr.head.sha} of the PR. This might not be a problem if you made additional local changes. Continue?", default=True).ask()
    if not do_continue:
        exit()
else:
    logger.debug("Local commit SHA matches PR.")
g.close()

##################################################
# Generate packageEnables.cmake

packageEnablesFile = Path(".")/f'packageEnables-PR-{pr_number}.cmake'
if packageEnablesFile.exists():
    overwrite = questionary.confirm(f"File {packageEnablesFile} already exists. Overwrite?", default=True).ask()
    if not overwrite:
        questionary.print("Not overwriting file. Aborting.")
        exit()

logger.debug("Fetching target branch")
subprocess.run(["git", "fetch", remote, f"{pr_base_ref}:{pr_number}-target"],
               capture_output=True, cwd=trilinos_source, universal_newlines=True)
logger.debug("Fetching merge branch")
subprocess.run(["git", "fetch", remote, f"refs/pull/{pr_number}/merge:{pr_number}-merge"],
               capture_output=True, cwd=trilinos_source, universal_newlines=True)

output = generate_package_enables(trilinos_source, packageEnablesFile, f"{pr_number}-target", f"{pr_number}-merge")
if not packageEnablesFile.exists():
    logger.info("Output from script that generates CMake fragment with package enables:\n"+output)
    assert False, "No package enables generated. Does the PR change any packages?"

with open(packageEnablesFile, 'r') as f:
    packageEnables = f.readlines()
logger.debug(f"{packageEnablesFile}:" + ''.join(packageEnables))

##################################################
# Summary

questionary.print(f"PR number:          {pr_number}")
questionary.print(f"Title:              {pr.title}")
questionary.print(f"Head ref:           {pr_head_ref}")
questionary.print(f"Base ref:           {pr_base_ref}")
questionary.print(f"Build:              {pr_build}")
questionary.print(f"Container image:    {image}:{tag}")
questionary.print(f"Genconfig build ID: {genconfig_build_id}")

##################################################
# Pull image

image_avail = image_available(image, tag)
logger.debug(f"image_avail = {image_avail}")
if not image_avail:
    pull_image(image, tag)

##################################################
# Launch container

launch_container(trilinos_source, packageEnablesFile, image, tag, genconfig_build_id)
