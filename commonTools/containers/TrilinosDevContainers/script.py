"""
Helper script for launching development containers

This script uses parses .github/workflows/AT2.yml for container images.
"""
from pathlib import Path
import questionary
import argparse
from shutil import which
import os
import sys
from pprint import pformat
from . tools import (findGitRoot,
                     is_valid_source_dir,
                     convert_to_valid_source_dir,
                     parse_workflow_images,
                     image_available,
                     pull_image,
                     launch_container,
                     ensure_podman_machine)
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

##################################################

def main():

    parser = argparse.ArgumentParser(description="""This script helps with launching development containers for the Trilinos project.
    Run the script without any arguments in interactive mode or specify required information using the command line flags.""")
    parser.add_argument("--source", type=convert_to_valid_source_dir, help="Trilinos source directory. Allows to skip the interactive prompt.")
    parser.add_argument("--build", type=Path, help="Trilinos build directory. Allows to skip the interactive prompt.")
    parser.add_argument("--debug", action="store_true", help="debug mode")
    parser.add_argument("--container-engine", choices=["docker", "podman", "cask"], default="podman",
                        help="Container engine to use (default: podman).")
    parser.add_argument("--container-image", default=None, help="Override container image")
    parser.add_argument("--container-tag", default=None, help="Override container tag")
    args, unparsed_args = parser.parse_known_args()
    container_engine = args.container_engine

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
    # Check that the requested container engine is available
    container_engine_cmd = which(container_engine)
    logger.debug(f"container_engine_cmd = {container_engine_cmd}")
    if container_engine_cmd is None:
        logger.warning(f"No command \"{container_engine}\" in path. Aborting.")
        exit(1)
    # On Mac: make sure that we can get a running machine
    if container_engine == "podman" and sys.platform == "darwin":
        ensure_podman_machine()

    ##################################################
    # Trilinos source directory

    if args.source is not None:
        trilinos_source = args.source
    else:
        trilinos_default_source = findGitRoot(Path(os.path.abspath(__file__)))
        logger.debug(f"trilinos_default_source = {trilinos_default_source}")
        trilinos_source = questionary.path("Where is the Trilinos source code checked out?",
                                           default=str(trilinos_default_source) if trilinos_default_source.exists() else "",
                                           validate=is_valid_source_dir).ask()
        trilinos_source = Path(trilinos_source)
    logger.debug(f"trilinos_source = {trilinos_source}")

    ##################################################
    # Trilinos build directory
    if args.build is not None:
        trilinos_build = args.build
    else:
        trilinos_build = questionary.path("Where should Trilinos be built?", default=str(Path.home()/"builds"/"trilinos")).ask()
    trilinos_build = Path(trilinos_build)
    if not trilinos_build.exists():
        if questionary.confirm(f"Build directory \"{trilinos_build}\" does not exist. Create?", default=True).ask():
            trilinos_build.mkdir(parents=True)
        else:
            exit()

    ##################################################
    # Available container images

    images = parse_workflow_images(trilinos_source, container_engine)
    logger.debug("images = \n" + pformat(images))
    assert images is not None and len(images) > 0

    ##################################################
    # Choose image

    if args.container_image is not None and args.container_tag is not None:
        image, tag = args.container_image, args.container_tag
    else:
        image, tag = questionary.select("Which image do you want to use?", choices=[questionary.Choice(":".join(key), key) for key in sorted(images)]).ask()

    ##################################################
    # Summary

    questionary.print(f"\nContainer image:         {image}:{tag}")
    questionary.print(f"Container engine:        {container_engine}")
    questionary.print(f"Host source directory:   {trilinos_source}")
    questionary.print(f"Host build directory:    {trilinos_build}\n")

    ##################################################
    # Pull image

    image_avail = image_available(container_engine, image, tag)
    logger.debug(f"image_avail = {image_avail}")
    if not image_avail:
        pull_image(container_engine, image, tag)

    ##################################################
    # Launch container

    launch_container(container_engine, trilinos_source, trilinos_build, image, tag, unparsed_args)
