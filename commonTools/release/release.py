#!/usr/bin/env python3
import argparse
import sys
import shutil
from utils import *
import logging
logger = logging.getLogger(__name__)


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
        "--debug",
        action='store_true',
        help="Debug outputs."
    )
    args = parser.parse_args()
    return args


def pre_checks():
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

    pre_checks()
    rel_version = parse_semver(args.rel_version)
    dev_version = parse_semver(args.dev_version)
    return 0


if __name__ == "__main__":
    sys.exit(main())
