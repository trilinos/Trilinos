#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys
import subprocess
import os

# Packages are snapshotted via install_reqs.sh or these are in Python's site-packages
try:                                                                                # pragma: no cover
    from .determinesystem import DetermineSystem

except ImportError:                                                                 # pragma: no cover
    try:  # Perhaps this file is being imported from another directory
        p = Path(__file__).parents[0]
        sys.path.insert(0, str(p))

        from determinesystem import DetermineSystem

    except ImportError:  # Perhaps LoadEnv was snapshotted and these packages lie up one dir.
        p = Path(__file__).parents[1]  # One dir up from the path to this file
        sys.path.insert(0, str(p))

        from determinesystem import DetermineSystem



def get_launch_env(system : str):
  """
  Gets the launch environment based on the detected system.
  This is an early environment that's required for running the driver.

  Returns:
      str: The environment used to launch the driver.
  """
  env = ""

  if env == "":
      return ""
  else:
      return "env" + env + " "


def get_launch_cmd(system : str):
  """
  Gets the launch command based on the detected system.

  Returns:
      str: The command used to launch the driver.
  """
  cmd = ""

  return cmd + " "


def get_driver_args(system : str):
  """
  Gets the driver arguments based on the detected system.

  Returns:
      str: The arguments passed to the driver.
  """
  return " " + "--on_" + system


def main(argv):
  """
  This python script determines what system it is running on and then launches
  the trilinos driver script appropriately.

  The script returns 0 upon success and non-zero otherwise.
  """
  parser = argparse.ArgumentParser(description='Launch a trilinos driver script on this system.')
  parser.add_argument('--build-name', required=True,
                      help='The name of the build being launched')
  parser.add_argument('--driver', required=False,
                      default='./PullRequestLinuxDriver.sh',
                      help='The driver script to launch')
  parser.add_argument('--supported-systems', required=False,
                      default='./LoadEnv/ini_files/supported-systems.ini',
                      help='The INI file containing supported systems')
  parser.add_argument('--in-container', default=False, action="store_true",
                      help="Build is happening in a container")
  parser.add_argument("--kokkos-develop", default=False, action="store_true",
                      help="Build is requiring to pull the current develop of kokkos and kokkos-kernels packages")
  parser.add_argument("--extra-configure-args",
                      help="Extra arguments that will be passed to CMake for configuring Trilinos.")
  args = parser.parse_args(argv)

  if os.getenv("TRILINOS_DIR") == None:
    print("LaunchDriver> ERROR: Please set TRILINOS_DIR.", flush=True)
    sys.exit(1)

  print("LaunchDriver> INFO: TRILINOS_DIR=\"" + os.environ["TRILINOS_DIR"] + "\"", flush=True)

  ds = DetermineSystem(args.build_name, args.supported_systems, force_build_name=True)

  launch_env = get_launch_env(ds.system_name)
  launch_cmd = get_launch_cmd(ds.system_name)
  driver_args = get_driver_args(ds.system_name)

  # Specify, and override the driver script for ATDM ATS2 builds. Note that
  # args.build_name is a required argument so it will be valid by the time it
  # reaches this check.
  if args.build_name.startswith("ats2_cuda"):
      args.driver = "./Trilinos/packages/framework/pr_tools/PullRequestLinuxCudaVortexDriver.sh"

  cmd = launch_env + launch_cmd + args.driver + driver_args

  if args.build_name.startswith("rhel8"):
    cmd += " --on_rhel8"

  if args.in_container:
     cmd += " --no-bootstrap"

  if args.kokkos_develop:
     cmd += " --kokkos-develop"

  if args.extra_configure_args:
     cmd += f" --extra-configure-args=\"{args.extra_configure_args}\""

  print("LaunchDriver> EXEC: " + cmd, flush=True)

  cmd_output = subprocess.run(cmd, shell=True)

  sys.exit(cmd_output.returncode)


if __name__ == "__main__":
    main(sys.argv[1 :])
