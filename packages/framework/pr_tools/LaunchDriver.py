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



def get_launch_env(build_name : str, system : str):
  """
  Gets the launch environment based on the detected system.
  This is an early environment that's required for running the driver.
  
  Returns:
      str: The environment used to launch the driver.
  """
  env = ""
  if "_rdc" in build_name:
      env += " TRILINOS_MAX_CORES=10"

  if system == "weaver" or system == "ats2":
      env += " Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE"

  if env == "":
      return ""
  else:
      return "env" + env + " "


def get_launch_cmd(build_name : str, system : str):
  """
  Gets the launch command based on the detected system.
  
  Returns:
      str: The command used to launch the driver.
  """
  if system == "weaver" or system == "ats2":
    cmd = "bsub -Is -J " + build_name + " -W 12:00"
  else:
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
  the trilinos driver script appropriatly.

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
  args = parser.parse_args(argv)

  if os.getenv("TRILINOS_DIR") == None:
    print("LaunchDriver> ERROR: Please set TRILINOS_DIR.", flush=True)
    sys.exit(1)

  print("LaunchDriver> INFO: TRILINOS_DIR=\"" + os.environ["TRILINOS_DIR"] + "\"", flush=True)

  ds = DetermineSystem(args.build_name, args.supported_systems)

  launch_env = get_launch_env(args.build_name, ds.system_name)
  launch_cmd = get_launch_cmd(args.build_name, ds.system_name)
  driver_args = get_driver_args(ds.system_name)

  cmd = launch_env + launch_cmd + args.driver + driver_args

  print("LaunchDriver> EXEC: " + cmd, flush=True)

  cmd_output = subprocess.run(cmd, shell=True)

  sys.exit(cmd_output.returncode)


if __name__ == "__main__":
    main(sys.argv[1 :])

