#!/usr/bin/env python3

import os
import subprocess
import shutil

import argparse


REGISTRY = "registry-ex.sandia.gov/trilinos-project/trilinos-containers"

DEPLOYS = [
    {"dockerfile": "gnu-openmpi",
    "build_args": {"compiler_version": "@12.1.0", "mpi_version": "@4.1.6"},
    "production": False,
    "image_name": "ubi8-gcc-12.1.0-openmpi-4.1.6"},
    {"dockerfile": "gnu-openmpi",
    "build_args": {"compiler_version": "@10.3.0", "mpi_version": "@4.1.6"},
    "production": False,
    "image_name": "ubi8-gcc-10.3.0-openmpi-4.1.6"},
    {"dockerfile": "cuda-gnu-openmpi",
    "build_args": {"compiler_version": "@10.3.0", "mpi_version": "@4.1.6", "cuda_version": "@11.4.2"},
    "production": False,
    "image_name": "ubi8-cuda-11.4.2-gcc-10.3.0-openmpi-4.1.6"},
    {"dockerfile": "cuda-gnu-openmpi",
    "build_args": {"compiler_version": "@10.3.0", "mpi_version": "@4.1.6", "cuda_version": "@12.4.1"},
    "production": False,
    "image_name": "ubi8-cuda-12.4.1-gcc-10.3.0-openmpi-4.1.6"},
    {"dockerfile": "gnu-serial",
    "build_args": {"compiler_version": "@8.3.0"},
    "production": False,
    "image_name": "ubi8-gcc-8.3.0-serial"},
    {"dockerfile": "python",
    "build_args": {},
    "production": False,
    "image_name": "ubi8-python-3.9"},
    {"dockerfile": "intel-intelmpi",
    "build_args": {},
    "production": True,
    "image_name": "ubi8-intel-intelmpi"}
]


parser = argparse.ArgumentParser()
parser.add_argument("IMAGES", nargs="*")
parser.add_argument("--skip-login", action="store_true", help="Skip registry login")
args = parser.parse_args()

repo_root = os.path.abspath(os.path.dirname(__file__))

if not args.skip_login:
    subprocess.check_call(['docker', 'login', '-u', 'gitlab+deploy-token-23', '-p', os.environ['AT2_BUILD_TOKEN'], 'registry-ex.sandia.gov'])

date_format = "%Y%m%d"
os.chdir(repo_root)

if args.IMAGES:
    deploys = [x for x in DEPLOYS if x["image_name"] in args.IMAGES]
else:
    deploys = DEPLOYS

for image in deploys:
    dockerfile = os.path.join(repo_root, "dockerfiles", image["dockerfile"])
    path = dockerfile + "/Dockerfile"
    # get the timestamp from the Dockerfile commit
    dockerfile_ts = subprocess.check_output(["git", "--no-pager", "log", "-1", "--format=%cd", f"--date=format:{date_format}", "--", path]).decode().strip()
    build_args = [k + "=" + v for k, v in image["build_args"].items()]
    tag = REGISTRY + ("/production/" if image["production"] else "/experimental/") + image["image_name"] + ":" + dockerfile_ts
    build_args.append(f"AT2_image_fullpath={tag}")
    build_args.append(f"AT2_image={image['image_name']}")
    print(f"Building Dockerfile '{dockerfile}' with args {build_args}")
    print(f"Tagging as {tag}")
    f = []
    for e in build_args:
        f.append("--build-arg")
        f.append(e)
    if os.path.exists(os.path.join(dockerfile, "GenConfig")):
        shutil.rmtree(os.path.join(dockerfile, "GenConfig"))
    if not os.path.exists(os.path.join(repo_root, "GenConfig", "ini_files")):
        raise SystemExit("The GenConfig dependency does not appear to exist!  Did you initialize git submodules?")
    shutil.copytree(os.path.join(repo_root, "GenConfig"), os.path.join(dockerfile, "GenConfig"), symlinks=True)
    try:
        subprocess.check_call(["podman", "build", "--tag", tag] + f + [dockerfile])
    except subprocess.CalledProcessError as e:
        print(f"check_call() returned {e.returncode}")
        continue
    subprocess.check_call(["podman", "push", tag])
    latest_tag = REGISTRY + ("/production/" if image["production"] else "/experimental/") + image["image_name"] + ":latest"
    subprocess.check_call(["podman", "push", tag, latest_tag])
