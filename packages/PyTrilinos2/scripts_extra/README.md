# Build scripts

This directory contains scripts that can be used to build PyTrilinos2.


## Docker image and container

The provided Dockerfile allows to build Trilinos and PyTrilinos2.

Run

```
podman build --build-arg REPO=https://github.com/trilinos/Trilinos --build-arg BRANCH=develop --build-arg BUILD_PARALLELISM=4 --format=Docker -t ghcr.io/trilinos/trilinos:latest -f Dockerfile

```
