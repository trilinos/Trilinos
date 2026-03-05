# Container image for PyTrilinos2

Container image: [`trilinos`](https://github.com/users/trilinos/packages/container/package/trilinos)

The image can also be pulled from the GitHub Container Registry or be built locally.


## Running locally

To use the image locally, please install `podman` or `docker`.

The Dockerfile is [Dockerfile](https://github.com/trilinos/Trilinos/blob/develop/packages/PyTrilinos2/scripts_extra/Dockerfile).

To build the image:
```
podman build --build-arg REPO=https://github.com/trilinos/Trilinos --build-arg BRANCH=develop --build-arg BUILD_PARALLELISM=4 --format=Docker -t ghcr.io/trilinos/trilinos:latest -f Dockerfile
```

The configuration for the build can be found in [trilinos-build.cmake](https://github.com/trilinos/Trilinos/blob/develop/packages/PyTrilinos2/scripts_extra/trilinos-build.cmake).
