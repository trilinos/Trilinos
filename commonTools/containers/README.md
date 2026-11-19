Trilinos development environment in a container
==========================================================

Requirements: Python (>= 3.10), podman or docker

Run
```bash
commonTools/containers/dev.sh
```
to launch a container image for Trilinos development. Source code 
and build directory get mounted in from the host systems.

To see the available commandline flags, run
```bash
commonTools/containers/dev.sh -h
```

