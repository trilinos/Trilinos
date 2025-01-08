#/bin/bash
export GCC_VERSION=8.3.0

nohup podman build \
        --build-arg=compiler_version="@${GCC_VERSION}" \
        -t gcc-${GCC_VERSION}-serial-trilinos-env:${USER}-test . \
  &> gcc-${GCC_VERSION}-serial-trilinos-env.output &

