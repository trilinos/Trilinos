Environments and configurations for Trilinos builds
================================================================

This folder contains:

- `get_dependencies.sh` script

  This script pulls in all dependencies for managing runtime environments and build configurations.


- `ini-files` folder

  The folder contains environments and configurations for Trilinos PR builds (in containers) and select open systems.


Usage
======

To initialize the scrips run
```
$TRILINOS_SRC/packages/framework/get_dependencies.sh
```
All scripts will be checked out in the Trilinos source directory under `packages/framework`.

Two scripts are particularly useful: LoadEnv and GenConfig. 
The scripts map a configuration string with keywords that are separated by underscores to an environment or a build configuration respectively.

For example:
```
source $TRILINOS_SRC/packages/framework/GenConfig/LoadEnv/load-env.sh frontier_mi250x
source $TRILINOS_SRC/packages/framework/GenConfig/gen-config.sh       frontier_mi250x_release_static_AMDgfx90a_Zen3_no-asan_no-complex_fpic_mpi_no-pt_no-rdc_no-uvm_deprecated-on_nightly-performance $TRILINOS_SRC
```

LoadEnv: Managing environments
--------------------------------------

LoadEnv manages runtime environments (environment variables, modules, etc).


Getting help:
```
source $TRILINOS_SRC/packages/framework/GenConfig/LoadEnv/load-env.sh --help
```

Listing the available environments on the current system:
```
source $TRILINOS_SRC/packages/framework/GenConfig/LoadEnv/load-env.sh --list-envs
```

LoadEnv tries to identify the system by its hostname. This mechanism can be overridden by
```
source $TRILINOS_SRC/packages/framework/GenConfig/LoadEnv/load-env.sh --list-envs --force <SYSTEM_NAME>
```

Loading an environment in a new subshell:
```
source $TRILINOS_SRC/packages/framework/GenConfig/LoadEnv/load-env.sh <ENV_STRING>
```

Loading an environment in the current shell:
```
source $TRILINOS_SRC/packages/framework/GenConfig/LoadEnv/load-env.sh --ci-mode <ENV_STRING>
```


GenConfig: Managing build configurations and environments
-----------------------------------------------------------------------

In addition to loading an environment with LoadEnv, GenConfig also handles build configurations.

Getting help:
```
source $TRILINOS_SRC/packages/framework/GenConfig/gen-config.sh --help
```

Listing the available configurations on the current system:
```
source $TRILINOS_SRC/packages/framework/GenConfig/gen-config.sh --list-configs
```

Listing  available configuration flags:
```
source $TRILINOS_SRC/packages/framework/GenConfig/gen-config.sh --list-config-flags
```

Loading an environment in a new subshell and configuring a build in the current directory:
```
source $TRILINOS_SRC/packages/framework/GenConfig/gen-config.sh <CONFIG_STRING> $TRILINOS_SRC
```

Loading an environment in the current shell and configuring a build in the current directory:
```
source $TRILINOS_SRC/packages/framework/GenConfig/gen-config.sh --ci-mode <CONFIG_STRING> $TRILINOS_SRC
```
