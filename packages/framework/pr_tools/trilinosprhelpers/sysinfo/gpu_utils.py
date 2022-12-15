import json
import os
import subprocess
from shutil import which


def has_nvidia_gpus():
    return bool(list_nvidia_gpus())


def _nvidia_smi():
    if which('nvidia-smi'):
        with open(os.devnull) as errout:
            return subprocess.check_output('nvidia-smi --list-gpus', stderr=errout, shell=True).splitlines()
    return []


def list_nvidia_gpus():
    gpu_ids = []
    try:
        gpu_ids = [str(x) for x in range(0, len(_nvidia_smi()))]
    except Exception as e:
        raise RuntimeError("Failed to acquire list of gpus: {0}".format(str(e)))
    return gpu_ids


def write_ctest_gpu_resource_file(filename, gpu_indices=list_nvidia_gpus(), slots_per_gpu=1):
    content = {
        "version": {
            "major": 1,
            "minor": 0
        },
        "local": [
            {
                "gpus": [{"id": x, "slots": slots_per_gpu} for x in gpu_indices]
            }
        ]
    }
    with open(filename, 'w') as outf:
        json.dump(content, outf)
