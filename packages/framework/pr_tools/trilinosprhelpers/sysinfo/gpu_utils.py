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
