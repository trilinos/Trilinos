#!/usr/bin/env python
import json
import os
import subprocess
from argparse import ArgumentParser
from shutil import which


def has_nvidia_gpus():
    gpu_ids = []
    if which('nvidia-smi'):
        try:
            with open(os.devnull) as errout:
                gpu_ids = [str(x) for x in range(0, len(subprocess.check_output('nvidia-smi --list-gpus', stderr=errout, shell=True).splitlines()))]
        except Exception as e:
            raise RuntimeError("Failed to acquire list of gpus: {0}".format(str(e)))
    return gpu_ids != []


def write_ctest_gpu_resource_file(filename, gpu_indices, slots_per_gpu=1):
    content = {
        "version": {
            "major": 1,
            "minor": 0
        },
        "local": [
            {
                "gpus": [{"id": x, "slots": 1} for x in gpu_indices]
            }
        ]
    }
    with open(filename, 'w') as outf:
        json.dump(content, outf)
