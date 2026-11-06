import argparse
from pathlib import Path
import yaml
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader
import json
import os
from shutil import which
import subprocess
import logging
logger = logging.getLogger(__name__)


def findGitRoot(p):
    if (p/".git").exists():
        return p
    else:
        if p.parent is None:
            return None
        else:
            return findGitRoot(p.parent)


def is_valid_pr_number(s):
    try:
        pr_number = int(s)
        return pr_number >= 1
    except:
        return False


def convert_to_pr_number(v):
    try:
        pr_number = int(v)
    except:
        raise argparse.ArgumentTypeError(f"Cannot convert \"{v}\" to a valid PR number.")
    if pr_number >= 1:
        return pr_number
    else:
        raise argparse.ArgumentTypeError(f"Cannot convert \"{v}\" to a valid PR number.")


def is_valid_source_dir(p):
    return Path(p).exists() and (Path(p)/'.github'/'workflows'/'AT2.yml').exists()


def convert_to_valid_source_dir(p):
    if is_valid_source_dir(p):
        return Path(p)
    else:
        raise argparse.ArgumentTypeError(f"Cannot convert \"{p}\" to a valid source directory.")


def get_github_upstream_of_local_branch(source_dir):
    try:
        git_upstream = subprocess.run(["git", "rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"],
                                      capture_output=True, cwd=source_dir, universal_newlines=True).stdout[:-1]
        git_upstream, git_upstream_branch = git_upstream.split('/')
        git_push_url = subprocess.run(["git", "remote", "get-url", "--push", git_upstream],
                                      capture_output=True, cwd=source_dir, universal_newlines=True).stdout[:-1]
        git_repo = git_push_url.split(":")[-1]
        git_owner = git_repo.split('/')[0]
        head = f"{git_owner}:{git_upstream_branch}"
        return head
    except:
        return None


def parse_workflow(workflow_file):
    images = []

    with open(workflow_file, 'r') as f:
        data = yaml.load(f, Loader=Loader)

    # loop over jobs
    for jobname in data['jobs']:
        logger.debug(f"Checking if \"{jobname}\" is a PR build job that can be reproduced")

        job = data['jobs'][jobname]
        if (('uses' in job) and (job['uses'] == "./.github/workflows/ci.yml")
              and ('with' in job) and ('image' in job['with'])):
            withBlock = job['with']
            image, tag = withBlock['image'].split(':')
        else:
            logger.debug(f"Could not parse job {jobname}")
            continue

        images.append((image, tag))
        logger.debug(f"\"{jobname}\" added")

    return set(images)


def parse_workflow_images(source_dir, container_engine):
    if container_engine in ("podman", "docker"):
        images = set()
        for workflow in [source_dir/'.github'/'workflows'/'AT2.yml',
                         source_dir/'.github'/'workflows'/'nightly.yml',
                         source_dir/'.github'/'workflows'/'weekly.yml']:
            assert workflow.exists()
            images |= parse_workflow(workflow)
        return images
    elif container_engine == "cask":
        images = set()
        # TODO: hook up cask ready images
        raise NotImplementedError()
        return images
    else:
        raise NotImplementedError(f"Unknown container_engine: \"{container_engine}\"")


def image_available(container_engine, image, tag):
    if container_engine == "podman":
        cmd = [container_engine, "image", "exists", f"{image}:{tag}"]
    else:
        cmd = [container_engine, "image", "inspect", f"{image}:{tag}"]
    ret = subprocess.run(cmd,
                         capture_output=True, universal_newlines=True)
    return ret.returncode == 0


def pull_image(container_engine, image, tag):
    cmd = [container_engine, "pull", f"{image}:{tag}"]
    subprocess.run(cmd,
                   capture_output=False, universal_newlines=True)


def launch_container(container_engine, source_dir, build_dir, image, tag, unparsed_args):
    """
    We map:
     source_dir -> /workspace/trilinos/source
     build_dir  -> /workspace/trilinos/build
    """
    script_location = Path(os.path.abspath(__file__)).parent

    if container_engine in ("podman", "docker"):
        # command
        cmd = [container_engine, "run"]
        # delete container after exit
        cmd += ["--rm"]
        # interactive
        cmd += ["-it"]
        # mount source directory from host
        cmd += ["-v", f"{source_dir}:/workspace/trilinos/source"]
        # mount build directory from host
        cmd += ["-v", f"{build_dir}:/workspace/trilinos/build"]
        # set environment variables for source, build and install directories
        cmd += ["-e" "TRILINOS_DIR=/workspace/trilinos/source",
                "-e" "TRILINOS_BUILD_DIR=/workspace/trilinos/build",
                "-e" "TRILINOS_INSTALL_DIR=/workspace/trilinos/install"]
        # workdir
        cmd += ["--workdir", "/workspace/trilinos"]
        # map GPU
        if which("nvidia-smi") is not None and image.find("cuda") >= 0:
            if container_engine == "podman":
                cmd += ["--device", "nvidia.com/gpu=all"]
            else:
                cmd += ["--gpus", "all"]

        # deal with certificates
        if Path("/etc/ssl/certs/ca-bundle.crt").exists():
            cmd += ["-v", "/etc/ssl/certs/ca-bundle.crt:/certs/ca-bundle.crt:ro",
                    "-e", "SSL_CERT_FILE=/certs/ca-bundle.crt",
                    "-e", "PIP_CERT=/certs/ca-bundle.crt",
                    "-e", "CURL_CA_BUNDLE=/certs/ca-bundle.crt",
                    "-e", "REQUESTS_CA_BUNDLE=/certs/ca-bundle.crt"]
        elif Path("/etc/ssl/certs/ca-certificates.crt").exists():
            cmd += ["-v", "/etc/ssl/certs/ca-certificates.crt:/certs/ca-certificates.crt:ro",
                    "-e", "SSL_CERT_FILE=/certs/ca-certificates.crt",
                    "-e", "PIP_CERT=/certs/ca-certificates.crt",
                    "-e", "CURL_CA_BUNDLE=/certs/ca-certificates.crt",
                    "-e", "REQUESTS_CA_BUNDLE=/certs/ca-certificates.crt"]
        else:
            if "SSL_CERT_FILE" in os.environ:
                ssl_cert_file = os.environ["SSL_CERT_FILE"]
                cmd += ["-v", f"{ssl_cert_file}:/certs/ssl-cert.crt:ro",
                        "-e", "SSL_CERT_FILE=/certs/ssl-cert.crt"]
            if "PIP_CERT" in os.environ:
                pip_cert = os.environ["PIP_CERT"]
                cmd += ["-v", f"{pip_cert}:/certs/pip-cert.crt:ro",
                        "-e", "PIP_CERT=/certs/pip-cert.crt"]
            if "CURL_CA_BUNDLE" in os.environ:
                curl_ca_bundle = os.environ["CURL_CA_BUNDLE"]
                cmd += ["-v", f"{curl_ca_bundle}:/certs/curl-cert.crt:ro",
                        "-e", "CURL_CA_BUNDLE=/certs/curl-cert.crt"]
            if "REQUESTS_CA_BUNDLE" in os.environ:
                requests_ca_bundle = os.environ["REQUESTS_CA_BUNDLE"]
                cmd += ["-v", f"{requests_ca_bundle}:/certs/requests-cert.crt:ro",
                        "-e", "REQUESTS_CA_BUNDLE=/certs/requests-cert.crt"]

        # image to run
        cmd += [f"{image}:{tag}"]
        # command to be executed
        cmd += ["bash", "-c", "echo \". /workspace/trilinos/source/sampleScripts/configurations/helpers.sh\" >> $HOME/.bashrc; bash"]
        logger.debug(f"{container_engine} command = {cmd}")
        ret = subprocess.run(cmd,
                             capture_output=False, universal_newlines=True)
        if ret.returncode != 0:
            message = f"""
            The {container_engine} container launch failed. Unfortunately there are lots of situations
            in which this can happen. Some suggestions for narrowing down the problem:

            - The container engine installation is simply broken. Try running
                {container_engine} run -it --rm rockylinux/rockylinux bash -c "echo \\"That worked\\""
              to run a minimal Linux container as a test.

            - Docker may require elevated privileges. Try running the command with sudo,
              or ensure your user belongs to the docker group. Podman might instead be
              provided by a module (for example, module load aue/podman).
            """
            print(message)
            exit(1)
    elif container_engine == "cask":
        cmd = [container_engine]
        cmd += unparsed_args
        # mount source directory from host
        cmd += ["--dir", f"{source_dir}:/workspace/trilinos/source"]
        # mount build directory from host
        # TODO: Check that CASK supports multiple "--dir" args. If not, print a warning and use a container-only build dir
        cmd += ["--dir", f"{build_dir}:/workspace/trilinos/build"]
        cmd += ["--agents_image", f"{image}:{tag}"]
        # map GPU
        if which("nvidia-smi") is not None and image.find("cuda") >= 0:
            cmd += ["--gpu", "nvidia"]
        logger.debug(f"{container_engine} command = {cmd}")
        ret = subprocess.run(cmd,
                             capture_output=False, universal_newlines=True)
    else:
        raise NotImplementedError(f"Unknown container engine: \"{container_engine}\"")


def ensure_podman_machine():
    proc = subprocess.run(["podman", "machine", "list", "--format", "json"], capture_output=True, text=True)
    assert proc.returncode == 0

    try:
        machines = json.loads(proc.stdout.strip())
    except json.JSONDecodeError:
        # Fallback if older podman version doesn't support json format perfectly
        machines = []

    # Handle initialization if no machines exist
    if not machines or len(machines) == 0:
        logger.debug("No Podman machine found. Initializing 'podman-machine-default'...")
        proc = subprocess.run(["podman", "machine", "init"], capture_output=False, text=True)
        assert proc.returncode == 0
        logger.debug("Initialization successful.")
        # Reload machine list to include the newly created machine
        proc = subprocess.run(["podman", "machine", "list", "--format", "json"], capture_output=True, text=True)
        assert proc.returncode == 0
        try:
            machines = json.loads(proc.stdout.strip())
        except json.JSONDecodeError:
            exit(1)

    # Check running state of the default/first machine
    target_machine = machines[0]
    machine_name = target_machine.get("Name", "podman-machine-default")
    is_running = target_machine.get("Running", False)

    if is_running:
        logger.debug(f"Podman machine '{machine_name}' is already running.")
    else:
        logger.debug(f"Podman machine '{machine_name}' is stopped. Starting it now.")
        proc = subprocess.run(["podman", "machine", "start"], capture_output=False, text=True)
        if proc.returncode != 0:
            # Handle edge case where it's already starting or active but misreported
            if "already running" in proc.stderr.strip().lower():
                logger.debug(f"Podman machine '{machine_name}' is now active.")
            else:
                logger.debug(f"Failed to start machine")
                exit(1)
        else:
            logger.debug(f"Podman machine '{machine_name}' started successfully.")
