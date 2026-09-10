import argparse
from pathlib import Path
import yaml
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader
import os
from shutil import which
import subprocess
from github import Github
import logging
logger = logging.getLogger(__name__)


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


def get_pr_number(repo, upstream):
    try:
        g = Github()
        repo = g.get_repo(repo)
        pulls = repo.get_pulls(state='open', sort='created', base='develop', head=upstream)
        pr_number = pulls[0].number
        g.close()
        return pr_number
    except:
        return None


def parse_workflow(workflow_file):
    pr_builds = {}

    with open(workflow_file, 'r') as f:
        data = yaml.load(f, Loader=Loader)

    # loop over jobs
    for jobname in data['jobs']:
        logger.debug(f"Checking if \"{jobname}\" is a PR build job that can be reproduced")

        job = data['jobs'][jobname]

        cmake_extra_args = ""
        if (('uses' in job) and (job['uses'] == "./.github/workflows/ci.yml")
              and ('with' in job) and ('image' in job['with']) and ('genconfig-string' in job['with'])):
            withBlock = job['with']
            image, tag = withBlock['image'].split(':')
            genconfig_build_id = withBlock['genconfig-string']
            if "extra-pr-driver-args" in withBlock:
                cmake_extra_args = withBlock["extra-pr-driver-args"].split("=")[-1].replace("${GITHUB_WORKSPACE}", "/workspace/trilinos/source").replace(";", " ")
        else:
            logger.debug("Could not parse job")
            continue

        pr_builds[jobname] = {"image": image,
                              "tag": tag,
                              "genconfig_build_id": genconfig_build_id,
                              "cmake_extra_args": cmake_extra_args}
        logger.debug(f"\"{jobname}\" added")

    return pr_builds


def parse_workflows(source_dir):
    AT2_workflow = source_dir/'.github'/'workflows'/'AT2.yml'
    assert AT2_workflow.exists()
    builds = parse_workflow(AT2_workflow)

    nightly_workflow = source_dir/'.github'/'workflows'/'nightly.yml'
    assert nightly_workflow.exists()
    builds.update(parse_workflow(nightly_workflow))

    return builds


def generate_package_enables(source_dir, packageEnablesFile, pr_base, pr_target, package_list_file=None):
    cmd = [str(source_dir/"commonTools"/"framework"/"get-changed-trilinos-packages.sh"),
           pr_base,
           pr_target,
           str(packageEnablesFile)]
    if package_list_file is not None:
        cmd += [str(package_list_file)]
    out = subprocess.run(cmd,
                         capture_output=True, cwd=packageEnablesFile.parent, universal_newlines=True).stdout[:-1]
    return out


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


def launch_container(container_engine, source_dir, packageEnablesFile, image, tag,
                     genconfig_build_id, extra_cmake_args=""):
    script_location = Path(os.path.abspath(__file__)).parent

    # command
    cmd = [container_engine, "run"]
    # delete container after exit
    cmd += ["--rm"]
    # interactive
    cmd += ["-it"]
    # mount source directory from host
    cmd += ["-v", f"{source_dir}:/workspace/trilinos/source"]
    # set environment variables for source, build and install directories
    cmd += ["-e" "TRILINOS_DIR=/workspace/trilinos/source",
            "-e" "TRILINOS_BUILD_DIR=/workspace/trilinos/build",
            "-e" "TRILINOS_INSTALL_DIR=/workspace/trilinos/install"]
    # mount package enables file
    if packageEnablesFile.exists():
        cmd += ["-v", f"{packageEnablesFile.absolute()}:/workspace/packageEnables.cmake"]
    # set environment variable for GenConfig
    cmd += ["-e", f"GENCONFIG_BUILD_ID={genconfig_build_id}"]
    # add cmake args env variable
    cmd += ["-e", f"CMAKE_EXTRA_ARGS={extra_cmake_args}"]
    # mount script with environment
    cmd += ["-v" f"{script_location}/container_commands.sh:/scripts/container_commands.sh:ro"]
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
    cmd += ["bash", "-c", "echo \". /scripts/container_commands.sh\" >> $HOME/.bashrc; bash"]
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
