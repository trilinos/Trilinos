from pathlib import Path
from shutil import which
import subprocess
from argparse import ArgumentParser

def main():

    ##################################################

    parser = ArgumentParser(description="A script for making an existing container image ready for CASK.")
    parser.add_argument("sourceImage", help="The name of the container image that should be CASK-ified")
    parser.add_argument("targetImage", help="The name of the CASK-ified image")
    parser.add_argument("--os", default="ubi", choices=["ubi", "ubuntu"], help="OS type of sourceImage")
    parser.add_argument("--os_version", default="", type=str, help="OS version. Currently not used by CASK images")
    parser.add_argument("--no-nexus", action="store_false", dest="useNexus", help="Skip usage of Nexus and instead download directly from the internet")
    args = parser.parse_args()

    ##################################################
    # Check CASK files

    caskExecutable = which("cask")
    if caskExecutable is None:
        common_cask_location = Path.home()/"cask"/"launchers"/"cask"
        if common_cask_location.exists():
            caskExecutable = common_cask_location
    assert caskExecutable is not None, "Cask executable could not be found in PATH"
    caskExecutable = Path(caskExecutable)
    caskPath = caskExecutable.parent.parent

    containerfileBaselineUbuntu = caskPath/"containers"/"baselines"/"Containerfile.ubuntu"
    containerfileBaselineUbi = caskPath/"containers"/"baselines"/"Containerfile.ubi"
    assert containerfileBaselineUbuntu.exists(), f"{containerfileBaselineUbuntu} does not exist"
    assert containerfileBaselineUbi.exists(), f"{containerfileBaselineUbi} does not exist"

    containerfileAgentBase = caskPath/"containers"/"agents"/"Containerfile.agent_base"
    assert containerfileAgentBase.exists(), f"{containerfileAgentBase} does not exist"

    containerfileAgentTuis = caskPath/"containers"/"agents"/"Containerfile.agent_tuis"
    assert containerfileAgentTuis.exists(), f"{containerfileAgentTuis} does not exist"

    ##################################################

    def buildImage(containerfile, targetImage, buildArgs):
        containerfile = containerfile.relative_to(caskPath)
        args = ["--http-proxy=false"]
        env = {"HTTP_PROXY": "http://proxy.sandia.gov:80",
               "HTTPS_PROXY": "http://proxy.sandia.gov:80",
               "http_proxy": "http://proxy.sandia.gov:80",
               "https_proxy": "http://proxy.sandia.gov:80"}

        args += [f"--build-arg {arg}={val}" for arg, val in buildArgs.items()]
        args += [f"--env {arg}=\"{val}\"" for arg, val in env.items()]
        args += [f"-f {containerfile}"]
        args += [f"-t {targetImage}"]
        args += ["."]
        cmd = "podman build " + " ".join(args)
        print(cmd)
        ret = subprocess.run(cmd, capture_output=False, text=True, shell=True, cwd=caskPath)
        assert ret.returncode == 0


    def rmi(image):
        cmd = f"podman rmi {image}"
        print(cmd)
        ret = subprocess.run(cmd, capture_output=False, text=True, shell=True)
        assert ret.returncode == 0


    if args.os == "ubuntu":
        containerfileBaseline = containerfileBaselineUbuntu
    elif args.os == "ubi":
        containerfileBaseline = containerfileBaselineUbi
    else:
        raise NotImplementedError(args.os)


    ##################################################
    # baseline image

    miniforgeVersion = "26.3.2-0"
    ripgrepVersion = "14.1.1"
    gosuVersion = "1.17"

    buildArgs = {"BASE_IMAGE": args.sourceImage,
                 "OS_VERSION": args.os_version}
    if not args.useNexus:
        buildArgs.update({"MINIFORGE_VERSION": miniforgeVersion,
                          "MINIFORGE_URL": f"https://github.com/conda-forge/miniforge/releases/download/{miniforgeVersion}/Miniforge3-{miniforgeVersion}-Linux-x86_64.sh",
                          "MINIFORGE_ARCH": "x86_64",
                          "RIPGREP_VERSION": ripgrepVersion,
                          "RIPGREP_URL": f"https://github.com/BurntSushi/ripgrep/releases/download/{ripgrepVersion}/ripgrep-{ripgrepVersion}-x86_64-unknown-linux-musl.tar.gz",
                          "GOSU_VERSION": gosuVersion,
                          "GOSU_URL": f"https://github.com/tianon/gosu/releases/download/{gosuVersion}/gosu-amd64"})
    baselineImage = f"{args.sourceImage}_temp0"
    buildImage(containerfileBaseline, baselineImage, buildArgs)

    ##################################################
    # Agent base image

    buildArgs = {"BASE_IMAGE": baselineImage}
    agentBaseImage = f"{args.sourceImage}_temp1"
    buildImage(containerfileAgentBase, agentBaseImage, buildArgs)

    ##################################################
    # Agent TUIs image

    buildArgs = {"AGENT_BASE_IMAGE": agentBaseImage}
    if not args.useNexus:
        buildArgs["NPM_REGISTRY"] = "https://registry.npmjs.org"
    buildImage(containerfileAgentTuis, args.targetImage, buildArgs)

    ##################################################
    # remove intermediate images

    rmi(baselineImage)
    rmi(agentBaseImage)
