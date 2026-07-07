import re
import git
import os
import logging
logger = logging.getLogger(__name__)

def get_git_root(path: str) -> str:
    """Get the root directory of a git repository.

    Args:
        path: Path to a file or directory within a git repository
    Returns:
        Absolute path to the git repository root
    """
    logger.debug(f"Finding git root for {path}")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Path does not exist: {path}")

    try:
        git_repo = git.Repo(path, search_parent_directories=True)
        git_root = git_repo.git.rev_parse("--show-toplevel")
        logger.debug(f"git root = {git_root}")
        return git_root
    except git.InvalidGitRepositoryError as e:
        raise ValueError(f"Path is not within a git repository: {path}")


def parse_semver(version: str) -> dict:
    """
    Parse a given semver string into its semver components.

    Args:
        version: Semantic version string
    Returns:
        Dictionary of semver regex matched groups.
    Raises:
        ValueError: If version string is invalid semver string.
    """
    logger.debug(f"Semver parsing {version}")
    # https://semver.org/#is-there-a-suggested-regular-expression-regex-to-check-a-semver-string
    SEMVER_REGEX = r"""^(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)(?:-(?P<prerelease>(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+(?P<buildmetadata>[0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$
"""
    SEMVER_RE = re.compile(SEMVER_REGEX, re.VERBOSE)
    m = SEMVER_RE.match(version)
    if not m:
        raise ValueError(f"Unable to parse semantic versioning from {version}")

    mgroups = m.groupdict()
    logger.debug(f"{version} parsed into {mgroups}")
    return mgroups


def fetch_branch(branch: str, repo_path: str, merge: bool=False) -> None:
    """Fetch a git branch from remote.

    Args:
        branch: Name of the branch to fetch
        repo_path: Path to the git repository (default: current directory)
        merge: If True, merge the fetched branch with local branch after fetching
    """
    try:
        git_root = get_git_root(repo_path)
        git_repo = git.Repo(git_root)

        logger.info(f"Fetching latest changes for repository at {git_root}")
        origin = git_repo.remotes.origin
        origin.fetch(branch)
        if merge:
            if branch in [b.name for b in git_repo.branches]:
                current_branch = git_repo.active_branch.name

                if current_branch != branch:
                    logger.info(f"Checking out branch {branch} before merging")
                    git_repo.git.checkout(branch)

                logger.info(f"Merging fetched branch {branch} with local branch")
                git_repo.git.merge(f"origin/{branch}")
            else:
                logger.info(f"Branch {branch} doesn't exist locally, cannot merge")
    except git.InvalidGitRepositoryError as e:
        raise ValueError(f"Invalid git repository at {repo_path}: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Failed to fetch branch {branch}: {e}") from e





def checkout_branch(branch: str, repo_path: str, remote: bool=False) -> None:
    """Checkout a git branch.

    Args:
        branch: Name of the branch to checkout
        repo_path: Path to the git repository (default: current directory)
        remote: If True, pull latest changes from remote after checkout
    """
    try:
        git_root = get_git_root(repo_path)
        git_repo = git.Repo(git_root)

        if branch in [b.name for b in git_repo.branches]:
            logger.info(f"Branch {branch} exists locally, checking it out")
            git_repo.git.checkout(branch)
        else:
            logger.info(f"Checking out new branch {branch}")
            git_repo.git.checkout(b=branch)

        if remote:
            logger.info(f"Remote specified, pulling latest remote changes of {branch}")
            git_repo.git.pull('origin', branch)

    except git.InvalidGitRepositoryError as e:
        raise ValueError(f"Invalid git repository at {repo_path}: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Failed to checkout branch {branch}: {e}") from e




