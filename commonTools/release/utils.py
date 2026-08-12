import re
import git
import os
from github import (Github,
                    Auth,
                    PullRequest,
                    BranchProtection,
                    GithubException)
import logging
logger = logging.getLogger(__name__)


ORG_REPO = "trilinos-cicd2/Trilinos-test"


def get_git_root(path: str) -> str:
    """Get the root directory of a git repository.

    Args:
        path: Path to a file or directory within a git repository
    Returns:
        Absolute path to the git repository root
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Path does not exist: {path}")

    try:
        git_repo = git.Repo(path, search_parent_directories=True)
        git_root = git_repo.git.rev_parse("--show-toplevel")
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

        origin = git_repo.remotes.origin
        origin.fetch(branch)
        if merge:
            if branch in [b.name for b in git_repo.branches]:
                current_branch = git_repo.active_branch.name

                if current_branch != branch:
                    git_repo.git.checkout(branch)

                git_repo.git.merge(f"origin/{branch}")
            else:
                git_repo.git.checkout("-b", branch, f"origin/{branch}")
    except git.InvalidGitRepositoryError as e:
        raise ValueError(f"Invalid git repository at {repo_path}: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Failed to fetch branch {branch}: {e}") from e


def verify_remote_branch_exists(branch: str, repo_path: str) -> None:
    """Verify that a branch exists on the origin remote.
    """
    try:
        git_root = get_git_root(repo_path)
        git_repo = git.Repo(git_root)

        refs = git_repo.git.ls_remote("--heads", "origin", branch)
        logger.debug(f"ls-remote for {branch}:\n{refs}")
        if not refs.strip():
            raise RuntimeError(
                f"Branch '{branch}' does not exist on the origin remote."
            )
    except git.InvalidGitRepositoryError as e:
        raise ValueError(f"Invalid git repository at {repo_path}: {e}") from e
    except RuntimeError:
        raise
    except Exception as e:
        raise RuntimeError(f"Failed to verify remote branch {branch}: {e}") from e


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
            git_repo.git.checkout(branch)
            if remote:
                git_repo.git.pull("origin", branch, "--rebase")
        else:
            if remote:
                git_repo.git.checkout('--track', f"origin/{branch}")
            else:
                git_repo.git.checkout(b=branch)

    except git.InvalidGitRepositoryError as e:
        raise ValueError(f"Invalid git repository at {repo_path}: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Failed to checkout branch {branch}: {e}") from e


def update_version_cmake(version: str, rel_branch: str, dev_mode: bool, repo_path: str):
    """
    Update the Version.cmake file in root of the Trilinos repo.

    Example Version.cmake in Trilinos:

    ```
    SET(Trilinos_VERSION 17.0.0)
    SET(Trilinos_MAJOR_VERSION 17)
    SET(Trilinos_MAJOR_MINOR_VERSION 170000)
    SET(Trilinos_VERSION_STRING "17.0.0-dev")
    SET(Trilinos_ENABLE_DEVELOPMENT_MODE_DEFAULT ON)
    ```

    when called to update for 17.1.1 in Release Mode,
    Version.cmake gets updated to:

    ```
    SET(Trilinos_VERSION 17.1.1)
    SET(Trilinos_MAJOR_VERSION 17)
    SET(Trilinos_MAJOR_MINOR_VERSION 170101)
    SET(Trilinos_VERSION_STRING "17.1.1")
    SET(Trilinos_ENABLE_DEVELOPMENT_MODE_DEFAULT OFF)
    ```
    """
    try:
        file_path = f"{repo_path}/Version.cmake"
        with open(file_path, 'r') as f:
            content = f.read()
    except Exception as e:
        raise RuntimeError(f"Failed to read Version.cmake: {e}") from e

    semver = parse_semver(version)
    major = semver['major']
    minor = semver['minor']
    patch = semver['patch']
    major_minor_patch = int(major) * 10000 + int(minor) * 100 + int(patch)
    base_version = version.split("-")[0]

    dev_mode_default = "ON" if dev_mode else "OFF"

    new_content = content
    new_content = re.sub(r'SET\(Trilinos_VERSION \d+\.\d+\.\d+\)', f'SET(Trilinos_VERSION {base_version})', new_content)
    new_content = re.sub(r'SET\(Trilinos_MAJOR_VERSION \d+\)', f'SET(Trilinos_MAJOR_VERSION {major})', new_content)
    new_content = re.sub(r'SET\(Trilinos_MAJOR_MINOR_VERSION \d+\)', f'SET(Trilinos_MAJOR_MINOR_VERSION {major_minor_patch})', new_content)
    new_content = re.sub(r'SET\(Trilinos_VERSION_STRING "[^"]+"\)', f'SET(Trilinos_VERSION_STRING "{version}")', new_content)
    new_content = re.sub(r'SET\(Trilinos_ENABLE_DEVELOPMENT_MODE_DEFAULT (?:ON|OFF)\)', f'SET(Trilinos_ENABLE_DEVELOPMENT_MODE_DEFAULT {dev_mode_default})', new_content)

    if not dev_mode:
        new_content = re.sub(r'SET\(Trilinos_REPOSITORY_BRANCH "[^"]+" CACHE INTERNAL ""\)', f'SET(Trilinos_REPOSITORY_BRANCH "{rel_branch}" CACHE INTERNAL "")', new_content)

    try:
        with open(file_path, 'w') as f:
            f.write(new_content)
    except Exception as e:
        raise RuntimeError(f"Failed to update Version.cmake: {e}") from e


def read_version_cmake(repo_path: str) -> dict:
    """
    Read the current version information from the Version.cmake file in the
    root of the Trilinos repo.

    Args:
        repo_path: Path to the git repository
    Returns:
        Dictionary with keys:
            - "version": value of SET(Trilinos_VERSION X.Y.Z)
            - "version_string": value of SET(Trilinos_VERSION_STRING "...")
            - "dev_mode": True if Trilinos_ENABLE_DEVELOPMENT_MODE_DEFAULT is ON,
              False if OFF
    """
    try:
        file_path = f"{repo_path}/Version.cmake"
        with open(file_path, 'r') as f:
            content = f.read()
    except Exception as e:
        raise RuntimeError(f"Failed to read Version.cmake: {e}") from e

    version_m = re.search(r'SET\(Trilinos_VERSION (\d+\.\d+\.\d+)\)', content)
    version_string_m = re.search(r'SET\(Trilinos_VERSION_STRING "([^"]+)"\)', content)
    dev_mode_m = re.search(r'SET\(Trilinos_ENABLE_DEVELOPMENT_MODE_DEFAULT (ON|OFF)\)', content)

    if not version_m or not version_string_m or not dev_mode_m:
        raise RuntimeError(
            f"Failed to parse expected version fields from {file_path}"
        )

    info = {
        "version": version_m.group(1),
        "version_string": version_string_m.group(1),
        "dev_mode": dev_mode_m.group(1) == "ON",
    }
    logger.debug(f"Read Version.cmake: {info}")
    return info


def is_release_mode(info: dict) -> bool:
    """
    Determine whether a Version.cmake info dict (see read_version_cmake) is in
    release mode: the version string carries no "-dev" (or other) prerelease
    suffix and development mode is disabled.
    """
    return info["dev_mode"] is False and "-" not in info["version_string"]


def commit_tracked(msg: str, repo_path: str):
    """
    Commit all modified tracked files in the working tree with the given message.
    """
    try:
        git_root = get_git_root(repo_path)
        git_repo = git.Repo(git_root)

        git_repo.git.add(update=True)

        git_status = git_repo.git.status()
        logger.debug(f"Git status after add:\n{git_status}")

        git_repo.git.commit('-s', '-m', msg)

        latest_log = git_repo.git.log('-1', '--pretty=fuller')
        logger.debug(f"Last commit:\n{latest_log}")

    except Exception as e:
        raise RuntimeError(f"Failed to commit: {e}") from e


def push(branch: str, repo_path: str, remote: str="origin", force: bool=False):
    try:
        git_root = get_git_root(repo_path)
        git_repo = git.Repo(git_root)

        if force:
            git_repo.git.push(origin, branch)
        else:
            git_repo.git.push(remote, '-f' ,branch)
    except Exception as e:
        raise RuntimeError(f"Failed to push: {e}") from e


def create_pull_request(base: str, head: str, title: str, body: str) -> PullRequest:
    """
    Open a pull request from head branch to base branch in the origin upstream repository.

    Args:
        base: target branch to merge head branch into
        head: source branch to merge into target branch. Use <fork:branch> form for branches from forks.
        title: title of opened pull request
        body: body text of opened pull request
    """
    token = os.getenv('GITHUB_TOKEN')
    if not token:
        raise RuntimeError(f"GITHUB_TOKEN environment variable not set.")

    gh = Github(auth=Auth.Token(os.environ['GITHUB_TOKEN']))

    if not base or not head or not title:
        raise ValueError("base, head, and title must be non-empty strings.")

    try:
        repo = gh.get_repo(ORG_REPO)
        pr = repo.create_pull(base=base, head=head, title=title, body=body)
        logger.debug(f"Open pull request at {pr.html_url}")
        return pr
    except GithubException as e:
        raise RuntimeError(f"GitHub API error when creating PR") from e
    finally:
        gh.close()


def set_release_branch_protection(rel_branch: str):
    """
    Set the GitHub branch protection rules for a release branch.

    This branch's required status checks are copied from the "develop" branch. All other branch protection rules are hard-coded to the ususual release branch rules.
    """
    token = os.getenv('GITHUB_TOKEN')
    if not token:
        raise RuntimeError(f"GITHUB_TOKEN environment variable not set.")

    try:
        gh = Github(auth=Auth.Token(os.environ['GITHUB_TOKEN']))
        repo = gh.get_repo(ORG_REPO)

        rel_branch_info = repo.get_branch(rel_branch)
        if rel_branch_info.protected:
            raise RuntimeError(f"Branch protection rules already exist for branch {rel_branch}.")

        # Get develop branch protection rules to extract current required status checks
        source_branch = "develop"
        source_branch_prots = repo.get_branch(source_branch).get_protection()

        required_status_checks = source_branch_prots.required_status_checks
        strict = required_status_checks.strict
        contexts = required_status_checks.contexts
        checks = [(c.context, c.app_id) for c in required_status_checks.checks]

        kwargs = {}

        kwargs["strict"] = strict
        kwargs["contexts"] = list(contexts)
        kwargs["checks"] = list(checks)
        kwargs["dismiss_stale_reviews"] = True
        kwargs["required_approving_review_count"] = 1
        kwargs["enforce_admins"] = True

        kwargs["team_push_restrictions"] = []
        kwargs["block_creations"] = True

        # Create branch protection rule for the release branch
        rel_branch_info.edit_protection(**kwargs)

    except GithubException as e:
        raise RuntimeError(f"GitHub API error when creating release branch protection rules.") from e
    finally:
        gh.close()


def create_tag(tag_name: str, tag_message: str, branch: str, repo_path: str):
    """
    Create an annotated git tag on the specified branch.

    Args:
        tag_name: Name of the tag to create
        tag_message: Message for the annotated tag
        branch: Branch to tag
        repo_path: Path to the git repository
    """
    try:
        git_root = get_git_root(repo_path)
        git_repo = git.Repo(git_root)

        git_repo.git.checkout(branch)

        git_repo.git.tag('-a', tag_name, '-m', tag_message)
        logger.debug(f"Created tag {tag_name} on branch {branch}")

    except Exception as e:
        raise RuntimeError(f"Failed to create tag {tag_name}: {e}") from e


def push_tag(tag_name: str, repo_path: str, remote: str="origin"):
    """Push a git tag to remote repository."""
    try:
        git_root = get_git_root(repo_path)
        git_repo = git.Repo(git_root)

        git_repo.git.push(remote, tag_name)
        logger.debug(f"Pushed tag {tag_name} to {remote}")
    except Exception as e:
        raise RuntimeError(f"Failed to push tag {tag_name}: {e}") from e


def create_release_label(label_name: str, color: str="658045"):
    """
    Create a GitHub issue label in the origin upstream repository.

    Used to track issues and PRs destined for a given release's release notes.
    If a label with the same name already exists, this is a no-op.

    Args:
        label_name: Name of the label (e.g. "17-1 release note")
        color: 6-character hex color code (without leading '#')
    """
    token = os.getenv('GITHUB_TOKEN')
    if not token:
        raise RuntimeError(f"GITHUB_TOKEN environment variable not set.")

    try:
        gh = Github(auth=Auth.Token(os.environ['GITHUB_TOKEN']))
        repo = gh.get_repo(ORG_REPO)

        if label_name in [l.name for l in repo.get_labels()]:
            logger.debug(f"Label '{label_name}' already exists; skipping creation.")
            return

        label = repo.create_label(name=label_name, color=color)
        logger.debug(f"Created label '{label.name}'")

    except GithubException as e:
        raise RuntimeError(f"GitHub API error when creating label '{label_name}'.") from e
    finally:
        gh.close()
