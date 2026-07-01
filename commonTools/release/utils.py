import re
import logging
logger = logging.getLogger(__name__)

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

