import subprocess
import pathlib

def find_files( root_path    : pathlib.Path,
                search_token : str,
                include      : list[str]=[],
                exclude      : list[str]=[]) -> list[pathlib.Path]:
    """
    Searches for files within a directory tree that contain a specified search token using
    the Unix/MacOS command line tool `grep`.

    This function wraps the Unix `grep` command to recursively search through files
    starting from a root directory. It returns a list of `pathlib.Path` objects for
    files that contain the specified search token. The search can be further refined
    by specifying patterns for files to include or exclude.

    Parameters:
    - root_path (pathlib.Path): The root directory from which the search will begin.
                                Must be a valid directory path.
    - search_token (str): The token to search for within files. This is passed directly
                          to `grep`, so regular expressions can be used.
    - includes (list[str], optional): A list of patterns to include in the search.
                                      Patterns should match the file names to include.
                                      For example, ['*.py'] to include only Python files.
                                      Defaults to an empty list, which includes all files.
    - excludes (list[str], optional): A list of patterns to exclude from the search.
                                      Patterns should match the file names to exclude.
                                      For example, ['*.txt'] to exclude all text files.
                                      Defaults to an empty list, which excludes no files.

    Returns:
    - list[pathlib.Path]: A list of `pathlib.Path` objects, each representing a file
                          that contains the search token. The list will be empty if
                          no matching files are found.

    Raises:
    - Exception: If the `grep` command fails for any reason (e.g., due to an invalid
                 root_path or issues executing `grep`), an exception is raised with
                 the error message from `grep`.

    Example:
    >>> find_files(pathlib.Path('/path/to/search'), 'def main', includes=['*.py'])
    [PosixPath('/path/to/search/script1.py'), PosixPath('/path/to/search/dir/script2.py')]

    Note:
    - This function relies on the Unix `grep` command and may not be portable to
      environments without `grep` (e.g., some Windows environments without Unix-like
      tools installed).
    """

    # Ensure the root path is an existant directory
    assert( root_path.exists() )
    assert( root_path.is_dir() )

    if isinstance(include,str):
        include=[include]
    if isinstance(exclude,str):
        exclude=[exclude]


    cmd = ['grep','-rl',search_token] + \
          [f'--include={inc}' for inc in include] + \
          [f'--exclude={exc}' for exc in exclude] + \
          [str(root_path)]

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Check if the command was successful
    if result.returncode != 0:
        raise Exception(f"Error executing grep: {result.stderr}")

    # Parse the output into a list of Path objects
    return [pathlib.Path(line.strip()) for line in result.stdout.splitlines()]


