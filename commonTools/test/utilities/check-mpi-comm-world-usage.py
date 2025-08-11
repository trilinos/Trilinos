import sys
import subprocess
import re
import argparse


def parse_diff_output(changed_files):
    # Regex to capture filename and the line numbers of the changes
    file_pattern = re.compile(r"^\+\+\+ b/(.*?)$", re.MULTILINE)
    line_pattern = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@", re.MULTILINE)

    files = {}
    for match in file_pattern.finditer(changed_files):
        file_name = match.group(1)

        # Filtering for C/C++ files and excluding certain directories
        if file_name.endswith((".c", ".cpp", ".h", ".hpp")) and all(
            excluded not in file_name
            for excluded in [
                "doc/",
                "test_utils/",
                "test/",
                "tests/",
                "unit_test",
                "perf_test",
                "example/",
                "examples/",
            ]
        ):
            # Find the lines that changed for this file
            lines_start_at = match.end()
            next_file_match = file_pattern.search(changed_files, pos=match.span(0)[1])

            # Slice out the part of the diff that pertains to this file
            file_diff = changed_files[
                lines_start_at : next_file_match.span(0)[0] if next_file_match else None
            ]

            # Extract line numbers of the changes
            changed_lines = []
            for line_match in line_pattern.finditer(file_diff):
                start_line = int(line_match.group(1))
                num_lines = int(line_match.group(2) or 1)

                # The start and end positions for this chunk of diff
                chunk_start = line_match.end()
                next_chunk = line_pattern.search(file_diff, pos=line_match.span(0)[1])
                chunk_diff = file_diff[
                    chunk_start : next_chunk.span(0)[0] if next_chunk else None
                ]

                lines = chunk_diff.splitlines()
                line_counter = 0
                for line in lines:
                    if line.startswith("+"):
                        if (
                            "MPI_COMM_WORLD" in line
                            and not "CHECK: ALLOW MPI_COMM_WORLD" in line
                        ):
                            # Only include lines where "MPI_COMM_WORLD" is added
                            # and "CHECK: ALLOW MPI_COMM_WORLD" is not present
                            changed_lines.append(start_line + line_counter)

                        line_counter += 1

            if changed_lines:
                files[file_name] = changed_lines

    return files


def get_common_ancestor(target_branch, feature_branch):
    cmd = ["git", "merge-base", target_branch, feature_branch]
    return subprocess.check_output(cmd).decode("utf-8").strip()


def get_changed_files(target_branch, feature_branch):
    """Get a dictionary of files and their changed lines between the common ancestor and feature_branch."""
    start_commit = get_common_ancestor(target_branch, feature_branch)
    cmd = [
        "git",
        "diff",
        "-U0",
        "--ignore-all-space",
        start_commit,
        feature_branch
    ]
    out = subprocess.check_output(cmd)
    result = []
    for line in out.splitlines():
        try:
            result.append(line.decode("utf-8"))
        except UnicodeDecodeError:
            print(f"WARNING: Line {line} contains non-UTF8 characters; ignoring it")

    return parse_diff_output("\n".join(result))


def print_occurences(changed_files, title):
    print(title)
    for file_name, lines in changed_files.items():
        print("-----")
        print(f"File: {file_name}")
        print("Changed Lines:", ", ".join(map(str, lines)))
        print("-----")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--base", default="origin/develop", help="BASE commit (default: %(default)s)"
    )
    parser.add_argument(
        "--head", default="HEAD", help="HEAD commit (default: %(default)s)"
    )

    start_commit = parser.parse_args().base
    print(f"Start commit: {start_commit}")

    end_commit = parser.parse_args().head
    print(f"End commit: {end_commit}")

    mpi_comm_world_detected = get_changed_files(start_commit, end_commit)

    if mpi_comm_world_detected:
        print_occurences(
            mpi_comm_world_detected, "Detected MPI_COMM_WORLD in the following files:"
        )

        sys.exit(1)  # Exit with an error code to fail the GitHub Action
    else:
        print("No addition of MPI_COMM_WORLD detected.")
        sys.exit(0)
