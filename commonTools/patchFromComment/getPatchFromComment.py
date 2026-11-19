from argparse import ArgumentParser
from github import Github
from gitlab import Gitlab
import re
import os
from pathlib import Path
from subprocess import Popen
import tempfile
import questionary


def getCommentGithub(repo, issue_number, comment_id):
    g = Github()
    r = g.get_repo(repo)
    issue = r.get_issue(issue_number)
    comment = issue.get_comment(comment_id)
    return comment


def getCommentGitlab(instance, repo, work_item_type, work_item_number, comment_id):
    token = None
    if "GITLAB_TOKEN" in os.environ:
        token = os.environ["GITLAB_TOKEN"]
    g = Gitlab(url=instance, private_token=token)
    project = g.projects.get(repo)
    if work_item_type == "merge_requests":
        workItem = project.mergerequests.get(work_item_number)
    elif work_item_type == "work_items":
        workItem = project.issues.get(work_item_number)
    note = workItem.notes.get(comment_id)
    return note


def extractPatch(comment):
    patch = ""
    addToPatch = False
    for line in comment.body.split("\n"):
        if not addToPatch and line[:7] == "```diff":
            addToPatch = True
        elif addToPatch and line[:3] == "```":
            addToPatch = False
            patch += "\n"
        else:
            if addToPatch:
                if len(patch) > 0:
                    patch += "\n"+line
                else:
                    patch = line
    assert not addToPatch
    return patch


def main():
    parser = ArgumentParser(description="This script retrieves and applies or saves patches from Github or Gitlab comments.")
    parser.add_argument("op", type=str, help="What to do with the patch.", choices=["save", "apply"])
    parser.add_argument("url", type=str, help="Comment url")
    args = parser.parse_args()

    github_re = re.compile(r"([a-zA-Z]+/[a-zA-Z]+)/(pull|issue)/([0-9]+)#issuecomment-([0-9]+)")
    match_github = github_re.search(args.url)
    gitlab_re = re.compile(r"(https://[a-zA-Z-\.]+)/([a-zA-Z]+/[a-zA-Z]+)/-/(work_items|merge_requests)/([0-9]+)#note_([0-9]+)")
    match_gitlab = gitlab_re.search(args.url)

    if match_github:
        repo = match_github.group(1)
        work_item_type = match_github.group(2)
        work_item_number = int(match_github.group(3))
        comment_id = int(match_github.group(4))
    elif match_gitlab:
        instance = match_gitlab.group(1)
        repo = match_gitlab.group(2)
        work_item_type = match_gitlab.group(3)
        work_item_number = int(match_gitlab.group(4))
        comment_id = int(match_gitlab.group(5))
    else:
        questionary.print("Could not parse the url.\n"
                          "Expected format of Github URLs:\n"
                          " trilinos/Trilinos/pull/XXXXX#issuecomment-XXXXXXXXXX\n"
                          " trilinos/Trilinos/issue/XXXXX#issuecomment-XXXXXXXXXX\n"
                          "Expected format of Gitlab URLs:\n"
                          " work_items/XXXXX#note_XXXXXXXXXX\n"
                          " merge_requests/XXXXX#note_XXXXXXXXXX\n"
                          f"Got:      {args.url}")
        exit(1)

    try:
        if match_github:
            comment = getCommentGithub(repo, work_item_number, comment_id)
        elif match_gitlab:
            comment = getCommentGitlab(instance, repo, work_item_type, work_item_number, comment_id)
    except Exception as e:
        questionary.print(f"Failed to retrieve comment. Error:\n {e}")
        exit(1)
    try:
        patch = extractPatch(comment)
    except:
        questionary.print(f"Could not extract patch from comment:\n\n{comment.body}")
        exit(1)

    if args.op == "save":
        patchFile = Path().cwd()/f"{work_item_type}{work_item_number}_comment{comment_id}.patch"
        if patchFile.exists():
            overwrite = questionary.confirm(f"File {patchFile} already exists. Overwrite?", default=True).ask()
            if not overwrite:
                questionary.print("Not overwriting file. Aborting.")
                exit(1)
        patchFile.write_text(patch)
        questionary.print(f"Wrote patch to {patchFile}")
    elif args.op == "apply":
        # relies on this file being in commonTools/patchFromComment/
        rootDir = Path(os.path.abspath(__file__)).parent.parent.parent
        assert (rootDir/"commonTools"/"patchFromComment"/"getPatchFromComment.py").exists()
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write(patch)
            patchFile = Path(f.name)
        proc = Popen(f"git apply {patchFile}", cwd=rootDir, shell=True)
        proc.wait()
        if proc.returncode == 0:
            patchFile.unlink()
            questionary.print(f"Applied patch to source tree in {rootDir}")
        else:
            questionary.print(f"Patch file {patchFile} did not apply to source tree in {rootDir}")
            exit(1)
    else:
        raise NotImplementedError(args.op)


if __name__ == "__main__":
    main()
