from argparse import ArgumentParser
from github import Github
import re
import os
from pathlib import Path
from subprocess import Popen
import tempfile
import questionary

repo = "trilinos/Trilinos"


def getComment(pr_number, comment_id):
    g = Github()
    r = g.get_repo(repo)
    pr = r.get_pull(pr_number)
    comment = pr.as_issue().get_comment(comment_id)
    return comment


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
    parser = ArgumentParser(description="This script retrieves and applies or saves patches from Github comments.")
    parser.add_argument("op", type=str, help="What to do with the patch.", choices=["save", "apply"])
    parser.add_argument("url", type=str, help="Comment url")
    args = parser.parse_args()

    try:
        m = re.compile(repo+r"/pull/([0-9]+)#issuecomment-([0-9]+)").search(args.url)
        pr_number = int(m.group(1))
        comment_id = int(m.group(2))
    except:
        questionary.print("Could not parse the url.\n"
                          "Expected: trilinos/Trilinos/pull/XXXXX#issuecomment-XXXXXXXXXX\n"
                          f"Got:      {args.url}")
        exit(1)

    try:
        comment = getComment(pr_number, comment_id)
    except Exception as e:
        questionary.print(f"Failed to retrieve comment. Error:\n {e}")
        exit(1)
    try:
        patch = extractPatch(comment)
    except:
        questionary.print(f"Could not extract patch from comment:\n\n{comment.body}")
        exit(1)

    if args.op == "save":
        patchFile = Path().cwd()/f"pr{pr_number}_comment{comment_id}.patch"
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
