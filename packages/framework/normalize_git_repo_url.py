#!/usr/bin/env python3

import sys

def removeStartingSubstr(inputStr, substrToRemove):
  if inputStr.startswith(substrToRemove):
    substrToRemoveLen = len(substrToRemove)
    return inputStr[substrToRemoveLen:]
  return inputStr

def removeTrailingSubstr(inputStr, substrToRemove):
  if inputStr.endswith(substrToRemove):
    substrToRemoveLen = len(substrToRemove)
    return inputStr[:-substrToRemoveLen]
  return inputStr

def normalizeGitRepoUrl(inputUrl):
  url = inputUrl
  url = removeStartingSubstr(url, "https://")
  url = removeStartingSubstr(url, "git@")
  url = removeTrailingSubstr(url, ".git")
  url = url.replace(":", "/")
  url = url.lower()
  return url

# Main

def main():
  if len(sys.argv) != 2:
    print("Usage: normalize_git_repo_url.py <string>")
    sys.exit(1)
  inputUrl = sys.argv[1]
  print(normalizeGitRepoUrl(inputUrl))

if __name__ == "__main__":
    main()
