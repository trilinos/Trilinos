#!/usr/bin/env python3
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
from .JenkinsEnv import JenkinsEnv



class TrilinosJenkinsEnv(JenkinsEnv):
    """
    This class extends the JenkinsEnv class to add Trilinos PR specific envvars.
    """
    def __init__(self):
        pass


    @property
    def trilinos_source_repo(self):
        return self.get_envvar_str("TRILINOS_SOURCE_REPO", error_if_missing=True)


    @property
    def trilinos_source_sha(self):
        return self.get_envvar_str("TRILINOS_SOURCE_SHA", error_if_missing=True)


    @property
    def trilinos_target_branch(self):
        return self.get_envvar_str("TRILINOS_TARGET_BRANCH", error_if_missing=True)


    @property
    def trilinos_target_repo(self):
        return self.get_envvar_str("TRILINOS_TARGET_REPO", error_if_missing=True)


    @property
    def trilinos_target_sha(self):
        return self.get_envvar_str("TRILINOS_TARGET_SHA", error_if_missing=True)


    @property
    def pullrequestnum(self):
        return self.get_envvar_str("PULLREQUESTNUM", error_if_missing=True)


    @property
    def pullrequest_cdash_track(self):
        return self.get_envvar_str("PULLREQUEST_CDASH_TRACK", error_if_missing=True)


    @property
    def force_clean(self):
        return self.get_envvar_str("FORCE_CLEAN", error_if_missing=True)
