#!/usr/bin/env python3
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
from .EnvvarHelper import EnvvarHelper



class JenkinsEnv(EnvvarHelper):
    """
    This is a helper class for Jenkins envvars that we might want to use
    in our scripts. This class contains properties that are generic to Jenkins.

    Note:
        Information on the standard set of Jenkins parameters can be found
        at the following website:
        - https://wiki.jenkins.io/display/JENKINS/Building+a+software+project
    """
    def __init__(self):
        pass


    @property
    def build_number(self):
        return self.get_envvar_str("BUILD_NUMBER")


    @property
    def build_id(self):
        return self.get_envvar_str("BUILD_ID")


    @property
    def build_url(self):
        return self.get_envvar_str("BUILD_URL")


    @property
    def node_name(self):
        return self.get_envvar_str("NODE_NAME")


    @property
    def job_name(self):
        return self.get_envvar_str("JOB_NAME")


    @property
    def job_base_name(self):
        return self.get_envvar_str("JOB_BASE_NAME")


    @property
    def build_tag(self):
        return self.get_envvar_str("BUILD_TAG")


    @property
    def jenkins_home(self):
        return self.get_envvar_str("JENKINS_HOME")


    @property
    def jenkins_url(self):
        return self.get_envvar_str("JENKINS_URL")


    @property
    def executor_number(self):
        return self.get_envvar_str("EXECUTOR_NUMBER")


    @property
    def node_labels(self):
        return self.get_envvar_str("NODE_LABELS")


    @property
    def java_home(self):
        return self.get_envvar_str("JAVA_HOME")


    @property
    def workspace(self):
        return self.get_envvar_str("WORKSPACE")


    @property
    def svn_revision(self):
        return self.get_envvar_str("SVN_REVISION")


    @property
    def cvs_branch(self):
        return self.get_envvar_str("CVS_BRANCH")


    @property
    def git_commit(self):
        return self.get_envvar_str("GIT_COMMIT")


    @property
    def git_url(self):
        return self.get_envvar_str("GIT_URL")


    @property
    def git_branch(self):
        return self.get_envvar_str("GIT_BRANCH")


    #
    # Other Jenkins Vars with unknown origin
    #
    @property
    def jenkins_job_weight(self):
        return self.get_envvar_str("JENKINS_JOB_WEIGHT", error_if_missing=False)


    @property
    def jenkins_test_weight(self):
        return self.get_envvar_str("JENKINS_TEST_WEIGHT", error_if_missing=False)




