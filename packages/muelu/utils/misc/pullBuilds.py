#!/usr/bin/env python3

# Pulls output for failing autotester runs from CDash and dumps them to file.
# Usage: ./pullBuilds.py THE_PR_NUMBER

import requests
from sys import argv

# The CDash instance
base = 'https://trilinos-cdash.sandia.gov/api/v1/'
testbase = 'https://trilinos-cdash.sandia.gov/'

# get the number of the PR from the commandline
prNum = argv[1]

# all builds for a particular PR are named like this
pr = 'PR-{}-test'.format(prNum)

outputFilename = 'PR_{}'.format(prNum)
print('writing to {}\n'.format(outputFilename))

with open(outputFilename, 'w') as f:
    # find all failing tests for a given PR
    responsePR = requests.get(base+'index.php?project=Trilinos&filtercombine=and&filtercombine=and&filtercombine=and&filtercombine=and&filtercombine=and&filtercombine=and&filtercount=2&showfilters=1&filtercombine=and&field1=buildname&compare1=63&value1={}&field2=buildstarttime&compare2=84&value2=NOW'.format(pr))
    assert responsePR.ok

    # get the names of all the builds, and reduce to the latest iteration of each build type
    latest_builds = {}
    for build in responsePR.json()['buildgroups'][0]['builds']:
        buildID = build['id']

        # Go through all the builds and build a list of the most recent iterations
        buildurl = base+'viewTest.php?onlyfailed&buildid={}'.format(buildID)
        responseBuild = requests.get(buildurl)
        assert responseBuild.ok

        build_name = responseBuild.json()['build']['name']
        build_name_base, build_iter = build_name.rsplit('-', 1)
        build_iter = int(build_iter)

        print("Build found:", build_name)

        if build_name_base in latest_builds:
            if build_iter > latest_builds[build_name_base]:
                latest_builds[build_name_base] = build_iter
        else:
            latest_builds[build_name_base] = build_iter

    latest_builds = set([build_name_base+'-'+str(latest_builds[build_name_base]) for build_name_base in latest_builds])
    print("Selected builds:", latest_builds)

    for build in responsePR.json()['buildgroups'][0]['builds']:
        buildID = build['id']

        # Go through all the builds and check the number of failed tests
        buildurl = base+'viewTest.php?onlyfailed&buildid={}'.format(buildID)
        #print('Grabbing test URL = '+testurl)
        responseBuild = requests.get(buildurl)
        assert responseBuild.ok
        build_name = responseBuild.json()['build']['name']

        if build_name not in latest_builds:
            continue

        num_failed = responseBuild.json()['numFailed']
        print('Build '+build_name+' has '+str(num_failed)+' failing tests')
        # print('- Test URL = '+buildurl)

        # Get details for each test
        for i in range(num_failed):
            test_id = responseBuild.json()['tests'][i]['buildtestid']

            testurl = base+'testDetails.php?buildtestid={}'.format(test_id)
            # print('- ' + testurl)

            responseTest = requests.get(testurl)
            assert responseTest.ok

            build_name = responseTest.json()['test']['build']
            test_name = responseTest.json()['test']['test']
            test_output = responseTest.json()['test']['output']

            print('\n' + '='*120, file=f)
            print('Build: {}'.format(build_name), file=f)
            print('Test: {}\n'.format(test_name), file=f)
            print(test_output, file=f)

            # print('- '+test_output)


print('done')
