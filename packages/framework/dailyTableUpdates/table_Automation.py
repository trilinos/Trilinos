#/usr/bin/env python3
"""
Trilinos Status Table Automation Script
      
      table_Automation.py <pr_status, pr_merged, failed_pr, 
                                waiting_pr, open_pr, mm_status, successful_mm[, jira_ticket]>
    
      This script prints out todays Trilinos status in MarkDown syntax so that it can be added to the Trilinos Wiki page.
    
** Dependencies needed to run script:
    
    pip3 install python-dateutil
    pip3 install pytablewriter
    
    
SYNOPSIS                                                                                                                                 
<python> table_Automation.py <[0-3] 1 2 3 4 [0-3] 5 6 7 8 9 [435]> [TrilFrame-435]

#STATUS: [0] == success, [1] == warning [3] == failure -> Pull Request(s) and Master Merge Status codes.

<User inputs twelve charter(s) from terminal, 2* is not used>
updated table string example: python table_Automation2.py 0 0 2* 22 21 0 5 3 4 57 1 000
"""
from time import sleep as pause
from datetime import timedelta

from datetime import date
import subprocess
import sys
import os

STOP = 0b11
ONE_OFF = -0b1

try:
  import pytablewriter as ptw
  from dateutil.relativedelta import relativedelta, FR
except ImportError:
  print("\n\nHere are a list of needed dependencies ... \n\n\n")
  pause(STOP)

  print("\t*** Package(s) that need to be installed ***\n \n\tpip3 install python-dateutil " 
        + "\n\tpip3 install pytablewriter \n\n\n\t*** Install and try running again ***\n\n\n")
  pause(STOP)  

def main():
    #STATUS: [0] == success, [1] == warning [3] == failure
    mm_status = [":white_check_mark:", ":warning:", ":x:"]
    pr_status = [":white_check_mark:", ":warning:", ":x:"]
    status_example = ":white_check_mark:"
    
    stat_container = sys.argv[1:]
   
    today = date.today()
    yesterday = today - timedelta(days = 1)
    try:
        last_Friday = today + relativedelta(weekday=FR(ONE_OFF))
    
    except:
        print("Veiw needed dependencies above\n\n")

    try:
        number_of_pr_merged = stat_container[1]
        number_of_failed_pr = stat_container[2]
        number_wip_prs = stat_container[3]
        number_reviewed_required = stat_container[4]
        number_change_requested = stat_container[6]
        number_review_approved = stat_container[7]
        number_of_waiting_pr = stat_container[8]
        number_open_pr = stat_container[9]
        number_of_successful_mm = stat_container[10]
        jira_ticket_number = stat_container[11]
    except:
        print("Requires more arguments ... Example: table_Automation.py 0 0 2* 22 21 0 5 3 4 57 1 000")
    try:
        NUMBER_OF_PRs_MERGED = "["+number_of_pr_merged+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=is%3Apr+merged%3A"+str(yesterday)+"T12%3A00%3A00-07%3A00.."+str(today)+"T12%3A00%3A00-07%3A00+base%3Adevelop)"
        NUMBER_OF_PRs_WAITING = "["+number_of_waiting_pr+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=is%3Apr+is%3Aopen+base%3Adevelop+review%3Aapproved+status%3Afailure+-label%3A%22AT%3A+WIP%22)"
        NUMBER_SUCCESSFUL_MM = "["+number_of_successful_mm+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=is%3Apr+merged%3A"+str(last_Friday)+"T12%3A00%3A00-07%3A00.."+str(today)+"T12%3A00%3A00-07%3A00+base%3Amaster+)"
        JIRA_TICKETS = "[TrilFrame-"+jira_ticket_number+"]"+"(https://sems-atlassian-son.sandia.gov/jira/browse/TRILFRAME-"+str(jira_ticket_number)+")"
        Open_PRs = "["+number_open_pr+"]"+"(https://github.com/trilinos/Trilinos/pulls)"
        WIP_PRs = "["+number_wip_prs+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=+is%3Apr+is%3Aopen+base%3Adevelop+label%3A%22AT%3A+WIP%22)"
        REVIEW_REQUIRED = "["+number_reviewed_required+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=is%3Apr+is%3Aopen+base%3Adevelop+review%3Arequired+-label%3A%22AT%3A+WIP%22)"
        CHANGE_REQUESTED = "["+number_change_requested+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=+is%3Apr+is%3Aopen+base%3Adevelop+review%3Achanges-requested+-label%3A%22AT%3A+WIP%22)"
        REVIEW_APROVED = "["+number_review_approved+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=+is%3Apr+is%3Aopen+base%3Adevelop+review%3Aapproved+-status%3Afailure+-label%3A%22AT%3A+WIP%22)"
    except:
        print("Missing required argument(s)... Example: table_Automation.py 0 0 2* 22 21 0 5 3 4 57 1 000")
    try:
        writer = ptw.MarkdownTableWriter(
              table_name="Trilinos Status Table",
              headers=["Date", "PR Status", "PRs Merged (Past 24 Hrs from 12pm)", "WIP PRs (@ 12pm)", "Review-Required PRs (@ 12pm) ", "Change-Requested PRs (@ 12pm) ", "Review-Approved PRs (@ 12pm)", " Failed PRs (@ 12pm)", "Total Open PRs (@ 12pm)", "MM Status", "Master Merges (Past 24 hrs from 12pm)", "Jira Ticket #"],
              value_matrix=[
                    [str(today), pr_status[int(stat_container[0])], NUMBER_OF_PRs_MERGED, WIP_PRs, REVIEW_REQUIRED, CHANGE_REQUESTED, REVIEW_APROVED, NUMBER_OF_PRs_WAITING, Open_PRs, mm_status[int(stat_container[5])], NUMBER_SUCCESSFUL_MM, JIRA_TICKETS]
              ],
          )
    except:
        print("Not enough arguments ... Example: table_Automation.py 0 0 2* 22 21 0 5 3 4 57 1 000")
    try:
        print(writer)
    except:
        print("Can not write to file, missing arguements ... Example: table_Automation.py 0 0 2* 22 21 0 5 3 4 57 1 000")
            
if __name__ == "__main__":
    main()