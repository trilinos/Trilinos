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
<python> table_Automation.py <0 4 6 4 72 0 4 [435]>

#STATUS: [0] == success, [1] == warning [3] == failure -> Pull Request(s) and Master Merge Status codes.

"""
from time import sleep as pause
from datetime import timedelta
from datetime import date
import sys

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
        print("View needed dependencies above\n\n")
        
    try:
        number_of_pr_merged = stat_container[1]
        number_of_failed_pr = stat_container[2]
        number_of_waiting_pr = stat_container[3]
        number_open_pr = stat_container[4]
        number_of_successful_mm = stat_container[6]
        jira_ticket_number = stat_container[7]
    except:
        print("Requires more arguments ... Example: table_Automation.py 0 4 6 4 72 0 4 435")
    try:
        NUMBER_OF_PRs_MERGED = "["+number_of_pr_merged+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=is%3Apr+merged%3A"+str(yesterday)+"T12%3A00%3A00-07%3A00.."+str(today)+"T12%3A00%3A00-07%3A00+base%3Adevelop)"
        NUMBER_OF_FAILED_PRs = "["+number_of_failed_pr+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=is%3Apr+updated%3A"+str(yesterday)+"T12%3A00%3A00-07%3A00.."+str(today)+"T12%3A00%3A00-07%3A00+base%3Adevelop+status%3Afailure+)"
        NUMBER_OF_PRs_WAITING = "["+number_of_waiting_pr+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=is%3Apr+is%3Aopen+updated%3A"+str(yesterday)+"T12%3A00%3A00-07%3A00.."+str(today)+"T12%3A00%3A00-07%3A00+base%3Adevelop+status%3Apending+)"
        NUMBER_SUCCESSFUL_MM = "["+number_of_successful_mm+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=is%3Apr+merged%3A"+str(last_Friday)+"T12%3A00%3A00-07%3A00.."+str(today)+"T12%3A00%3A00-07%3A00+base%3Amaster+)"
        JIRA_TICKETS = "[TrilFrame-"+jira_ticket_number+"]"+"(https://sems-atlassian-son.sandia.gov/jira/browse/TRILFRAME-"+str(jira_ticket_number)+")"
    except:
        print("Missing required argument(s)... Example: table_Automation.py 0 4 6 4 72 0 4 435")
    try:
        writer = ptw.MarkdownTableWriter(
              table_name="Trilinos Status Table",
              headers=["Date", "PR Status", "Number of PRs Merged", "Number of Failed PRs", "Number of PRs Waiting for Build & Test", "Total Number of Open PRs", "MM Status", "Number of Successful Master Merges","Jira Ticket #"],
              value_matrix=[
                    [str(today), pr_status[int(stat_container[0])], NUMBER_OF_PRs_MERGED, NUMBER_OF_FAILED_PRs, NUMBER_OF_PRs_WAITING, number_open_pr , mm_status[int(stat_container[5])], NUMBER_SUCCESSFUL_MM, JIRA_TICKETS]
              ],
          )
    except:
        print("Not enough arguments ... Example: table_Automation.py 0 4 6 4 72 0 4 435")
    try:
        print(writer)
    except:
        print("Can not write to file, missing arguements ... Example: table_Automation.py 0 4 6 4 72 0 4 435")
    
if __name__ == "__main__":
    main()