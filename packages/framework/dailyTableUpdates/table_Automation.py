#/usr/bin/env python3
"""
Trilinos Status Table Automation Script
      
      <python> table_Automation.py --prstatus <int>[0-2] --merged <int> --wip <int> --reviewed <int> 
                                        --changerequest <int> --reviewaproved <int> --failedprs <int> --totalopen <int> --mmstatus <int>[0-2] --mmmerges <int> --jira <int>

      python table_Automation.py --prstatus [0-2] --merged 2 --wip 3 --reviewed 4 --changerequest 5 --reviewaproved 6 --failedprs 7 --totalopen 8 --mmstatus [0-2] --mmmerges 10 --jira 11

      This script prints out todays Trilinos status in MarkDown syntax so that it can be added to the Trilinos Wiki page.
    
** Dependencies needed to run script:
    
    pip3 install python-dateutil
    pip3 install pytablewriter
    
    
SYNOPSIS                                                                                                                                 

<python> table_Automation.py --prstatus <int>[0-2] --merged <int> --wip <int> --reviewed <int> --changerequest <int> --reviewaproved <int> --failedprs <int> --totalopen <int> --mmstatus <int>[0-2] --mmmerges <int> --jira <int>

#STATUS --prstatus --mmstatus : [0] == success, [1] == warning [2] == failure -> Pull Request(s) and Master Merge Status codes.

<User inputs eleven int value(s) from terminal>

python table_Automation.py --prstatus [0-2] --merged 2 --wip 3 --reviewed 4 --changerequest 5 --reviewaproved 6 --failedprs 7 --totalopen 8 --mmstatus [0-2] --mmmerges 10 --jira 11
"""
from time import sleep as pause
from datetime import timedelta

from datetime import date
import argparse
import sys

STOP = 0b11
ONE = 0b1
ONE_OFF = -ONE
DUNDER_MAIN = '__main__'

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
    #STATUS: [0] == success, [1] == warning [2] == failure

    pr_status = [":white_check_mark:", ":warning:", ":x:"] 

    mm_status = [":white_check_mark:", ":warning:", ":x:"]
    
    today = date.today()

    yesterday = today - timedelta(days = ONE)

    try:
        last_Friday = today + relativedelta(weekday=FR(ONE_OFF))
    
    except:
        print("Veiw needed dependencies above\n\n")

    try:
        parser = argparse.ArgumentParser()
        
        parser.add_argument("--prstatus", type=int, required=True, help="Pull Request Status")
        
        parser.add_argument("--merged", type=int, required=True, help="Pull Request Merged Past 24 Hrs")
        
        parser.add_argument("--wip", type=int, required=True, help="WIP Pull Request @ 12pm")
        
        parser.add_argument("--reviewed", type=int, required=True, help="Review Required, Pull Request @ 12pm")
        
        parser.add_argument("--changerequest", type=int, required=True, help="Change Requested Pull Request @ 12pm")
        
        parser.add_argument("--reviewaproved", type=int, required=True, help="Review Approved, Pull Request @ 12pm")
        
        parser.add_argument("--failedprs", type=int, required=True, help="Failed Pull Request @ 12pm")

        parser.add_argument("--totalopen", type=int, required=True, help="Total Open Pull Request @ 12pm")
        
        parser.add_argument("--mmstatus", type=int, required=True, help="Master Merge Status")
        
        parser.add_argument("--mmmerges", type=int, required=True, help="Master Merge(s) last 24Hr's")
        
        parser.add_argument("--jira", type=int, required=True, help="Jira tickets, if needed")
        
        args = parser.parse_args()
       
        pull_request_status = args.prstatus

        number_of_pr_merged = args.merged 
        
        number_wip_prs = args.wip
        
        number_reviewed_required = args.reviewed 
        
        number_change_requested = args.changerequest 
        
        number_review_approved = args.reviewaproved 
        
        number_of_waiting_pr = args.failedprs 
        
        number_of_successful_mm = args.mmmerges 
        
        master_merge_status = args.mmstatus 
        
        number_open_pr = args.totalopen 
        
        jira_ticket_number = args.jira 

    except:
        print("Error in variable initialization ....")
        
        print("\n ") 
        
        print("table_Automation.py --prstatus <int>[0-2] --merged <int> --wip <int> --reviewed <int> --changerequest <int> --reviewaproved <int> --failedprs <int> --totalopen <int> --mmstatus <int>[0-2]")
        
        print("    ")
        
        print(sys.stderr)

    try:
        NUMBER_OF_PRs_MERGED = "["+str(number_of_pr_merged)+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=is%3Apr+merged%3A"+str(yesterday)+"T12%3A00%3A00-07%3A00.."+str(today)+"T12%3A00%3A00-07%3A00+base%3Adevelop)"
        
        WIP_PRs = "["+str(number_wip_prs)+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=+is%3Apr+is%3Aopen+base%3Adevelop+label%3A%22AT%3A+WIP%22)"
        
        REVIEW_REQUIRED = "["+str(number_reviewed_required)+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=is%3Apr+is%3Aopen+base%3Adevelop+review%3Arequired+-label%3A%22AT%3A+WIP%22)"
        
        CHANGE_REQUESTED = "["+str(number_change_requested)+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=+is%3Apr+is%3Aopen+base%3Adevelop+review%3Achanges-requested+-label%3A%22AT%3A+WIP%22)"
        
        REVIEW_APROVED = "["+str(number_review_approved)+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=+is%3Apr+is%3Aopen+base%3Adevelop+review%3Aapproved+-status%3Afailure+-label%3A%22AT%3A+WIP%22)"
        
        NUMBER_OF_PRs_WAITING = "["+str(number_of_waiting_pr)+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=is%3Apr+is%3Aopen+base%3Adevelop+review%3Aapproved+status%3Afailure+-label%3A%22AT%3A+WIP%22)"
        
        Open_PRs = "["+str(number_open_pr)+"]"+"(https://github.com/trilinos/Trilinos/pulls)"
        
        NUMBER_SUCCESSFUL_MM = "["+str(number_of_successful_mm)+"]"+"(https://github.com/trilinos/Trilinos/pulls?q=is%3Apr+merged%3A"+str(last_Friday)+"T12%3A00%3A00-07%3A00.."+str(today)+"T12%3A00%3A00-07%3A00+base%3Amaster+)"
        
        JIRA_TICKETS = "[TrilFrame-"+str(jira_ticket_number)+"]"+"(https://sems-atlassian-son.sandia.gov/jira/browse/TRILFRAME-"+str(jira_ticket_number)+")"    
        
    except:
        print("Error in Link initialization ... ")

    try:
        print("\n  ")

        writer = ptw.MarkdownTableWriter(
              table_name="Trilinos Status Table",
              headers=["Date", "PR Status", "PRs Merged (Past 24 Hrs from 12pm)", "WIP PRs (@ 12pm)", "Review-Required PRs (@ 12pm) ", "Change-Requested PRs (@ 12pm) ", "Review-Approved PRs (@ 12pm)", " Failed PRs (@ 12pm)", "Total Open PRs (@ 12pm)", "MM Status", "Master Merges (Past 24 hrs from 12pm)", "Jira Ticket #"],
              value_matrix=[
                    [str(today), pr_status[pull_request_status], NUMBER_OF_PRs_MERGED, WIP_PRs, REVIEW_REQUIRED, CHANGE_REQUESTED, REVIEW_APROVED, NUMBER_OF_PRs_WAITING, Open_PRs, mm_status[master_merge_status], NUMBER_SUCCESSFUL_MM, JIRA_TICKETS]
              ],
          )

        print("\n  ") 

    except:
        print("Error in table creation")
        
        print("\n ") 
        
        print("table_Automation.py --prstatus <int>[0-2] --merged <int> --wip <int> --reviewed <int> --changerequest <int> --reviewaproved <int> --failedprs <int> --totalopen <int> --mmstatus <int>[0-2]")
        
        print("    ")

    try:
        print(writer)

        print("\n  ")

    except:
        print("Error in generating table to terminal ....")
        
        print("\n ") 
        
        print("table_Automation.py --prstatus <int>[0-2] --merged <int> --wip <int> --reviewed <int> --changerequest <int> --reviewaproved <int> --failedprs <int> --totalopen <int> --mmstatus <int>[0-2]")
        
        print("    ")
            
if __name__ == DUNDER_MAIN:
    main()