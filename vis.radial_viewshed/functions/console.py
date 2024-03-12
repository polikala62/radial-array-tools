'''
Created on Dec 13, 2023

@author: Karl
'''

import datetime

#------------------------------------------------------------------------------ 

# Prints message to console.
def console(message_text, leading_spaces=0):
    
    time_string = datetime.datetime.now().strftime("%H:%M:%S")
    
    leading_space_string = "".join([" " for i in range(0,leading_spaces)]) #@UnusedVariable
    
    print("{} | {}{}".format(time_string, leading_space_string, message_text))

#------------------------------------------------------------------------------ 

def console_breakline(break_char, div_char, line_length):
    
    print("{}{}{}".format("".join([break_char for i in range(0,9)]), #@UnusedVariable
                          div_char,
                         "".join([break_char for i in range(0,line_length)]))) #@UnusedVariable

#------------------------------------------------------------------------------ 

# Prints 'percent complete' messages to console at specified intervals.
# Intervals are in percentage points: '5' will produce messages at each 5% interval.
def prcnt_complete(pr_count, pr_total, prcnt_inc, start_time, leading_spaces=0, leading_text=""):
    
    # Create list of breakpoints.
    breakpoints = [(pr_total / (100/prcnt_inc)) * i  for i in range(1, int((100/prcnt_inc))+1)]
    
    # Create dictionary to hold breakpoints.
    breakpoint_dict = {}
    
    # Round values in breakpoint list to integers, so that they correspond to pr_count values.
    for i, j in enumerate(breakpoints):
    
        breakpoint_dict[int(j) + bool(j%1)] = ((i+1)*prcnt_inc)
    
    # Estimate time remaining.
    avg_time = ((datetime.datetime.now() - start_time).total_seconds())/pr_count
    est_secs = (pr_total - pr_count) * avg_time
    est_hours = str(int(est_secs // 3600)).zfill(2)
    est_minutes = str(int((est_secs % 3600)) // 60).zfill(2)
    est_seconds = str(int((est_secs % 3600)) % 60).zfill(2)
    est_string = "{}:{}:{}".format(est_hours, est_minutes, est_seconds)
    
    # If pr_count is in the breakpoint dict, print message to console.
    if pr_count in breakpoint_dict.keys():
        
        console("{}{}% complete, approx. {} left...".format(leading_text, breakpoint_dict[pr_count], est_string), leading_spaces)