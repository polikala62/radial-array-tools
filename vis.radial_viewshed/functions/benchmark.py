'''
Created on Mar 12, 2024

@author: Karl
'''

import datetime

# Define timing function for the end of target function, and write output to CSV.
def benchmark(start_time, dictionary, function_name, ignore_sub_millisecond=True):
    
    # Get timedelta by subtracting current time from start time.
    timedelta = datetime.datetime.now() - start_time
    
    # Ignore value if less than a millisecond.
    if timedelta.total_seconds() > 0.0001 or ignore_sub_millisecond:
        
        # Add time elapsed to list in dictionary, create new list if key does not exist.
        if str(function_name) not in dictionary.keys(): # Create a list of elapsed times for each unique 'function name'.
            dictionary[function_name] = [timedelta]
        else:
            dictionary[function_name].append(timedelta)
    
    # Return updated dictionary.
    return dictionary

#------------------------------------------------------------------------------ 
# Define function for returning the finished benchmarking dictionary.
def print_benchmark_to_console(start_time, dictionary):
    
    sort_dict = {}
    total_time = 0
    
    # Calculate runtime.
    runtime = float((datetime.datetime.now() - start_time).total_seconds())
    
    # Calculate the total time for all functions.
    for func in dictionary.keys():
        
        # Calculate sum of time measured.
        func_time = float(sum([i.total_seconds() for i in dictionary[func]]))
        
        # Add sum to overall total.
        total_time += func_time
        
        # Add values to new dictionary, for sorting.
        sort_dict[func_time] = func
    
    # Define template for row.
    row_template = "{:<40} | {:<10} | {:<15} | {:<10} | {:<10} | {:<10} | {:<10} | {:<10}"
    
    # Create array with headers.
    print(row_template.format("FUNCTION_NAME", "%_RUNTIME", "%_MEASURED_TIME", "CALL_COUNT", "TOTAL_TIME", "MIN_TIME", "MAX_TIME", "AVG_TIME"))
    
    # Calculate number of times fired, average time, time as proportion of total.
    for sort_value in sorted(sort_dict.keys(), reverse=True):
        
        # Calculate row values.
        func_name = sort_dict[sort_value]
        percent_runtime = round((sort_value / runtime) * 100,2)
        percent_functiontime = round((sort_value / total_time) * 100,2)
        call_count = len(dictionary[func_name])
        total_time_used = round(sort_value,4)
        min_time = round(float(min([i.total_seconds() for i in dictionary[func_name]])),4)
        max_time = round(float(max([i.total_seconds() for i in dictionary[func_name]])),4)
        time_average = round(sort_value / call_count,4)
        
        # Print row values to console.
        print(row_template.format(func_name, percent_runtime, percent_functiontime, call_count, total_time_used, min_time, max_time, time_average))