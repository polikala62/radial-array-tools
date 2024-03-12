'''
Created on Mar 12, 2024

@author: Karl
'''

import time

# Define timing function for the end of target function, and write output to CSV.
def benchmark(start_time, dictionary, function_name):
    
    time_diff = round((time.time() - start_time), 4)
    if time_diff > 0.00001: # Ignore value if less than 1/10,000 of a second.
        
        # Add time elapsed to list in dictionary, create new list if key does not exist.
        if str(function_name) not in dictionary.keys(): # Create a list of elapsed times for each unique 'function name'.
            dictionary[function_name] = [time_diff]
        else:
            dictionary[function_name].append(time_diff)
            
    return dictionary
                
#------------------------------------------------------------------------------ 
# Define function for returning the finished benchmarking dictionary.
def print_benchmark_to_console(runtime, dictionary):
    
    sort_dict = {}
    total_time = 0
    # Calculate the total time for all functions.
    for key in dictionary.keys():
        total_time += sum(dictionary[key])
        # Add values to new dictionary, for sorting.
        sort_dict[sum(dictionary[key])] = key
    
    # Define template for row.
    row_template = "{:<50} ; {:<10} ; {:<10} ; {:<10} ; {:<10} ; {:<10}"
    
    # Create array with headers.
    print(row_template.format("FUNC_NAME", "%_RUNTIME", "%_FUNCTIME", "CALL_COUNT", "TOTAL_TIME", "AVG_TIME"))
    
    # Calculate number of times fired, average time, time as proportion of total.
    for sort_value in sorted(sort_dict.keys(), reverse=True):
        
        # Add values to row. MUST BE CALLED IN SAME ORDER AS HEADER IN ROW LIST.
        key = sort_dict[sort_value]
        percent_runtime = round((float(sum(dictionary[key])) / float(runtime)) * 100,4)
        percent_functiontime = round((float(sum(dictionary[key])) / float(total_time)) * 100,2)
        call_count = len(dictionary[key])
        total_time_used = round(sum(dictionary[key]),2)
        time_average = round(float(sum(dictionary[key])) / float(call_count),2)
        
        print(row_template.format(key, percent_runtime, percent_functiontime, call_count, total_time_used, time_average))