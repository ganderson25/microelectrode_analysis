#Program to split one CA file containing mulptiple conditions into a separate file with corresponding trials for each condition 

from pathlib import Path
import os
import pandas as pd
import re
import numpy as np

#Inputs
path = Path("C:/Users/Grace/Documents/Data/data_processing/convert/split/")
path_done = Path("C:/Users/Grace/Documents/Data/data_processing/convert/reformated")

condition_type = 'RH'

num_trials = 3
manually_assign_loops = False
manual_loop_assignment = {
    "60%" : 3,
    "70%" : 3,
    "80%" : 3
}


# Functions

def record_loop_start(text_lines, full_df):
    old_value = 0
    loop_counter = 0
    i = 0
    line_recorder = []
    ns_line_list = text_lines[0].split()
    ns_max = ns_line_list[-1]
    for item in full_df['Ns']:
        new_value = item
        if ((old_value == 0)& (new_value == 1)):
            loop_counter += 1
        if i == 0:
            line_recorder.append(i)
        elif ((old_value == ns_max)& (new_value == 0)):
            line_recorder.append(i)
        old_value = item
        i += 1
    return line_recorder

def extract_conditions(filename, condition_type):
    experimental_conditions = experimental_conditions_to_dictionary(filename)
    if condition_type in experimental_conditions.keys():
        conditions_combined = experimental_conditions[condition_type]
        conditions = conditions_combined.split(',')
    else:
        conditions = 'Error - condition type not found in conditions'
    return conditions


def experimental_conditions_to_dictionary(filename):
    experimental_inputs = extract_experimental_conditions(filename)
    condition_to_value = dict()
    for element in experimental_inputs:
        condition = element.split("=")[0]
        value = element.split("=")[1]
        condition_to_value[condition] = value
    return condition_to_value


def extract_experimental_conditions(filename):
    experimental_inputs = []
    filename_details = filename.split()
    for item in filename_details:
        if "=" in item:
            experimental_inputs.append(item)
    return experimental_inputs

def assign_num_loops(manually_assign_loops, manual_loop_assignment, conditions, num_trials):
    if manually_assign_loops:
        loop_assignments = manual_loop_assignment
    else:
        loop_assignments = automatic_num_loops(conditions, num_trials)
    return loop_assignments        


def automatic_num_loops(conditions, num_trials):
    i=0
    num_trials_list = []
    while i < len(conditions):
        num_trials_list.append(num_trials)
        i += 1
    auto_loop_assignments = dict(zip(conditions, num_trials_list))
    return auto_loop_assignments


def loop_df(conditions, loop_number_assignments):
    i=0
    loop_assignments = []
    while i < len(conditions):
        num_loops = loop_number_assignments[conditions[i]]
        num_loops_previous = loop_number_assignments[conditions[i-1]]
        if i == 0:
            loops = [0, num_loops - 1]
        else:
            loops = [num_loops_previous + 1, num_loops_previous +  num_loops]    
        loop_assignments.append(loops)
        i += 1
        loops = dict(zip(conditions, loop_assignments))
        loop_df = pd.DataFrame(data=loops)
    return loop_df
    

# Data extraction

dfs = []
fileList = os.listdir(path)
for filename in fileList:
    Ns = 0
    Ewe = 0
    filepath = os.path.join(path, filename)

    with open(filepath, "r") as f:
        loops = []
        text_lines = []
        while True:
            line = f.readline()

            if bool(re.match("Ns", line)):
                if not bool(re.search('mA', line)):
                    text_lines.append(line)
            if "Ei (V)" in line:
                text_lines.append(line)
            #if "Number of loops" in line:
            #    num_loops = line.split()[-1]
            #if "from point number" in line:
            #    loop_details = re.findall(r'\d+', line)
            #    loops.append(loop_details)

            # 'Ewe/V' key to start reading data
            if (bool(re.match("Ns", line)) & bool(re.search('mA', line))):
                text_lines.append(line)
                columns = line.split("\t")[:-1]
                full_df = pd.read_csv(f, delimiter="\t", names=columns)
                break

        #loop_details = pd.DataFrame(np.array(loops), columns=['Loop Number', 'Starting Index', 'Ending Index'])

        print(text_lines)

        #get loops - sometimes ec-lab doesn't write loop summaries correctly
        loop_starting_point = record_loop_start(text_lines, full_df)
        
        

        #Assign loops to conditions
        conditions = extract_conditions(filename, condition_type)
        if conditions == 'Error - condition type not found in conditions':
            print("Error conditions not found in filename for " + filename)
            break

        else:
            loop_number_assignments = assign_num_loops(manually_assign_loops, manual_loop_assignment, conditions, num_trials)
            loop_assignments = loop_df(conditions, loop_number_assignments)
            
            for condition in conditions:
                loops = loop_assignments[condition]



        print("yay")


