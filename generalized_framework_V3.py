# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 14:17:42 2023

@author: Computer
"""

###Imports
import pandas as pd
import numpy as np
import os
import copy
from itertools import product

dirname = os.getcwd()
filename = os.path.join(dirname, 'Acid_Base_Data')

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, 'Acid_Base_Data')
print(filename)

##########  File loading  ##########

##Import acid data
Alanine_df = pd.read_excel(f"{filename}/Alanine.xlsx")
Arginine_df = pd.read_excel(f"{filename}/Arginine.xlsx")
Asparagine_df = pd.read_excel(f"{filename}/Asparagine.xlsx")
Aspartic_acid_df = pd.read_excel(f"{filename}/Aspartic Acid.xlsx")
Cysteine_df = pd.read_excel(f"{filename}/Cysteine.xlsx")
Glutamic_acid_df = pd.read_excel(f"{filename}/Glutamic Acid.xlsx")
Glutamine_df = pd.read_excel(f"{filename}/Glutamine.xlsx")
Glycine_df = pd.read_excel(f"{filename}/Glycine.xlsx")
Isoleucine_df = pd.read_excel(f"{filename}/Isoleucine.xlsx")
Leucine_df = pd.read_excel(f"{filename}/Leucine.xlsx")
Lysine_df = pd.read_excel(f"{filename}/Lysine.xlsx")
Methionine_df = pd.read_excel(f"{filename}/Methionine.xlsx")
Selenocysteine_df = pd.read_excel(f"{filename}/Selenocysteine.xlsx")
Serine_df = pd.read_excel(f"{filename}/Serine.xlsx")
Threonine_df = pd.read_excel(f"{filename}/Threonine.xlsx")
Valine_df = pd.read_excel(f"{filename}/Valine.xlsx")

##Acids that were excluded from the first wave of analysis (due to ring structures)
Histidine_df = pd.read_excel(f"{filename}/New/Histidine.xlsx")
Phenylalanine_df = pd.read_excel(f"{filename}/New/Phenylalanine.xlsx")
Proline_df = pd.read_excel(f"{filename}/New/Proline.xlsx")
Tryptophan_df = pd.read_excel(f"{filename}/New/Tryptophan.xlsx")
Tyrosine_df = pd.read_excel(f"{filename}/New/Tyrosine.xlsx")

##########  File cleaning  ##########
# Dictionary containing DFs before they are transformed to correct data representation
# Used to provide easy access to dfs on command.
acid_processing_df_dict = {}
acid_processing_df_dict = {"Alanine" : Alanine_df, "Arginine" : Arginine_df,
                "Asparagine" : Asparagine_df, "Aspartic Acid" : Aspartic_acid_df,
                "Cysteine" : Cysteine_df, "Glutamic Acid" : Glutamic_acid_df,
                "Glutamine" : Glutamine_df, "Glycine" : Glycine_df,
                "Histidine" : Histidine_df, 
                "Isoleucine" : Isoleucine_df,
                "Leucine" : Leucine_df, "Lysine" : Lysine_df,
                "Methionine" : Methionine_df, "Phenylalanine" : Phenylalanine_df,
                "Proline" : Proline_df,
                #"Selenocysteine": Selenocysteine_df,
                "Serine" : Serine_df, "Threonine" : Threonine_df,
                "Tryptophan" : Tryptophan_df, "Tyrosine" : Tyrosine_df,
                "Valine" : Valine_df}


#This is what we reference while our code runs
acid_df_dict = {}

##Original data not in a format that's easy to manipulate during normalization
##Turns DF into new format (with one entry per row) and stores transformed result in acid_df_dict
for df in acid_processing_df_dict:
    #Shows us where we are in the process
    print(df)
    #Starts from first node and breaks data into rows so that each edge (graph formalism) and accompanying nodes are a single data point
    conversion_dict = {}
    node_value = 1
    for i in range(len(acid_processing_df_dict[df])):
        split_cell = acid_processing_df_dict[df]['Edges'][i].split()
        cleaned_cell = []
        for j in split_cell:
            cleaned_cell.append(j.strip(','))
        conversion_dict[str(node_value)] = cleaned_cell
        node_value += 1
    ##At this point conversion_dict is full of the cleaned data in dict format
    new_df = pd.DataFrame({'Node_1': [], 'Node_2': [], 'Node_1_type': [], 'Node_2_type': []})
    ##This adds rows to the end of the dataframe
    counter = 1
    ##There should be a better way to do this with maps (or something)...
    for entry in conversion_dict:
        list_of_values = conversion_dict[str(entry)]
        for node2 in list_of_values:
            ##Set like this to easily convert last two columns to node type
            new_df.loc[len(new_df.index)] = [entry, node2, entry, node2]
    #We now have paired values for the acid in one value per row.
    #And now need to build node types (atomic element) for Node_1 and Node_2
    #Make a dict using node as key and node type as value from df columns, then create new columns from that
    type_dict ={}
    for i in range(len(acid_processing_df_dict[df])):
        type_dict[f'{i + 1}'] = acid_processing_df_dict[df]['Type'].loc[acid_processing_df_dict[df].index[i]]
    new_df = new_df.replace({'Node_1_type': type_dict, 'Node_2_type': type_dict})
    #new_df = pd.DataFrame.from_dict(conversion_dict)
    ##Now change node and edge columns back to integers so we can sort later
    new_df['Node_1'] = new_df['Node_1'].astype(int)
    new_df['Node_2'] = new_df['Node_2'].astype(int)
    ##Add columns for normalization
    new_df['Normalized_node_1'] = new_df['Node_1']
    new_df['Normalized_node_2'] = new_df['Node_2']
    #Now store the new df
    acid_df_dict[df] = copy.deepcopy(new_df)

##########  Solution set derivation  ##########

###Can be adjusted based on model parameters (w/ or w/out Glycine and Proline -- currently w/out).
###Taken from "solution_set_derivation.py"
#lower and upper bound, respectively. Upper bound chosen based on the number of atoms in the amino
#acids that come from codons at the physiological pH. (Crucial assumption, see discussion in doc about environmental transformations) 
l = 0
m = 27
counter = 1

##soluton_set_derivation
##This section creates the node groupings that are tested within an Interpretation Framework
## For a five node graph we could see: [[1, 2], [3], [4, 5]] or [[1, 2, 3, 4, 5], [], []]
solution_set_derivation = {}
acid_length_dict = {}
for i in range(l, m + 1):
    #Setting up the stuff we change each loop
    hold_list = []
    for j in range(i):
        hold_list.append(j + 1)
        
    #Where we are going to put the numbers to get our answers
    bucket_1 = []

    bucket_2 = []
    
    bucket_3 = []

    #Set up ticker to track progress, need to use negative sign in front of ticker to preserve order
    ticker = 0
    solution_set_list = []
    while ticker <= i:
        print()
        if ticker == 0:
            bucket_1 = hold_list
            #store results as a solution set
            solution_set = [bucket_1, bucket_2, bucket_3]
            print(solution_set)
            solution_set_list.append(copy.deepcopy(solution_set))
            ticker += 1
            continue
        
        else:
            #for loop?
            #Get the numbers for bucket 2 and for bucket 1
            bucket_2 = hold_list[-ticker:]
            #Don't modify hold_list directly, just place it into the correct bucket
            bucket_1 = hold_list[:-ticker]
            #Other buckets get assigned correct values, but we need to manually clear bucket_3
            bucket_3 = []
            #Now store results as a solution set
            solution_set = [bucket_1, bucket_2, bucket_3]
            print(solution_set)
            solution_set_list.append(copy.deepcopy(solution_set))
            #Moves largest element from bucket_2 to bucket_3, ends when bucket_2 is empty
            while len(bucket_2) != 0:
                #Add largest element from bucket_2 to bucket_3
                bucket_3.insert(0, bucket_2[-1])
                #Remove largest element from bucket_2
                ### This is probably where the error is occuring
                bucket_2.pop()
                #Store solution set
                solution_set = [bucket_1, bucket_2, bucket_3]
                solution_set_list.append(copy.deepcopy(solution_set))
                print(solution_set)
                ##I think our use of continue here could be problematic
                #continue
            ticker += 1
            if ticker > i:
                solution_set_derivation[f"list_n{i}_k{len(solution_set_list)}"] = solution_set_list
                acid_length_dict[i] = f"list_n{i}_k{len(solution_set_list)}"
            #continue


##########  Variable declarations (Info)  ##########

# Dictionary containing acids and their sequences
sequence_dict = {"Alanine" : ["GCA", "GCC", "GCG", "GCU"],
           "Arginine" : ["CGA", "CGC", "CGG", "CGU"],
           "Asparagine" : ["AAC", "AAU"],
           "Aspartic Acid" : ["GAC", "GAU"],
           "Cysteine" : ["UGC", "UGU"],
           "Glutamic Acid" : ["GAA", "GAG"],
           "Glutamine" : ["CAA", "CAG"],
           "Glycine" : ["GGA", "GGC", "GGG", "GGU"],
           "Histidine" : ["CAC", "CAU"],
           "Isoleucine" : ["AUA", "AUC", "AUU"],
           "Leucine" : ["CUA", "CUC", "CUG", "CUU", "UUA", "UUG"],
           "Lysine" : ["AAA", "AAG"],
           "Methionine" : ["AUG"],
           "Phenylalanine" : ["UUC", "UUU"],
           "Proline" : ["CCA", "CCC", "CCG", "CCU"],
           "Serine" : ["AGC", "AGU", "UCA", "UCC", "UCG", "UCU"],
           #"Selenocysteine": [],
           "Threonine" : ["ACA", "ACC", "ACG", "ACU"],
           "Tryptophan" : ["UGG"],
           "Tyrosine" : ["UAC", "UAU"],
           "Valine" : ["GUA", "GUC", "GUG", "GUU"]
           }

#List of the amino acids we are concerned with
full_acid_list = ['Alanine', "Arginine", 'Aspartic Acid', 'Asparagine', 'Cysteine', 'Glutamic Acid',
                  'Glutamine', 'Glycine', 'Histidine', 'Isoleucine', 'Leucine',
                  'Lysine', 'Methionine', 'Phenylalanine', 'Proline', 
                  #"Selenocysteine",
                  'Serine',
                  'Threonine', 'Tryptophan', 'Tyrosine', 'Valine']

#This is the version where we do not include Glycine and Proline
#Argument for excluding Glycine is that since it is included in all other graphs
#it might be a precondition e.g. it is the base acid constructed once an acid starts the build process
#This would then mean that ALL rules for Glycine base pairs ["GGA", "GGC", "GGG", "GGU"] are empty steps.
#This is a possibility that ultimately deserves to be considered in isolation (and alongside start/stop sequences)
solution_set_dict_alt1 = {
    "Alanine" : solution_set_derivation['list_n4_k15'], 
    "Arginine" : solution_set_derivation['list_n18_k190'],
    "Asparagine" : solution_set_derivation['list_n8_k45'], "Aspartic Acid" : solution_set_derivation['list_n6_k28'],
    "Cysteine" : solution_set_derivation['list_n5_k21'], "Glutamic Acid" : solution_set_derivation['list_n9_k55'], "Glutamine" : solution_set_derivation['list_n11_k78'],
    #"Glycine" : <doesn't work in this model>
    "Histidine": solution_set_derivation['list_n11_k78'], "Isoleucine" : solution_set_derivation['list_n13_k105'], "Leucine" : solution_set_derivation['list_n13_k105'],
    "Lysine" : solution_set_derivation['list_n16_k153'], "Methionine" : solution_set_derivation['list_n11_k78'], "Phenylalanine" : solution_set_derivation['list_n14_k120'],
    #"Proline" : <doesn't work in this model>
    #"Selenocysteine" : solution_set_derivation['list_n5_k21'], 
    "Serine" : solution_set_derivation['list_n5_k21'],
    "Threonine" : solution_set_derivation['list_n8_k45'], "Tryptophan" : solution_set_derivation['list_n18_k190'],
    "Tyrosine" : solution_set_derivation['list_n15_k136'], "Valine" : solution_set_derivation['list_n10_k66'] 
    }

#This version includes Glycine and Proline
solution_set_dictalt2 = {
    "Alanine" : solution_set_derivation['list_n4_k15'], 
    "Arginine" : solution_set_derivation['list_n18_k190'],
    "Asparagine" : solution_set_derivation['list_n8_k45'], "Aspartic Acid" : solution_set_derivation['list_n6_k28'],
    "Cysteine" : solution_set_derivation['list_n5_k21'], "Glutamic Acid" : solution_set_derivation['list_n9_k55'], "Glutamine" : solution_set_derivation['list_n11_k78'],
    #"Glycine" : <doesn't work in this model>
    "Histidine": solution_set_derivation['list_n11_k78'], "Isoleucine" : solution_set_derivation['list_n13_k105'], "Leucine" : solution_set_derivation['list_n13_k105'],
    "Lysine" : solution_set_derivation['list_n16_k153'], "Methionine" : solution_set_derivation['list_n11_k78'], "Phenylalanine" : solution_set_derivation['list_n14_k120'],
    #"Proline" : <doesn't work in this model>
    #"Selenocysteine" : solution_set_derivation['list_n5_k21'], 
    "Serine" : solution_set_derivation['list_n5_k21'],
    "Threonine" : solution_set_derivation['list_n8_k45'], "Tryptophan" : solution_set_derivation['list_n18_k190'],
    "Tyrosine" : solution_set_derivation['list_n15_k136'], "Valine" : solution_set_derivation['list_n10_k66'] 
    }

#This assigns the various node combinations to the amino acids. (And then changing them whenever we adjust criteria)
solution_set_dict = {}
for acid in full_acid_list:
    print(acid)
    print(len(acid_processing_df_dict[acid]))
    #Before normalization
    len_current_acid = len(acid_processing_df_dict[acid])
    #print(type (len_current_acid))
    solution_set_name = acid_length_dict[len_current_acid]
    solution_set_dict[acid] = solution_set_derivation[solution_set_name]

##########  IF selection (edit the values to select IF to test)  ##########
## This lets us generates all possible (assuming no references outside the codon) IFs to test en masse
permutation_results = []
for permutation in product(range(3), repeat=6):
    ## Requires that one of: [y > x, r > q, z > w] is true, and that none of these pairs has x, q, w exceeding y,r, z
    if (permutation[1] > permutation[0] or permutation[3] > permutation[2] or permutation[5] > permutation[4]) and permutation[1] >= permutation[0] and permutation[3] >= permutation[2] and permutation[5] >= permutation[4]:
        permutation_results.append(permutation)
        
#This stores results
answer_dict = {}
total_counter = 0
results_df = pd.DataFrame(['IF', 'Solution_array', 'Counter'])
for possible_IF in permutation_results:
    ###Set selection criteria for t_d_x here. First variable is lower bound, second is upper. (x,y), (q, r), (w, z)
    x = possible_IF[0]
    y = possible_IF[1]
    q = possible_IF[2]
    r = possible_IF[3]
    w = possible_IF[4]
    z = possible_IF[5]
    print('Current IF: ', x, y, q, r, w, z)
    IF_name = f'{x, y, q, r, w, z}'
    
    ##Original implementation broke all acids into three groups called overlap groups
    ##We're going to revisit this
    gate_order = []
    ##########  Time_dict_x generation  ##########
    time_dict_1 = {}
    time_dict_2 = {}
    time_dict_3 = {}
    
    #Builds the dictionaries we use to store values for rules.
    #Rules are generated based on our values selected for (x, y, q, r, w, z)
    #The combination of [(x, y), (q, r), (w, z)] represents an Interpretation framework
    #We use the variables to select parts of a base to assign to a rule.
    #For example, if acid_base = 'AGC', then acid_base[x:y] uses (x, y) to select how
    #the first time step of this problem is treated (we assume three time steps -- hence the number of time_dicts -- and allow empty steps)
    for acid in full_acid_list:
        acid_of_interest = sequence_dict[acid]
        
        # This layer then goes through the sequences of the selected acid
        # For example, Alanine has four bases, so we check each ["GCA", "GCC", "GCG", "GCU"]
        for sequence in acid_of_interest:
            #time_dict_1[sequence[x:y]] = pd.DataFrame()
            time_dict_1[sequence[x:y]] = pd.DataFrame()
            
        for sequence in acid_of_interest:
            #time_dict_2[sequence[q:r]] = pd.DataFrame()
            time_dict_2[sequence[q:r]] = pd.DataFrame()
            
        for sequence in acid_of_interest:
            #time_dict_3[sequence[w:z]] = pd.DataFrame()
            time_dict_3[sequence[w:z]] = pd.DataFrame()
            
    ## save_state_dict generation
    ### This seems to work, but should be rechecked ###
    ##This might not work because we later modify the time dicts. (But that might be the actual intent...)
    #Save state is supposed to allow us to create a reference point to revert to when running the main code (since conflicts will occur...)
    save_state_dict = {}
    for acid in full_acid_list:
        print(acid)
        save_state_dict[acid] = [copy.deepcopy(time_dict_1), copy.deepcopy(time_dict_2), copy.deepcopy(time_dict_3)]
    
    #Used to create the template for the IF's time_dict_Xs
    ###Should stay exactly the same
    empty_save_state_dict = copy.deepcopy(save_state_dict)
    
    save_state_dict = {
        "Alanine" : [time_dict_1, time_dict_2, time_dict_3], "Arginine" : [time_dict_1, time_dict_2, time_dict_3],
        "Asparagine" : [time_dict_1, time_dict_2, time_dict_3], "Aspartic Acid" : [time_dict_1, time_dict_2, time_dict_3],
        "Glycine" : [time_dict_1, time_dict_2, time_dict_3], "Proline": [time_dict_1, time_dict_2, time_dict_3],
        "Cysteine" : [time_dict_1, time_dict_2, time_dict_3], "Glutamic Acid" : [time_dict_1, time_dict_2, time_dict_3],
        "Glutamine" : [time_dict_1, time_dict_2, time_dict_3], "Histidine" : [time_dict_1, time_dict_2, time_dict_3],
        "Isoleucine" : [time_dict_1, time_dict_2, time_dict_3], "Leucine" : [time_dict_1, time_dict_2, time_dict_3],
        "Lysine" : [time_dict_1, time_dict_2, time_dict_3], "Methionine" : [time_dict_1, time_dict_2, time_dict_3],
        "Phenylalanine" : [time_dict_1, time_dict_2, time_dict_3], 
        #"Selenocysteine" : [time_dict_1, time_dict_2, time_dict_3],
        "Serine" : [time_dict_1, time_dict_2, time_dict_3], "Threonine" : [time_dict_1, time_dict_2, time_dict_3],
        "Tryptophan" : [time_dict_1, time_dict_2, time_dict_3], "Tyrosine": [time_dict_1, time_dict_2, time_dict_3],
        "Valine" : [time_dict_1, time_dict_2, time_dict_3]    
         
        }
    
    ##########  Function definition  ##########
    #(used to be k) tracks gate position
    gate_position = 0
    #o tracks position within solution_sets
    o = 0
    #Working values for the current run
    acid_solution_set = []
    #Where we store the values for all solutions
    answer_list = []
    #The order we want to test the acids in. Arranged loosely by eye test, so likely a high degree of optimization possible here
    #Trick is to arrange the acids so as to maximize initial overlap, so computational time is not wasted checking acids that don't overlap
    #More effort not put into optimization here because empirically this works well enough (most IFs finish in ~10-15K runs out of a possible ~10^15)
    #The "reason" this works is that we no longer need to keep checking branches once incompatability is determined between acids.
    #For example, if Glycine and Proline are incompatible with a given solution set, there is no reason to look at what the remaining acids do.
    gate_order = ['Glycine', 'Proline', 'Alanine', "Arginine", 'Aspartic Acid', 'Asparagine', 'Cysteine', 'Glutamic Acid',
                      'Glutamine', 'Histidine', 'Isoleucine', 'Leucine',
                      'Lysine', 'Methionine', 'Phenylalanine', 
                      #"Selenocysteine",
                      'Serine',
                      'Threonine', 'Tryptophan', 'Tyrosine', 'Valine']
    
    ###### Global part of Gate_Change, required for initialization...?
    ##The acid being used based on gate_order
    current_acid = gate_order[gate_position]
    #That acid's dataframe
    acid_df = acid_df_dict[current_acid]
    #The codons for that acid
    bases_current_acid = sequence_dict[current_acid]
    #The [[...], [...], [...]] used to group nodes and assign them to a time step
    #If the [] is empty, then no change is made at that step.
    testing_sets = solution_set_dict[current_acid]
    ########
    # Global part of testing set update
    testing_set = testing_sets[o]
    #We call this function when switching between acids, in turn updating key values
    def gate_change() :
        # We then use k to select the gate from our chosen gate_order
        # current_acid is the acid name as a string
        global current_acid
        #Does changing this to k - 1 screw us up?
        current_acid = gate_order[gate_position]
        
        # Returns a df of the acid for the current gate
        global acid_df
        acid_df = acid_df_dict[current_acid]
        
        # Pulls the list of bases for the acid from sequence_dict
        # bases_current_acid is a list of strings, each string is a base
        global bases_current_acid
        bases_current_acid = sequence_dict[current_acid]
        
        # Pulls appropriate list_n_x from solution_set_dict, which will be iterated over
        # testing_sets is a list of lists
        global testing_sets
        testing_sets = solution_set_dict[current_acid]
        
    # Changes the testing_set being used in the current run
    def testing_set_update() :
        # testing_set[o] gives the list we are using to determine time step composition
        global testing_set
        testing_set = testing_sets[o]
        
    #After a failed run, this reverts us to the last set of working values  
    def load_save_state():
        #When k = 0 we need special behavior to avoid going out of bounds on initial run
        #And then to clear save states when k = 0 is reached again as part of the process
        global time_dict_1
        global time_dict_2
        global time_dict_3
        global save_state_dict
        #No lower gate_position. 
        if gate_position == 0 and o > 0:
            #Should these all be copy.deepcopy() ?
            save_state_dict = copy.deepcopy(empty_save_state_dict)
            time_dict_1 = save_state_dict[gate_order[0]][0]
            time_dict_2 = save_state_dict[gate_order[0]][1]
            time_dict_3 = save_state_dict[gate_order[0]][2]
        else:
            # Use k - 1 and not k because we are going to retest k, and so k - 1 has last full success
            loading_acid = gate_order[gate_position - 1]
            # This load structure works because we always save the dicts in the same order
    
            time_dict_1 = save_state_dict[loading_acid][0]
            time_dict_2 = save_state_dict[loading_acid][1]
            time_dict_3 = save_state_dict[loading_acid][2]
    #The, ahem, "actual run" that does the heavy listing to check rules and their values        
    def actual_run():
        #These should all be working as intended but if the current round of simulations does not yield a notable result
        #I'll be going through it again for the n^fifth (?) time.
        global o
        global gate_position
        global acid_solution_set
        
        #First of three for loops. Each for loop checks to see if all bases for the acid work/match in the current IF 
        #(since rules for the base pairs of an acid might conflict with themselves -- recall that there is usually more than one "base pair -> amino acid" relationship per acid)
        #Uses (x, y)
        for base in bases_current_acid:
            
            #Checks to see if the first position of the testing_set is an empty step
            if testing_set[0] == [] :
                #Then checks for conflicts at specified location
                if time_dict_1[base[x:y]].equals(pd.DataFrame()) == False :
                    print(f"Overlap error in t_d_1 for {current_acid} at o = {o}")
                    o = o + 1
                    return
                
            else :
                #Start by normalizing the selection. This makes it so that the nodes from one graph can be compared to the other safely without
                #worrying about differing graph sizes. That is, two graph patterns might be the same, however without normalization
                #we wouldn't know that unless they had exactly the same size graph before the selection.
                #Technically normalization isn't necessary for the first time step as all graph subsection start at node 1 here.
                norm_dict = {}
                norm = 1
                for i in range(len(testing_set[0])):
                    #Creates a dictionary entry for the selected node number and assigns it a norm value
                    #For t_d_1 this seems unnecessary as all t_d_1 selections start with min node = 1.
                    norm_dict[testing_set[0][i]] = norm
                    norm += 1
                ##Now build norm_dict to handle exceptions above and below bounds; test length to ensure not running on an empty array []
                #Key assumption here: we treat the node type of a node outside the time_dict to be irrelevant,
                #meaning we only count the number of incoming and outgoing edges.
                #This is done because the information should (still testing this assumption) be contained in the other time dicts. (additionally, testing a found solution should provide clarity if needed)
                if (len(testing_set[0])) > 0:
                    #Lower bound
                    for i in range(1, min(norm_dict)):
                        norm_dict[i] = "Out of bounds lower"
                    #Upper bound
                    for i in range(max(norm_dict) + 1, max(testing_set)[-1] + 1):
                        norm_dict[i] = "Out of bounds upper"
                
                #This is the range of nodes we want to select for the given time step (note: recheck zero indexing in analysis)
                #Get acid info, selecting rows that appear in above expression
                ##This should overwrite previous df_selection value, so I don't think we need to declare it empty (will check)
                df_selection = copy.deepcopy(acid_df_dict[current_acid][acid_df_dict[current_acid].Node_1.isin(range(min(testing_set[0]), max(testing_set[0]) + 1))])
                #Drop Node_1, Node_2 values to avoid comparing against them
                df_selection = df_selection.drop('Node_1', axis=1)
                df_selection = df_selection.drop('Node_2', axis=1)
                #Normalize the normalization columns
                df_selection = df_selection.replace({'Normalized_node_1': norm_dict, 'Normalized_node_2': norm_dict})
                #The value in t_d_1 is empty, so take normalized info and store appropriately.
                #We only store data when the associated rule is empty, otherwise we just check for equivalence.
                if time_dict_1[base[x:y]].equals(pd.DataFrame()) == True :
                    
                    #Store in t_d_1
                    time_dict_1[base[x:y]] = copy.deepcopy(df_selection)
                #The value in t_d_1 is not empty, so we need to normalize and check info
                elif time_dict_1[base[x:y]].equals(pd.DataFrame()) == False :
                    #Compare normalized df_selection to stored value
                    if time_dict_1[base[x:y]].equals(df_selection) == False:
                        print(f"Overlap error in t_d_1 for {current_acid} at o = {o}")
                        o = o + 1
                        return
        
        #Loop for time_dict_2 (uses (q, r))     
        for base in bases_current_acid:
            
            #Checks to see if the first position of the testing_set is an empty step
            if testing_set[1] == [] :
                #Then checks for conflicts at specified location
                if time_dict_2[base[q:r]].equals(pd.DataFrame()) == False :
                    print(f"Overlap error in t_d_2 for {current_acid} at o = {o}")
                    o = o + 1
                    return
            
            #If not an empty time step we execute this code block   
            else :
                #See normalization explanation in time_dict_1 section
                norm_dict = {}
                norm = 1
                for i in range(len(testing_set[1])):
                    #Creates a dictionary entry for the selected node number and assigns it a norm value
                    #For t_d_1 this seems unnecessary as all t_d_1 selections start with min node = 1.
                    norm_dict[testing_set[1][i]] = norm
                    norm += 1
                ##Now build norm_dict to handle exceptions above and below bounds; test length to ensure not running on an empty array []
                if (len(testing_set[1])) > 0:
                    #Lower bound
                    for i in range(1, min(norm_dict)):
                        norm_dict[i] = "Out of bounds lower"
                    #Upper bound
                    for i in range(max(norm_dict) + 1, max(testing_set)[-1] + 1):
                        norm_dict[i] = "Out of bounds upper"
                #This is the range of nodes we want to select for the given time step
                #careful with zero-indexing...
                #Get acid info, selecting rows that appear in above expression
                ##This should overwrite previous df_selection value, so I don't think we need to declare it empty (will check)
                df_selection = copy.deepcopy(acid_df_dict[current_acid][acid_df_dict[current_acid].Node_1.isin(range(min(testing_set[1]), max(testing_set[1]) + 1))])
                #Drop Node_1, Node_2 values
                df_selection = df_selection.drop('Node_1', axis=1)
                df_selection = df_selection.drop('Node_2', axis=1)
                #Normalize the normalization columns
                df_selection = df_selection.replace({'Normalized_node_1': norm_dict, 'Normalized_node_2': norm_dict})
                #The value in t_d_1 is empty, so take normalized info and store appropriately
                if time_dict_2[base[q:r]].equals(pd.DataFrame()) == True :
                    #Store in t_d_1
                    time_dict_2[base[q:r]] = copy.deepcopy(df_selection)
                #The value in t_d_1 is not empty, so we need to normalize and check info
                elif time_dict_2[base[q:r]].equals(pd.DataFrame()) == False :
                    #Normalized df_selection and stored value are not equivalent, so we need to move on to the next o.
                    if time_dict_2[base[q:r]].equals(df_selection) == False:
                        print(f"Overlap error in t_d_2 for {current_acid} at o = {o}")
                        o = o + 1
                        return
        
        #Loop for time_dict_3 (same as time_dict_2, but uses (w, z))
        for base in bases_current_acid:
            
            #Checks to see if the first position of the testing_set is an empty step
            if testing_set[2] == [] :
                #Then checks for conflicts at specified location
                if time_dict_3[base[w:z]].equals(pd.DataFrame()) == False :
                    print(f"Overlap error in t_d_3 for {current_acid} at o = {o}")
                    o = o + 1
                    return
            
            #If not an empty time step we execute this code block   
            else :
                #See normalization explanation in time_dict_1 section
                norm_dict = {}
                norm = 1
                for i in range(len(testing_set[2])):
                    #Creates a dictionary entry for the selected node number and assigns it a norm value
                    #For t_d_1 this seems unnecessary as all t_d_1 selections start with min node = 1.
                    norm_dict[testing_set[2][i]] = norm
                    norm += 1
                ##Now build norm_dict to handle exceptions above and below bounds
                if (len(testing_set[2])) > 0:
                    #Lower bound
                    for i in range(1, min(norm_dict)):
                        norm_dict[i] = "Out of bounds lower"
                    #Upper bound
                    for i in range(max(norm_dict) + 1, max(testing_set)[-1] + 1):
                        norm_dict[i] = "Out of bounds upper"
                #This is the ran
                
                #This is the range of nodes we want to select for the given time step
                #careful with zero-indexing...
                #Get acid info, selecting rows that appear in above expression
                ##This should overwrite previous df_selection value, so I don't think we need to declare it empty (will check)
                df_selection = copy.deepcopy(acid_df_dict[current_acid][acid_df_dict[current_acid].Node_1.isin(range(min(testing_set[2]), max(testing_set[2]) + 1))])
                #Drop Node_1, Node_2 values
                df_selection = df_selection.drop('Node_1', axis=1)
                df_selection = df_selection.drop('Node_2', axis=1)
                #Normalize the normalization columns
                df_selection = df_selection.replace({'Normalized_node_1': norm_dict, 'Normalized_node_2': norm_dict})
                #The value in t_d_1 is empty, so take normalized info and store appropriately
                if time_dict_3[base[w:z]].equals(pd.DataFrame()) == True :
                    
                    #Store in t_d_1
                    time_dict_3[base[w:z]] = copy.deepcopy(df_selection)
                #The value in t_d_1 is not empty, so we need to normalize and check info
                elif time_dict_3[base[w:z]].equals(pd.DataFrame()) == False :
                    #Normalized df_selection and stored value are not equivalent, so we need to move on to the next o.
                    if time_dict_3[base[w:z]].equals(df_selection) == False:
                        print(f"Overlap error in t_d_3 for {current_acid} at o = {o}")
                        o = o + 1
                        return
        
        #If we clear the preceding three for loops, then there were no conflicts and we've found a potential solution
        #We store the current o value.
        acid_solution_set.append(o)
        print(f"Value found for {current_acid} at gate_position = {gate_position}; o = {o}")
        #Captures current values for each t_d_X. These represent current working answers and will be reverted
        #to if subsequent iterations in the next step reveal conflicts.
        save_state_dict[current_acid] = copy.deepcopy([time_dict_1, time_dict_2, time_dict_3])
        #Values for the current acid work, move on to the next gate
        gate_position = gate_position + 1
        #Not sure about this one
        #Since the above continues until a solution is found (or exit functions when impossible), we need to trigger gate_change() to move on
        if gate_position < len(gate_order):
            gate_change()
    
    ########## The code that executes the search ##########
    counter = 0
    #Our running condition
    while gate_position >= 0:
        print(f'Counter iteration: {counter}')
        #This is the first condition run, it is only run once per Interpretation Framework.
        if gate_position == 0 and o == 0:
            actual_run()
            counter += 1
            print()
            
        #This is the condition that runs the most. 
        elif gate_position <= len(gate_order) - 1 and (o < len(testing_sets) and o > 0):
            #Likely modified a t_d_x in previous run, so need to load the last working values
            load_save_state()
            testing_set_update()
            print(f"acid: {current_acid} gate_position: {gate_position} o: {o}")
            actual_run()
            counter += 1 
            
        #This condition occurs right after a gate_change() and so there is no need to call load_save_state
        elif gate_position <= len(gate_order) - 1 and o == 0:
            testing_set_update()
            print(f"acid: {current_acid} gate_position: {gate_position} o: {o}")
            actual_run()
            counter += 1 
            
        #If this condition is met, we've found a solution (tested every gate and found no conflicts, so the last gate finding a working value pushes us past len(gate_order) - 1)    
        elif gate_position > len(gate_order) - 1:
            print("An answer has been found using acid_solution_set = ", acid_solution_set)
            answer = copy.deepcopy(acid_solution_set)
            answer_list.append(answer)
            print("Number of answers is now: ", len(answer_list))
            #Return to previous gate to check for more solutions
            gate_position = gate_position - 1
            #We've changed gates, so we update our values
            gate_change()
            #Pick up from where we left off
            o = acid_solution_set[-1]
            #Remove the last value to search for another solution
            acid_solution_set.pop()
            #Increment by 1, so we're checking a new combination
            o = o + 1
            #Now we load the save_state
            load_save_state()
            continue
        
        #This is our ending condition for the current Interpretation framework.
        #It stores information for later analysis.
        #Should this be len(testing_sets) - 1? (Pretty sure it shouldn't be, but will check again)
        elif gate_position == 0 and o >= len(testing_sets):
            print(f"Program has finished running. There were {len(answer_list)} answers found. Counter: {counter}")
            #This is where we put code to move between overlap_groups if we're testing acid combinations separately.
            answer_dict[IF_name] = answer_list
            results_df.loc[len(results_df.index)] = [IF_name, answer_list, counter]
            total_counter += counter
            break
        #All testing_set values for given gate have been exhausted. Move back to previous gate and increment by 1.
        elif o >= len(testing_sets) :
            print(f"No solution found for {current_acid}. Moving to previous gate.")
            gate_position = gate_position - 1
            #The gate has changed, so we must once again update
            gate_change()
            #Set value for testing set
            o = acid_solution_set[-1]
            #Remove old value
            acid_solution_set.pop()
            #print()
            o = o + 1
            #Don't need to load_save_state() since that will be done in the next cycle of the while loop.
          
        elif counter >= 1000000:
            print(f"Program has exceeded acceptable computational bounds, moving on to next framework. There were {len(answer_list)} answers found. Counter: {counter}")
            #This is where we put code to move between overlap_groups if we're testing acid combinations separately.
            answer_dict[IF_name] = answer_list
            results_df.loc[len(results_df.index)] = [IF_name, answer_list, counter]
            break
            
        else:
            print("An unaccounted for scenario occurred. Exiting program.")
            print(f"current_acid: {current_acid} gate_position: {gate_position} o: {o}")
            break