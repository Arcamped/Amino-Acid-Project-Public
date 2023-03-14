# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 16:28:22 2021

@author: Michael Campbell
"""
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

# Load the CSV file containing simple Amino Acid data
df2 = pd.read_csv("C:/Users/Computer/Documents/GitHub/Amino-Acid-Project/Fixed Acid Data.csv")

# Data is still a collection of strings, so it needs to be split on position
# into three new columns
df_mod = df2
df_mod['AminoAcid'] = df_mod['AminoAcid'].astype(str)
df_mod['BasePairs'] = df_mod['BasePairs'].astype(str)

# Create a new column and populate it with the timestep of the first position
## We should be able to change this to fit framework needs
df_mod['pos_1'] = df_mod['BasePairs']
df_mod['pos_1'] = df_mod['pos_1'].str[0:1]

# Same process but for the second timestep position
df_mod['pos_2'] = df_mod['BasePairs']
df_mod['pos_2'] = df_mod['pos_2'].str[1:2]

# Third timestep position 
df_mod['pos_3'] = df_mod['BasePairs']
df_mod['pos_3'] = df_mod['pos_3'].str[2:3]


# Interpretation Framework mapping
# We now neeed to assign rules to each acid based on their pos_x value

df_mod['r_pos_1'] = ""
df_mod['r_pos_2'] = ""
df_mod['r_pos_3'] = ""
# for the elif statement of rules we are just going to use r_i with i = x
# and specify the x


# For the first position
# The current issue is with mistreating the column as a series
# The for and elif loop needs to be adjusted to put one
# rule for each row...
for i in range(len(df_mod['pos_1'])) :
    if df_mod['pos_1'][i] in ['A'] :
        df_mod['r_pos_1'][i] = '1'
    elif df_mod['pos_1'][i] in ['U'] :
        df_mod['r_pos_1'][i] = '2'
    elif df_mod['pos_1'][i] in ['C'] :
        df_mod['r_pos_1'][i] = '3'
    elif df_mod['pos_1'][i] in ['G']:
        df_mod['r_pos_1'][i] = '4'
     
# For the second position

for i in range(len(df_mod['pos_2'])) :
    if df_mod['pos_2'][i] in ['A'] :
        df_mod['r_pos_2'][i] = '5'
    elif df_mod['pos_2'][i] in ['U'] :
        df_mod['r_pos_2'][i] = '6'
    elif df_mod['pos_2'][i] in ['C'] :
        df_mod['r_pos_2'][i] = '7'
    elif df_mod['pos_2'][i] in ['G']:
        df_mod['r_pos_2'][i] = '8'


# For the third position

for i in range(len(df_mod['pos_3'])) :
    if df_mod['pos_3'][i] in ['A'] :
        df_mod['r_pos_3'][i] = '9'
    elif df_mod['pos_3'][i] in ['U'] :
        df_mod['r_pos_3'][i] = '10'
    elif df_mod['pos_3'][i] in ['C'] :
        df_mod['r_pos_3'][i] = '11'
    elif df_mod['pos_3'][i] in ['G']:
        df_mod['r_pos_3'][i] = '12'

# Interpretation Framework results/analysis

# Build a table with acids on both axes, and show # of rule overlap

# First we need to extract the names of the acids 

acid_list = df_mod['AminoAcid'].unique()

#And take just the columns with rules
r_cols = ['r_pos_1', 'r_pos_2', 'r_pos_3']
# This gives unique rules by position
grouped_by = df_mod.groupby('AminoAcid')[r_cols].nunique()

# This lets us extract only the rows of interest
t_cols = ['AminoAcid', 'r_pos_1', 'r_pos_2', 'r_pos_3']
rule_df = df_mod[t_cols]
# This currently gives us unique values in the columns and we still need
# to group by AminoAcid... (We might want to write a function...?)
np.unique(df_mod[r_cols].values)

#This currently gives us list objects of unique rules for each acid
# We need to make them usable to compare to each other for overlap

acid_dict = {}
for acid in acid_list :
    work_df = rule_df[rule_df['AminoAcid'] == acid]
    acid_dict[acid] = np.unique(work_df[r_cols].values)
    print(np.unique(work_df[r_cols].values))

# Creating the df to store our counts
final_df = pd.DataFrame(index=acid_list, columns=acid_list)
    
# Using the dictionary that has unique values for each acid
# The next line lists how we can examine which rules overlap
#list(set(acid_dict['Valine']) & set(acid_dict['Alanine']))
#final_df.loc['Valine', 'Alanine'] = len(list(set(acid_dict['Valine']) & set(acid_dict['Alanine'])))
# if we wrap the expression with len() then we get the count

# First for loop set the row, then iterates
# Second for loop moves through each column filling in values
for acid in acid_list :
    
    for acid_2 in acid_list :
        total = len(list(set(acid_dict[acid]) & set(acid_dict[acid_2])))
        final_df.loc[acid, acid_2] = total

# Acid(s) with the largest rule set

# Acid(s) with the smallest rule set (n = 3)

# Acids with no overlapping rules?

# Any acids that subsume another(one contains all rules of the other)
# How would we write a formula for that?
