# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 17:00:52 2022

@author: Computer
"""
import copy
import csv

#lower and upper bound, respectively
l = 0
m = 18
counter = 1
#soluton_set_derivation
solution_set_derivation = {}
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
            #continue