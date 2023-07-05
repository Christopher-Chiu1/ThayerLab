# -*- coding: utf-8 -*-
"""
Created on Mon Feburary 27th, 19:32:00 2023
Working on Thursday April 13th, 11:21:00 2023
mat_rearrange complete on Monday June 5th, 06:05:00 2023
switch_sector_values complete on Thursday June 8th, 02:50:00 2023

@author: Christopher Chiu, cchiu01@wesleyan.edu
"""

import pandas as pd
import seaborn
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#from hmmlearn import hmm

"""
This is the importing of the datasets section. You can theoretically import any size matrix,
but the data has to be run through the CC_MD_SECTORS.sh function in the 
/zfshomes/cchiu01/UpdatedMD_SECTORS/ directory.
"""

#base_dir4 = "C:/Work Stuff/4-10-23outputlarge.csv" #Import CSV File to Python
#largefixed = pd.read_csv(base_dir4) #Make the CSV File into a readable Dataframe

#seaborn.heatmap(largefixed, cmap="rocket_r")

#base_dir5 = "C:/Work Stuff/4-13-231tupstrip.csv"
#largefixed2 = pd.read_csv(base_dir5) #The 1tup Stripped Dataset

#base_dir6 = "C:/Work Stuff/4-17-23stripoutput.csv"
#datazeroed = pd.read_csv(base_dir6)

#base_dir7 = "C:/Work Stuff/20_y220c_1ms_protein_only_output2.csv"
#y220c_1ms = pd.read_csv(base_dir7)

#base_dir8 = "C:/Work Stuff/wt_1ms_protien_only_output.csv"
#wt_1ms = pd.read_csv(base_dir8)

#base_dir9 = "C:/Work Stuff/6sm_1ms_protein_only_output.csv"
#_6sm_1ms = pd.read_csv(base_dir9)

# Zero Bonded Residues 

def zero_bonded_residues(dataset): #calls for a dataframe variable through pandas
    """
    This function takes a dataset from the pd.read_csv, and finds any one in the dataset
    and converts it into a zero, along with any residue above or below it. This removes
    bonded residues, which are obviously correlated already, and not important to our
    calculations.
    
    Parameters
    ----------
    dataset : DataFrame
        The Input is a Pandas DataFrame, and relies on a loaded in CSV file from Excel.
    
    Returns
    -------
    There is no real return, but, In the Variable Table, the dataset's ones and 
    neighbors are converted to zeros
    """
    z = 0 #counter for the entire function
    i = 0 #header iterator
    k = 0 #column iterator
    while z <= len(dataset)-2: 
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."
        
        y = dataset[str(accheader)] #Selects each column of the datast
        for c in y: #Iterates through column
            if (y == 1).bool: #Anything next to a 1 is 0'd. Doesn't work on edge cases yet
                y[k] = 0
                y[k-1] = 0
                y[k+1] = 0
        del y #Deletes Column Variables as when you create y, a new column variable is created
        k += 1
        
        z += 1
    
    #This is the edge case of the first column of the inputted dataset
    y = dataset["c0001"]
    y[1] = 0
    del y
    
    #This is the edge case for the last column of the inputted dataset
    dataset.iloc[len(dataset)-1, -1] = 0
    dataset.iloc[len(dataset)-2, -1] = 0
    
    #Creates a heatmap with the new data
    #seaborn.heatmap(dataset, cmap="rocket_r")
    return(dataset)

def calc_pica_values(dataset):
    """
    From each column, it adds up all of the values into one, and assigns each column
    a total vpica value, ex ['c0001', .5666]
    
    Parameters
    ----------
    dataset : DataFrame
        This function takes a zeroed dataset from the zero_bonded_residues function.

    Returns
    -------
    sorted_vpica_values: a list of tuples that show each column's overall vpica value
    """
    z = 0 #counter for the entire function
    i = 0 #header iterator
    r = 0 #Vpica Value Iterator
    counter = 0 #iterates for the Vpica values. Simply makes new variables in the column_pica_values list
    column_pica_values = {} #Create a Dictionary to Contain the Data Set's Vpica values
    while z <= len(dataset)-1: 
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."
        
        y = dataset[str(accheader)] #Selects each column of the datast
        placeholdervpica = 0.0
        for data in y: #Iterates through column
            placeholdervpica = placeholdervpica + data
        del data
        placeholdervpica = placeholdervpica/len(dataset)
        column_pica_values[str(accheader)] = placeholdervpica
        del y #Deletes Column Variables as when you create y, a new column variable is created
        
        z += 1
        counter += 1
    
    Vpica_Value = 0.0
    i = 0
    while r < len(dataset): #Loop that adds all of the columns together into a single Vpica Value
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."
        Vpica_Value = Vpica_Value + column_pica_values[str(accheader)]
        r += 1
    #Vpica_Value = Vpica_Value/len(dataset) #Takes Vpica Value and Divides it by the dataset's len
    
    
    sorted_vpica_values = sorted(column_pica_values.items(), key=lambda x:x[1], reverse=True) #Sorts Vpica Values of the Residues from Highest to Lowest
    #print("The Residues Sorted by Highest to Lowest Vpica Values are as Follows:",sorted_vpica_values)
    print("The Overall Vpica Value of Your Dateset is:",Vpica_Value)
    
    highest_20_percent = []
    i = 0 #Reset "i" iterator to zero
    d = int(len(dataset)/5) #Top 20% of Residues
    while i < d: #Takes the top 20% of residues and 
        highest_20_percent.append(sorted_vpica_values[i])
        i += 1
        
    #print("The Sector, The Top 20% of Residues, is:", highest_20_percent)
    return(sorted_vpica_values)


def calc_sector_residues(vpica_values, dataset):
    """
    This function take the vpica values list, and shortens it to only the residues in
    the sector. So the top 20% of the residues are the only ones included into a new list

    Parameters
    ----------
    vpica_values : list
        output from the calc_pica_values function
    dataset : DataFrame
        same DataFrame that has been used in the calc_pica_values function

    Returns
    -------
    highest_20_percent: a new list that only includes the sector residues

    """
    highest_20_percent = []
    i = 0 #Reset "i" iterator to zero
    d = int(len(dataset)/5) #Top 20% of Residues
    while i < d: #Takes the top 20% of residues and 
        highest_20_percent.append(vpica_values[i])
        i += 1
    
    return(highest_20_percent)


def mat_rearrange(vpica_values, dataset):
    """
    The dataset inputted is rearranged by column according to the highest vpica values
    from the calc_pica_values function. 

    Parameters
    ----------
    vpica_values : list (of tuples)
        A list of tuples is read reassign each column according to highest vpica values
    dataset : DataFrame
        This is the same input for the calc_pica_values function.

    Returns
    -------
    Rearranged Dataset

    """
    out_mat = dataset.copy()
    z = 0 #counter for the entire function
    i = 0 #header iterator
    r = 0 #iterator that copies values from the original matrix
    t = 0 #iterator that copies through the new matrix
    list_of_rearranged_dataframes = []
    while z <= len(dataset)-1: 
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."
        
        x = out_mat[str(accheader)] #Selects each column of the datast
        y = dataset[str(vpica_values[z][0])] #Selects each column of the dataset
        r = 0 #Have to reset r and t each time to count through each column
        t = 0
        while t <= len(dataset)-1: #Iterates through columns
            #So, we have to jump from point to point in the column, so, we have to go to
            #point 5 first, which we can achieve by making a new variable, m that represents
            #the row we want, we essentially need to redefine it in rank order of the 
            #highest residue vpica values.
            m = vpica_values[r][0]
            if m[1] != str(0): #This section of ifs takes the 'c0001, c0002' and narrows it down to just '0001, 0002' etc so we can choose the row more precisely.
                m = m[1:5]
                #print("Success 1")
                #print(m)
                
            elif m[2] != str(0):
                m = m[2:5]
                #print("Success 2")
                #print(m)
                
            elif m[3] != str(0):
                m = m[3:5]
                #print("Success 3")
                #print(m)
                
            elif m[4] != str(0):
                m = m[4:5]
                #print("Success 4")
                #print(m)
              
            x[t] = y[int(m)-1] #Makes the original data set point equal to the new dataset point
            
            t += 1 #iterates counter up
            r += 1
            
        p = z % 5
        if p == 0:
            list_of_rearranged_dataframes.append(out_mat.copy())
        
        z += 1
        
    return(out_mat, list_of_rearranged_dataframes)


def calc_sector_values(dataset, sectorvalues):
    """
    This function does the same thing as the calc_pica_value function, however, it only
    includes residues that are in the sector, so therefore, if there are 10 residues in
    the sector, this function will only take the first 10 by 10 box of the DataFrame.
    
    Parameters
    ----------
    vpica_values : list
        This is the output of the calc_vpica_value function, and is input here to find the
        same list, but only with sector values
    dataset : DataFrame
        This should be the rearranged dataset as the output from mat_rearrange

    Returns
    -------
    sorted_vpica_values: This is the output vpica values, but only for the residues within
    the sector.
    """
    z = 0 #counter for the entire function
    i = 0 #header iterator
    r = 0 #Vpica Value Iterator
    counter = 0 #iterates for the Vpica values. Simply makes new variables in the column_pica_values list
    column_pica_values = {} #Create a Dictionary to Contain the Data Set's Vpica values
    while z <= len(sectorvalues)-1: 
        r = 0 #Reset the Counter Variable
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."
        
        y = dataset[str(accheader)] #Selects each column of the datast
        placeholdervpica = 0.0
        while r <= len(sectorvalues)-1: #Iterates through column
            placeholdervpica = placeholdervpica + y[r]
            r += 1
        placeholdervpica = placeholdervpica/len(sectorvalues)
        column_pica_values[str(accheader)] = placeholdervpica
        del y #Deletes Column Variables as when you create y, a new column variable is created
        
        z += 1
        counter += 1
    
    Vpica_Value = 0.0
    i = 0 #Reset Counter Variable
    r = 0 #Reset Counter Variable
    while r < len(sectorvalues): #Loop that adds all of the columns together into a single Vpica Value
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."
        Vpica_Value = Vpica_Value + column_pica_values[str(accheader)]
        r += 1
    
    
    sorted_vpica_values = sorted(column_pica_values.items(), key=lambda x:x[1], reverse=True) #Sorts Vpica Values of the Residues from Highest to Lowest
    print("The Overall Vpica Value of Your Sector is:",Vpica_Value)
    
    return(sorted_vpica_values)


def calc_outside_sector_values(sorteddataset, sector_values):
    """
    This function is used when you want to calculate the value of the matrix excluding the 
    sector. This can be applied to see how the sector's value changes when you run other 
    functions on it. This function takes a sorted dataset, and the output of the 
    calc_sector_values function. 
    
    Parameters
    ----------
    vpica_values : list
        This is the output of the calc_vpica_value function, and is input here to find the
        same list, but only without sector values
    dataset : DataFrame
        This should be the rearranged dataset as the output from mat_rearrange

    Returns
    -------
    sorted_vpica_values: This is the output vpica values, but only for the residues out of
    the sector.
    """
    z = 0 #counter for the entire function
    i = 0 #header iterator
    r = 0 #Local Iterator
    counter = 0 #iterates for the Vpica values. Simply makes new variables in the column_pica_values list
    column_pica_values = {} #Create a Dictionary to Contain the Data Set's Vpica values
    while z <= len(sorteddataset)-1: 
        r = 0 #Reset the Counter Variable
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."

        y = sorteddataset[str(accheader)] #Selects each column of the datast
        placeholdervpica = 0.0
        if z <= len(sector_values)-1:
            r = len(sector_values)
            while r <= len(sorteddataset)-1: #Iterates through column
                placeholdervpica = placeholdervpica + y[r]
                #print(r)
                r += 1
            placeholdervpica = placeholdervpica/len(sorteddataset)
            column_pica_values[str(accheader)] = placeholdervpica
            #print(placeholdervpica)
            del y #Deletes Column Variables as when you create y, a new column variable is created
        
        else:
            while r <= len(sorteddataset)-1: #Iterates through column
                placeholdervpica = placeholdervpica + y[r]
                r += 1
            placeholdervpica = placeholdervpica/len(sorteddataset)
            column_pica_values[str(accheader)] = placeholdervpica
            del y #Deletes Column Variables as when you create y, a new column variable is created
            
        
        z += 1
        counter += 1
    
    Vpica_Value = 0.0
    i = (len(sorteddataset)/5)-1 #Set Counter Variable To Start After the Sector
    r = len(sorteddataset)/5 #Set Counter Variable To Start After the Sector
    r = int(r)
    i = int(i)
    while r < len(sorteddataset): #Loop that adds all of the columns together into a single Vpica Value
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."
        Vpica_Value = Vpica_Value + column_pica_values[str(accheader)]
        r += 1
    #Vpica_Value = Vpica_Value/len(dataset) #Takes Vpica Value and Divides it by the dataset's len
    
    
    sorted_vpica_values = sorted(column_pica_values.items(), key=lambda x:x[1], reverse=True) #Sorts Vpica Values of the Residues from Highest to Lowest

    return(sorted_vpica_values)


def tuple_to_list(tuplelist):
    """
    This function is just a quick way to change tuples made inside of functions into lists
    which are callable. The function allows for easier manipulation of the data, and is
    only run inside other functions.
    
    Parameters
    ----------
    tuplelist : list
        A list of tuples. Such as an output of calc_pica_values
    
    Returns
    -------
    A list of lists, instead of a list of tuples.
    """
    y = 0
    while y < len(tuplelist):
        o = tuplelist[y]
        o = list(o)
        tuplelist[y] = o
       #print("success", y)
        y += 1
        
        
def list_to_tuple(listlist):
    
    y = 0
    while y < len(listlist):
        o = listlist[y]
        o = tuple(o)
        listlist[y] = o
       #print("success", y)
        y += 1

def switch_sector_values(sorteddataset, sector_values):
    """
    This function requires a sorted dataset that has had the function, "mat_rearrange,"
    used on it. After that, the input of sector_values is the output for the, 
    "calc_sector_values" function. Once those two inputs are satisfied, run the function.
    
    Parameters
    ----------
    sorteddataset : DataFrame
        The output from the mat_rearrange function
    sector_values : TYPE
        Output from the calc_sector_values function

    Returns
    -------
    A list that is rearranged to check the last five residues in a sector if they need a 
    switch. If one of the next five residues out of the sector is found to be higher than one
    of the last five residues in the sector, it is swapped out. 
    """
    z = 0 #counter for the entire function
    i = 0 #header iterator
    r = 0 #Vpica Value Iterator
    counter = 0 #iterates for the Vpica values. Simply makes new variables in the column_pica_values list
    column_pica_values = {} #Create a Dictionary to Contain the Data Set's Vpica values
    test_vpica_values = {} #Create a Dictionary to Contain the Test Vpica Values
    while z <= len(sector_values)-1: 
        r = 0 #Reset the Counter Variable
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."
        
        y = sorteddataset[str(accheader)] #Selects each column of the datast
        placeholdervpica = 0.0
        while r <= len(sector_values)-1: #Iterates through column
            placeholdervpica = placeholdervpica + y[r]
            r += 1
        placeholdervpica = placeholdervpica/len(sector_values)
        column_pica_values[str(accheader)] = placeholdervpica
        del y #Deletes Column Variables as when you create y, a new column variable is created
        
        z += 1
        counter += 1

    
    while z <= len(sorteddataset)-1:
        r = 0 #Reset the Counter Variable
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."
        
        y = sorteddataset[str(accheader)] #Selects each column of the datast
    
        placeholdervpica = 0.0
        while r <= len(sector_values)-1: #Iterates through column
            placeholdervpica = placeholdervpica + y[r]
            r += 1
        placeholdervpica = placeholdervpica/len(sector_values)
        test_vpica_values[str(accheader)] = placeholdervpica
        del y #Deletes Column Variables as when you create y, a new column variable is created
        
        z += 1
        counter += 1
    
    sorted_vpica_values = sorted(column_pica_values.items(), key=lambda x:x[1], reverse=True)
    sorted_test_vpica = sorted(test_vpica_values.items(), key=lambda x:x[1], reverse=True)
    
    tuple_to_list(sorted_vpica_values)
    tuple_to_list(sorted_test_vpica)
    
    
    Vpica_Value = 0.0
    i = 0
    r = 0
    while r < len(sector_values)-1: #Loop that adds all of the columns together into a single Vpica Value
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."
        Vpica_Value = Vpica_Value + column_pica_values[str(accheader)]
        r += 1  
    
    
    r = 0

    all_vpica_values = sorted_vpica_values + sorted_test_vpica
    
    list_to_tuple(all_vpica_values)
    
    sorted_all_pica = sorted(all_vpica_values, key=lambda x:x[1], reverse=True)
    
    splicedvpica = sorted_all_pica[:len(sector_values)]
    
    Vpica_Value_New = 0.0
    i = 0
    r = 0
    while r < len(sector_values)-1: #Loop that adds all of the columns together into a single Vpica Value
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."
        Vpica_Value_New = Vpica_Value_New + splicedvpica[r][1]
        r += 1  
    
    print("Your Old Vpica Sector Value is:", Vpica_Value) 
    print("Your New Vpica Sector Value is:", Vpica_Value_New)
    
    return(sorted_all_pica)

def cut_out_second_sector(dataset, sector):
    sec_mat = dataset.copy()
    
    z = 0 #counter for the entire function
    i = 0 #header iterator
    o = len(dataset) #header backwards iterator
    r = 0 #iterator that copies values from the original matrix
    t = 0 #iterator that copies through the new matrix
    #list_of_rearranged_dataframes = []
    while z <= len(dataset)-1: 
        i += 1
        accheader = ""
        acc = str(i)
        while len(acc) < 4:
            acc = "0" + acc
        
        accheader = "c"+acc # "c0001, c0002, c0003....."
        
        x = sec_mat[str(accheader)] #Selects each column of the datast
        
        occheader = ""
        occ = str(o)
        while len(occ) < 4:
            occ = "0" + occ
            
        occheader = "c"+occ # "c0001, c0002, c0003....."
        
        y = dataset[str(occheader)] #Selects each column of the dataset
        r = 0 #Have to reset r and t each time to count through each column
        t = len(dataset)-1
        while t >= 0: #Iterates through columns
            
            x[r] = y[t] #Makes the original data set point equal to the new dataset point
            
            t -= 1 #iterates counter up
            r += 1
            
        o -= 1
        z += 1
        
    z = len(dataset)-len(sector)
    o = len(dataset) #header backwards iterator
    t = len(dataset)-1
    while z < o:
        occheader = ""
        occ = str(o)
        while len(occ) < 4:
            occ = "0" + occ
            
        occheader = "c"+occ # "c0001, c0002, c0003....."
        
        sec_mat = sec_mat.drop(t)
        del sec_mat[str(occheader)]
    
        o -= 1   
        t -= 1
    
    seaborn.heatmap(sec_mat, cmap="rocket_r")
    return(sec_mat)


def auto_MD_SECTORS_Calc(dataset):
    """
    

    Parameters
    ----------
    dataset : DataFrame
        A Pandas Dataframe. As long as it is a square matrix, outputted by the 
        CC_MD_SECTORS.sh function in the Thayer Lab's HPCC, this function can be
        ran, and outputs will appear in the python console

    Returns
    -------
    None.

    """
    print("How many times would you like to repeat recalculation? Input Must Be An Integer. \n")
    z = int(input())
    i = 0
    rearrangedlists = []
    vpica_values_list = []
    datasetzeroed = zero_bonded_residues(dataset)
    print("A Heatmap of the original inputted dataset has been made in the 'Plots' area in the Variable Table \n")
    generate_heatmap(datasetzeroed)
    print("Zeroed Success!")
    pica_values = calc_pica_values(datasetzeroed)
    print("Calculate Vpica Values Sucess!")
    sector_values = calc_sector_values(datasetzeroed, pica_values)
    #print(len(sector_values))
    print("Calculate Sector Values Sucess!")
    sector_residues = calc_sector_residues(sector_values, datasetzeroed)
    #print(len(sector_residues))
    print("Calculate Sector Residues Sucess!")
    rearranged, rearrangesteps = mat_rearrange(pica_values, datasetzeroed)
    #print("Rearrange Success!")
    sector_processed = calc_sector_values(rearranged, sector_residues)
    print("Sector Processed Success!")
    vpica_switch = switch_sector_values(rearranged, sector_processed)
    print("Vpica Switch Success! \n")
    
    while (i <= z):
        rearrangedlists.append(rearranged.copy())
        vpica_values_list.append(sector_processed.copy())
        rearranged, rearrangesteps = mat_rearrange(vpica_switch, rearranged)
        #seaborn.heatmap(after_swap_large, cmap="rocket_r")
        
        print("Rearrange #",i)
        
        pica_values = calc_pica_values(rearranged)
        print("Calculate Vpica Values Sucess!")
        sector_values = calc_sector_values(rearranged, pica_values)
        print("Calculate Sector Values Sucess!")
        sector_residues = calc_sector_residues(sector_values, rearranged)
        print("Calculate Sector Residues Sucess!")
        sector_processed = calc_sector_values(rearranged, sector_residues)
        print("Sector Processed Success!")  
        vpica_switch = switch_sector_values(rearranged, sector_processed)
        print("Vpica Switch Success! \n")
        i += 1
        
    seaborn.heatmap(rearranged, cmap="rocket_r") 
    return(rearrangedlists, sector_values, vpica_values_list, datasetzeroed)


def generate_heatmap(data):
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(data, cmap="rocket_r", ax=ax)
    plt.show()

# Example usage
#rearranged = pd.DataFrame(...)  # Your input DataFrame

def generate_heatmap_steps(dataframelist):
    i = 0
    while i <= len(dataframelist)-1:
        generate_heatmap(dataframelist[i])
        i += 1
    
def run_everything():
    print("This function runs the entire sectors process. Simply input the .csv file, and the process will run automatically.")
    print("Please Enter The Directory of the .csv File: (Example: C:/WorkFiles/4-10-23outputlarge.csv)")
    base_dir_dataset = str(input())
    test_dataset = pd.read_csv(base_dir_dataset)
    #print("Please Enter the Name of Your Output Variables")
    #placeholdervariablename = input()
    
    test_listoflists, test_sectorvalues, test_vpica_values_list, datasetzeroed = auto_MD_SECTORS_Calc(test_dataset)
    counter = 0
    print("These are the sector values being changed by each loop. Look at how the residue numbers change after each print, and you can see how the sector changes. Do this in pair with the seaborn heatmaps to see the sector visually change at the top left of the function.")
    for i in test_vpica_values_list:
        counter += 1
        print("Sector Swap #",counter)
        print(i)
        print("\n")    
        
    print("Now, to visualize the steps of the changing vpica values, use the test_lists of lists variable.")
    print("This is done using the generate_heatmap_steps function, the test_listoflists variable will be ran into it and the heatmaps will be visualized in the 'plots' section in the variable explorer")
    print("\n")
    
    #print("Would you like to cut out a potential second sector in the bottom right of the final dataset? You can find it modelled in the 'plots' menu in the variable explorer. 1 = yes, 0 = no")
    #cutout = int(input())
    cutout = 1
    if cutout == 1:
        print("Cutout Variables are stored as the same variable names, but with 'cut' in front of them.")
        final_dataset = test_listoflists[len(test_listoflists)-1]
        test_cutout = cut_out_second_sector(final_dataset, test_vpica_values_list[0])
        
        cut_test_listoflists, cut_test_sectorvalues, cut_test_vpica_values_list, datasetzeroed = auto_MD_SECTORS_Calc(test_cutout)
        
        counter = 0
        print("These are the sector values being changed by each loop. Look at how the residue numbers change after each print, and you can see how the sector changes. Do this in pair with the seaborn heatmaps to see the sector visually change at the top left of the function.")
        for i in cut_test_vpica_values_list:
            print("Cutout Sector Swap #",counter)
            print(i)
            print("\n")
            counter += 1
    
    print("Heatmaps Are Now Being Generated")
    generate_heatmap_steps(test_listoflists)
    print("Original Heatmap to Seperate Full Matrix from Cut Matrix")
    generate_heatmap(datasetzeroed)
    print("The Cutout Heatmaps Will Now Be Generated")
    generate_heatmap_steps(cut_test_listoflists)
    print("MD Sectors Process Complete!")
    print("After each switch is essentially a rerun of a rearrange function that tests if the residues in the sector actually belong inside of the sector. So in order to see how the sector changes, just look at the residues after each sector swap, and there potentially can be a change in the sector residues.")
    return(test_listoflists, test_sectorvalues, test_vpica_values_list, cut_test_listoflists, cut_test_sectorvalues, cut_test_vpica_values_list)
    
    
test_listoflists, test_sectorvalues, test_vpica_values_list, cut_test_listoflists, cut_test_sectorvalues, cut_test_vpica_values_list = run_everything()

