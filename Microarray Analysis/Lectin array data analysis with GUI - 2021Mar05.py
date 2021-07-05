#############################################################################
############################ Import Modules #################################
#############################################################################

import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
import numpy as np
import pandas as pd

#############################################################################
########################## Initialize Parameters ############################
#############################################################################

root = tk.Tk()
root.title("Lectin Array Data Analysis  - Version 2021.Mar.05")

#The types of following parameters only denote their types when they are eventually passed to data processing
analysis_mode = 0
raw_file_path_list = []
sample_file_path_list = []
lectin_file_path = "0"
output_dir = "0"
n_reps = 0
Grubbs_co = 0
SNR_co_532 = 0
SNR_co_635 = 0
min_sample_co = 0
#min_sample_co --- minimum value of [number of samples in which a lectin passed the SNR cut-offs / number of all samples] to keep the lectin (between 0-1)
raw_value_type = 0 
#raw_value_type --- 1: median; 2: mean
array_norm_method = 0 
#array_norm_method ---
#1: for each array, normalize intensities of individual spots to the median intensity of the corresponding channel of all spots
#2: for each array, normalize intensities of individual spots to the mean intensity  of the corresponding channel of all spots
#3: for each array, normalize the ratios of individual spots to the median ratio of all spots
#4: for each array, normalize the ratios of individual spots to the mean ratio of all spots
#5: no normalization
initiate_analysis = False


### Variables in tkinter
raw_data_path = ()
sample_path = ()
lectin_path = ()

raw_data_number = tk.StringVar() # This is the only way to change the text of an existing label
sample_number = tk.StringVar()
lectin_file_number = tk.StringVar()
output_path = tk.StringVar()

raw_data_number.set('No file chosen...')
sample_number.set('No file chosen...')
lectin_file_number.set('No file chosen...')
output_path.set('No directory chosen...')

analysis_mode_tk = tk.IntVar()
analysis_mode_tk.set(2)
raw_value_type_tk = tk.IntVar()
raw_value_type_tk.set(1)
array_norm_method_tk = tk.IntVar()
array_norm_method_tk.set(1)

#############################################################################
############################ Define Functions ###############################
#############################################################################

###########
########### tkinter functions
###########
def entrystate():
    if analysis_mode_tk.get() == 1:
        SNR_co_532_entry.config(state = 'disabled')
        SNR_co_532_entry.update()
        SNR_co_635_entry.config(state = 'disabled')
        SNR_co_635_entry.update()
        minsample_cutoff_entry.config(state = 'disabled')
        minsample_cutoff_entry.update()
    else:
        SNR_co_532_entry.config(state = 'normal')
        SNR_co_532_entry.update()
        SNR_co_635_entry.config(state = 'normal')
        SNR_co_635_entry.update()
        minsample_cutoff_entry.config(state = 'normal')
        minsample_cutoff_entry.update()

def choosefile_raw_data():
    global raw_data_path
    raw_data_path = tk.filedialog.askopenfilename(filetypes = [("TXT", ".txt")], multiple = True)
    raw_data_number.set(str(len(raw_data_path)) + " file(s) chosen")

def choosefile_lectin():
    global lectin_path
    lectin_path = tk.filedialog.askopenfilename(filetypes = [("TXT", ".txt")], multiple = True)
    lectin_file_number.set(str(len(lectin_path)) + " file(s) chosen")

def choosefile_sample():
    global sample_path
    sample_path = tk.filedialog.askopenfilename(filetypes = [("TXT", ".txt")], multiple = True)
    sample_number.set(str(len(sample_path)) + " file(s) chosen")

def choosedir_output():
    output_tk = tk.filedialog.askdirectory()
    output_path.set(output_tk)


def checkfilenumbers():
    global raw_file_path_list, sample_file_path_list

    if type(raw_data_path) == type((0,0)) and type(sample_path) == type((0,0)):
        raw_file_path_list = list(raw_data_path)
        raw_file_path_list.sort(reverse = False)
        sample_file_path_list = list(sample_path)
        sample_file_path_list.sort(reverse = False)
        if len(raw_file_path_list) == len(sample_file_path_list):
            passparameters()
        else:
            tk.messagebox.showwarning(title = "Error", message = "The number of raw data files and the number of sample list files must be the same!")
    else:
        tk.messagebox.showwarning(title = "Error", message = "An error occured when reading the files. Please double check your inputs or contact the author for assistance.")
        root.destroy()

def passparameters():
        global raw_file_path_list, sample_file_path_list, analysis_mode, lectin_file_path, output_dir, n_reps, Grubbs_co, min_sample_co, raw_value_type, array_norm_method, SNR_co_635, SNR_co_532, initiate_analysis

        raw_file_path_list = list(raw_data_path)
        raw_file_path_list.sort(reverse = False)
        sample_file_path_list = list(sample_path)
        sample_file_path_list.sort(reverse = False)

        analysis_mode = int(analysis_mode_tk.get())
        lectin_file_path = lectin_path[0]
        output_dir = output_path.get()
        n_reps = int(n_reps_entry.get())
        Grubbs_co = float(Grubbs_cutoff_entry.get())
        raw_value_type = int(raw_value_type_tk.get())
        initiate_analysis = True

        if analysis_mode == 2:
            SNR_co_532 = float(SNR_co_532_entry.get())
            SNR_co_635 = float(SNR_co_635_entry.get())
            min_sample_co = float(minsample_cutoff_entry.get())
            array_norm_method = int(array_norm_method_tk.get())

        tk.messagebox.showinfo(title = "Message", message = "The analysis process can take up to a few minutes. You will get a message about the status of the analysis later. Press the OK button to proceed.")
        root.destroy() ###close the main input window



###########
########### functions for data processing
###########

def replist(x, y):
    return [item for item in x for i in range(y)]
    #x: list; y: times each element in x is repeated

def Gval_532(x):
    if x["STD 532"] == 0:
        return 0
    else:
        return abs(x["F532"] - x["Mean 532"])/x["STD 532"]

def Gval_635(x):
    if x["STD 635"] == 0:
        return 0
    else:
        return abs(x["F635"] - x["Mean 635"])/x["STD 635"]

def Grubbs_check_dc(x):
    global Grubbs_co
    if (x[13] == x[15] and x[14] == x[16] and x[13] > Grubbs_co and x[14] > Grubbs_co):
        return False
    else:
        return True

def Grubbs_check_sc(x):
    global Grubbs_co
    if (x["Gval 635 max"] == x["Gval 635"] and x["Gval 635 max"] > Grubbs_co):
        return False
    else:
        return True

#############################################################################
############################## User Interface ###############################
#############################################################################

### row 0: info
textLabel = tk.Label(root, text = "Author: Rui (Ric) Qin; rq2@ualberta.ca", padx = 10, pady = 10, justify = tk.LEFT) #justify: alignment
textLabel.grid(row = 0)

### row 1-2: Choose analysis mode
tk.Label(root, text = "Choose analysis mode").grid(row = 1, column = 0, sticky = tk.W)
tk.Radiobutton(root, text = "Single color", value = 1, variable = analysis_mode_tk, command = entrystate).grid(row = 2, column = 0, sticky = tk.W)
tk.Radiobutton(root, text = "Dual color", value = 2, variable = analysis_mode_tk, command = entrystate).grid(row = 2, column = 1, sticky = tk.W)

###row 3: divider
tk.Label(root, text = "*************************************************************************************************************************").grid(row = 3, column = 0, columnspan = 2, sticky = tk.W)

### row 4: Choose raw data file
openrawfilebotton = tk.Button(root, text = "Choose Genepix Raw Data File(s)", command = choosefile_raw_data)
openrawfilebotton.grid(row = 4, column = 0, sticky = tk.NSEW)
tk.Label(root, textvariable = raw_data_number, width = 40).grid(row = 4, column = 1, sticky = tk.W)

### row 5: Choose sample list file
opensamplebotton = tk.Button(root, text = "Choose Sample List File(s)", command = choosefile_sample)
opensamplebotton.grid(row = 5, column = 0, sticky = tk.NSEW)
tk.Label(root, textvariable = sample_number, width = 40).grid(row = 5, column = 1, sticky = tk.W)

###row 6: Choose lectin list file
openlectinbotton = tk.Button(root, text = "Choose Lectin List File", command = choosefile_lectin)
openlectinbotton.grid(row = 6, column = 0, sticky = tk.NSEW)
tk.Label(root, textvariable = lectin_file_number, width = 40).grid(row = 6, column = 1, sticky = tk.W)

###row 7: Choose output directory
openlectinbotton = tk.Button(root, text = "Choose Output Folder", command = choosedir_output)
openlectinbotton.grid(row = 7, column = 0, sticky = tk.NSEW)
tk.Label(root, textvariable = output_path, width = 40).grid(row = 7, column = 1, sticky = tk.W)

###row 8: divider
tk.Label(root, text = "*************************************************************************************************************************").grid(row = 8, column = 0, columnspan = 2, sticky = tk.W)

###row 9: Number of replicates
tk.Label(root, text = "Number of replicates of lectins per array").grid(row = 9, column = 0, sticky = tk.W)
n_reps_entry = tk.Entry(root)
n_reps_entry.insert(0, 3)
n_reps_entry.grid(row = 9, column = 1, sticky = tk.W)

###row 10: Cut-off for Grubbs' test
tk.Label(root, text = "Cut-off for Grubbs' test").grid(row = 10, column = 0, sticky = tk.W)
Grubbs_cutoff_entry = tk.Entry(root)
Grubbs_cutoff_entry.insert(0, 1.15)
Grubbs_cutoff_entry.grid(row = 10, column = 1, sticky = tk.W)

###row 11: Cut-off of the 532nm channel
tk.Label(root, text = "SNR cut-off (532nm channel)").grid(row = 11, column = 0, sticky = tk.W)
SNR_co_532_entry = tk.Entry(root)
SNR_co_532_entry.insert(0, 3)
SNR_co_532_entry.grid(row = 11, column = 1, sticky = tk.W)

###row 12: Cut-off of the 635nm channel
tk.Label(root, text = "SNR cut-off (635nm channel)").grid(row = 12, column = 0, sticky = tk.W)
SNR_co_635_entry = tk.Entry(root)
SNR_co_635_entry.insert(0, 3)
SNR_co_635_entry.grid(row = 12, column = 1, sticky = tk.W)

###row 13: Sample number cut-off
tk.Label(root, text = "Percentage samples cutoff (between 0-1)").grid(row = 13, column = 0, sticky = tk.W)
minsample_cutoff_entry = tk.Entry(root)
minsample_cutoff_entry.insert(0, 0.5)
minsample_cutoff_entry.grid(row = 13, column = 1, sticky = tk.W)

###row 14: divider
tk.Label(root, text = "*************************************************************************************************************************").grid(row = 14, column = 0, columnspan = 2, sticky = tk.W)

###row 15-16: Raw fluorescence value type
tk.Label(root, text = "Choose type of intensity data for analysis").grid(row = 15, column = 0, sticky = tk.W)
tk.Radiobutton(root, text = "Median intensities of spots", value = 1, variable = raw_value_type_tk).grid(row = 16, column = 0, sticky = tk.W)
tk.Radiobutton(root, text = "Mean intensities of spots", value = 2, variable = raw_value_type_tk).grid(row = 16, column = 1, sticky = tk.W)

###row 17: divider
tk.Label(root, text = "*************************************************************************************************************************").grid(row = 17, column = 0, columnspan = 2, sticky = tk.W)

###row 18-23: Array normalization method
tk.Label(root, text = "Choose array normalization method (DUAL COLOR ONLY)").grid(row = 18, column = 0, sticky = tk.W)
tk.Radiobutton(root, text = "For each array, normalize intensities of individual spots to the median intensity of all spots in the corresponding channel", value = 1, variable = array_norm_method_tk).grid(row = 19, column = 0, columnspan = 2, sticky = tk.W)
tk.Radiobutton(root, text = "For each array, normalize intensities of individual spots to the mean intensity of all spots in the corresponding channel", value = 2, variable = array_norm_method_tk).grid(row = 20, column = 0, columnspan = 2, sticky = tk.W)
tk.Radiobutton(root, text = "For each array, normalize the ratios of individual spots to the median ratio of all spots", value = 3, variable = array_norm_method_tk).grid(row = 21, column = 0, columnspan = 2, sticky = tk.W)
tk.Radiobutton(root, text = "For each array, normalize the ratios of individual spots to the mean ratio of all spots", value = 4, variable = array_norm_method_tk).grid(row = 22, column = 0, columnspan = 2, sticky = tk.W)
tk.Radiobutton(root, text = "Do not normalize", value = 5, variable = array_norm_method_tk).grid(row = 23, column = 0, columnspan = 2, sticky = tk.W)

###row 24: run analysis
tk.Button(root, text = "Run Analysis", command = passparameters).grid(row = 24, sticky = tk.NSEW)

#mainloop
root.mainloop()

#############################################################################
###################### Single Color Data Processing #########################
#############################################################################

###########
###########Step 1: annotate each raw input file and combine them
###########

if initiate_analysis and analysis_mode == 1:
    try:
        lectin_list = pd.read_csv(lectin_file_path)
        lectin_list_list = lectin_list.iloc[:,0].tolist()

        ##Global parameters
        n_lectins = lectin_list.iloc[:,0].count()
        n_slides = len(raw_file_path_list)

        for i in range(n_slides):
            ##Read raw input file and sample list, one slide at a time
            raw_file_path = raw_file_path_list[i]
            sample_file_path = sample_file_path_list[i]

            raw_data_input = pd.read_csv(raw_file_path, delimiter = "\t", skiprows = 32, na_values = "Error")
            sample_list = pd.read_csv(sample_file_path)

            ##Slide-specific parameters
            n_rows = raw_data_input["Row"].max()
            n_cols = raw_data_input["Column"].max()
            n_samples_i = sample_list.iloc[:,0].count()

            #select the columns of interest from the raw_input dataframe
            if raw_value_type == 1:
                raw_input_coi = raw_data_input.loc[:,["Block", "Column", "Flags", "F635 Median - B635", "SNR 635"]]
            else:
                raw_input_coi = raw_data_input.loc[:,["Block", "Column", "Flags", "F635 Mean - B635", "SNR 635"]]

            #generate a list of lectins, each lectin repeated n_reps times in each block
            lectin_list_rep = (replist(lectin_list_list, n_reps))*n_samples_i
            lectin_list_rep_df = pd.DataFrame(lectin_list_rep)
            lectin_list_rep_df.columns = ["Lectin"] #Another way: lectin_list_rep_df.rename(columns = {0:'Lectin'}, inplace = True)

            #generate a list of samples, each lectin repeated n_rows*n_cols times in each block
            sample_list_list = sample_list.iloc[:,0].tolist()
            sample_list_rep = replist(sample_list_list,n_rows*n_cols)
            sample_list_rep_df = pd.DataFrame(sample_list_rep)
            sample_list_rep_df.columns = ["Sample"]
            
            if i == 0:
                sample_list_all = sample_list_list
            else:
                sample_list_all.extend(sample_list_list) #extend the master sample list

            #merge the above three dataframes and write to file
            raw_data_anno = pd.concat([raw_input_coi, lectin_list_rep_df, sample_list_rep_df], axis = 1) #axis = 1: join horizontally
            raw_data_anno.columns = ["Block", "Column", "Flags", "F635", "SNR 635", "Lectin", "Sample"]

            if i == 0:
                raw_data_anno_combined = raw_data_anno
            else:
                raw_data_anno_combined = pd.concat([raw_data_anno_combined, raw_data_anno], axis = 0)

        n_samples = len(sample_list_all)
        raw_data_anno_combined.reset_index(drop = True, inplace = True)
        raw_data_anno_combined.to_csv(output_dir + "/data_raw_annotated.txt", index = False, sep='\t')


        ###########
        ###########Step 2: Grubbs' test
        ###########
        data_step2_grouped = raw_data_anno_combined.groupby(by = ["Sample", "Lectin"], sort = False)

        Grubbs_mean_635 = data_step2_grouped.mean()[["F635"]] 
        Grubbs_mean_635_list = Grubbs_mean_635.iloc[:,0].tolist()
        Grubbs_mean_635_rep = replist(Grubbs_mean_635_list,3)
        Grubbs_mean_635_rep_df = pd.DataFrame(Grubbs_mean_635_rep)
        Grubbs_mean_635_rep_df.columns = ["Mean 635"]

        Grubbs_std_635 = data_step2_grouped.std()[["F635"]] 
        Grubbs_std_635_list = Grubbs_std_635.iloc[:,0].tolist()
        Grubbs_std_635_rep = replist(Grubbs_std_635_list,3)
        Grubbs_std_635_rep_df = pd.DataFrame(Grubbs_std_635_rep)
        Grubbs_std_635_rep_df.columns = ["STD 635"]

        raw_data_anno_combined_Grubbs = pd.concat([raw_data_anno_combined, Grubbs_mean_635_rep_df, Grubbs_std_635_rep_df], axis = 1)

        raw_data_anno_combined_Grubbs["Gval 635"] = raw_data_anno_combined_Grubbs.apply(Gval_635, axis = 1)

        data_step2_grouped_2 = raw_data_anno_combined_Grubbs.groupby(by = ["Sample", "Lectin"], sort = False)

        Grubbs_max_635 = data_step2_grouped_2.max()[["Gval 635"]]
        Grubbs_max_635_list = Grubbs_max_635.iloc[:,0].tolist()
        Grubbs_max_635_rep = replist(Grubbs_max_635_list,3)
        Grubbs_max_635_rep_df = pd.DataFrame(Grubbs_max_635_rep)
        Grubbs_max_635_rep_df.columns = ["Gval 635 max"]

        raw_data_anno_combined_Grubbs = pd.concat([raw_data_anno_combined_Grubbs, Grubbs_max_635_rep_df], axis = 1)
        raw_data_anno_combined_Grubbs["Result"] = raw_data_anno_combined_Grubbs.apply(Grubbs_check_sc, axis = 1)
        data_step2 = raw_data_anno_combined_Grubbs[raw_data_anno_combined_Grubbs["Result"]]
        data_step2.reset_index(drop = True, inplace = True)

        ###########
        ###########Step 3: Construct a final dataframe and write to file
        ###########
        data_step3_grouped = data_step2.groupby(by = ["Sample", "Lectin"], sort = False)
        data_step3 = data_step3_grouped.mean()["F635"]
        data_step3 = pd.DataFrame(np.array(data_step3).reshape((n_samples, n_lectins)))
        data_step3_int = pd.DataFrame(np.rint(data_step3))

        data_step3_int.columns = lectin_list_list
        data_step3_int.index = sample_list_all
        data_step3_int = data_step3_int.transpose()
        data_step3_int.to_csv(output_dir + "/data_final.txt", index = True, sep = "\t")
        tk.messagebox.showinfo(title = "Message", message = "Congratulations! The analysis is finished.\nPlease find the output files in the designated folder.")
    except:
        tk.messagebox.showwarning(title = "Error", message = "An error occured during the analysis process. However, output files may have been generated if the error is not fatal. In case of a fatal error, Please double check your inputs or contact the author for assistance.") 


#############################################################################
####################### Dual Color Data Processing ##########################
#############################################################################

if initiate_analysis and analysis_mode == 2:
    try:
        lectin_list = pd.read_csv(lectin_file_path)
        lectin_list_list = lectin_list.iloc[:,0].tolist()

        ##Global parameters
        n_lectins = lectin_list.iloc[:,0].count()
        n_slides = len(raw_file_path_list)

        ###########
        ###########Step 1: annotate each raw input file and combine them
        ###########
        for i in range(n_slides):
            ##Read raw input file and sample list, one slide at a time
            raw_file_path = raw_file_path_list[i]
            sample_file_path = sample_file_path_list[i]

            raw_data_input = pd.read_csv(raw_file_path, delimiter = "\t", skiprows = 32, na_values = "Error")
            sample_list = pd.read_csv(sample_file_path)

            ##Slide-specific parameters
            n_rows = raw_data_input["Row"].max()
            n_cols = raw_data_input["Column"].max()
            n_samples_i = sample_list.iloc[:,0].count()

            #select the columns of interest from the raw_input dataframe
            if raw_value_type == 1:
                raw_input_coi = raw_data_input.loc[:,["Block", "Column", "Flags", "F532 Median - B532", "F635 Median - B635", "SNR 532", "SNR 635"]]
            else:
                raw_input_coi = raw_data_input.loc[:,["Block", "Column", "Flags", "F532 Mean - B532", "F635 Mean - B635", "SNR 532", "SNR 635"]]

            #generate a list of lectins, each lectin repeated n_reps times in each block
            lectin_list_rep = (replist(lectin_list_list, n_reps))*n_samples_i
            lectin_list_rep_df = pd.DataFrame(lectin_list_rep)
            lectin_list_rep_df.columns = ["Lectin"] #Another way: lectin_list_rep_df.rename(columns = {0:'Lectin'}, inplace = True)

            #generate a list of samples, each lectin repeated n_rows*n_cols times in each block
            sample_list_list = sample_list.iloc[:,0].tolist()
            sample_list_rep = replist(sample_list_list,n_rows*n_cols)
            sample_list_rep_df = pd.DataFrame(sample_list_rep)
            sample_list_rep_df.columns = ["Sample"]
            
            if i == 0:
                sample_list_all = sample_list_list
            else:
                sample_list_all.extend(sample_list_list) #extend the master sample list

            #merge the above three dataframes and write to file
            raw_data_anno = pd.concat([raw_input_coi, lectin_list_rep_df, sample_list_rep_df], axis = 1) #axis = 1: join horizontally
            raw_data_anno.columns = ["Block", "Column", "Flags", "F532", "F635", "SNR 532", "SNR 635", "Lectin", "Sample"]

            if i == 0:
                raw_data_anno_combined = raw_data_anno
            else:
                raw_data_anno_combined = pd.concat([raw_data_anno_combined, raw_data_anno], axis = 0)

        n_samples = len(sample_list_all)
        raw_data_anno_combined.reset_index(drop = True, inplace = True)
        raw_data_anno_combined.to_csv(output_dir + "/data_raw_annotated.txt", index = False, sep='\t')


        ###########
        ###########Step 2: Grubbs' test
        ###########
        data_step2_grouped = raw_data_anno_combined.groupby(by = ["Sample", "Lectin"], sort = False)

        Grubbs_mean_532 = data_step2_grouped.mean()[["F532"]] #if use [] this will be a Series. [[]] makes it a dataframe.
        Grubbs_mean_532_list = Grubbs_mean_532.iloc[:,0].tolist()
        Grubbs_mean_532_rep = replist(Grubbs_mean_532_list,3)
        Grubbs_mean_532_rep_df = pd.DataFrame(Grubbs_mean_532_rep)
        Grubbs_mean_532_rep_df.columns = ["Mean 532"]

        Grubbs_std_532 = data_step2_grouped.std()[["F532"]]
        Grubbs_std_532_list = Grubbs_std_532.iloc[:,0].tolist()
        Grubbs_std_532_rep = replist(Grubbs_std_532_list,3)
        Grubbs_std_532_rep_df = pd.DataFrame(Grubbs_std_532_rep)
        Grubbs_std_532_rep_df.columns = ["STD 532"]

        Grubbs_mean_635 = data_step2_grouped.mean()[["F635"]] 
        Grubbs_mean_635_list = Grubbs_mean_635.iloc[:,0].tolist()
        Grubbs_mean_635_rep = replist(Grubbs_mean_635_list,3)
        Grubbs_mean_635_rep_df = pd.DataFrame(Grubbs_mean_635_rep)
        Grubbs_mean_635_rep_df.columns = ["Mean 635"]

        Grubbs_std_635 = data_step2_grouped.std()[["F635"]] 
        Grubbs_std_635_list = Grubbs_std_635.iloc[:,0].tolist()
        Grubbs_std_635_rep = replist(Grubbs_std_635_list,3)
        Grubbs_std_635_rep_df = pd.DataFrame(Grubbs_std_635_rep)
        Grubbs_std_635_rep_df.columns = ["STD 635"]

        raw_data_anno_combined_Grubbs = pd.concat([raw_data_anno_combined, Grubbs_mean_532_rep_df, Grubbs_std_532_rep_df, Grubbs_mean_635_rep_df, Grubbs_std_635_rep_df], axis = 1)

        raw_data_anno_combined_Grubbs["Gval 532"] = raw_data_anno_combined_Grubbs.apply(Gval_532, axis = 1)
        raw_data_anno_combined_Grubbs["Gval 635"] = raw_data_anno_combined_Grubbs.apply(Gval_635, axis = 1)
        #spot pass this test if at least one channel has a grubbs statistics less than cutoff

        data_step2_grouped_2 = raw_data_anno_combined_Grubbs.groupby(by = ["Sample", "Lectin"], sort = False)

        Grubbs_max_532 = data_step2_grouped_2.max()[["Gval 532"]]
        Grubbs_max_532_list = Grubbs_max_532.iloc[:,0].tolist()
        Grubbs_max_532_rep = replist(Grubbs_max_532_list,3)
        Grubbs_max_532_rep_df = pd.DataFrame(Grubbs_max_532_rep)
        Grubbs_max_532_rep_df.columns = ["Gval 532 max"]

        Grubbs_max_635 = data_step2_grouped_2.max()[["Gval 635"]]
        Grubbs_max_635_list = Grubbs_max_635.iloc[:,0].tolist()
        Grubbs_max_635_rep = replist(Grubbs_max_635_list,3)
        Grubbs_max_635_rep_df = pd.DataFrame(Grubbs_max_635_rep)
        Grubbs_max_635_rep_df.columns = ["Gval 635 max"]

        raw_data_anno_combined_Grubbs = pd.concat([raw_data_anno_combined_Grubbs, Grubbs_max_532_rep_df, Grubbs_max_635_rep_df], axis = 1)
        raw_data_anno_combined_Grubbs["Result"] = raw_data_anno_combined_Grubbs.apply(Grubbs_check_dc, axis = 1)
        data_step2 = raw_data_anno_combined_Grubbs[raw_data_anno_combined_Grubbs["Result"]]
        data_step2.reset_index(drop = True, inplace = True)

        ###########
        ###########Step 3: SNR annotation
        ###########
        data_step3 = data_step2.iloc[:,list(range(13))]
        data_step3_grouped = data_step3.groupby(by = ["Sample", "Lectin"], sort = False)

        SNR_mean_532 = data_step3_grouped.mean()["SNR 532"].tolist()
        SNR_count_532 = data_step3_grouped.count()["SNR 532"].tolist()
        SNR_mean_635 = data_step3_grouped.mean()["SNR 635"].tolist()
        SNR_count_635 = data_step3_grouped.count()["SNR 635"].tolist()

        SNR_mean_532_rep = []
        for i in range(len(SNR_count_532)):
            SNR_mean_532_temp = [SNR_mean_532[i]]*int(SNR_count_532[i])
            SNR_mean_532_rep.extend(SNR_mean_532_temp)
        SNR_mean_532_rep_df = pd.DataFrame(SNR_mean_532_rep)
        SNR_mean_532_rep_df.columns = ["SNR 532 Mean"]

        SNR_mean_635_rep = []
        for i in range(len(SNR_count_635)):
            SNR_mean_635_temp = [SNR_mean_635[i]]*int(SNR_count_635[i])
            SNR_mean_635_rep.extend(SNR_mean_635_temp)
        SNR_mean_635_rep_df = pd.DataFrame(SNR_mean_635_rep)
        SNR_mean_635_rep_df.columns = ["SNR 635 Mean"]

        SNR_count_rep = []
        for i in range(len(SNR_count_635)):
            SNR_count_temp = [SNR_count_635[i]]*int(SNR_count_635[i])
            SNR_count_rep.extend(SNR_count_temp)
        SNR_count_rep_df = pd.DataFrame(SNR_count_rep)
        SNR_count_rep_df.columns = ["Count"]

        data_step3 = pd.concat([data_step3, SNR_mean_532_rep_df, SNR_mean_635_rep_df, SNR_count_rep_df], axis = 1)
        data_step3["SNR Check Result"] = data_step3.apply(lambda x: 1/x["Count"] if (x["SNR 532 Mean"] >= SNR_co_532 and x["SNR 635 Mean"] >= SNR_co_635) else 0, axis = 1)


        ###########
        ###########Step 4: Remove lectins
        ###########
        data_step4_grouped = data_step3.groupby(by = ["Sample", "Lectin"], sort = False)
        SNR_result_by_sample = np.array(data_step4_grouped.sum()["SNR Check Result"].tolist()).reshape((n_samples, n_lectins))
        SNR_final_count = SNR_result_by_sample.sum(axis = 0).flatten().tolist()

        lectin_failed = [lectin_list_list[i] for i in range(n_lectins) if SNR_final_count[i] < (min_sample_co * n_samples)] #List comprehension
        lectin_passed = [i for i in lectin_list_list if i not in lectin_failed]

        SNR_filter = ~data_step3.Lectin.isin(lectin_failed) #

        data_step4 = data_step3[SNR_filter]
        data_step4.reset_index(drop = True, inplace = True)

        ###########
        ###########Step 5: Centering and compute ratios
        ###########

        data_step5 = data_step4.iloc[:,[8,7,3,4]] #Generate a simpified dataframe
        data_step5_ratio = pd.DataFrame(data_step5["F532"] / data_step5["F635"])
        data_step5_ratio.columns = ["Ratio 532/635"]
        data_step5 = pd.concat([data_step5, data_step5_ratio], axis = 1)

        data_step5_grouped = data_step5.groupby(by = "Sample", sort = False)
        lectin_pass_count = data_step5_grouped.count()["Lectin"].tolist()

        calibrator_by_sample_532 = [] 
        calibrator_by_sample_635 = []
        def array_norm_1():
            global calibrator_by_sample_532, calibrator_by_sample_635
            calibrator_by_sample_532 = data_step5_grouped.median()["F532"].tolist()
            calibrator_by_sample_635 = data_step5_grouped.median()["F635"].tolist()

        def array_norm_2():
            global calibrator_by_sample_532, calibrator_by_sample_635
            calibrator_by_sample_532 = data_step5_grouped.mean()["F532"].tolist()
            calibrator_by_sample_635 = data_step5_grouped.mean()["F635"].tolist()

        def array_norm_3():
            global calibrator_by_sample_532, calibrator_by_sample_635
            calibrator_by_sample_532 = data_step5_grouped.median()["Ratio 532/635"].tolist()
            calibrator_by_sample_635 = data_step5_grouped.median()["Ratio 532/635"].tolist()

        def array_norm_4():
            global calibrator_by_sample_532, calibrator_by_sample_635
            calibrator_by_sample_532 = data_step5_grouped.mean()["Ratio 532/635"].tolist()
            calibrator_by_sample_635 = data_step5_grouped.mean()["Ratio 532/635"].tolist()

        def array_norm_5():
            global calibrator_by_sample_532, calibrator_by_sample_635
            calibrator_by_sample_532 = [1]*n_samples
            calibrator_by_sample_635 = [1]*n_samples

        def choose_array_norm_method(x):
            return {1:array_norm_1, 2:array_norm_2, 3:array_norm_3, 4:array_norm_4, 5:array_norm_5}.get(x).__call__()
            #when the function to be called in the dictionary has no parameters, one must: (1) NOT include "()" in the dictionary. For example, in this case, "array_norm_1()" in the dictionary will cause error. (2) include .__call__() in the end
        choose_array_norm_method(array_norm_method) #Choose function through a dictionary

        calibrator_by_sample_532_rep = []
        for i in range(n_samples):
            calibrator_by_sample_532_temp = [calibrator_by_sample_532[i]]*int(lectin_pass_count[i])
            calibrator_by_sample_532_rep.extend(calibrator_by_sample_532_temp)
        calibrator_by_sample_532_rep_df = pd.DataFrame(calibrator_by_sample_532_rep)
        calibrator_by_sample_532_rep_df.columns = ["Calibrator 532"]

        calibrator_by_sample_635_rep = []
        for i in range(n_samples):
            calibrator_by_sample_635_temp = [calibrator_by_sample_635[i]]*int(lectin_pass_count[i])
            calibrator_by_sample_635_rep.extend(calibrator_by_sample_635_temp)
        calibrator_by_sample_635_rep_df = pd.DataFrame(calibrator_by_sample_635_rep)
        calibrator_by_sample_635_rep_df.columns = ["Calibrator 635"]

        data_step5 = pd.concat([data_step5, calibrator_by_sample_532_rep_df, calibrator_by_sample_635_rep_df], axis = 1)

        if array_norm_method <=2:
            data_step5["Normalized Ratio 532/635"] = (data_step5["F532"] / data_step5["Calibrator 532"]) / (data_step5["F635"] / data_step5["Calibrator 635"])
        else:
            data_step5["Normalized Ratio 532/635"] = data_step5["Ratio 532/635"] / data_step5["Calibrator 532"]


        ###########
        ###########Step 6: Construct a final dataframe and write to file
        ###########
        data_step6_grouped = data_step5.groupby(by = ["Sample", "Lectin"], sort = False)
        data_step6 = data_step6_grouped.mean()["Normalized Ratio 532/635"]
        data_step6 = pd.DataFrame(np.array(data_step6).reshape((n_samples, len(lectin_passed))))
        data_step6_log2 = pd.DataFrame(np.log2(data_step6))

        data_step6_log2.columns = lectin_passed
        data_step6_log2.index = sample_list_all
        data_step6_log2 = data_step6_log2.transpose()
        data_step6_log2.to_csv(output_dir + "/data_final.txt", index = True, sep = "\t")
        
        tk.messagebox.showinfo(title = "Message", message = "Congratulations! The analysis is finished.\nPlease find the output files in the designated folder.")
    except:
        tk.messagebox.showwarning(title = "Error", message = "An error occured during the analysis process. However, output files may have been generated if the error is not fatal. In case of a fatal error, Please double check your inputs or contact the author for assistance.") 
