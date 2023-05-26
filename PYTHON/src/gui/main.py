import numpy as np
import os
import pandas as pd
import statsmodels.api as sm
import sys

from pymap3d.vincenty import vdist
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from tkinter import messagebox
from tkinter import *
from tkinter import filedialog
from tkinter import ttk

# Create the main window
window = Tk()
window.title("GapMET")

# Create a frame to hold the widgets
frame = Frame(window)
frame.pack()

# Creating frames
## Input Parameters frame
input_parameter_frame = LabelFrame(frame, text="Input Parameters")
input_parameter_frame.grid(row=0, column=0, padx=(20,0), pady=10)

## Input quality control frame
input_quality_cntrl_frame = LabelFrame(frame, text="Input Quality control parameters")
input_quality_cntrl_frame.grid(row=0, column=2, padx=(0,20), pady=(10,115))

## Input data frame
input_data_frame = LabelFrame(frame, text="Input Data")
input_data_frame.grid(row=0, column=1, padx=(10,10), pady=10)

## Output frame
output_frame = LabelFrame(frame, text="Output")
output_frame.grid(row=0, column=2, padx=(0,20), pady=(180,10))

# Creating Labels
## Input Parameter Labels
### Label and Combobox for selecting the Reference Series Type
ref_type_label = Label(input_parameter_frame, text="Reference Serie Type")
ref_type_combobox = ttk.Combobox(input_parameter_frame, values=["nearby", "external"])
ref_type_label.grid(row=0, column=0, pady=(20,0))
ref_type_combobox.grid(row=1, column=0)

#### Label and Combobox for selecting the Gap-filling Method
list_all_methods = ["RLS: Simple linear regression", "RLM: Multiple linear regression", "MPR: Regional weighting method", "MUK: UK traditional method", "IDM: Inverse distance method ", "MAS: Simple arithmetic mean"]
list_ext_methods = ["RLS: Simple linear regression", "MUK: UK traditional method"]
pree_meth_label = Label(input_parameter_frame, text="Gap-filling Method")
pree_meth_combobox = ttk.Combobox(input_parameter_frame, values=list_all_methods)
pree_meth_label.grid(row=3, column=0, pady=(20,0))
pree_meth_combobox.grid(row=4, column=0)

### Label and Entry for setting the Maximum Distance for Gap Filling
max_dis_label = Label(input_parameter_frame, text="Maximum Distance \nfor Gap Filling (km)")
max_dis_label.grid(row=5, column=0, pady=(20,0))
max_dis_entry = Entry(input_parameter_frame)
max_dis_entry.grid(row=6, column=0)

### Label and Entry for setting the Minimum Nearby Reference Series
min_est_label = Label(input_parameter_frame, text="Minimum Nearby \nReference Series")
min_est_label.grid(row=7, column=0, pady=(20,0))
min_est_entry = Entry(input_parameter_frame)
min_est_entry.grid(row=8, column=0, pady=(0,25))

## Input Data Labels
### Label, Entry, and Browse Button for selecting the CSV file with series coordinates/ID
input_coordinate_label = Label(input_data_frame, text="Input CSV File: \nseries coordinates/ID")
input_coordinate_label.grid(row=1, column=1, pady=(20, 0))

### Label, Entry, and Browse Button for selecting the CSV file with gap-filling series
input_data_label = Label(input_data_frame, text="Input CSV File: \ngap-filling series")
input_data_label.grid(row=3, column=1, pady=(20, 0))

### Label, Entry, and Browse Button for selecting the CSV file with external series
input_external_label = Label(input_data_frame, text="Input CSV File: \nexternal series")
input_external_label.grid(row=5, column=1, pady=(20,0))

### Input and empty label to ajust the layout
input_empty_label = Label(input_data_frame, text="")
input_empty_label.grid(row=8, column=1, pady=(20,12))

## Quality Control Labels
### Label and Entry for setting type of quality control
limit_type_label = Label(input_quality_cntrl_frame, text="Type of \nQuality control")
limit_type_combobox = ttk.Combobox(input_quality_cntrl_frame, values=["range", "percentile"])
limit_type_label.grid(row=0, column=0, padx=53, pady=(20,0))
limit_type_combobox.grid(row=1, column=0, padx=20)

### Label and Entry for setting interval of quality control
interval_label = Label(input_quality_cntrl_frame, text="Interval")
interval_label.grid(row=2, column=0, padx=0, pady=(20,0))

### Label and Combobox for setting max and min intervals
percentiles_interval = [int(i) for i in range(0, 101)]

interval_min_label = Label(input_quality_cntrl_frame, text="\nMin")
interval_min_combobox = ttk.Combobox(input_quality_cntrl_frame, values=[''], width=5)
interval_min_label.grid(row=2, column=0, padx=(0, 80), pady=(40,0),ipadx=0)
interval_min_combobox.grid(row=3, column=0, padx=(0, 80), pady=(0, 20), ipadx=0)

interval_max_label = Label(input_quality_cntrl_frame, text="\nMax")
interval_max_combobox = ttk.Combobox(input_quality_cntrl_frame, values=[''], width=5)
interval_max_label.grid(row=2, column=0, padx=(80, 0), pady=(40, 0), ipadx=0)
interval_max_combobox.grid(row=3, column=0, padx=(80, 0), pady=(0, 20), ipadx=0)

## Output Labels
output_path_label = Label(output_frame, text="Output Folder")
output_path_label.grid(row=0, column=0, padx=0, pady=5)

## Apply padding to the widgets inside the input_parameter_frame
for widget in input_parameter_frame.winfo_children():
    widget.grid_configure(padx=20, ipadx=22)

# Handle combobox and widgets:
def handle_combobox_selection(event):
    selected_ref_type = ref_type_combobox.get()
    selected_method = pree_meth_combobox.get()
    
    if selected_ref_type == "external":
        pree_meth_combobox.configure(state="readonly", values=list_ext_methods)
        if selected_method != "RLS: Simple linear regression" and selected_method != "MUK: UK traditional method":
            pree_meth_combobox.current(0)
            selected_method = pree_meth_combobox.get()
        min_est_entry.delete(0, 'end')
        min_est_entry.configure(state="disabled")
        max_dis_entry.delete(0, 'end') 
        max_dis_entry.configure(state="disabled")
        input_external_entry.delete(0, 'end')
        input_external_entry.configure(state="normal")
        input_external_browse_button.configure(state="normal")
    elif selected_ref_type == "nearby" and selected_method in ["RLS", "MUK"]:
        min_est_entry.configure(state="disabled")
        max_dis_entry.configure(state="normal")
        pree_meth_combobox.configure(state="readonly", values=list_all_methods)
    else:
        min_est_entry.configure(state="normal")
        max_dis_entry.configure(state="normal")
        pree_meth_combobox.configure(state="readonly", values=list_all_methods)
        input_external_browse_button.configure(state="disabled")
        input_external_entry.configure(state="disabled")

def handle_combobox_selection_quality_contrl(event):
    selected_limit_type = limit_type_combobox.get()
    selected_interval_min = interval_min_combobox.get()
    selected_interval_max = interval_max_combobox.get()
    
    if limit_type_combobox.get() == "percentile":
        interval_min_combobox.configure(state="normal", values=percentiles_interval)
        interval_max_combobox.configure(state="normal", values=percentiles_interval)    
    else:
        previous_limit_type = 0
        interval_min_combobox.configure(state="normal", values=["Min"])
        interval_max_combobox.configure(state="normal", values=["Max"])

## Bind the combobox selection events to the handle_combobox_selection_quality_contrl and handle_combobox_selection function
limit_type_combobox.bind("<<ComboboxSelected>>", handle_combobox_selection_quality_contrl)
interval_min_combobox.bind("<<ComboboxSelected>>", handle_combobox_selection_quality_contrl)
interval_max_combobox.bind("<<ComboboxSelected>>", handle_combobox_selection_quality_contrl)
ref_type_combobox.bind("<<ComboboxSelected>>", handle_combobox_selection)
pree_meth_combobox.bind("<<ComboboxSelected>>", handle_combobox_selection)

# browse csv file functions and buttons
df_coord =pd.DataFrame()
def browse_csv_file_coordinate():
    global df_coord
    file_path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
    if file_path:
        file_name = os.path.basename(file_path)
        input_coordinate_entry.delete(0, END)  # Clear previous entry
        input_coordinate_entry.insert(0, file_name)  # Set selected file name
        df_coord = pd.read_csv(file_path)
        print("Selected CSV file:", file_path)
input_coordinate_browse_button = Button(input_data_frame, text="Browse", command=browse_csv_file_coordinate)
input_coordinate_browse_button.grid(row=2, column=2, padx=(0, 20))  # Remove padx
input_coordinate_entry = Entry(input_data_frame)
input_coordinate_entry.grid(row=2, column=1, padx=(20, 0))

df_data =pd.DataFrame()
def browse_csv_file_data():
    global df_data
    file_path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
    if file_path:
        file_name = os.path.basename(file_path)
        input_data_entry.delete(0, END)  # Clear previous entry
        input_data_entry.insert(0, file_name)  # Set selected file name
        df_data = pd.read_csv(file_path)
        print("Selected CSV file:", file_path)
input_data_browse_button = Button(input_data_frame, text="Browse", command=browse_csv_file_data)
input_data_browse_button.grid(row=4, column=2, padx=(0, 20))  # Remove padx
input_data_entry = Entry(input_data_frame)
input_data_entry.grid(row=4, column=1, padx=(20, 0))

df_exter =pd.DataFrame()
def browse_csv_file_external():
    global df_exter
    file_path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
    if file_path:
        file_name = os.path.basename(file_path)
        input_external_entry.delete(0, END)  # Clear previous entry
        input_external_entry.insert(0, file_name)  # Set selected file name
        df_exter = pd.read_csv(file_path)
        print("Selected CSV file:", file_path)
input_external_browse_button = Button(input_data_frame, text="Browse", command=browse_csv_file_external)
input_external_browse_button.grid(row=6, column=2, padx=(0, 20))  # Remove padx
input_external_entry = Entry(input_data_frame)
input_external_entry.grid(row=6, column=1, padx=(20, 0))

def browse_output_folder():
    folder_path = filedialog.askdirectory()
    if folder_path:
        output_path_entry.delete(0, END)  # Clear previous entry
        output_path_entry.insert(0, folder_path)  # Set selected folder path
        print("Selected output folder:", folder_path)
output_path_browse_button = Button(output_frame, text="Browse", command=browse_output_folder)
output_path_browse_button.grid(row=0, column=1, padx=(0, 10))  # Remove padx
output_path_entry = Entry(output_frame)
output_path_entry.grid(row=0, column=0, padx=(10, 0))

# Class
class OutputRedirector(object):
    def __init__(self, text_widget, tag):
        self.text_widget = text_widget
        self.tag = tag

    def write(self, message):
        self.text_widget.insert(END, message, self.tag)

# Additional Functions
## Function to clear all fields and data
def clear_all():
    # Clear input fields and data
    ref_type_combobox.set('')
    pree_meth_combobox.set('')
    max_dis_entry.delete(0, 'end')
    min_est_entry.delete(0, 'end')
    input_coordinate_entry.delete(0, 'end')
    input_data_entry.delete(0, 'end')
    input_external_entry.delete(0, 'end')
    limit_type_combobox.set('')
    interval_min_combobox.set('')
    interval_max_combobox.set('')
    output_path_entry.delete(0, 'end')

    # Clear data variables
    global df_coord, df_data, df_exter
    df_coord = None
    df_data = None
    df_exter = None

## Validate Parameter functions:
def validate_ref_type(ref_type):
    if ref_type == "":
        ref_type = 1
    elif "nearby" in ref_type:
        ref_type = 1
    elif "external" in ref_type:
        ref_type = 2
    
    return ref_type

def validate_pree_meth(pree_meth):
    if len(pree_meth) != 3 or pree_meth not in ["RLS", "RLM", "MPR", "MUK", "IDM", "MAS"]:
        raise ValueError("'pree_meth' must be a string containing one of the gap-filling methods: (RLS, RLM, MPR, MUK, IDM, or MAS)")

def validate_min_est(min_est):
    if min_est == "":
        min_est = 3
    elif not isinstance(min_est, int) and any(val < 0 for val in min_est):
        raise ValueError('"min_est" must be a positive integer')

def validate_max_dis(max_dis):
    if max_dis == "":
        max_dis = 100
    elif not isinstance(max_dis, int) and any(val < 0 for val in max_dis):
        raise ValueError('"max_dis" must be a positive integer')

def validate_limit_type(limit_type):
    if limit_type == "":
        limit_type = "range"
    elif "range" not in limit_type:
        if "percentile" not in limit_type:
            raise ValueError('"limit_type" must be "range" or "percentile"')

def validate_interval(interval, limit_type, dt_data):
    if "percentile" in limit_type:
        if interval == "":
            interval = [5, 95]
        elif not isinstance(interval, list) or len(interval) != 2 or not all(isinstance(x, (int, float)) for x in interval):
            raise ValueError('The "percentile" "interval" must be an array with min and max percentile')
        elif interval[0] >= interval[1]:
            raise ValueError('The "percentile" min interval must be lower than max interval')
        elif any(x < 0 or x > 100 for x in interval):
            raise ValueError('The "percentile" min and max interval must be between 0 and 100')
    elif "range" in limit_type:
        if interval == "":
            interval = [np.nanmin(dt_data.iloc[:, 4:].values), np.nanmax(dt_data.iloc[:, 4:].values)]
        if interval[0] == -9999:
             interval[0] = np.nanmin(dt_data.iloc[:, 4:].values)
        if interval[1] == -9999:
             interval[1] = np.nanmax(dt_data.iloc[:, 4:].values)        
        elif not isinstance(interval, list) or len(interval) != 2 or not all(isinstance(x, (int, float)) for x in interval):
            raise ValueError('The "range" "interval" must be an array with min and max values')
        elif interval[0] >= interval[1]:
            raise ValueError('The min interval must be lower than max interval')

## Quality control gap-filled data functions:
def LimCheck(dt_ori, gap, gap_filled, limit_type, interval):
    gap_filled_checked = gap_filled.copy()

    if limit_type != "range" and limit_type != "percentile":
        raise ValueError('"limit_type" must be "range" or "percentile"')

    if limit_type == "percentile":
        if interval is None:
            interval = [5, 95]
        elif not isinstance(interval, list) or len(interval) != 2:
            raise ValueError('The "percentile" "interval" must be an array with min and max percentile')
        elif interval[0] >= interval[1]:
            raise ValueError('The "percentile" min interval must be lower than max interval')
        elif any(x < 0 or x > 100 for x in interval):
            raise ValueError('The "percentile" min and max interval must be between 0 and 100')

        min_lim = np.nanpercentile(dt_ori.values.flatten(), interval[0])
        max_lim = np.nanpercentile(dt_ori.values.flatten(), interval[1])

    elif limit_type == "range":
        if interval is None:
            interval = [dt_ori.min().min(), dt_ori.max().max()]
        elif not isinstance(interval, list) or len(interval) != 2:
            raise ValueError('The "range" "interval" must be an array with min and max values')
        elif interval[0] >= interval[1]:
            raise ValueError('The min interval must be lower than max interval')

        min_lim = interval[0]
        max_lim = interval[1]

    for line_check in range(len(dt_ori)):
        if pd.isna(gap[line_check]):
            if gap_filled[line_check] <= min_lim or gap_filled[line_check] >= max_lim:
                gap_filled_checked[line_check] = np.nan

    return gap_filled_checked

## gap-filling functions:
def RLS(dt_ori, limit_type, interval, ref_type, index_est, dt_ext):
        
    dt = dt_ori.copy()
    dt_flag = np.zeros((dt_ori.shape[0], dt_ori.shape[1]))
    #dt_flag = pd.DataFrame(np.zeros((dt_ori.shape[0], dt_ori.shape[1])))
    coef_lin = np.full((dt_ori.shape[1], dt_ori.shape[1]), np.nan)
    coef_ang = np.full((dt_ori.shape[1], dt_ori.shape[1]), np.nan)
    R2 = np.full((dt_ori.shape[1], dt_ori.shape[1]), np.nan)

    # RLS: Using nearby stations to gap-filling (ref_type==1)
    if ref_type == 1:
        for p in range(dt_ori.shape[1]):
            kf = index_est.iloc[p]
            kf = kf.dropna()
            for k in range(len(kf)):
                df = pd.DataFrame()
                df[0] = dt.iloc[:, p]
                df[1] = dt.iloc[:, int(index_est.iloc[p, k])]
                df = df.dropna(subset=[0, 1])
                gap = df[0]
                ref = df[1]
                ref_array = ref.values.reshape(-1, 1)
                gap_array = gap.values.reshape(-1, 1)
                mdl = LinearRegression().fit(ref_array, gap_array)
                coef_lin[p][k] = mdl.intercept_[0]
                coef_ang[p][k] = mdl.coef_[0][0]
                R2[p][k] = mdl.score(ref_array, gap_array)

        # Gap-filling
        numNans_fim = np.count_nonzero(np.isnan(dt))
        iter = 1

        if numNans_fim > 0:
            numNans_inicio = np.count_nonzero(np.isnan(dt))
            numNans_fim = 0
            while numNans_fim < numNans_inicio:
                numNans_inicio = np.count_nonzero(np.isnan(dt))
                if numNans_inicio == 0:
                    break
                print(f'Running Iteration {iter}. Number of Gaps {numNans_inicio}')
                for p in range(dt_ori.shape[1]):
                    print(f'Gapfilling station {p+1} of {dt_ori.shape[1]}')
                    kf = index_est.iloc[p]
                    kf = kf.dropna()
                    for k in range(len(kf)):
                        if iter == 1:
                            ref = dt_ori.iloc[:, int(index_est.iloc[p, k])]
                        else:
                            ref = dt.iloc[:, int(index_est.iloc[p, k])]
                        gap = dt.iloc[:, p]
                        gap_filled = gap.copy()
                        for i in range(len(gap)-1):
                            if np.isnan(gap[i]):
                                y = coef_lin[p, k] + (coef_ang[p, k] * ref[i])
                                gap_filled[i] = y
                                dt_flag[i, p] = iter
                        gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
                        dt.iloc[:, p] = gap_filled_checked
                numNans_fim = np.count_nonzero(np.isnan(dt))
                iter += 1

    # RLS: Using external time series (ref_type==2)
    elif ref_type == 2:
        dt = dt_ori.copy()
        dt = pd.concat([dt, dt_ext], axis=1)
        dt_flag = np.zeros((dt.shape[0], dt.shape[1]))
        #dt_flag = pd.DataFrame(np.zeros((dt_ori.shape[0], dt_ori.shape[1])))
        coef_lin = np.full((dt.shape[1], dt.shape[1]), np.nan)
        coef_ang = np.full((dt.shape[1], dt.shape[1]), np.nan)
        R2 = np.full((dt.shape[1], dt.shape[1]), np.nan)

        for p in range(dt_ori.shape[1]):
            for k in range(dt_ext.shape[1]):
                df = pd.DataFrame()
                df[0] = dt.iloc[:, p]
                df[1] = dt.iloc[:, dt_ori.shape[1] + k]
                df = df.dropna(subset=[0, 1])
                gap = df[0]
                ref = df[1]
                ref_array = ref.values.reshape(-1, 1)
                gap_array = gap.values.reshape(-1, 1)
                mdl = LinearRegression().fit(ref_array, gap_array)
                coef_lin[p][k] = mdl.intercept_[0]
                coef_ang[p][k] = mdl.coef_[0][0]
                R2[p][k] = mdl.score(ref_array, gap_array)

        # Gap-filling
        numNans_fim = np.count_nonzero(np.isnan(dt))
        iter = 1

        if numNans_fim > 0:
            numNans_inicio = np.count_nonzero(np.isnan(dt))
            numNans_fim = 0
            while numNans_fim < numNans_inicio:
                numNans_inicio = np.count_nonzero(np.isnan(dt))
                if numNans_inicio == 0:
                    break
                print(f'Running Iteration {iter}. Number of Gaps {numNans_inicio}')
                for p in range(dt_ori.shape[1]):
                    print(f'Gapfilling station {p+1} of {dt_ori.shape[1]}')
                    for k in range(dt_ext.shape[1]):
                        if iter == 1:
                            ref = dt_ext.iloc[:, k]
                        else:
                            ref = dt.iloc[:, dt_ori.shape[1] + k]
                        gap = dt.iloc[:, p]
                        gap_filled = gap.copy()
                        for i in range(len(gap)):
                            if np.isnan(gap[i]):
                                y = coef_lin[p, k] + (coef_ang[p, k] * ref[i])
                                gap_filled[i] = y
                                dt_flag[i, p] = iter
                        gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
                        dt.iloc[:, p] = gap_filled_checked
                numNans_fim = np.count_nonzero(np.isnan(dt))
                iter += 1

    # Finishing gap-filled dataset
    dt_pree = dt.copy()
    dt_flag[np.isnan(dt_pree)] = np.nan
    #dt_flag = dt_flag.mask(pd.isnull(dt_pree))

    if numNans_fim != 0:
        print(f'{numNans_fim} Gaps left on the dataset')

    return dt_pree, dt_flag

def RLM(dt_ori, limit_type, interval, min_est, index_est):

    dt = np.copy(dt_ori)
    dt_flag = np.zeros((dt_ori.shape[0], dt_ori.shape[1]))

    numNans_inicio = np.sum(np.isnan(dt))
    numNans_fim = 0
    iter = 1

    while numNans_fim < numNans_inicio:
        numNans_inicio = np.sum(np.isnan(dt))
        if numNans_inicio == 0:
            break

        print(f'Running Iteration {iter}. Number of Gaps {numNans_inicio}')
        for p in range(dt_ori.shape[1]):
            print(f'Gapfilling station {p + 1} of {dt_ori.shape[1]}')
            gap = dt[:, p]
            gap_filled = np.copy(gap)
            kf_all = index_est.iloc[p]
            kf = kf_all.dropna()
            X = pd.DataFrame(np.full((dt_ori.shape[0], len(kf)), np.nan))
            for j in range(len(kf)):
                if iter == 1:
                    X.iloc[:,j] = dt_ori.iloc[:, kf[j]]
                else:
                    X.iloc[:,j] = dt[:, kf[j]]
            for i in range(len(X)):
                if np.isnan(dt[i, p]):
                    X1 = X.iloc[i, :]
                    ind = np.where(~np.isnan(X1))[0]
                    if len(ind) >= min_est:
                        X0 = X
                        X0['y'] = gap
                        X0 = X0.dropna()
                        X0 = sm.add_constant(X0)
                        X2 = X0.iloc[:,0:-1]
                        y = X0['y']
                        model = sm.OLS(y, X2)
                        results = model.fit()
                        coefficients = results.params[1:]
                        constant = results.params['const']
                        gap_filled[i] = np.nansum(X1 * coefficients) + constant
                        dt_flag[i, p] = iter

            gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
            dt[:, p] = gap_filled_checked

        numNans_fim = np.sum(np.isnan(dt))
        iter += 1

    dt_pree = np.copy(dt)
    dt_flag[np.isnan(dt_pree)] = np.nan

    if numNans_fim != 0:
        print(f'{numNans_fim} Gaps left on the dataset')

    return dt_pree, dt_flag

def MUK(dt_ori, dt_mean_month, dt_months, limit_type, interval, ref_type, index_est, dt_ext, dt_mean_ext):

    dt = dt_ori
    dt_flag = np.zeros((dt_ori.shape[0], dt_ori.shape[1]))

    numNans_inicio = np.count_nonzero(np.isnan(dt))
    numNans_fim = 0
    iter = 1

    if ref_type == 1:
        while numNans_fim < numNans_inicio:
            numNans_inicio = np.count_nonzero(np.isnan(dt))
            if numNans_inicio == 0:
                break
            print(f'Running Iteration {iter}. Number of Gaps {numNans_inicio}')
            for p in range(dt_ori.shape[1]):
                print(f'Gapfilling station {p+1} of {dt_ori.shape[1]}')
                kf = index_est.iloc[p]
                kf = kf.dropna()
                for k in range(len(kf)):
                    gap = dt.iloc[:, p]
                    gap_filled = gap.copy()
                    if iter == 1:
                        ref = dt_ori.iloc[:, kf[k]]
                    else:
                        ref = dt.iloc[:, kf[k]]
                    for i in range(len(gap)):
                        if np.isnan(gap[i]):
                            gap_filled[i] = ref[i] + (dt_mean_month.iloc[dt_months[i]-1, p] - dt_mean_month.iloc[dt_months[i]-1, kf[k]])
                            dt_flag[i, p] = iter
                    gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
                    dt.iloc[:, p] = gap_filled_checked
            numNans_fim = np.count_nonzero(np.isnan(dt))
            iter += 1

    if ref_type == 2:
        while numNans_fim < numNans_inicio:
            numNans_inicio = np.count_nonzero(np.isnan(dt))
            if numNans_inicio == 0:
                break
            for p in range(dt_ori.shape[1]):
                print(f'Gapfilling station {p+1} of {dt_ori.shape[1]}')
                gap = dt.iloc[:, p]
                gap_filled = gap.copy()
                ref = dt_ext.iloc[:, p]
                for i in range(len(gap)):
                    if np.isnan(gap[i]):
                        gap_filled[i] = ref[i] + (dt_mean_month.iloc[dt_months[i]-1, p] - dt_mean_ext.iloc[dt_months[i]-1, p])
                        dt_flag[i, p] = iter
                gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
                dt.iloc[:, p] = gap_filled_checked
            numNans_fim = np.count_nonzero(np.isnan(dt))

    dt_pree = dt.copy()
    dt_flag[np.isnan(dt_pree)] = np.nan

    if numNans_fim != 0:
        print(f'{numNans_fim} Gaps left on the dataset')
    
    return dt_pree, dt_flag

def MPR(dt_ori, dt_mean_month, dt_months, limit_type, interval, min_est, index_est):
    
    dt = dt_ori
    dt_flag = np.zeros((np.size(dt_ori, 0), np.size(dt_ori, 1)))
    
    # Gap filling
    numNans_inicio = np.count_nonzero(np.isnan(dt))
    numNans_fim = 0
    iter = 1
    
    while numNans_fim < numNans_inicio:
        numNans_inicio = np.count_nonzero(np.isnan(dt))
        if numNans_inicio == 0:
            break
        print('Running Iteration', iter, '. Number of Gaps', numNans_inicio)
        
        for p in range(np.size(dt_ori, 1)):
            print('Gapfilling station', p + 1, 'of', np.size(dt_ori, 1))
            gap = dt.iloc[:, p]
            gap_filled = np.copy(gap)
            kf_all = index_est.iloc[p]
            kf = kf_all.dropna()
            X = np.full((dt_ori.shape[0],len(kf)), np.nan)           
            for j in range(len(kf)):
                if iter == 1:
                    X[:,j] = dt_ori.iloc[:, kf[j]]
                else:
                    X[:,j] = dt.iloc[:, kf[j]]            
            for i in range(np.size(X, 0)):
                if np.isnan(dt.iloc[i, p]):
                    X1 = X[i, :]
                    ind = np.where(~np.isnan(X1))[0]
                    if len(ind) >= min_est:
                        m_ref = dt_months[i]
                        Xm = dt_mean_month.iloc[m_ref - 1,:].values
                        Xm = Xm[~pd.isna(kf_all)].astype(float)
                        Xm = Xm[~np.isnan(X1)]
                        X2 = X1[~np.isnan(X1)]
                        gap_filled[i] = (1 / len(Xm)) * (np.sum(X2 / Xm)) * dt_mean_month.iloc[m_ref - 1, p]
                        dt_flag[i, p] = iter
                        del X1, X2, Xm            
            gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
            dt.iloc[:, p] = gap_filled_checked
        
        numNans_fim = np.count_nonzero(np.isnan(dt))
        iter += 1
    
    dt_pree = np.copy(dt)
    dt_flag[np.isnan(dt_pree)] = np.nan
    
    if numNans_fim != 0:
        print(numNans_fim, 'Gaps left in the dataset')
    
    return dt_pree, dt_flag

def MAS(dt_ori, limit_type, interval, min_est, index_est):

    dt = dt_ori
    dt_flag = np.zeros((dt_ori.shape[0], dt_ori.shape[1]))

    numNans_inicio = np.count_nonzero(np.isnan(dt))
    numNans_fim = 0
    iter = 1

    while numNans_fim < numNans_inicio:
        numNans_inicio = np.count_nonzero(np.isnan(dt))
        if numNans_inicio == 0:
            break
        print(f'Running Iteration {iter}. Number of Gaps {numNans_inicio}')
        for p in range(dt_ori.shape[1]):
            print(f'Gapfilling station {p+1} of {dt_ori.shape[1]}')
            gap = dt.iloc[:, p]
            gap_filled = gap.copy()
            kf = index_est.iloc[p]
            kf = kf.dropna()
            X = np.full((dt_ori.shape[0],len(kf)), np.nan)
            for j in range(len(kf)):  
                if iter == 1:
                    X[:,j] = dt_ori.iloc[:, kf[j]]
                else:
                    X[:,j] = dt.iloc[:, kf[j]]
                for i in range(len(gap)):
                    if np.isnan(dt.iloc[i, p]):
                        X1 = X[i]
                        ind = np.where(~np.isnan(X1))[0]
                        if len(ind) >= min_est:
                            if np.isnan(X1).all():
                                gap_filled[i] = np.nan
                            else:
                                gap_filled[i] = np.nanmean(X1)
                            dt_flag[i, p] = iter
                            del X1
                gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
                dt.iloc[:, p] = gap_filled_checked
        numNans_fim = np.count_nonzero(np.isnan(dt))
        iter += 1

    dt_pree = dt.copy()
    dt_flag[np.isnan(dt_pree)] = np.nan

    if numNans_fim != 0:
        print(f'{numNans_fim} Gaps left on the dataset')
    return dt_pree, dt_flag

def IDM(dt_ori, dt_dist, limit_type='range', interval=None, min_est=3, index_est=None):

    dist = np.array(dt_dist)
    dt = np.array(dt_ori)
    dt_flag = np.zeros((np.size(dt_ori, 0), np.size(dt_ori, 1)))

    # Gap filling
    num_nans_inicio = np.count_nonzero(np.isnan(dt))
    num_nans_fim = 0
    iter = 1

    while num_nans_fim < num_nans_inicio:
        num_nans_inicio = np.count_nonzero(np.isnan(dt))
        if num_nans_inicio == 0:
            break
        print(f'Running Iteration {iter}. Number of Gaps {num_nans_inicio}')
        for p in range(np.size(dt_ori, 1)):
            print(f'Gapfilling station {p + 1} of {np.size(dt_ori, 1)}')
            gap = dt[:, p]
            gap_filled = gap.copy()
            kf_all = index_est.iloc[p]
            kf = kf_all.dropna()
            X = np.full((dt_ori.shape[0],len(kf)), np.nan)   
            for j in range(len(kf)):
                if iter == 1:
                    X[:,j] = dt_ori.iloc[:, kf[j]]
                else:
                    X[:,j] = dt[:, kf[j]]                
            for i in range(len(X)):
                if np.isnan(dt[i, p]):
                    X1 = X[i]
                    ind = np.where(~np.isnan(X1))[0]
                    if len(ind) >= min_est:
                        Xd = np.full((len(kf)), np.nan)  
                        for n in range(len(X1)):
                            Xd[n] = dist[p, kf[n]]
                        X2 = X1[~np.isnan(X1)]
                        Xd = Xd[~np.isnan(X1)]
                        gap_filled[i] = np.sum(X2 / Xd) / np.sum(1 / Xd)
                        dt_flag[i, p] = iter
            gap_filled_checked = LimCheck(dt_ori, gap, gap_filled, limit_type, interval)
            dt[:, p] = gap_filled_checked
        num_nans_fim = np.count_nonzero(np.isnan(dt))
        iter += 1
    
    # Finishing gap-filled dataset
    dt_pree = dt
    dt_flag[np.isnan(dt_pree)] = np.nan

    if num_nans_fim != 0:
        print(f'{num_nans_fim} Gaps left on the dataset')

    return dt_pree, dt_flag

## GapMET functions:
def GapMet(dt_coord, dt_data, pree_meth, min_est, max_dis, limit_type, interval, ref_type, dt_extern):
    
    # SECTION 1: CHECK PARAMETERS
    if 'ref_type' not in locals():
        ref_type = 1
    else: 
        ref_type = validate_ref_type(ref_type)
    
    if 'min_est' not in locals():
        min_est = 3
    else:
        validate_min_est(min_est)

    validate_pree_meth(pree_meth)
    
    if 'max_dis' not in locals():
        max_dis = 100
    else:
        validate_max_dis(max_dis)

    if 'limit_type' not in locals():
        limit_type = "range"
    else:
        validate_limit_type(limit_type)

    if 'interval' not in locals():
        if "percentile" in limit_type:
            interval = [5, 95]
        elif "range" in limit_type:
            interval = [np.nanmin(dt_data.iloc[:, 4:].values), np.nanmax(dt_data.iloc[:, 4:].values)]
    else:
        validate_interval(interval, limit_type, dt_data)
    
    # SECTION 2: CHECK DATASETS
    if dt_data.empty:
            raise ValueError("The gapfilling data serie was not select. Select the CSV file dataset to be gap-filled")     
    if dt_data.shape[1] < (min_est + 4):
        raise ValueError("The number of columns (n) in 'dt_data' must be equal to or higher than n = min_est + 4")

    if dt_coord.shape[1] != 3:
        raise ValueError("The number of columns (m) in 'dt_coord' must be equal to 3 (station_id, latitude, longitude)")
    if dt_coord.shape[0] != (dt_data.shape[1] - 3):
        raise ValueError("The number of stations (rows) in 'dt_coord' is different than the number of stations (columns) in 'dt_data'")
    if dt_coord.empty:
            raise ValueError("The gapfilling coordinates/ID file was not select. Select the CSV file with the time series ID and coordinates")  
    
    if ref_type == 2:
        if df_exter.empty:
            raise ValueError("The external data serie was not select. Select an external data series CSV file")    
        if dt_data.shape[1] != dt_extern.shape[1]:
            raise ValueError("The number of stations (columns) in the external dataset is different than the number of stations (columns) in 'dt_data'")
        if dt_data.shape[0] != dt_extern.shape[0]:
            raise ValueError("The timeseries size (time steps) in the external dataset do not match the timeseries size (time steps) in 'dt_data'")

    dt_ori = dt_data.iloc[:, 3:]
    if ref_type == 2:
        dt_ext = dt_extern.iloc[:, 3:]

    dt_cod = dt_coord.iloc[:, 0].values
    dt_lat = dt_coord.iloc[:, 1].values
    dt_lon = dt_coord.iloc[:, 2].values

    if np.any(dt_lat > 90) and np.any(dt_lat < -90) and np.any(dt_lon > 180) and np.any(dt_lon < -180):
        lat_inf = np.unique(dt_cod[dt_lat < -90])
        lat_sup = np.unique(dt_cod[dt_lat > 90])
        lon_inf = np.unique(dt_cod[dt_lon < -180])
        lon_sup = np.unique(dt_cod[dt_lon > 180])
        est_err = np.unique(np.concatenate((lat_inf, lat_sup, lon_inf, lon_sup)))

        if est_err.size > 0:
            est_err_str = ",".join(est_err.astype(str))
            raise ValueError(f"The geographic coordinates in the station {est_err_str} must be between -90<lat<90 and -180<lon<180")

    # SECTION 3: CALCULATE THE DISTANCE BETWEEN DATA SERIES
    # Calculate the distances 
    dt_dist = np.zeros((len(dt_lat), len(dt_lat)))
    for i in range(len(dt_lat)):
        for j in range(len(dt_lat)):
            if i != j :
                lat1 = dt_lat[i]
                lon1 = dt_lon[i]
                lat2 = dt_lat[j]
                lon2 = dt_lon[j]
                dist_km, _ = vdist(lat1, lon1, lat2, lon2)   # Calculate distance in kilometers
                dt_dist[i, j] = dist_km/1000
            else:
                dt_dist[i, j] = 0         
    dt_dist = pd.DataFrame(dt_dist, index=dt_cod, columns=dt_cod)

    # Validate the distances with seted max_dis and min_est parameters
    if ref_type == 1:
        dt = dt_dist.copy()
        dt[dt == 0] = np.nan

        min_dis = np.round(np.nanmin(dt))
        if max_dis < min_dis:
            raise ValueError(f"The shortest distance to all stations has at least 1 reference station within the 'max_dis' is {min_dis} km")

        dt1 = dt.copy()
        dt1[dt1 <= max_dis] = 1
        dt1[dt1 > max_dis] = 0
        min_est_in_max_dis = np.nanmin(np.nansum(dt1, axis=1))
        dt1 = np.sort(dt, axis=1)
        max_dis_in_min_est = np.round(np.nanmax(dt1[:, min_est]))

        if max_dis < max_dis_in_min_est and min_est > min_est_in_max_dis:
            raise ValueError("The number of reference stations within the 'max_dis' radius is lower than set in 'min_est'. Decrease the 'min_est' parameter to "
                            f"{min_est_in_max_dis} stations or increase 'max_dis' to {max_dis_in_min_est} km")

    # SECTION 4: CALCULATE THE CORRELATION BETWEEN DATA SERIES
    if ref_type == 1:
        dt_corr = np.full((len(dt_cod), len(dt_cod)), np.nan)
        
        for i in range(dt_ori.shape[1]):
            for j in range(dt_ori.shape[1]):
                if i != j :
                    ponts = pd.DataFrame()
                    ponts[0] = dt_ori.iloc[:, i]
                    ponts[1] = dt_ori.iloc[:, j]
                    ponts = ponts.dropna(subset=[0, 1])
                    rho, pval = pearsonr(ponts[0], ponts[1])
                    # Check if the correlations are significant at 0.05
                    if pval <= 0.05:
                        dt_corr[i, j] = rho
                    else:
                        dt_corr[i, j] = np.nan
                else:
                    dt_corr[i, j] = np.nan        
        # Create correlations matrix table
        dt_corr = pd.DataFrame(dt_corr, index=dt_cod, columns=dt_cod)

    elif ref_type == 2:
        dt_corr = np.full((len(dt_cod), 1), np.nan)

        for i in range(dt_ori.shape[1]):
            ponts = pd.DataFrame()
            ponts[0] = dt_ori.iloc[:, i]
            ponts[1] = dt_ext.iloc[:, j]
            ponts = ponts.dropna(subset=[0, 1])
            rho, pval = pearsonr(ponts[0], ponts[1])
            # Check if the correlation is significant at 0.05
            if pval <= 0.05:
                dt_corr[i, 0] = rho
            else:
                dt_corr[i, 0] = np.nan
        # Create correlations matrix table
        dt_corr = pd.DataFrame(dt_corr, index=dt_cod)

    # SECTION 5: DEFINE REFERENCE SERIES BASED ON THE DISTANCE AND CORRELATION
    if ref_type == 1:
        dt = dt_dist.replace(0, np.nan).to_numpy()
        dt1 = dt_corr.replace(np.nan, 0).to_numpy()
        dt1[dt1 == 1] = 0
        dt1[dt > max_dis] = 0
        dt1 = np.abs(dt1)
        M = np.sort(dt1, axis=1)[:, ::-1]
        I = np.argsort(dt1, axis=1)[:, ::-1]
        I = pd.DataFrame(I, dtype=object)
        M = pd.DataFrame(M, dtype=object)
        I[M == 0] = np.nan
        index_est = I

    if pree_meth in ("MPR", "MUK"):
        dt_months = dt_data.iloc[:, 1].to_numpy()
        table_mes=dt_data.iloc[:, 1:]
        table_mes.drop(columns=['day'], inplace=True)
        dt_mean_month = table_mes.groupby(['month']).mean().to_numpy()
        dt_mean_month = pd.DataFrame(dt_mean_month, index=list(range(1, 13)), columns=dt_cod)
        
        if ref_type == 2:
            table_mes=dt_extern.iloc[:, 1:]
            table_mes.drop(columns=['day'], inplace=True)
            dt_mean_month_ext = table_mes.groupby(['month']).mean().to_numpy()
            dt_mean_ext = pd.DataFrame(dt_mean_month_ext , index=list(range(1, 13)), columns=dt_cod)

    # SECTION 6: GAP-FILLING
    # RLS: Simple linear regression
    if pree_meth == "RLS":
        if ref_type == 1:
            dt_pree, dt_flag = RLS(dt_ori, limit_type, interval, ref_type, index_est, [])
        else:
            dt_pree, dt_flag = RLS(dt_ori, limit_type, interval, ref_type, [], dt_ext)

    # RLM: Multiple linear regression
    if pree_meth == "RLM":
        if ref_type == 1:
            dt_pree, dt_flag = RLM(dt_ori, limit_type, interval, min_est, index_est)
        else:
            raise ValueError('This methodology (RLM) can only be used with nearby reference stations (Ref_type="nearby")')
    
    # MUK: UK traditional method
    if pree_meth == "MUK":
        if ref_type == 1:
            dt_pree, dt_flag = MUK(dt_ori, dt_mean_month, dt_months, limit_type, interval, ref_type, index_est, '','')
        else:
            dt_pree, dt_flag = MUK(dt_ori, dt_mean_month, dt_months, limit_type, interval, ref_type, '', dt_ext, dt_mean_ext)

    # MAS: Simple arithmetic mean   
    if pree_meth == "MAS":
        if ref_type == 1:
            dt_pree, dt_flag = MAS(dt_ori, limit_type, interval, min_est, index_est)
        else:
            raise ValueError("This methodology (MAS) can only be used with nearby reference stations (Ref_type='nearby')")

    #  MPR: Regional weighting method
    if pree_meth == "MPR":
        if ref_type == 1:
            dt_pree, dt_flag = MPR(dt_ori, dt_mean_month, dt_months, limit_type, interval, min_est, index_est)
        else:
            raise ValueError('This methodology (MPR) can only be used with nearby reference stations (Ref_type="nearby")')

    # IDM: Inverse distance method
    if pree_meth == "IDM":
        if ref_type == 1:
            dt_pree, dt_flag = IDM(dt_ori, dt_dist, limit_type, interval, min_est, index_est)
        else:
            raise ValueError('This methodology (IDM) can only be used with nearby reference stations (Ref_type="nearby")')

    return dt_pree, dt_flag, dt_dist, dt_corr

## run GapMET function
def run_gapmet():

    # Create a new window for displaying the output
    output_window = Toplevel(window)
    output_window.title("Output")
    output_window.geometry("500x300")

    # Create a Text widget to display the printed output and raised errors
    output_text = Text(output_window, width=60, height=20)
    output_text.pack()

    # Configure the Text widget tags for different styles
    output_text.tag_configure("normal", foreground="black", font=("Times", 10))
    output_text.tag_configure("error", foreground="red", font=("Arial", 12, "bold"))

    # Redirect the standard output to the Text widget
    sys.stdout = OutputRedirector(output_text, "normal")
    sys.stderr = OutputRedirector(output_text, "error")

    # Validate parameters
    try:
        ref_type = ref_type_combobox.get()
        list_methods = ["RLS","RLM","MPR","MUK","IDM","MAS"]
        index = list_all_methods.index(pree_meth_combobox.get())
        pree_meth=list_methods[index]

        max_dis = int(max_dis_entry.get()) if max_dis_entry['state'] != 'disabled' and max_dis_entry.get() != '' else 100
        min_est = int(min_est_entry.get()) if min_est_entry['state'] != 'disabled' and min_est_entry.get() != '' else 3

        if df_coord.empty:
            raise ValueError("The gapfilling coordinates/ID file was not select. Select the CSV file with the time series ID and coordinates")  
        if df_data.empty:
            raise ValueError("The gapfilling data serie was not select. Select the CSV file dataset to be gap-filled") 
        dt_coord = df_coord
        dt_data = df_data

        dt_extern = df_exter
        limit_type = limit_type_combobox.get()
    
        if interval_min_combobox.get() == 'Min' or interval_min_combobox.get() == 'min' or interval_min_combobox.get() == '':
            min_interval = np.nanmin(dt_data.iloc[:, 4:].values)
        else:
            min_interval = float(interval_min_combobox.get())

        if interval_max_combobox.get() == 'Max' or interval_max_combobox.get() == 'max' or interval_max_combobox.get() == '':
            max_interval = np.nanmax(dt_data.iloc[:, 4:].values)
        else:
            max_interval = float(interval_max_combobox.get())
        
        interval = [min_interval, max_interval]
        output_folder = output_path_entry.get()

        # Run the GapMet algorithm with the specified parameters
        dt_pree, dt_flag, dt_dist, dt_corr = GapMet(dt_coord, dt_data, pree_meth, min_est, max_dis, limit_type, interval, ref_type, dt_extern)

        # Gapfilled Dataset
        dt_pree = dt_pree.round(3)
        dt_fill = dt_data.copy()
        dt_fill.iloc[:, 3:] = dt_pree

        # Gapfilling Flags
        dt_flag = pd.DataFrame(dt_flag)
        dt_flags = dt_data.copy()
        dt_flags.iloc[:, 3:] = dt_flag

        # Save the gap-filled data, flags, distance, and correlation to CSV files
        dt_fill.to_csv(os.path.join(output_folder, 'gap_filled_data.csv'), index=False, na_rep=np.nan)
        dt_flags.to_csv(os.path.join(output_folder, 'gap_flags.csv'), index=False, na_rep=np.nan)
        dt_dist.to_csv(os.path.join(output_folder, 'gap_distance.csv'), index=True)
        dt_corr.to_csv(os.path.join(output_folder, 'gap_correlation.csv'), index=True)

        print("Gap-filling completed.")
        print("Results saved to:", output_folder)
    
    except Exception as e:
        print("Error occurred during gap-filling:")
        print(e)

    # Reset the standard output
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

# Run and Clear Buttons
gapfilling_button = Button(frame, text="Gap-filling", command=run_gapmet)
gapfilling_button.grid(row=0, column=2, padx=(0,80), pady=(285,10),ipadx=30)

clear_button = Button(frame, text="Clear", command=clear_all)
clear_button.grid(row=0, column=2, padx=(120,0), pady=(285,10),ipadx=1)

window.mainloop()