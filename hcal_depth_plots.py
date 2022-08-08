#!/usr/bin/env python
import matplotlib

matplotlib.use("Agg")

#import ROOT as ROOT
import uproot as up
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np
import os

def search_files(directory = None):
    assert os.path.isdir(directory)
    file_list = []
    for cur_path, directories, files in os.walk(directory):
        for file in files:
            if ".root" in file:
                file_list.append(os.path.join(cur_path, file))
    return file_list

def merge_arrays(arr_list, branch):
    arrs = []
    for file in arr_list:
        arrs.append(up.open(file)[branch].array())
    return ak.concatenate(arrs, axis=0)



ct900_dir = "crab_hcal_testLLP-MH125_MFF_12_CT900_to4mu_def/results"
ct10000_dir ="crab_hcal_testLLP-MH350_MFF160_CT10000_with_new_MC_def/results"
nu_dir = "crab_hcal_testLLP-Nu_Pt20_with_reco_def/results"

ct900_files = search_files(ct900_dir)
ct10000_files = search_files(ct10000_dir)
nu_files = search_files(nu_dir)[:5]

ct900_counts = []
ct10000_counts = []
nu_counts = []
# hcal depth loop
for i in range(1,8):
    hcal_str = "l1CaloTowerEmuTree/L1CaloTowerTree/CaloTP/hcalTPDepth{}".format(i)
    print(hcal_str)
    ct900_arr = merge_arrays(ct900_files, hcal_str)
    ct10000_arr = merge_arrays(ct10000_files, hcal_str)
    nu_arr = merge_arrays(nu_files, hcal_str)

    # fill total counts for each arr
    ct900_counts.append(ak.sum(ct900_arr) / len(ct900_arr))
    ct10000_counts.append(ak.sum(ct10000_arr)/ len(ct10000_arr))
    nu_counts.append(ak.sum(nu_arr)/ len(nu_arr))



ct900_time = []
ct10000_time = []
nu_time = []
# hcal timing loop
for i in range(1,8):
    hcal_str = "l1CaloTowerEmuTree/L1CaloTowerTree/CaloTP/hcalTPtiming{}".format(i)
    print(hcal_str)
    ct900_arr = merge_arrays(ct900_files, hcal_str)
    ct10000_arr = merge_arrays(ct10000_files, hcal_str)
    nu_arr = merge_arrays(nu_files, hcal_str)

    # fill total counts for each arr, over number of events
    ct900_time.append(ak.sum(ct900_arr[ct900_arr > 0]) / len(ct900_arr))
    ct10000_time.append(ak.sum(ct10000_arr[ct10000_arr > 0])/ len(ct10000_arr))
    nu_time.append(ak.sum(nu_arr[nu_arr > 0])/ len(nu_arr))

# convert to arrays
ct900_counts = np.array(ct900_counts)
ct10000_counts = np.array(ct10000_counts)
nu_counts = np.array(nu_counts)

x_ticks = np.arange(1, 8, 1)

fig1 = plt.figure()
plt.scatter(x_ticks, ct900_counts, label="ct 900", marker="x")
plt.scatter(x_ticks, ct10000_counts, label="ct 10000", marker="x")
plt.scatter(x_ticks, nu_counts, label="Nu baseline", marker="x")
plt.xlabel("HCal depth")
plt.ylabel("Accumulated counts")
plt.title("Counts per HCal depth")
plt.legend()
#plt.yscale("log")
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_depth_plot.png")

# plot timing
fig2 = plt.figure()
plt.scatter(x_ticks, ct900_time, label="ct 900", marker="x")
plt.scatter(x_ticks, ct10000_time, label="ct 10000", marker="x")
plt.scatter(x_ticks, nu_time, label="Nu baseline", marker="x")
plt.xlabel("HCal TP timing")
plt.ylabel("Accumulated counts over number of events")
plt.title("Counts per HCal timing")
plt.legend()
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_timing_plot.png")

# time vs depth
fig3 = plt.figure()
plt.scatter(ct900_counts, ct900_time, label="ct 900", marker="x")
plt.scatter(ct10000_counts, ct10000_time, label="ct 10000", marker="x")
plt.scatter(nu_counts, nu_time, label="Nu baseline", marker="x")
plt.ylabel("HCal TP timing")
plt.xlabel("HCAL depth")
plt.title("HCal depth vs. timing")
plt.legend()
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_depth_timing_ct900.png")

"""
fig2 = plt.figure()
plt.scatter(x_ticks, ct900_counts/np.sum(ct900_counts), label="ct 900", marker="x")
plt.scatter(x_ticks, ct10000_counts/np.sum(ct10000_counts), label="ct 10000", marker="x")
plt.scatter(x_ticks, nu_counts/np.sum(nu_counts), label="Nu baseline", marker="x")
plt.xlabel("HCal depth")
plt.ylabel("Accumulated counts [Normed]")
plt.title("Normed counts per HCal depth")
plt.legend()
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_depth_plot_normed.png")
"""
