#!/usr/bin/env python

import uproot as up
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import pickle

with open("key_names_1.pkl", "rb") as f:
    key_names = pickle.load(f)

with open("125_CT900_with_reco_def.pkl", "rb") as f:
    arr_125 = pickle.load(f)

with open("350_CT10000_with_reco_def.pkl", "rb") as f:
    arr_350 = pickle.load(f)

with open("Nu_Pt20_with_reco_def.pkl", "rb") as f:
    arr_nu = pickle.load(f)

# with open('arr_QCD_1.pkl', 'rb') as f:
#    arr_QCD = pickle.load(f)
with open("QCD_Pt3000_reduced.pkl", "rb") as f:
    arr_QCD = pickle.load(f)


print("#####Plotting#####")
for i, k in enumerate(key_names):
    print(k)
    # define multi axes variales
    # if any(string in k for string in ["hcalDetIdIEta","hcalDetIdIPhi","hcalQIESample","egEt", "eg", "Et"]):

    # broken arrays
    if k in ["orbit", "hlt", "jetM", "partId"]:
        continue

    elif type(arr_125[i][1]) == ak.highlevel.Array:
        # k in ["hcalDetIdIEta","hcalDetIdIPhi","hcalQIESample","egEt", "egEt"]:  #
        #
        # from IPython import embed;embed()
        # embed()
        # continue
        h_125 = ak.mean(arr_125[i], axis=1)
        h_350 = ak.mean(arr_350[i], axis=1)
        h_nu = ak.mean(arr_nu[i], axis=1)
        h_QCD = ak.mean(arr_QCD[i], axis=1)

    # define regualar arrays
    elif type(arr_125[i][1]) == int or type(arr_125[i][1]) == float:
        h_125 = ak.to_numpy(arr_125[i])
        h_350 = ak.to_numpy(arr_350[i])
        h_nu = ak.to_numpy(arr_nu[i])
        h_QCD = ak.to_numpy(arr_QCD[i])

    # everything else
    else:
        from IPython import embed

        embed()
        continue
    # everything else
    fig = plt.figure()  # arr_125[i]
    plt.hist(
        h_125,
        label="MH 125, total yields={}".format(np.round(ak.sum(h_125), 2)),
        density=True,
        histtype="step",
        linewidth=2,
    )
    plt.hist(
        h_350,
        label="MH 350, total yields={}".format(np.round(ak.sum(h_350), 2)),
        density=True,
        histtype="step",
        linewidth=2,
    )
    plt.hist(
        h_QCD,
        label="QCD, total yields={}".format(np.round(ak.sum(h_QCD), 2)),
        density=True,
        histtype="step",
        linewidth=2,
    )
    if type(h_nu[0]) != type(None):
        # from IPython import embed
        # embed()
        plt.hist(
            h_nu,
            label="nu Pt2_20, total yields={}".format(np.round(ak.sum(h_nu), 2)),
            density=True,
            histtype="step",
            linewidth=2,
        )
    # plt.title('Time per events')
    plt.xlabel(k)
    plt.ylabel("Normalised")
    plt.legend(loc="upper center", bbox_to_anchor=(0.5, 1.1), ncol=2)
    plt.savefig("reco_plots/reco_comp_-" + k + ".png")
    plt.close()
