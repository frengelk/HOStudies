#!/usr/bin/env python

import uproot as up
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import pickle

with open("old_key_names.pkl", "rb") as f:
    key_names = pickle.load(f)

with open("old_arr_125.pkl", "rb") as f:
    arr_125 = pickle.load(f)

with open("old_arr_350.pkl", "rb") as f:
    arr_350 = pickle.load(f)

# with open('arr_nu.pkl', 'rb') as f:
# arr_nu = pickle.load(f)
#
# with open('arr_QCD.pkl', 'rb') as f:
# arr_QCD = pickle.load(f)

# from IPython import embed;embed()

# hist_pairs = [('nJet', 'time'),
#            ('jetPt', 'sumEt'),
#            ('muonEta','muonEt')
# ]

hist_pairs = [
    (pair1, pair2)
    for pair1 in key_names
    for pair2 in key_names
    if key_names.index(pair2) > key_names.index(pair1)
]

from IPython import embed

embed()

for pair in hist_pairs:
    # print(pair)

    if pair[0] in ["orbit", "hlt", "jetM"] or pair[1] in ["orbit", "hlt", "jetM"]:
        continue

    if type(arr_125[key_names.index(pair[0])][1]) == ak.highlevel.Array:
        h_125_1 = ak.to_numpy(ak.mean(arr_125[key_names.index(pair[0])], axis=1))
        h_350_1 = ak.to_numpy(ak.mean(arr_350[key_names.index(pair[0])], axis=1))

    # define regualar arrays
    elif type(arr_125[key_names.index(pair[0])][1]) == int:
        h_125_1 = ak.to_numpy(arr_125[key_names.index(pair[0])])
        h_350_1 = ak.to_numpy(arr_350[key_names.index(pair[0])])

    if type(arr_125[key_names.index(pair[1])][1]) == ak.highlevel.Array:
        h_125_2 = ak.to_numpy(ak.mean(arr_125[key_names.index(pair[1])], axis=1))
        h_350_2 = ak.to_numpy(ak.mean(arr_350[key_names.index(pair[1])], axis=1))

    elif type(arr_125[key_names.index(pair[1])][1]) == int:
        h_125_2 = ak.to_numpy(arr_125[key_names.index(pair[1])])
        h_350_2 = ak.to_numpy(arr_350[key_names.index(pair[1])])

    # from IPython import embed;embed()

    fig = plt.figure()
    plt.hist2d(h_125_1, h_125_2, bins=20, label="125")
    plt.xlabel(pair[0])
    plt.ylabel(pair[1])
    # plt.legend()
    plt.savefig("reco_plots/plots_2D/2D_reco_125_{}_{}.png".format(pair[0], pair[1]))
    plt.clf()
    plt.close()

    fig = plt.figure()
    plt.hist2d(h_350_1, h_350_2, bins=20, label="350")
    plt.xlabel(pair[0])
    plt.ylabel(pair[1])
    # plt.legend()
    plt.savefig("reco_plots/plots_2D/2D_reco_350_{}_{}.png".format(pair[0], pair[1]))
    plt.clf()
    plt.close()


"""
fig = plt.figure()
plt.hist(
    h_125,
    label="MH 125, total yields={}".format(ak.sum(h_125)),
    density=True,
    histtype="step",
)
plt.hist(
    h_350,
    label="MH 350, total yields={}".format(ak.sum(h_350)),
    density=True,
    histtype="step",
)
plt.hist(
    h_nu,
    label="nu Pt2_20, total yields={}".format(ak.sum(h_nu)),
    density=True,
    histtype="step",
)
plt.hist(
    h_QCD,
    label="QCD, total yields={}".format(ak.sum(h_QCD)),
    density=True,
    histtype="step",
)
# plt.title('Time per events')
plt.xlabel(k)
plt.ylabel("Normalised")
plt.legend()
plt.savefig("reco_plots/reco_comp_-" + k + ".png")
plt.close()
"""
