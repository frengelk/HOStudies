#!/usr/bin/env python

import uproot as up
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import pickle

if __name__ == "__main__":

    loc = "/nfs/dust/cms/user/frengelk/Code/cmssw/CMSSW_11_0_2/src/HcalTrigger/Validation/scripts/"

    # with open(loc + "key_names_1.pkl", "rb") as f:
        # key_names = pickle.load(f)
#
    # with open(loc + "125_CT900_with_reco_def.pkl", "rb") as f:
        # arr_125 = pickle.load(f)
#
    # with open(loc + "350_CT10000_with_reco_def.pkl", "rb") as f:
        # arr_350 = pickle.load(f)

    with open("muon_reco.rootdict_real.pkl", "rb") as f:
        arr_muon = pickle.load(f)

    # with open('arr_nu.pkl', 'rb') as f:
    # arr_nu = pickle.load(f)
    #
    # with open('arr_QCD.pkl', 'rb') as f:
    # arr_QCD = pickle.load(f)

    # hist_pairs = [('nJet', 'time'),
    #            ('jetPt', 'sumEt'),
    #            ('muonEta','muonEt')
    # ]

    key_names = arr_muon.keys()

    hist_pairs = [
        (pair1, pair2)
        for i,pair1 in enumerate(key_names)
        for j,pair2 in enumerate(key_names)
        #if i>j
        #if key_names.index(pair2) > key_names.index(pair1)
        #  list(a.keys()).index('L1UpgradeKBmtfMuon/tfMuonWh') a=collections.OrderedDict(arr_muon)
    ]

    #from IPython import embed;embed()

    for pair in hist_pairs:
        #print(pair)

        if pair[0] in ["orbit", "hlt", "jetM"] or pair[1] in ["orbit", "hlt", "jetM"]:
            continue

        if not ("L1HO" in pair[0]):
            continue
        if not ("L1HO" in pair[1]):
            continue

        if type(arr_muon[pair[0]][0]) == ak.highlevel.Array:
            #from IPython import embed;embed()

            h_muon_1 = ak.to_numpy(ak.mean(arr_muon[pair[0]], axis=1))
            #h_350_1 = ak.to_numpy(ak.mean(arr_350[key_names.index(pair[0])], axis=1))

        # define regualar arrays
        elif type(arr_muon[pair[0]][0]) == int:
            h_muon_1 = ak.to_numpy(arr_muon[pair[0]])
            #h_350_1 = ak.to_numpy(arr_350[key_names.index(pair[0])])

        if type(arr_muon[pair[1]][0]) == ak.highlevel.Array:
            h_muon_2 = ak.to_numpy(ak.mean(arr_muon[pair[1]], axis=1))
            #h_350_2 = ak.to_numpy(ak.mean(arr_350[key_names.index(pair[1])], axis=1))

        elif type(arr_muon[pair[1]][0]) == int:
            h_muon_2 = ak.to_numpy(arr_muon[pair[1]])
            #h_350_2 = ak.to_numpy(arr_350[key_names.index(pair[1])])

        # from IPython import embed;embed()

        #if np.sum(h_muon_1) < 1 or np.sum(h_muon_2) < 1:
        #    continue

        print(pair)


        if (
            type(h_muon_1) == memoryview
            or type(h_muon_2) == memoryview
            #or type(h_350_1) == np.ma.core.MaskedArray
            #or type(h_350_2) == np.ma.core.MaskedArray
        ):
            print("skipped")
            continue

        if (
            type(h_muon_1) == np.ma.core.MaskedArray
            or type(h_muon_2) == np.ma.core.MaskedArray
            #or type(h_350_1) == np.ma.core.MaskedArray
            #or type(h_350_2) == np.ma.core.MaskedArray
        ):
            h_muon_1=h_muon_1.data
            h_muon_2=h_muon_2.data

        if (
            type(h_muon_1) == memoryview
            or type(h_muon_2) == memoryview
            #or type(h_350_1) == np.ma.core.MaskedArray
            #or type(h_350_2) == np.ma.core.MaskedArray
        ):
            print("skipped")
            continue

        if pair == ('L1HO/nHcalDetIds', 'L1UpgradeBmtfMuon/tfMuonProcessor') or pair == ('L1HO/nHcalDetIds', 'L1Upgrade/egTowerIEta'):
            from IPython import embed

            #embed()

        fig = plt.figure()
        plt.hist2d(h_muon_1, h_muon_2, bins=20, label="reco muons")
        plt.xlabel(pair[0])
        plt.ylabel(pair[1])
        # plt.legend()
        plt.savefig(
            loc + "reco_plots/muon_plots/2D_reco_muon_{}_{}.png".format(pair[0].replace("/", "_"), pair[1].replace("/", "_"))
        )
        plt.clf()
        plt.close()

        # fig = plt.figure()
        # plt.hist2d(h_350_1, h_350_2, bins=40, label="350")
        # plt.xlabel(pair[0])
        # plt.ylabel(pair[1])
        # # plt.legend()
        # plt.savefig(
            # loc + "reco_plots/plots_2D/2D_reco_350_{}_{}.png".format(pair[0], pair[1])
        # )
        # plt.clf()
        # plt.close()


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
