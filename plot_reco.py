#!/nfs/dust/cms/user/frengelk/Anaconda/envs/EPR_env/bin/python

import uproot as up
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import mplhep as hep
import pickle

# with open("key_names_HO.pkl", "rb") as f:
# key_names = pickle.load(f)
#
# with open("125_CT900_with_reco_def_HO.pkl", "rb") as f:
# arr_125 = pickle.load(f)
#
# with open("350_CT10000_with_reco_def_HO.pkl", "rb") as f:
# arr_350 = pickle.load(f)

# with open("key_names_1.pkl", "rb") as f:
with open("key_names_HO_real.pkl", "rb") as f:
    key_names = pickle.load(f)
# with open("125_CT900_with_old_reco_real_new_conddict_real.pkl", "rb") as f:
# with open("125_CT900_with_new_reco_real_def_defdict_real.pkl", "rb") as f:
# arr_125 = pickle.load(f)
# with open("350_CT10000_with_new_reco_real_new_conddict_real.pkl", "rb") as f:
# with open("350_CT10000_with_new_reco_real_def_defdict_real.pkl", "rb") as f:
# arr_350 = pickle.load(f)
# #with open("Nu_Pt20_with_new_reco_real_new_cond_reduceddict_real.pkl", "rb") as f:
# with open("Nu_Pt20_with_new_reco_real_def_reduceddict_real.pkl", "rb") as f:
# arr_nu = pickle.load(f)
# with open("testMuons_MH125_CT900_defdict_real.pkl", "rb") as f:
# arr_muon = pickle.load(f)

with open("muon_reco.rootdict_real.pkl", "rb") as f:
    arr_muon = pickle.load(f)


# from IPython import embed;embed()
"""
with open('arr_QCD_1.pkl', 'rb') as f:
arr_QCD = pickle.load(f)
with open("QCD_Pt3000_reduced.pkl", "rb") as f:
arr_QCD = pickle.load(f)
"""

print("#####Plotting def#####")
for i, k in enumerate(arr_muon.keys()):  # key_names
    # if not "hcal" in k:

    # define multi axes variales
    # if any(string in k for string in ["hcalDetIdIEta","hcalDetIdIPhi","hcalQIESample","egEt", "eg", "Et"]):

    # broken arrays
    if any(
        check in k for check in ["orbit", "hlt", "jetM", "partId", "tfMuonDecodedTrAdd"]
    ):
        print("skipped", k)
        continue

    #    if k=="muonEta":
    #        from IPython import embed;embed()

    elif type(arr_muon[k][1]) == ak.highlevel.Array:
        # k in ["hcalDetIdIEta","hcalDetIdIPhi","hcalQIESample","egEt", "egEt"]:  #
        #

        # continue
        # h_125 = ak.to_numpy(ak.fill_none(ak.max(arr_125[k], axis=1),0))
        # h_350 = ak.to_numpy(ak.fill_none(ak.flatten(arr_350[k]), 0))
        # h_nu = ak.to_numpy(ak.fill_none(ak.max(arr_nu[k], axis=1),0))
        h_muon = ak.to_numpy(ak.fill_none(ak.flatten(arr_muon[k]), 0))
        # h_125 = ak.mean(arr_125[i], axis=1)
        # h_350 = ak.mean(arr_350[i], axis=1)
        # h_nu = ak.mean(arr_nu[i], axis=1)
        # h_QCD = ak.mean(arr_QCD[i], axis=1)
        # from IPython import embed;embed()

    # define regualar arrays
    elif type(arr_muon[k][1]) == int or type(arr_muon[k][1]) == float:
        # h_125 = ak.to_numpy(ak.fill_none(arr_125[k],0))
        # h_350 = ak.to_numpy(ak.fill_none(arr_350[k], 0))
        # h_nu = ak.to_numpy(ak.fill_none(arr_nu[k],0))
        h_muon = ak.to_numpy(ak.fill_none(arr_muon[k], 0))
        # h_QCD = ak.to_numpy(arr_QCD[i])

    # everything else
    else:
        print(k, "went wrong")
        # embed()
        continue
    # everything else

    print(k)

    bins = np.histogram(h_muon, bins=30)[
        1
    ]  # np.concatenate((h_350, h_muon) np.hstack((h_125,h_350, h_nu)), bins=30)[1]

    plt.figure(figsize=(18, 9))
    # plt.hist(
    # h_125,
    # label="MH 125, ",#total yields={}".format(np.round(ak.sum(h_125), 2)),
    # density=True,
    # bins=bins,
    # histtype="step",
    # linewidth=1.5,
    # #linestyle="--",
    # )
    # plt.hist(
    # h_350,
    # label="MH 350 into 4b, ",  # total yields={}".format(np.round(ak.sum(h_350), 2)),
    # density=True,
    # bins=bins,
    # histtype="step",
    # linewidth=1.5,
    # )
    plt.hist(
        h_muon,
        label="MH 125 into 4 muon, ",  # total yields={}".format(np.round(ak.sum(h_350), 2)),
        density=True,
        bins=bins,
        histtype="step",
        linewidth=1.5,
    )
    # plt.hist(
    # h_QCD,
    # label="QCD, total yields={}".format(np.round(ak.sum(h_QCD), 2)),
    # density=True,
    # histtype="step",
    # linewidth=1.5,
    # linestyle=":",
    # )
    #
    # if type(h_nu[0]) != type(None):
    # # from IPython import embed
    # # embed()
    # plt.hist(
    # h_nu,
    # label="nu Pt2_20, ",#total yields={}".format(np.round(ak.sum(h_nu), 2)),
    # density=True,
    # histtype="step",
    # linewidth=1.5,
    # bins=bins,
    # linestyle="--",
    # #linestyle="-.",
    # )
    # plt.title('Time per events')
    hep.style.use("CMS")
    hep.cms.label(llabel="Work in progress", loc=0)
    plt.xlabel(k)
    plt.ylabel("Normalised")
    # plt.legend(loc="upper right", bbox_to_anchor=(1.01, 0.99), ncol=1)
    plt.legend(loc="best", ncol=1)
    plt.tight_layout()
    plt.savefig(
        "reco_plots/muon_plots/muon_def-" + k.replace("/", "_") + ".png",
        bbox_inches="tight",
    )
    plt.yscale("log")
    plt.savefig(
        "reco_plots/muon_plots/muon_def_log-" + k.replace("/", "_") + ".png",
        bbox_inches="tight",
    )
    plt.close()

    print("plotted", k)
