#!/usr/bin/env python
import matplotlib

matplotlib.use("Agg")

import ROOT as ROOT
import uproot as up
import matplotlib.pyplot as plt
import mplhep as hep
import awkward as ak
import numpy as np
import os


def search_files(directory=None):
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


def create_TChain(treename, file_list):
    chain = ROOT.TChain(treename)
    for file in file_list:
        chain.AddFile(file)
    return chain


def IPhiToCMSPhi(iphi):
    dPhi = 2.0 * np.pi / 72.0
    phi = (abs(iphi) - 1 + abs(iphi)) / 2 * dPhi + np.pi
    return ak.where(phi > np.pi, phi - 2 * np.pi, phi)


def IEtaToCMSEta(ieta):
    dEta = 2.0 * np.pi / 72.0
    eta = (abs(ieta) - 1 + abs(ieta)) / 2 * dEta
    return ak.where(ieta > 0, eta, -1 * eta)


def dR(eta, phi, eta_H=0, phi_H=0):
    delta_R = np.sqrt((eta - eta_H) ** 2 + (phi - phi_H) ** 2)
    # cut on jets in close proximity
    return ak.where(delta_R < 0.4, delta_R, 999)

ct900_dir = "crab_hcal_testLLP-MH125_MFF_12_CT900_to4mu_def/results"
ct10000_dir = "crab_hcal_testLLP-MH350_MFF160_CT10000_with_new_MC_def/results"
nu_dir = "crab_hcal_testLLP-Nu_Pt20_with_reco_def/results"

ct900_files = search_files(ct900_dir)
ct10000_files = search_files(ct10000_dir)[:1]
nu_files = search_files(nu_dir)[:1]

# H0/H02 pdgid 35
genId_str = "l1GeneratorTree/L1GenTree/Generator/partId"
ct900_H0mask = merge_arrays(ct900_files, genId_str) == 35
ct10000_H0mask = merge_arrays(ct10000_files, genId_str) == 35
nu_H0mask = merge_arrays(nu_files, genId_str) == 35

geneta_str = "l1GeneratorTree/L1GenTree/Generator/partEta"

ct900_geneta = merge_arrays(ct900_files, geneta_str)  # [ct900_H0mask]
ct10000_geneta = merge_arrays(ct10000_files, geneta_str)[ct10000_H0mask]
nu_geneta = merge_arrays(nu_files, geneta_str)[nu_H0mask]

genphi_str = "l1GeneratorTree/L1GenTree/Generator/partPhi"

ct900_genphi = merge_arrays(ct900_files, genphi_str)[ct900_H0mask]
ct10000_genphi = merge_arrays(ct10000_files, genphi_str)[ct10000_H0mask]
nu_genphi = merge_arrays(nu_files, genphi_str)[nu_H0mask]

gen_parent_str = "l1GeneratorTree/L1GenTree/Generator/partParent"

ct900_gen_parentmask = merge_arrays(ct900_files, gen_parent_str) == 35
ct900_geneta_parent = merge_arrays(ct900_files, geneta_str)[
    ct900_gen_parentmask
]  # [:,0]
ct900_genphi_parent = merge_arrays(ct900_files, genphi_str)[
    ct900_gen_parentmask
]  # [:,0]

# unpacking hcal eta, calculating distance to gen-level daughters
# finding the TP wth the smallest distance, generating a mask to find hcal TPs
# ho
ho_str = "l1HOTree/L1HOTree/L1HO/"
ct900_ho_eta = IEtaToCMSEta(merge_arrays(ct900_files, ho_str + "hcalDetIdIEta"))
ct900_ho_phi = IPhiToCMSPhi(merge_arrays(ct900_files, ho_str + "hcalDetIdIPhi"))

# ho cut eta < 1.4
ct900_geneta_parent = ak.where(abs(ct900_geneta_parent) < 1.4, ct900_geneta_parent, -999)

ho_parent_mask_0 = ak.argmin(
    dR(
        ct900_ho_eta, ct900_ho_phi, ct900_geneta_parent[:, 0], ct900_genphi_parent[:, 0]
    ),
    axis=1,
    keepdims=True,
)
ho_parent_mask_1 = ak.argmin(
    dR(
        ct900_ho_eta, ct900_ho_phi, ct900_geneta_parent[:, 1], ct900_genphi_parent[:, 1]
    ),
    axis=1,
    keepdims=True,
)
ho_parent_mask = ak.concatenate((ho_parent_mask_0, ho_parent_mask_1), axis=1)
ho_parent_matched_eta = ct900_ho_eta[ho_parent_mask]

#redo ct900 with ct10000
gen_parent_str = "l1GeneratorTree/L1GenTree/Generator/partParent"

ct10000_gen_parentmask = merge_arrays(ct10000_files, gen_parent_str) == 35
ct10000_geneta_parent = merge_arrays(ct10000_files, geneta_str)[
    ct10000_gen_parentmask
]  # [:,0]
ct10000_genphi_parent = merge_arrays(ct10000_files, genphi_str)[
    ct10000_gen_parentmask
]  # [:,0]

# unpacking hcal eta, calculating distance to gen-level daughters
# finding the TP wth the smallest distance, generating a mask to find hcal TPs
# ho
ho_str = "l1HOTree/L1HOTree/L1HO/"
ct10000_ho_eta = IEtaToCMSEta(merge_arrays(ct10000_files, ho_str + "hcalDetIdIEta"))
ct10000_ho_phi = IPhiToCMSPhi(merge_arrays(ct10000_files, ho_str + "hcalDetIdIPhi"))

# ho cut eta < 1.4
ct10000_geneta_parent = ak.where(abs(ct10000_geneta_parent) < 1.4, ct10000_geneta_parent, -999)

ho_10000_parent_mask_0 = ak.argmin(
    dR(
        ct10000_ho_eta, ct10000_ho_phi, ct10000_geneta_parent[:, 0], ct10000_genphi_parent[:, 0]
    ),
    axis=1,
    keepdims=True,
)
ho_10000_parent_mask_1 = ak.argmin(
    dR(
        ct10000_ho_eta, ct10000_ho_phi, ct10000_geneta_parent[:, 1], ct10000_genphi_parent[:, 1]
    ),
    axis=1,
    keepdims=True,
)
ho_10000_parent_mask = ak.concatenate((ho_10000_parent_mask_0, ho_10000_parent_mask_1), axis=1)
ho_10000_parent_matched_eta = ct10000_ho_eta[ho_10000_parent_mask]


# comparing hist of H0 and matched
fig_H0_matched = plt.figure(figsize=(12, 9))
hep.style.use("CMS")
hep.cms.label(
    llabel="Work in progress",
    loc=0,
    fontsize=20,
)
#plt.hist(ak.flatten(ct900_ho_eta), label="HO 900", histtype="step", density=True, bins=25)  # , density=True)
#plt.hist(ak.flatten(ct10000_ho_eta), label="HO 10000", histtype="step", density=True, bins=25)  # , density=True)
plt.hist(
    ak.flatten(ho_parent_matched_eta),
    label="ct900: HO matched to daughter from H0",
    histtype="step", density=True, bins=72)  # , color="black") #
plt.hist(
    ak.flatten(ho_10000_parent_matched_eta),
    label="ct10000: HO matched to daughter from H0",
    histtype="step", density=True, bins=72)
plt.xlabel("Eta")
plt.ylabel("Counts")
# plt.title("Matching between hcal TP and H0 daughters")
plt.legend(
    ncol=1,
    loc="upper right",
    bbox_to_anchor=(1, 0.95),
    borderaxespad=0,
)  # loc='upper center', bbox_to_anchor=(.5,1.15))#, fancybox=True, shadow=True)
plt.yscale("log")
plt.savefig(
    "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/Eta_matching_daughter_only_LLP.png"
)

from IPython import embed; embed()
