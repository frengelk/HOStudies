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

"""
# let's do root drawing
width = 1200
#height = width
#height = 1500
#height = int(width / 1.618030)
height = int(width * 3/4.)
effName = "delta R"
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadLeftMargin(0.115)
ROOT.gStyle.SetPadRightMargin(0.03)
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetErrorX(0)

# defining root canvas
canvas = ROOT.TCanvas(effName, effName, width, height)
#pad = ROOT.gPad

# pad for drawing on
#pad1 = ROOT.TPad('pad1','This is pad1',0.02,0.52,0.48,0.98,21)
#pad1.cd()
#pad1.Range(0, 0, 1000, 40000)
#pad1.RangeAxis(0, 0, 1000, 40000)

# create legend to be shown laterleg = ROOT.TLegend(0.1, 0.7, 0.48, 0.9)
leg = ROOT.TLegend(0.8, 0.8, 0.95, 0.95)

# creating TChain objects to draw from several files at once
chain_ct10000 = create_TChain("l1GeneratorTree/L1GenTree", ct10000_files)
chain_ct10000.SetMarkerColor(2)
chain_ct10000.SetLineColor(2)
chain_ct10000.Draw("(jetEta**2 + jetPhi**2)**(1/2)", "", "E1")
leg.AddEntry(chain_ct10000, "ct10000", "l")


chain_ct900 = create_TChain("l1GeneratorTree/L1GenTree", ct900_files)
chain_ct900.SetMarkerColor(1)
chain_ct900.SetLineColor(1)
#chain_ct900.SetLineStyle(1)
#chain_ct900.SetFillColor(1)
chain_ct900.Draw("(jetEta**2 + jetPhi**2)**(1/2)", "", "SAME E1")
leg.AddEntry(chain_ct900, "ct900", "l")

# accessing first hist and altering range
#th1 = canvas.GetPrimitive("htemp")
#th1.GetXaxis().SetRange(0,10)
#th1.GetYaxis().SetRange(0,40000)

chain_nu = create_TChain("l1GeneratorTree/L1GenTree", nu_files)
chain_nu.SetMarkerColor(3)
chain_nu.SetLineColor(3)
chain_nu.Draw("(jetEta**2 + jetPhi**2)**(1/2)", "", "SAME E1")
#th1 = canvas.GetPrimitive("htemp")
#th1.GetXaxis().SetRange(0,10)
#th1.GetYaxis().SetRange(0,40000)


leg.AddEntry(chain_nu, "NuGun", "l")
leg.SetHeader("Datasets")
leg.SetTextSize(0.03)
leg.Draw()

#canvas.RangeAxis(0, 0, 1000, 40000)
#canvas.Range(0, 0, 1000, 40000)
#canvas.Update()
#canvas.Draw()
#pad.Range(0,0,10,10000)
#pad.RangeAxis(0,0,10,10000)
#canvas.RangeAxisChanged()
#canvas.RangeChanged()

canvas.SaveAs("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/hcal_TChain_dR.png")
canvas.SaveAs("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/hcal_TChain_dR.pdf")

a = input("wait")
print(a)

canvas.Close()
"""

ct900_counts = []
ct10000_counts = []
nu_counts = []

# creating eta masks for loop
eta_str = "l1GeneratorTree/L1GenTree/Generator/jetEta"
ct900_mask = ak.any(abs(merge_arrays(ct900_files, eta_str)) < 1.262, axis=1)
ct10000_mask = ak.any(abs(merge_arrays(ct10000_files, eta_str)) < 1.262, axis=1)
nu_mask = ak.any(abs(merge_arrays(nu_files, eta_str)) < 1.262, axis=1)

# eta phi dR
phi_str = "l1GeneratorTree/L1GenTree/Generator/jetPhi"
eta_900 = merge_arrays(ct900_files, eta_str)
phi_900 = merge_arrays(ct900_files, phi_str)
dR_900 = (phi_900**2 + eta_900) ** (1 / 2)

eta_10000 = merge_arrays(ct10000_files, eta_str)
phi_10000 = merge_arrays(ct10000_files, phi_str)
dR_10000 = (phi_10000**2 + eta_10000) ** (1 / 2)

eta_nu = merge_arrays(nu_files, eta_str)
phi_nu = merge_arrays(nu_files, phi_str)
dR_nu = (phi_nu**2 + eta_nu) ** (1 / 2)

a = up.open(ct900_files[0])
iphi = a["l1CaloTowerEmuTree/L1CaloTowerTree/CaloTP/hcalTPiphi"].array()
phi = IPhiToCMSPhi(iphi)

ieta = a["l1CaloTowerEmuTree/L1CaloTowerTree/CaloTP/hcalTPieta"].array()
eta = IPhiToCMSPhi(ieta)

hcal_ieta_str = "l1CaloTowerEmuTree/L1CaloTowerTree/CaloTP/hcalTPieta"
ct900_ieta = IEtaToCMSEta(merge_arrays(ct900_files, hcal_ieta_str))

hcal_iphi_str = "l1CaloTowerEmuTree/L1CaloTowerTree/CaloTP/hcalTPiphi"
ct900_iphi = IPhiToCMSPhi(merge_arrays(ct900_files, hcal_iphi_str))


delta_R = dR(eta, phi)
genEta = a["l1GeneratorTree/L1GenTree/Generator/partEta"].array()

genPhi = a["l1GeneratorTree/L1GenTree/Generator/partPhi"].array()

gendelta_R = dR(genEta, genPhi)

max_delta_r = ak.flatten(delta_R[ak.argmax(delta_R, axis=1, keepdims=True)])
delta_r_argmax = ak.argmax(delta_R, axis=1, keepdims=True)


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

parent_mask_0 = ak.argmin(
    dR(ct900_ieta, ct900_iphi, ct900_geneta_parent[:, 0], ct900_genphi_parent[:, 0]),
    axis=1,
    keepdims=True,
)
parent_mask_1 = ak.argmin(
    dR(ct900_ieta, ct900_iphi, ct900_geneta_parent[:, 1], ct900_genphi_parent[:, 1]),
    axis=1,
    keepdims=True,
)

# parent_matched_eta_0 = ct900_ieta[parent_distance]
# unpacking hcal eta, calculating distance to gen-level daughters
# finding the TP wth the smallest distance, generating a mask to find hcal TPs
parent_mask = ak.concatenate((parent_mask_0, parent_mask_1), axis=1)
parent_matched_eta = ct900_ieta[parent_mask]


# ho
ho_energy_str = "l1HOTree/L1HOTree/L1HO/SampleEnergy"
ct900_ho = merge_arrays(ct900_files, ho_energy_str)
ct10000_ho = merge_arrays(ct10000_files, ho_energy_str)
nu_ho = merge_arrays(nu_files, ho_energy_str)
ho_str = "l1HOTree/L1HOTree/L1HO/"
ho_eta = a[ho_str]["hcalDetIdIEta"].array()

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


TPet = a["l1CaloTowerEmuTree/L1CaloTowerTree;4"]["CaloTP/hcalTPet"].array()[
    delta_r_argmax
]
dat = ak.to_numpy(ak.flatten(TPet)).data


"""
# macthing H0 and hcal TP
b = ct900_H0mask[:10000]
H0_eta = ak.flatten(genEta[b])
H0_phi = ak.flatten(genPhi[b])
distance = dR(eta, phi, H0_eta, H0_phi)
matched_eta = eta[ak.argmin(distance, axis=1, keepdims=True)]
"""

# hcal depth loop
for i in range(1, 8):
    hcal_str = "l1CaloTowerEmuTree/L1CaloTowerTree/CaloTP/hcalTPDepth{}".format(i)
    print(hcal_str)
    ct900_arr = merge_arrays(ct900_files, hcal_str)
    ct10000_arr = merge_arrays(ct10000_files, hcal_str)
    nu_arr = merge_arrays(nu_files, hcal_str)

    # fill total counts for each arr
    ct900_counts.append(ak.sum(ct900_arr[ct900_mask]) / len(ct900_arr))
    ct10000_counts.append(ak.sum(ct10000_arr[ct10000_mask]) / len(ct10000_arr))
    nu_counts.append(ak.sum(nu_arr[nu_mask]) / len(nu_arr))


ct900_time = []
ct10000_time = []
nu_time = []
# hcal timing loop
for i in range(1, 8):
    hcal_str = "l1CaloTowerEmuTree/L1CaloTowerTree/CaloTP/hcalTPtiming{}".format(i)
    print(hcal_str)
    ct900_arr = merge_arrays(ct900_files, hcal_str)
    ct10000_arr = merge_arrays(ct10000_files, hcal_str)
    nu_arr = merge_arrays(nu_files, hcal_str)

    # fill total counts for each arr, over number of events
    ct900_time.append(ak.sum(ct900_arr[ct900_arr > 0]) / len(ct900_arr))
    ct10000_time.append(ak.sum(ct10000_arr[ct10000_arr > 0]) / len(ct10000_arr))
    nu_time.append(ak.sum(nu_arr[nu_arr > 0]) / len(nu_arr))

# convert to arrays
ct900_counts = np.array(ct900_counts)
ct10000_counts = np.array(ct10000_counts)
nu_counts = np.array(nu_counts)

# comparing hist of H0 and matched
fig_H0_matched = plt.figure(figsize=(12, 9))
hep.style.use("CMS")
hep.cms.label(
    llabel="Work in progress",
    loc=0,
    fontsize=20,
)
# plt.hist(ak.flatten(ct900_ieta), label="HCAL TP", histtype="step", density=True)
plt.hist(ak.flatten(ct900_ho_eta), label="HO", histtype="step")  # , density=True)
# plt.hist(H0_eta, label="Pdgid 35", histtype="step")
# plt.hist(ak.flatten(parent_matched_eta), label="HCAL TP matched to daughter from H0", histtype="step", density=True)
plt.hist(
    ak.flatten(ho_parent_matched_eta),
    label="HO matched to daughter from H0",
    histtype="step",
)  # , color="black") #
plt.xlabel("Eta")
plt.ylabel("Counts")
# plt.title("Matching between hcal TP and H0 daughters")
plt.legend(
    ncol=1,
    loc="upper right",
    bbox_to_anchor=(1, 0.5),
    borderaxespad=0,
)  # loc='upper center', bbox_to_anchor=(.5,1.15))#, fancybox=True, shadow=True)
plt.yscale("log")
plt.savefig(
    "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/Eta_matching_daughter_comp.png"
)

# eta comp
fig_eta_comp = plt.figure()
plt.hist(
    ak.flatten(ct900_geneta[abs(ct900_geneta) < 10]),
    bins=100,
    label="Generator particles",
    histtype="step",
)
# plt.hist(H0_eta, label="Pdgid 35", histtype="step")
plt.hist(ak.flatten(ct900_ieta), label="HCAl ieta converted", histtype="step")
plt.hist(ak.flatten(ieta), bins=50, label="ieta", histtype="step")
plt.xlabel("Eta")
plt.ylabel("Counts")
# plt.title("Matching between hcal TP and H0")
plt.legend()
plt.yscale("log")
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/Eta__comp.png")

# phi comp
fig_phi_comp = plt.figure()
plt.hist(
    ak.flatten(ct900_genphi), bins=100, label="Generator particles", histtype="step"
)
plt.hist(ak.flatten(ct900_iphi), label="HCAl ieta converted", histtype="step")
plt.hist(ak.flatten(iphi), bins=50, label="ieta", histtype="step")
plt.xlabel("Phi")
plt.ylabel("Counts")
# plt.title("Matching between hcal TP and H0")
plt.legend()
plt.yscale("log")
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/Phi__comp.png")

# ho energy
fig_ho_energy = plt.figure()
plt.hist(
    ak.flatten(ct900_ho[ct900_ho < 40]),
    label="HO",
    histtype="step",
    density=True,
    bins=20,
)
plt.hist(
    ak.flatten(ct900_ho[ho_parent_mask]),
    label="HO matched",
    histtype="step",
    density=True,
    bins=20,
)
plt.xlabel("HO E")
plt.ylabel("Counts")
plt.title("Matching between HO and H0 daughters")
plt.legend()
plt.xlim(0, 10)
plt.savefig(
    "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/ho_energy_matching_daughter_comp.png"
)
plt.yscale("log")
plt.savefig(
    "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/ho_energy_matching_daughter_comp_log.png"
)


# TPet
# fig_TPet = plt.figure()
# plt.hist(dat)
# plt.xlabel("HcalTPet")
# plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/TPet.png")

# eta Higgs2/H0 plot-> problem?

# hcal TPet
fig_ho = plt.figure()
plt.scatter(
    [1.0, 2],
    [ak.sum(ct900_ho) / len(ct900_ho), ak.sum(ct900_ho[ct900_mask]) / len(ct900_ho)],
    label="ct900",
)
plt.scatter(
    [1.0, 2],
    [
        ak.sum(ct10000_ho) / len(ct10000_ho),
        ak.sum(ct10000_ho[ct10000_mask]) / len(ct10000_ho),
    ],
    label="ct10000",
)
plt.scatter(
    [1.0, 2],
    [ak.sum(nu_ho) / len(nu_ho), ak.sum(nu_ho[nu_mask]) / len(nu_ho)],
    label="NuGun",
)
plt.xlabel("HO Sample Energy")
plt.ylabel("Accumulated counts/total events")
plt.title("No Eta mask vs Eta mask")
plt.legend()
# plt.yscale("log")
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/ho_plot.png")

# dR
fig_dR = plt.figure()
plt.hist(ak.flatten(dR_900), label="ct900", histtype="step")
plt.hist(ak.flatten(dR_10000), label="ct10000", histtype="step")
plt.hist(ak.flatten(dR_nu), label="nu", histtype="step")
plt.xlabel("dR")
plt.ylabel("Accumulated counts")
plt.title("(eta**2 + phi**2)**(1/2)")
plt.legend()
# plt.yscale("log")
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/dR_plot.png")

"""
# GenEta
fig_dR = plt.figure()
plt.hist(ct900_geneta, label="ct900", histtype="step")
plt.hist(ct10000_geneta, label="ct10000", histtype="step")
plt.hist(nu_geneta, label="nu", histtype="step")
plt.xlabel("GenEta")
plt.ylabel("Accumulated counts")
plt.title("Generator particles matched to id 35 (H0)")
plt.legend()
plt.yscale("log")
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/GenEta_35_plot.png")
"""

# GenPhi
fig_dR = plt.figure()
plt.hist(ct900_genphi, label="ct900", histtype="step")
plt.hist(ct10000_genphi, label="ct10000", histtype="step")
plt.hist(nu_genphi, label="nu", histtype="step")
plt.xlabel("Genphi")
plt.ylabel("Accumulated counts")
plt.title("Generator particles matched to id 35 (H0)")
plt.legend()
plt.yscale("log")
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/GenPhi_35_plot.png")


# leading jets are sorted:
# b = a["l1GeneratorTree/L1GenTree/Generator/jetPt"].array()
# In [42]: ak.sum(ak.any(abs(nu_jeta) < 1.4, axis=1)) / len(nu_jeta)
# Out[42]: 0.4421003717472119
#
# In [43]: ak.sum(ak.any(abs(jet_eta) < 1.4, axis=1)) / len(jet_eta)
# Out[43]: 0.9841
# mask = ak.any(abs(nu_jeta) < 1.4, axis=1)


x_ticks = np.arange(1, 8, 1)

fig1 = plt.figure()
plt.scatter(x_ticks, ct900_counts, label="ct 900", marker="x")
plt.scatter(x_ticks, ct10000_counts, label="ct 10000", marker="x")
plt.scatter(x_ticks, nu_counts, label="Nu baseline", marker="x")
plt.xlabel("HCal depth")
plt.ylabel("Accumulated counts/total events")
plt.title("Counts per HCal depth")
plt.legend()
# plt.yscale("log")
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/hcal_depth_plot.png")

# plot timing
fig2 = plt.figure()
plt.scatter(x_ticks, ct900_time, label="ct 900", marker="x")
plt.scatter(x_ticks, ct10000_time, label="ct 10000", marker="x")
plt.scatter(x_ticks, nu_time, label="Nu baseline", marker="x")
plt.xlabel("HCal TP timing")
plt.ylabel("Accumulated counts/total events")
plt.title("Counts per HCal timing")
plt.legend()
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/hcal_timing_plot.png")

# time vs depth
fig3 = plt.figure()
plt.scatter(ct900_counts, ct900_time, label="ct 900", marker="x")
plt.scatter(ct10000_counts, ct10000_time, label="ct 10000", marker="x")
plt.scatter(nu_counts, nu_time, label="Nu baseline", marker="x")
plt.ylabel("HCal TP timing")
plt.xlabel("HCAL depth")
plt.title("HCal depth vs. timing")
plt.legend()
plt.savefig(
    "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/hcal_depth_timing_ct900.png"
)

fig2 = plt.figure()
plt.scatter(x_ticks, ct900_counts / np.sum(ct900_counts), label="ct 900", marker="x")
plt.scatter(
    x_ticks, ct10000_counts / np.sum(ct10000_counts), label="ct 10000", marker="x"
)
plt.scatter(x_ticks, nu_counts / np.sum(nu_counts), label="Nu baseline", marker="x")
plt.xlabel("HCal depth")
plt.ylabel("Accumulated counts [Normed]")
plt.title("Normed counts per HCal depth")
plt.legend()
plt.savefig(
    "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/hcal_depth_plot_normed.png"
)

from IPython import embed

embed()
a = 1
