#!/usr/bin/env python
import matplotlib as mpl

mpl.use("Agg")

#import ROOT as ROOT
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
    return delta_R  # ak.where(delta_R < 0.4, delta_R, 999)


# function for dR matching with multiple particles
def HO_matching(eta, phi, eta_H=0, phi_H=0):
    # assume consistent structure
    num_gen = len(eta_H[0])
    # now do a matching for every interesting particle per event
    masks = []
    min_vals = []
    for i in range(num_gen):
        dR_HO_b = dR(eta, phi, eta_H[:, i], phi_H[:, i])
        min_val = ak.min(dR_HO_b, axis=1)  # , keepdims=True)
        # min_mask = min_val_04 == dR_HO_b
        # FIXME something is going wrong here, ma[1] is just a plain False, not a sub array
        # min_mask = ak.where(min_val < 0.4, min_val == dR_HO_b, False)
        min_mask = ak.where(dR_HO_b, min_val == dR_HO_b, False)
        masks.append(min_mask)
        # min_vals.append(min_val_04)
    ma = masks[0]
    for m in masks[1:]:
        ma = ma + m
    return ma


def weighted_fc(ho_sample):

    """
    for j, QIEsample in enumerate(ho_sample):
        #for i, sam in enumerate(QIEsample):

        #if i == 0 or i == len(QIEsample)-1:
        #    print("skipped", i)
        #    continue

        from IPython import embed; embed()

        ind_max = np.argmax(adc[j])

        weighted_adc = (
            adc[j][ind_max - 1] * QIEsample[ind_max - 1]
            + adc[j][ind_max] * QIEsample[ind_max]
            + adc[j][ind_max + 1] * QIEsample[ind_max + 1]
        ) / (adc[j][ind_max - 1] + adc[j][ind_max] + adc[j][ind_max + 1])

        arr_adc[QIEsample[ind_max]] += weighted_adc
        #arr_fc[sam] += fc[j][i]
        # arr_counts[sam]+=1

    total_counts_adc += arr_adc
    total_counts_fc += arr_fc

    fc_ws = []

    for fc in ho_sample.QIESampleFc:
        for i in range(1, len(fc) - 1):
             fc_w = (fc[i - 1]/2 + fc[i] + fc[i + 1]/2) / 2
             fc_ws.append(fc_w)
             from IPython import embed; embed()
    return fc_ws
    """
    fc_ws = []

    for fc in ho_sample.QIESampleFc:
        slices = int(len(fc)/10)
        #print(slices)
        for i in range(slices):
            #print(i * 10 + 4)
            fc_w = (fc[i * 10 + 4] + fc[i * 10 + 5] + fc[i * 10 + 6] ) / 3
            fc_ws.append(fc_w)
        #print(slices, i * 10 + 4)
    return fc_ws


def weighted_time(ho_sample):
    # maybe check
    # https://indico.cern.ch/event/799070/contributions/3320191/attachments/1799399/2934434/DN-17-047_auto.pdf,
    # eq. (1)

    # sum over all time sclices in product with fc devided by fc sum

    time_ws = []

    for index, timefc in enumerate(ho_sample.QIESampleFc):
        # always TS from to 9 -> 10er split
        # num_tp = int(len(timefc)/10)
        time_qie = ho_sample.hcalQIESample[index]
        num_tp = int(len(timefc))
        for i in range(num_tp - 1):  # to avoid failure at last i

            numen, denumen = 0, 0
            it_counter = 0

            # for j in range(i*10, (i+1)*10, 1):
            for j in range(i - 1, (i + 2), 1):
                # TS +1 so you eject 0
                it_counter += 1

                # numen += it_counter*timefc[j]
                numen += time_qie[j] * timefc[j]
                denumen += timefc[j]

            time_ws.append(numen / denumen)

    return time_ws


ct900_dir = "crab_hcal_testLLP-MH125_MFF_12_CT900_to4mu_def/results"
ct10000_dir = "crab_hcal_testLLP-MH350_MFF160_CT10000_with_new_MC_def/results"
nu_dir = "crab_hcal_testLLP-Nu_Pt20_with_reco_def/results"
qcd_dir = "crab_hcal_testLLP-QCD_Pt3000_with_reco__redone/results"
# "crab_hcal_testLLP-QCD_Pt3000_with_reco_AODEMUGEN_MC"

# ct900_files = search_files(ct900_dir)[:1]
ct10000_files = search_files(ct10000_dir)
qcd_files = search_files(qcd_dir)[:3]

# let's go to gen level
# H0/H02 pdgid 35
gen_str = "l1GeneratorTree/L1GenTree/Generator/"
genId_str = "l1GeneratorTree/L1GenTree/Generator/partId"
# ct900_H0mask = abs(merge_arrays(ct900_files, genId_str)) == 13
ct10000_H0mask = merge_arrays(ct10000_files, genId_str) == 35

geneta_str = "l1GeneratorTree/L1GenTree/Generator/partEta"

# ct900_geneta = merge_arrays(ct900_files, geneta_str)  # [ct900_H0mask]
ct10000_geneta = merge_arrays(ct10000_files, geneta_str)  # [ct10000_H0mask]

genphi_str = "l1GeneratorTree/L1GenTree/Generator/partPhi"

# ct900_genphi = merge_arrays(ct900_files, genphi_str)[ct900_H0mask]
ct10000_genphi = merge_arrays(ct10000_files, genphi_str)  # [ct10000_H0mask]

gen_parent_str = "l1GeneratorTree/L1GenTree/Generator/partParent"


# unpacking hcal eta, calculating distance to gen-level daughters
# finding the TP wth the smallest distance, generating a mask to find hcal TPs
# ho
ho_str = "l1HOTree/L1HOTree/L1HO/"
# ct900_ho_eta = IEtaToCMSEta(merge_arrays(ct900_files, ho_str + "hcalDetIdIEta"))
# ct900_ho_phi = IPhiToCMSPhi(merge_arrays(ct900_files, ho_str + "hcalDetIdIPhi"))

# finding the b daughters in the dataset
gen_id = merge_arrays(ct10000_files, genId_str)
par = merge_arrays(ct10000_files, gen_parent_str)
daughter_id = ak.flatten(gen_id[par == 35])[0]
mother_id = 35
b = gen_id[par == daughter_id]

LLP_daughters = par == daughter_id

"""
dR(ct900_ho_eta, ct900_ho_phi, ct900_geneta_parent, ct900_genphi_parent)
fails
Maybe need to do some cross product...
"""

# ho_parent_mask_0 = ak.argmin(
# dR(
# ct900_ho_eta, ct900_ho_phi, ct900_geneta_parent[:, 0], ct900_genphi_parent[:, 0]
# ),
# axis=1,
# keepdims=True,
# )
# ho_parent_mask_1 = ak.argmin(
# dR(
# ct900_ho_eta, ct900_ho_phi, ct900_geneta_parent[:, 1], ct900_genphi_parent[:, 1]
# ),
# axis=1,
# keepdims=True,
# )
# ho_parent_mask = ak.concatenate((ho_parent_mask_0, ho_parent_mask_1), axis=1)
# ho_parent_matched_eta = ct900_ho_eta[ho_parent_mask]

# redo ct900 with ct10000
gen_parent_str = "l1GeneratorTree/L1GenTree/Generator/partParent"

ct10000_gen_parentmask = merge_arrays(ct10000_files, gen_parent_str) == 35
ct10000_geneta = merge_arrays(
    ct10000_files, geneta_str
)  # [    ct10000_gen_parentmask]  # [:,0]
ct10000_genphi = merge_arrays(
    ct10000_files, genphi_str
)  # [    ct10000_gen_parentmask]  # [:,0]

# unpacking hcal eta, calculating distance to gen-level daughters
# finding the TP wth the smallest distance, generating a mask to find hcal TPs
# ho
ho_str = "l1HOTree/L1HOTree/L1HO/"

# QIE Sample = 4!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ct10000_ho = merge_arrays(ct10000_files, ho_str)

qie_sample_4_10000 = ct10000_ho.hcalQIESample == 4  # FIXME

ct10000_ho_eta = IEtaToCMSEta(ct10000_ho.hcalDetIdIEta[qie_sample_4_10000])
ct10000_ho_phi = IPhiToCMSPhi(ct10000_ho.hcalDetIdIPhi[qie_sample_4_10000])

# ho cut eta < 1.4
ct10000_geneta = ak.where(abs(ct10000_geneta) < 1.4, ct10000_geneta, -999)
b_eta = ct10000_geneta[LLP_daughters]
b_phi = ct10000_genphi[LLP_daughters]

print("HO matching")
mask = HO_matching(ct10000_ho_eta, ct10000_ho_phi, b_eta, b_phi)

ho_10000_parent_matched_eta = ct10000_ho_eta[mask]
# ak.where(mask == True, ct10000_ho_eta, np.nan)
# ho_10000_parent_matched_eta = ak.mask(ct10000_ho_eta, mask) #ct10000_ho_eta[mask]

# qcd
# somehow I am able to create an array with columns
qcd_gen = merge_arrays(qcd_files, gen_str)
qcd_ho = merge_arrays(qcd_files, ho_str)

# eta 1.4 for HO volume
qcd_gen_eta_cut = abs(qcd_gen.jetEta) < 1.4
qcd_eta_gencut = ak.where(qcd_gen_eta_cut, qcd_gen.jetEta, 999)

# ho TS 4 cut
qcd_ho_TS4_cut = qcd_ho.hcalQIESample == 4  # FIXME
qcd_ho_eta = IEtaToCMSEta(qcd_ho.hcalDetIdIEta[qcd_ho_TS4_cut])
qcd_ho_phi = IPhiToCMSPhi(qcd_ho.hcalDetIdIPhi[qcd_ho_TS4_cut])

# doing matching to leading jet
qcd_dR = dR(
    qcd_ho_eta,
    qcd_ho_phi,
    qcd_eta_gencut[:, 0],
    qcd_gen.jetPhi[:, 0],
)

qcd_min_dR = ak.argmin(qcd_dR, axis=1, keepdims=True)

qcd_dR_04 = ak.any(qcd_dR < 0.4, axis=1)

qcd_HOTP_eta = qcd_ho_eta[qcd_min_dR][qcd_dR_04]

# comparing quantities
qcd_sample_energy = qcd_ho.SampleEnergy[qcd_min_dR][qcd_dR_04]
qcd_sample_energy = ak.where(qcd_sample_energy < 5, qcd_sample_energy < 5, np.nan)
ct10000_sample_energy = ct10000_ho.SampleEnergy[qie_sample_4_10000][mask]
# ak.where(mask, ct10000_ho.SampleEnergy[qie_sample_4_10000], np.nan)

ct10000_sample_energy = ct10000_sample_energy[ct10000_sample_energy < 5]
# ak.where(ct10000_sample_energy < 5, ct10000_sample_energy, np.nan)

# hcalQIESampleAdc
qcd_qie_adc = qcd_ho.hcalQIESampleAdc[qcd_min_dR][qcd_dR_04]
ct10000_qie_adc = ct10000_ho.hcalQIESampleAdc[qie_sample_4_10000][mask]
# ak.where(mask, ct10000_ho.hcalQIESampleAdc[qie_sample_4_10000], np.nan)

# hcalQIESampleFC
qcd_qie_fc = qcd_ho.QIESampleFc[qcd_min_dR][qcd_dR_04]
ct10000_qie_fc = ct10000_ho.QIESampleFc[qie_sample_4_10000][mask]
# ak.where(mask, ct10000_ho.QIESampleFc[qie_sample_4_10000], np.nan)


# sumQ
qcd_sumQ = qcd_ho.sumQ[ak.any(qcd_min_dR, axis=1)]
ct10000_sumQ = ct10000_ho.sumQ


# from IPython import embed; embed()

# weighted averages considering time/fc
# takes a looooong time
#ct10000_time_fc = weighted_time(ct10000_ho)
#qcd_time_fc = weighted_time(qcd_ho)
qcd_w_fc = weighted_fc(qcd_ho)
ct10000_w_fc = weighted_fc(ct10000_ho)


# here I should do a log plot and somehow fix the binning to maybe 0, 40, 10 or something
bins_adc = np.arange(0, 100, 5)

# plottig in loop:
for (qcd, ct10000, var, x_label, log, bins) in [
    #(qcd_sumQ, ct10000_sumQ, "sumQ", "HO sumQ", "log", np.arange(0, 2000, 100)),
    #(qcd_qie_adc, ct10000_qie_adc, "qieadc", "QIE Sample Adc", "log", bins_adc),
    #(qcd_qie_fc, ct10000_qie_fc, "qiefc", "QIE Fc", "log", np.arange(0, 4000, 100)),
    (
        qcd_w_fc,
        ct10000_w_fc,
        "Test_Errors_456_weighted_fc",
        "Weighted Fc of TS 4|5|6",
        "log",
        np.arange(0, 4000, 100),
    ),
    #(qcd_time_fc, ct10000_time_fc, "weighted_time_fc", "Weighted TS by Fc", "log", np.arange(0, 11, 0.5)),
]:

    plt.figure(figsize=(12, 9))
    hep.style.use("CMS")
    hep.cms.label(
        llabel="Work in progress",
        loc=0,
        fontsize=20,
    )
    plt.hist(
        ak.flatten(qcd, axis=None),
        label="QCD jet matched ",
        histtype="step",
        bins=bins,
        density=True,
        linewidth=2,
    )
    qcd_hist = np.histogram(ak.flatten(qcd, axis=None), bins=bins, density=True)
    plt.hist(
        ak.flatten(ct10000, axis=None),
        # ho_10000_parent_matched_eta,
        label="ct10000: LLP -> b",
        histtype="step",
        bins=bins,
        density=True,
        linewidth=2,
    ),
    # plt.errorbar(
        # qcd_hist[1][:-1] + 50,
        # qcd_hist[0],
        # yerr = np.sqrt(qcd_hist[0]),
        # marker = '.',
        # drawstyle = 'steps-mid'
    # )
    plt.xlabel(x_label)
    plt.ylabel("a.u.")
    # plt.title("Matching between hcal TP and H0 daughters")
    plt.legend(
        ncol=1,
        loc="upper right",
        bbox_to_anchor=(1, 0.95),
        borderaxespad=0,
    )  # loc='upper center', bbox_to_anchor=(.5,1.15))#, fancybox=True, shadow=True)
    if log == "log":
        plt.yscale("log")
        log = "_" + log
    plt.savefig(
        "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/TS4_Qcd_vs_LLP4b_{}_v2.png".format(
            var + log
        )
    )
    plt.savefig(
            "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/TS4_Qcd_vs_LLP4b_{}_v2.pdf".format(
                var + log
            )
        )
    plt.close()

"""
# 2D hist qie vs adc
plt.figure(figsize=(12, 9))
hep.style.use("CMS")
hep.cms.label(
    llabel="Work in progress",
    loc=0,
    fontsize=20,
)
ct_2d = plt.hist2d(
    ak.flatten(ct10000_ho.hcalQIESample),
    #ak.flatten(ct10000_ho.hcalQIESampleAdc),
    ak.flatten(ct10000_ho.QIESampleFc),
    bins=[np.arange(-0.5, 10.5, 1), 14],#np.arange(0,150,10)],
    norm=mpl.colors.LogNorm(),
    density=True,
)
plt.xlabel("TS(QIE Sample)")
plt.xticks(np.arange(0, 10, 1), np.arange(0, 10, 1))
# plt.xlim([-0.5, 10.5])
plt.ylabel("Fc")
# plt.yscale("log")
plt.title("LLP", fontsize=28)
plt.colorbar()
# plt.legend()
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/2D_QIE_Fc_LLP4b.png")

# same for qcd
plt.figure(figsize=(12, 9))
hep.style.use("CMS")
hep.cms.label(
    llabel="Work in progress",
    loc=0,
    fontsize=20,
)
qcd_2d = plt.hist2d(
    ak.flatten(qcd_ho.hcalQIESample),
    #ak.flatten(qcd_ho.hcalQIESampleAdc),
    ak.flatten(qcd_ho.QIESampleFc),
    bins=[np.arange(-0.5, 10.5, 1), 14],#, np.arange(0,150,10)],
    norm=mpl.colors.LogNorm(),
    density=True,
)
plt.xlabel("TS(QIE Sample)")
plt.xticks(np.arange(0, 10, 1), np.arange(0, 10, 1))
plt.ylabel("Fc")
# plt.yscale("log")
plt.colorbar()
plt.title("QCD", fontsize=28)
# plt.legend()
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/2D_QIE_Fc_QCD.png")

# 2d ratio plot
new_2d = ct_2d[0] / qcd_2d[0]
# get rid of inf and nan and trasform to log
new_2d = np.log(np.where(new_2d < 2000, new_2d, 1))
plt.figure(figsize=(16, 9))
hep.style.use("CMS")
hep.cms.label(
    llabel="Work in progress",
    loc=0,
    fontsize=12,
)
im = plt.imshow(
    np.transpose(new_2d),
    origin="lower",
)
plt.title("Ratio QCD/LLP", fontsize=12)
plt.xlabel("TS(QIE Sample)")
plt.xticks(np.arange(0, 10, 1), np.arange(0, 10, 1))
plt.ylabel("Fc")
plt.yticks(np.arange(0, 15, 1), np.arange(0, 150, 10))
plt.colorbar(label="log(ratio)")
plt.savefig(
    "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/test_2D_QIE_Fc_LLPvsQCD.png"
)


checking weighted time fc
In [2]: np.mean(qcd_time_fc)
Out[2]: 4.488208411706821

In [3]: np.std(qcd_time_fc)
Out[3]: 2.1517970651507516

In [4]: np.std(ct10000_time_fc)
Out[4]: 2.138523730699388

In [5]: np.mean(ct10000_time_fc)
Out[5]: 4.493978481503885


from IPython import embed

embed()




# qcd plots
fig_qcd_llp = plt.figure(figsize=(12, 9))
hep.style.use("CMS")
hep.cms.label(
    llabel="Work in progress",
    loc=0,
    fontsize=20,
)
plt.hist(
    ak.flatten(qcd_HOTP_eta, axis=None),
    label="QCD Jets to HO TP",
    histtype="step",
    bins=28,
    density=True,
)
plt.hist(
    ak.flatten(ho_10000_parent_matched_eta, axis=None),
    # ho_10000_parent_matched_eta,
    label="ct10000: HO matched to b quarks from LLP",
    histtype="step",
    bins=28,
    density=True,
)  # , density=True
plt.xlabel("Eta")
plt.ylabel("Counts")
# plt.title("Matching between hcal TP and H0 daughters")
plt.legend(
    ncol=1,
    loc="upper right",
    bbox_to_anchor=(1, 0.95),
    borderaxespad=0,
)  # loc='upper center', bbox_to_anchor=(.5,1.15))#, fancybox=True, shadow=True)
# plt.yscale("log")
plt.savefig("/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/Qcd_vs_LLP4b_.png")

plt.figure(figsize=(12, 9))
hep.style.use("CMS")
hep.cms.label(
    llabel="Work in progress",
    loc=0,
    fontsize=20,
)
plt.hist(
    ak.flatten(qcd_sample_energy, axis=None),
    label="QCD Jets to HO TP",
    histtype="step",
    bins=20,
    density=True,
)
plt.hist(
    ak.flatten(ct10000_sample_energy, axis=None),
    # ho_10000_parent_matched_eta,
    label="ct10000: HO matched to b quarks from LLP",
    histtype="step",
    bins=20,
    density=True,
)  # , density=True
plt.xlabel("Sample Energy")
plt.ylabel("Counts")
# plt.title("Matching between hcal TP and H0 daughters")
plt.legend(
    ncol=1,
    loc="upper right",
    bbox_to_anchor=(1, 0.95),
    borderaxespad=0,
)
plt.savefig(
    "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/Sample_Energy_Qcd_vs_LLP4b_.png"
)

# Qie ADC
plt.figure(figsize=(12, 9))
hep.style.use("CMS")
hep.cms.label(
    llabel="Work in progress",
    loc=0,
    fontsize=20,
)
plt.hist(
    ak.flatten(qcd_qie_adc, axis=None),
    label="QCD Jets to HO TP",
    histtype="step",
    bins=20,
    density=True,
)
plt.hist(
    ak.flatten(ct10000_qie_adc, axis=None),
    # ho_10000_parent_matched_eta,
    label="ct10000: HO matched to b quarks from LLP",
    histtype="step",
    bins=20,
    density=True,
)  # , density=True
plt.xlabel("QIE ADC")
plt.ylabel("Counts")
# plt.title("Matching between hcal TP and H0 daughters")
plt.legend(
    ncol=1,
    loc="upper right",
    bbox_to_anchor=(1, 0.95),
    borderaxespad=0,
)
plt.savefig(
    "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/QIE_ADC_Qcd_vs_LLP4b_.png"
)# comparing hist of H0 and matched
fig_H0_matched = plt.figure(figsize=(12, 9))
hep.style.use("CMS")
hep.cms.label(
    llabel="Work in progress",
    loc=0,
    fontsize=20,
)

#plt.hist(ak.flatten(ct900_ho_eta), label="HO 900", histtype="step", density=True, bins=25)  # , density=True)
hist = plt.hist(ak.flatten(ct10000_ho_eta), label="HO TP ct10000", histtype="step", bins=28, density=True)  # , density=True)
hist_matched = plt.hist(
    ak.flatten(ho_10000_parent_matched_eta, axis=None),
    #ho_10000_parent_matched_eta,
    label="ct10000: HO matched to b quarks from LLP",
    histtype="step", bins=28, density=True) # , density=True
plt.xlabel("Eta")
plt.ylabel("Counts")
# plt.title("Matching between hcal TP and H0 daughters")
plt.legend(
    ncol=1,
    loc="upper right",
    bbox_to_anchor=(1, 0.95),
    borderaxespad=0,
)  # loc='upper center', bbox_to_anchor=(.5,1.15))#, fancybox=True, shadow=True)
#plt.yscale("log")
plt.savefig(
    "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/Eta_matching_4b_15cut_density.png"
)

for key in up.open(ct10000_files[0])[ho_str].keys():
    ho_base = (merge_arrays(ct10000_files, ho_str + key))
    plt.hist(ak.flatten(ho_base, axis=None), label="HO TP ct10000", histtype="step")
    ho_matched = ak.where(mask == True, ho_base, np.nan)
    plt.hist(
        ak.flatten(ho_10000_parent_matched_eta, axis=None),
        label="ct10000: HO matched to b quarks from LLP",
        histtype="step")
    plt.xlabel(key)
    plt.ylabel("Counts")
    plt.legend()
    plt.savefig(
        "/nfs/dust/cms/user/frengelk/examples/EPR/hcal_plots/Eta_matching_4b_{}.png".format(key)
    )
"""
