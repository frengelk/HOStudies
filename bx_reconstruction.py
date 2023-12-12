#!/nfs/dust/cms/user/frengelk/Anaconda/envs/EPR_env/bin/python

import os
import uproot as up
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import math
import pickle

# from tqdm import tqdm

plt.rcParams.update({"font.size": 16})
import mplhep as hep


def file_list(direc, filelist):
    for root, dirs, files in os.walk(direc):
        for filename in files:
            if filename.endswith(".root"):
                filelist.append(direc + "/results/" + filename)
    # return filelist


# DTTP conversion functions
def DttpPhiToCmsPhi(phi, dttpSection):
    # dttpSection must be [0, 11]
    cmsPhi = phi / 4096.0 + math.pi / 6 * dttpSection
    # from IPython import embed;embed()
    # cmsphi if  cmsPhi > math.pi else cmsPhi - 2 * math.pi
    return ak.where((cmsPhi > math.pi), cmsPhi, cmsPhi - 2 * math.pi)


def CmsPhiToHoIPhi(cmsPhi):
    # [-pi, pi] to [1, 72]
    cmsPhi = ak.where(cmsPhi <= 0, cmsPhi, cmsPhi + 2 * math.pi)
    # cmsPhi = (cmsPhi <= 0) ? cmsPhi + 2 * math.pi : cmsPhi;
    iPhi = cmsPhi * 72 / 2 / math.pi + 0.5
    return ak.values_astype(iPhi, "int64")


if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option(
        "--dir",
        dest="crabDir",
        default="crab_hcal_testLLP-350_CT10000_with_new_reco_real_def",
        help="crabDir to process",
    )
    (options, args) = parser.parse_args()
    crabDir = options.crabDir

    loc = "/nfs/dust/cms/user/frengelk/Code/cmssw/CMSSW_11_0_2/src/HcalTrigger/Validation/scripts/"

    files_crabDir = []
    # file_list(loc + crabDir, files_crabDir)

    # print(files_crabDir)

    # fname = files_crabDir[0]

    fname_luc = "Lucas_L1Ntuple_112-muon_reco.root"

    file = up.open(fname_luc)

    with open("muon_reco.rootdict_real.pkl", "rb") as f:
        arr_muon = pickle.load(f)

    total_counts_adc = np.zeros(10)
    total_counts_fc = np.zeros(10)

    adc = arr_muon["L1HO/hcalQIESampleAdc"]
    sample = arr_muon["L1HO/hcalQIESample"]
    fc = arr_muon["L1HO/QIESampleFc"]

    print(len(sample))

    from IPython import embed

    # embed()

    arr_adc = np.zeros(10)
    arr_fc = np.zeros(10)

    for j, QIEsample in enumerate(sample):
        #for i, sam in enumerate(QIEsample):

        #if i == 0 or i == len(QIEsample)-1:
        #    print("skipped", i)
        #    continue

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

    #embed()

    plt.figure(figsize=(18, 9))
    plt.hist(
        np.arange(-0.5, 10.5, 1.1),
        weights=total_counts_adc,
        bins=10,
        label="ADC per BX",
        histtype="step",
        linewidth=1.5,
    )
    # plt.hist(
        # np.arange(-0.5, 9.5, 1),
        # weights=total_counts_fc,
        # bins=10,
        # label="FC per BX",
        # histtype="step",
        # linewidth=1.5,
    # )
    hep.style.use("CMS")
    hep.cms.label(llabel="Work in progress", loc=0)
    plt.xticks(np.arange(0, 11, 1))
    plt.xlabel("BX ID")
    plt.ylabel("Weighted Sum per time slice")
    plt.legend(
        loc="upper right",
        #bbox_to_anchor=(1, 1),
        #borderaxespad=0,
    )
    plt.savefig("reco_plots/muon_plots/weighted_ADC_FC_QIEsampleHist.png")
    plt.close()
    """
    index = [
        ("l1EventTree", "L1EventTree", "Event"),
        ("l1UpgradeTree", "L1UpgradeTree", "L1Upgrade"),
        ("l1UpgradeEmuTree", "L1UpgradeTree", "L1Upgrade"),
        ("l1UpgradeTfMuonEmuTree", "L1UpgradeTfMuonTree", "L1UpgradeBmtfMuon"),
        ("l1UpgradeTfMuonEmuTree", "L1UpgradeTfMuonTree", "L1UpgradeBmtfInputs"),
        ("l1GeneratorTree", "L1GenTree", "Generator"),
        ("l1HOTree", "L1HOTree", "L1HO"),
        ("l1CaloTowerEmuTree", "L1CaloTowerTree", "CaloTP"),
        ("l1CaloTowerEmuTree", "L1CaloTowerTree", "L1CaloCluster"),
        ("l1CaloTowerEmuTree", "L1CaloTowerTree", "L1CaloTower"),
    ]

    # bmtf_input = arr_125['L1UpgradeBmtfInputs/thBx']
    # bmtf_size = arr_125['L1UpgradeBmtfInputs/phSize']

    HO_tree = file["l1HOTree/L1HOTree/L1HO"]
    bmtf_muons = file["l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree/L1UpgradeBmtfMuon"]
    bmtf_input = file["l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree/L1UpgradeBmtfInputs"]
    muon_reco = file["l1MuonRecoTree/Muon2RecoTree/Muon"]
    bmtf_muon_bx = ak.flatten(bmtf_muons["tfMuonBx"].array())
    wrong_bmtf_muons = bmtf_muons["tfMuonBx"].array() != 0
    event_time = file["l1EventTree/L1EventTree/Event/time"].array()

    counts = np.histogram(bmtf_muon_bx, bins=5)[0]
    tot = np.sum(counts)
    correct_bx = (counts / tot)[2]

    # dttp prefix -> The BMTF input variables, this is the TwinMux output
    dttpBx = bmtf_input["phBx"]
    dttpWheel = bmtf_input["phWh"]
    dttpSection = bmtf_input["phSe"]
    dttpStation = bmtf_input["phSt"]
    dttpQualityCode = bmtf_input["phCode"]
    dttpTs2Tag = bmtf_input["phTs2Tag"]
    dttpPhi = bmtf_input["phAng"]
    dttpPhiB = bmtf_input["phBandAng"]

    cmsPhi = DttpPhiToCmsPhi(dttpPhi.array(), dttpSection.array())
    iPhi = CmsPhiToHoIPhi(cmsPhi)

    fig = plt.figure()
    plt.hist2d(adc[0], sample[0], bins=20, label="reco muons")
    plt.xlabel("adc")
    plt.ylabel("sample")
    # plt.legend()
    plt.savefig(loc + "reco_plots/muon_plots/QIEsample_2D_0.png")

    plt.clf()
    plt.close()

    plt.figure(figsize=(12, 9))
    plt.hist(adc[0], bins=20, label="adc")
    hep.style.use("CMS")
    hep.cms.label(
    llabel="Work in progress",
    loc=0)
    plt.xlabel("Adc")
    plt.ylabel("Counts")
    plt.legend(shadow=True)
    plt.savefig("reco_plots/muon_plots/QIEsampleADC_0.png")
    plt.yscale("log")
    plt.savefig("reco_plots/muon_plots/QIEsampleADC_0_log.png")
    plt.close()


    plt.figure(figsize=(12, 9))
    plt.hist(adc[0], bins=20, label="Simple")
    hep.style.use("CMS")
    hep.cms.label(
    llabel="Work in progress",
    loc=0)
    plt.xlabel("Simple")
    plt.ylabel("Counts")
    plt.legend(shadow=True)
    plt.savefig("reco_plots/muon_plots/QIEsampleSimple_0.png")
    plt.yscale("log")
    plt.savefig("reco_plots/muon_plots/QIEsampleSimple_0_log.png")
    plt.close()
    """
    """
    plt.figure(figsize=(12, 9))
    plt.hist(bmtf_muon_bx, bins=5, label="correctly classified:{}".format(np.round(correct_bx,3)))
    hep.style.use("CMS")
    hep.cms.label(
    llabel="Work in progress",
    loc=0)
    plt.xlabel("BX ID")
    plt.ylabel("Counts")
    plt.legend(shadow=True)
    plt.savefig("reco_plots/bx_plots/bx_hist-reco_muons18.png")
    plt.yscale("log")
    plt.savefig("reco_plots/bx_plots/bx_hist-reco_muons18_log.png")
    plt.close()
    """

    # egEt_cut=(ak.mean(file[index[2][0]][index[2][1]][index[2][2]]['egEt'].array(), axis=1)>4)
    # hcal_TP1_timing_cut=(ak.mean(file['l1CaloTowerEmuTree']['L1CaloTowerTree']['CaloTP']['hcalTPtiming1'].array(),axis=1)>0)
    # hcal_TP1_depth_cut=(ak.mean(file['l1CaloTowerEmuTree']['L1CaloTowerTree']['CaloTP']['hcalTPDepth1'].array(),axis=1)>0.1)
    # eg_hoe_cut=(ak.sum(file['l1UpgradeEmuTree/L1UpgradeTree/L1Upgrade']['egTowerHoE'].array(), axis=1)>10)
    # hcal_TP_cut=(file['l1CaloTowerEmuTree/L1CaloTowerTree/CaloTP']['nHCALTP'].array()>900)
#
# run=file['l1EventTree/L1EventTree/Event']['run'].array()
#
# cutflow=[
# ak.sum(run),
# ak.sum(run[egEt_cut]),
# ak.sum(run[egEt_cut&hcal_TP1_timing_cut]),
# ak.sum(run[egEt_cut&hcal_TP1_timing_cut&hcal_TP1_depth_cut]),
# ak.sum(run[egEt_cut&hcal_TP1_timing_cut&hcal_TP1_depth_cut&eg_hoe_cut]),
# ak.sum(run[egEt_cut&hcal_TP1_timing_cut&hcal_TP1_depth_cut&eg_hoe_cut&hcal_TP_cut]),
# #ak.sum(run),
# ]
#
# cuts=['nominal', 'egEt', 'hcalTPtiming', 'hcalTPDepth1', 'egTowerHoE','nHCALTP']
#
# plt.figure(figsize=(12, 9))
# plt.plot(np.arange(0,6,1), cutflow,drawstyle='steps-mid')
# plt.xticks(ticks=np.arange(0,6,1), labels=cuts, )
# hep.style.use("CMS")
# hep.cms.label(
# llabel="Work in progress",
# loc=0)
# plt.xlabel("Cut")
# plt.ylabel("Counts")
# #plt.legend(loc="upper center", bbox_to_anchor=(1.1, 0.9), ncol=1)
# plt.savefig("reco_plots/cutflow/cutflow-" + crabDir + ".png")
# plt.close()
#
# print("\n plotted \n")
