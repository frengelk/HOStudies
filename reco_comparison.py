#!/nfs/dust/cms/user/frengelk/Anaconda/envs/EPR_env/bin/python

import os
import uproot as up
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import pickle


def file_list(direc, filelist):
    for root, dirs, files in os.walk(direc):
        for filename in files:
            if filename.endswith(".root"):
                filelist.append(direc + "/results/" + filename)
    # return filelist


if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option(
        "--dir",
        dest="crabDir",
        default="crab_hcal_testLLP-125_CT900_with_reco_def",
        help="crabDir to process",
    )
    (options, args) = parser.parse_args()
    crabDir = options.crabDir

    loc = "/nfs/dust/cms/user/frengelk/Code/cmssw/CMSSW_11_0_2/src/HcalTrigger/Validation/scripts/"

    # str_350 = "crab_hcal_testLLP-350_CT10000_with_reco_def"
    # str_125 = "crab_hcal_testLLP-125_CT900_with_reco_def"
    # str_nu = "crab_hcal_testLLP-Nu_Pt20_with_reco_def"
    # str_QCD = "crab_hcal_testLLP-QCD_Pt3000_with_reco_def"
    # files_350, files_125, files_nu, files_QCD = [], [], [], []
    # file_list(loc + str_350, files_350)
    # file_list(loc + str_125, files_125)
    # file_list(loc + str_nu, files_nu)
    # file_list(loc + str_QCD, files_QCD)

    files_crabDir = []
    file_list(loc + crabDir, files_crabDir)

    # files_350 = [
    # "crab_hcal_testLLP-350_CT10000_with_reco_def/results/L1Ntuple_1.root",
    # "crab_hcal_testLLP-350_CT10000_with_reco_def/results/L1Ntuple_2.root",
    # "crab_hcal_testLLP-350_CT10000_with_reco_def/results/L1Ntuple_3.root",
    # "crab_hcal_testLLP-350_CT10000_with_reco_def/results/L1Ntuple_4.root",
    # ]

    # files_125 = [
    # "crab_hcal_testLLP-125_CT900_with_reco_def/results/L1Ntuple_1.root",
    # "crab_hcal_testLLP-125_CT900_with_reco_def/results/L1Ntuple_2.root",
    # "crab_hcal_testLLP-125_CT900_with_reco_def/results/L1Ntuple_3.root",
    # ]

    # hacky way for getting keys
    """
    fname = files_crabDir[0]

    file = up.open(fname)
    """
    file = up.open(options.crabDir)

    # ev = file["l1EventTree"]["L1EventTree"]["Event"].keys()
    # upg = file["l1UpgradeTree"]["L1UpgradeTree"]["L1Upgrade"].keys()
    # Gen = file["l1GeneratorTree/L1GenTree"]["Generator"].keys()

    index = [
        ("l1EventTree", "L1EventTree", "Event"),
        ("l1UpgradeTree", "L1UpgradeTree", "L1Upgrade"),
        ("l1UpgradeEmuTree", "L1UpgradeTree", "L1Upgrade"),
        ("l1UpgradeTfMuonEmuTree", "L1UpgradeTfMuonTree", "L1UpgradeBmtfMuon"),
        ("l1UpgradeTfMuonEmuTree", "L1UpgradeTfMuonTree", "L1UpgradeOmtfMuon"),
        ("l1UpgradeTfMuonEmuTree", "L1UpgradeTfMuonTree", "L1UpgradeEmtfMuon"),
        ("l1UpgradeTfMuonEmuTree", "L1UpgradeTfMuonTree", "L1UpgradeBmtfInputs"),
        ("l1UpgradeTfMuonEmuTree", "L1UpgradeTfMuonTree", "L1UpgradeKBmtfMuon"),
        # ("l1GeneratorTree", "L1GenTree", "Generator"),
        ("l1HOTree", "L1HOTree", "L1HO"),
        ("l1CaloTowerEmuTree", "L1CaloTowerTree", "CaloTP"),
        ("l1CaloTowerEmuTree", "L1CaloTowerTree", "L1CaloCluster"),
        ("l1CaloTowerEmuTree", "L1CaloTowerTree", "L1CaloTower"),
    ]

    # from IPython import embed; embed()
    # for files in [files_350, files_125]:
    # arr_350, arr_125, arr_nu, arr_QCD = [], [], [], []
    arr_crabDict = {}
    key_names = []

    # print("Generator Tree")
    # for coll, files in (arr_350, files_350), (arr_125, files_125):
    # for k in file["l1GeneratorTree/L1GenTree"]["Generator"].keys():
    # if "weight" in k or "pthat" in k:
    # continue
    # for fname in files:
    #
    # file = up.open(fname)
    # print(fname, k)
    #
    # # if "hcalDetIdIEta" in k or "hcalDetIdIPhi" in k or "hcalQIESample" in k:
    #
    # # ar = ak.mean(file["l1HOTree"]["L1HOTree"]["L1HO"][k].array(), axis=1)
    # # else:
    # # ar = file["l1HOTree"]["L1HOTree"]["L1HO"][k].array()
    #
    # if "L1Ntuple_1" in fname:
    # arr = file["l1GeneratorTree/L1GenTree"]["Generator"][k].array()
    # # ak.zeros_like(ar)
    #
    # else:
    # # from IPython import embed;embed()
    # arr = ak.concatenate(
    # (
    # arr,
    # file["l1UpgradeTree"]["L1UpgradeTree"]["L1Upgrade"][k].array(),
    # ),
    # axis=0,
    # )
    # # ak.sum(arr, file["l1HOTree;1"]["L1HOTree;1"]["L1HO"][k].array())
    # coll.append(arr)
    # key_names.extend(file["l1GeneratorTree/L1GenTree"]["Generator"].keys())

    # from IPython import embed;embed()

    for ind in index:
        print(ind[0])
        # for coll, files in (
        # (arr_crabDir, files_crabDir)
        # (arr_350, files_350),
        # (arr_125, files_125),
        # (arr_nu, files_nu),
        # (arr_QCD, files_QCD),
        # ):
        for k in file[ind[0]][ind[1]][ind[2]].keys():
            if k == "tfMuonDecodedTrAdd":
                continue
            print(k)
            for fname in [options.crabDir]:
                # for fname in files_crabDir:

                # with up.open(fname) as file:

                file = up.open(fname)
                # print(fname, k)

                # if "hcalDetIdIEta" in k or "hcalDetIdIPhi" in k or "hcalQIESample" in k:

                # ar = ak.mean(file["l1HOTree"]["L1HOTree"]["L1HO"][k].array(), axis=1)
                # else:
                # ar = file["l1HOTree"]["L1HOTree"]["L1HO"][k].array()

                # if "L1Ntuple_1" in fname:
                arr = file[ind[0]][ind[1]][ind[2]][k].array()
                """
                if fname == files_crabDir[0]:
                    arr = file[ind[0]][ind[1]][ind[2]][k].array()
                    # ak.zeros_like(ar)

                else:
                    # from IPython import embed;embed()
                    arr = ak.concatenate(
                        (arr, file[ind[0]][ind[1]][ind[2]][k].array()), axis=0
                    )
                # ak.sum(arr, file["l1HOTree;1"]["L1HOTree;1"]["L1HO"][k].array())
                """
            comb_key = ind[2] + "/" + k
            arr_crabDict.update({comb_key: arr})
            del arr
            # key_names.extend(file[ind[0]][ind[1]][ind[2]].keys())
            key_names.append(comb_key)

    # with open(loc + "key_names_HO_real.pkl", "wb") as f:
    #    pickle.dump(key_names, f)

    with open(loc + crabDir.split("-")[1] + "dict_real.pkl", "wb") as f:
        pickle.dump(arr_crabDict, f)

    # with open(loc + "arr_125_1.pkl", "wb") as f:
    # pickle.dump(arr_125, f)

    # with open(loc + "arr_350_1.pkl", "wb") as f:
    # pickle.dump(arr_350, f)

    # with open(loc + "arr_QCD_1.pkl", "wb") as f:
    # pickle.dump(arr_QCD, f)

    # with open(loc + "arr_nu_1.pkl", "wb") as f:
    # pickle.dump(arr_nu, f)
"""

print("Events")
for coll, files in (arr_350, files_350), (arr_125, files_125):
    for k in file["l1EventTree"]["L1EventTree"]["Event"].keys():
        # if k=="nHcalDetIds":
        #    continue
        for fname in files:

            file = up.open(fname)
            print(fname, k)

            # if "hcalDetIdIEta" in k or "hcalDetIdIPhi" in k or "hcalQIESample" in k:

            # ar = ak.mean(file["l1HOTree"]["L1HOTree"]["L1HO"][k].array(), axis=1)
            # else:
            # ar = file["l1HOTree"]["L1HOTree"]["L1HO"][k].array()

            if "L1Ntuple_1" in fname:
                arr = file["l1EventTree"]["L1EventTree"]["Event"][k].array()
                # ak.zeros_like(ar)

            else:
                # from IPython import embed;embed()
                arr = ak.concatenate(
                    (arr, file["l1EventTree"]["L1EventTree"]["Event"][k].array()),
                    axis=0,
                )
            # ak.sum(arr, file["l1HOTree;1"]["L1HOTree;1"]["L1HO"][k].array())
        coll.append(arr)
key_names.extend(file["l1EventTree"]["L1EventTree"]["Event"].keys())


print("Upgrade Tree")
for coll, files in (arr_350, files_350), (arr_125, files_125):
    for k in file["l1UpgradeTree"]["L1UpgradeTree"]["L1Upgrade"].keys():
        # if k=="nHcalDetIds":
        #    continue
        for fname in files:

            file = up.open(fname)
            print(fname, k)

            # if "hcalDetIdIEta" in k or "hcalDetIdIPhi" in k or "hcalQIESample" in k:

            # ar = ak.mean(file["l1HOTree"]["L1HOTree"]["L1HO"][k].array(), axis=1)
            # else:
            # ar = file["l1HOTree"]["L1HOTree"]["L1HO"][k].array()

            if "L1Ntuple_1" in fname:
                arr = file["l1UpgradeTree"]["L1UpgradeTree"]["L1Upgrade"][k].array()
                # ak.zeros_like(ar)

            else:
                # from IPython import embed;embed()
                arr = ak.concatenate(
                    (
                        arr,
                        file["l1UpgradeTree"]["L1UpgradeTree"]["L1Upgrade"][k].array(),
                    ),
                    axis=0,
                )
            # ak.sum(arr, file["l1HOTree;1"]["L1HOTree;1"]["L1HO"][k].array())
        coll.append(arr)
key_names.extend(file["l1UpgradeTree"]["L1UpgradeTree"]["L1Upgrade"].keys())



print("#####Plotting#####")
for i, k in enumerate(key_names):
    print(k)
    # define multi axes variales
    if "hcalDetIdIEta" in k or "hcalDetIdIPhi" in k or "hcalQIESample" in k:  #
        #
        # from IPython import embed;embed()
        # # embed()
        #continue
        h_125 = ak.mean(arr_125[i], axis=1)
        h_350 = ak.mean(arr_350[i], axis=1)
    # define skips
    elif "orbit" in k:
        continue
    #everything else
    else:
        h_125 = ak.to_numpy(arr_125[i])
        h_350 = ak.to_numpy(arr_350[i])

    fig = plt.figure()  # arr_125[i]
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
    # plt.title('Time per events')
    plt.xlabel(k)
    plt.ylabel("Number of events")
    plt.legend()
    plt.savefig("reco_plots/mass_comp_-" + k + ".png")
    plt.clf()

from IPython import embed

# embed()

"""
