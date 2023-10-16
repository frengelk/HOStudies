#!/usr/bin/env python
import matplotlib

matplotlib.use("Agg")

#import ROOT as ROOT
import uproot as up
import matplotlib.pyplot as plt
import awkward as ak

mass_point = "125"

# fname_xrd = 'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-mass_point_MFF-12_CTau-900mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/181F649B-EDBE-2046-8CC5-AA23D78CE377.root'
# fname_xrd = "root://xrootd-cms.infn.it//store/user/lwiens/Run2018D_SingleMuon_ZMu_PromptReco_v1_L1Ntuple/SingleMuon/Run2018D_ZMu_PromptReco_v1_L1Ntuple/210614_124005/0006/L1Ntuple_6327.root"
fname = "./crab_hcal_testLLP_def/results/L1Ntuple_1.root"
fname_mass_point = "./crab_hcal_testLLP-{}_def/results/L1Ntuple_1.root".format(
    mass_point
)
fname_nu = "crab_hcal_testLLP_nu_def/results/L1Ntuple_10.root"
fname_test = "/afs/cern.ch/user/f/frengelk/Code/CMSSW_11_0_2/src/L1Ntuple.root"
fname_mu = "crab_projects/crab_SingleMuon/results/L1Ntuple_465.root"
fname_350 = "crab_hcal_testLLP-350_CT10000_def/results/L1Ntuple_2.root"
fname_125 = "crab_hcal_testLLP-125_CT900_def/results/L1Ntuple_2.root"
# "crab_hcal_LLP125_retry_def/results/L1Ntuple_1.root"
fname_new = "crab_hcal_testLLP-MH350_MFF160_CT10000_with_new_MC_def/results/L1Ntuple_1.root"
fname_local = "/nfs/dust/cms/user/frengelk/examples/EPR/MH350_MFF160_CT10000_L1Ntuple_1.root"

file = up.open(fname_new)
file_local = up.open(fname_local)
# file_xrd = up.open(fname_xrd)
# file_mu = up.open(fname_mu)
# file_test = up.open(fname_test)
#file_350 = up.open(fname_350)
#file_125 = up.open(fname_125)

#muon = file["l1MuonRecoTree"]["Muon2RecoTree"]["Muon"]
#time = file["l1EventTree"]["L1EventTree"]["Event"]["time"].array()

branch_local = file_local["l1CaloTowerEmuTree"]["L1CaloTowerTree"]["CaloTP"]
hcalTPdepth7 = branch_local['hcalTPDepth7'].array()

from IPython import embed;embed()

print(file_125["l1GeneratorTree;1"]["L1GenTree"]["Generator"].keys())
for k in file_350["l1GeneratorTree"]["L1GenTree"]["Generator"].keys():
    print(k)
    if "part" in k:
        continue
    if "jet" in k:
        h_125 = file_125["l1GeneratorTree;1"]["L1GenTree"]["Generator"][k].array()[:, 0]
        h_350 = file_350["l1GeneratorTree;1"]["L1GenTree"]["Generator"][k].array()[:, 0]
    else:
        h_125 = file_125["l1GeneratorTree;1"]["L1GenTree"]["Generator"][k].array()
        h_350 = file_350["l1GeneratorTree;1"]["L1GenTree"]["Generator"][k].array()

    fig = plt.figure()
    plt.hist(
        h_125,
        label="MH 125, total yields={}".format(h_125.sum()),
        density=True,
        histtype="step",
    )
    plt.hist(
        h_350,
        label="MH 350, total yields={}".format(h_350.sum()),
        density=True,
        histtype="step",
    )
    # plt.title('Time per events')
    plt.xlabel(k)
    plt.ylabel("Number of events")
    plt.legend()
    plt.savefig("mass_comp_-" + k + ".png")
    plt.clf()
