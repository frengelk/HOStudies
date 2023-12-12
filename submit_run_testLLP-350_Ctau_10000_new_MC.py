from CRABClient.UserUtilities import config# RUN = "testLLP-350_CT10000_with_new_MC"
#RUN = "testLLP-MH350_MFF160_CT10000_with_new_MC"
#RUN = "testLLP-MH-125_MFF-12_CT900_with_new_MC"
#RUN = "testLLP-Nu_Pt20_with_reco"
RUN = "testLLP-QCD_Pt3000_with_reco"
#RUN = "testLLP-MH125_MFF_12_CT900_to4mu"
NEWCONDITIONS = False
OUTPUTSITE = "T2_DE_DESY"
#DATASET = "/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-900mm_TuneCP5_14TeV_pythia8/Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/GEN-SIM-RAW"
#DATASET = "/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-10000mm_TuneCP5_14TeV_pythia8/Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/GEN-SIM-RAW"
#DATASET = "/Neutrino_Pt-2to20_gun/Run3Winter20DRPremixMiniAOD-SNB_110X_mcRun3_2021_realistic_v6-v1/GEN-SIM-RAW"
DATASET = "/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8_HCAL/Run3Winter20DRPremixMiniAOD-packHBTDC_110X_mcRun3_2021_realistic_v6-v1/GEN-SIM-DIGI-RAW"
#DATASET = "/HTo2LongLivedTo4mu_MH-125_MFF-12_CTau-900mm_TuneCP5_14TeV_pythia8/Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/GEN-SIM-RAW"

# template for crab submission.
# submit_jobs.py will insert definitions above
conditionType = "_redone"  # default
# if new L1TriggerObjects conditons have been specified with either
# a new tag or file
if NEWCONDITIONS:
    conditionType = "new_cond"

config = config()

config.General.requestName = "hcal_" + str(RUN) + "_" + conditionType
config.General.transferLogs = True
config.General.transferOutputs = True

# Name of the CMSSW configuration file
# config.JobType.psetName = 'ntuple_maker_' + conditionType + '.py'
config.JobType.psetName = "ntuple_maker.py" #_def
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = "Analysis"
# config.JobType.maxMemoryMB = 2500
config.JobType.outputFiles = ["L1Ntuple.root"]

config.Data.inputDataset = DATASET
config.Data.ignoreLocality = True
config.Data.inputDBS = "global"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 5
config.Data.useParent = False

# This string is used to construct the output dataset name
config.Data.outputDatasetTag = "Hcal" + str(RUN) + "_" + conditionType

# These values only make sense for processing data
#    Select input data based on a lumi mask
# config.Data.lumiMask = LUMIMASK

#    Select input data based on run-ranges
# config.Data.runRange = str(RUN)

# Where the output files will be transmitted to
config.Site.storageSite = OUTPUTSITE
config.Site.whitelist = [
    "T3_UK_London_QMUL",
    "T3_BG_UNI_SOFIA",
    "T2_CN_Beijing",
    "T1_RU_JINR",
    "T2_CH_CSCS",
    "T2_DE_RWTH",
    "T3_IT_Bologna",
    "T1_IT_CNAF",
    "T3_FR_IPNL",
]
