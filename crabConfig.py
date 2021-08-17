import CRABClient
from WMCore.Configuration import Configuration

config = Configuration()
config.section_("General")
config.General.transferOutputs = True
# config.General.requestName = 'ntuple_MH-350_directory'
config.General.requestName = "ntuple_MH-125_directory"
config.section_("JobType")
# config.JobType.psetName = 'my_CMSSW_config.py'
config.JobType.pluginName = "Analysis"
# config.JobType.psetName = 'pset_tutorial_analysis.py'
config.JobType.psetName = "ntuple_maker_350_def.py"
config.JobType.outputFiles = ["L1ntuple.root"]
config.section_("Data")
# config.Data.inputDataset = '/Neutrino_Pt-2to20_gun/Run3Winter20DRPremixMiniAOD-SNB_110X_mcRun3_2021_realistic_v6-v1/GEN-SIM-RAW'
# config.Data.inputDataset = '/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-10000mm_TuneCP5_14TeV_pythia8/Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/GEN-SIM-RAW'
config.Data.inputDataset = "/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-900mm_TuneCP5_14TeV_pythia8/Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v2/GEN-SIM-RAW"
config.Data.publication = True
config.Data.unitsPerJob = 180  # 20
config.Data.splitting = "Automatic"
# config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
# config.Data.splitting = 'FileBased'
config.Data.inputDBS = "global"
config.Data.outputDatasetTag = "testLLP"
config.section_("Site")
# config.Site.blacklist = ['T2_IT_Legnaro']
# config.Site.whitelist = ['T2_IT_Bari']
config.Site.storageSite = "T2_DE_DESY"
