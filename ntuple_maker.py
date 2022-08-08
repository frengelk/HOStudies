# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: l1Ntuple -s RAW2DIGI --python_filename=ntuple_maker_def.py -n 10 --no_output --era=Run3 --conditions=110X_mcRun3_2021_realistic_v6 --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAWsimHcalTP --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleAODRAWEMU --customise_commands=process.HcalTPGCoderULUT.LUTGenerationMode=cms.bool(False) --filein=/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-900mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/181F649B-EDBE-2046-8CC5-AA23D78CE377.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process("RAW2DIGI", Run3)

# import of standard configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(10),
    output=cms.optional.untracked.allowed(cms.int32, cms.PSet),
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(
        #"/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-900mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/181F649B-EDBE-2046-8CC5-AA23D78CE377.root"
        "/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4mu_MH-125_MFF-12_CTau-900mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/240000/40F6EEC0-86E7-6F44-A6DF-CF042AE287D6.root"
    ),
    secondaryFileNames=cms.untracked.vstring(),
)

process.options = cms.untracked.PSet(
    FailPath=cms.untracked.vstring(),
    IgnoreCompletely=cms.untracked.vstring(),
    Rethrow=cms.untracked.vstring(),
    SkipEvent=cms.untracked.vstring("ProductNotFound"),
    allowUnscheduled=cms.obsolete.untracked.bool,
    canDeleteEarly=cms.untracked.vstring(),
    emptyRunLumiMode=cms.obsolete.untracked.string,
    eventSetup=cms.untracked.PSet(
        forceNumberOfConcurrentIOVs=cms.untracked.PSet(),
        numberOfConcurrentIOVs=cms.untracked.uint32(1),
    ),
    fileMode=cms.untracked.string("FULLMERGE"),
    forceEventSetupCacheClearOnNewRun=cms.untracked.bool(False),
    makeTriggerResults=cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
    numberOfConcurrentRuns=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0),
    numberOfThreads=cms.untracked.uint32(1),
    printDependencies=cms.untracked.bool(False),
    sizeOfStackForThreadsInKB=cms.optional.untracked.uint32,
    throwIfIllegalParameter=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(False),
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation=cms.untracked.string("l1Ntuple nevts:10"),
    name=cms.untracked.string("Applications"),
    version=cms.untracked.string("$Revision: 1.19 $"),
)

# Output definition

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, "110X_mcRun3_2021_realistic_v6", "")

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step, process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask

associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAWsimHcalTP

# call to customisation function L1TReEmulFromRAWsimHcalTP imported from L1Trigger.Configuration.customiseReEmul
process = L1TReEmulFromRAWsimHcalTP(process)

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleAODRAWEMU
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleRAW
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleGEN
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleEMU

from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleAODEMU
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleAOD_MC
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleAOD

# call to customisation function L1NtupleAODRAWEMU imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleRAW(process)
process = L1NtupleGEN(process)
#process = L1NtupleRAWEMUGEN(process)
#process = L1NtupleAOD(process)
#process = L1NtupleAODRAWEMU(process)
process = L1NtupleEMU(process)
#process = L1NtupleAODEMU(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseSettings
from L1Trigger.Configuration.customiseSettings import L1TSettingsToCaloParams_2018_v1_4_1

process = L1TSettingsToCaloParams_2018_v1_4_1(process)

# End of customisation functions

# Customisation from command line

# process.HcalTPGCoderULUT.LUTGenerationMode=cms.bool(False)
# Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import (
    customiseLogErrorHarvesterUsingOutputCommands,
)

process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete

process = customiseEarlyDelete(process)
# End adding early deletion
