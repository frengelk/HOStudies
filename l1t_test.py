# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: --conditions 80X_mcRun2_asymptotic_v14 -s L1REPACK:FullMC,RAW2DIGI -n 10 --era=Run2_2016 --filein root://eoscms.cern.ch//eos/cms/store/mc/RunIISpring16DR80/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RAWAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v2/30000/040AFBF1-8E0B-E611-9C22-001E67E95A58.root --processName reL1T --python_filename l1t_test.py --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAW
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016

process = cms.Process("reL1T", Run2_2016)

# import of standard configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.SimL1EmulatorRepack_FullMC_cff")
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
        "root://eoscms.cern.ch//eos/cms/store/mc/RunIISpring16DR80/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RAWAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v2/30000/040AFBF1-8E0B-E611-9C22-001E67E95A58.root"
    ),
    secondaryFileNames=cms.untracked.vstring(),
)

process.options = cms.untracked.PSet(
    FailPath=cms.untracked.vstring(),
    IgnoreCompletely=cms.untracked.vstring(),
    Rethrow=cms.untracked.vstring(),
    SkipEvent=cms.untracked.vstring(),
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
    annotation=cms.untracked.string("--conditions nevts:10"),
    name=cms.untracked.string("Applications"),
    version=cms.untracked.string("$Revision: 1.19 $"),
)

# Output definition

process.RECOSIMoutput = cms.OutputModule(
    "PoolOutputModule",
    dataset=cms.untracked.PSet(
        dataTier=cms.untracked.string(""), filterName=cms.untracked.string("")
    ),
    fileName=cms.untracked.string("--conditions_L1REPACK_RAW2DIGI.root"),
    outputCommands=process.RECOSIMEventContent.outputCommands,
    splitLevel=cms.untracked.int32(0),
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, "80X_mcRun2_asymptotic_v14", "")

# Path and EndPath definitions
process.L1RePack_step = cms.Path(process.SimL1Emulator)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(
    process.L1RePack_step,
    process.raw2digi_step,
    process.endjob_step,
    process.RECOSIMoutput_step,
)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask

associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleRAW

# call to customisation function L1NtupleRAW imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleRAW(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete

process = customiseEarlyDelete(process)
# End adding early deletion
