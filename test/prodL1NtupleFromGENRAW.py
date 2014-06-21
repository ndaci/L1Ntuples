# Produce L1Trees containing L1Ntuple, L1Extra, L1Gen

# Imports
import FWCore.ParameterSet.Config as cms
from L1TriggerDPG.L1Menu.customL1NtupleGEN_cfg import *

process.p.remove(process.l1RecoTreeProducer)
process.p.remove(process.l1MuonRecoTreeProducer)

if options.useStage1Layer2:
    process.p *= process.Layer2
    process.p *= process.l1NtupleProducerStage1Layer2

# L1GenTreeProducer parameters
process.l1GenTreeProducer.maxJet = cms.uint32(20)
process.l1GenTreeProducer.excludePdgId = cms.vuint32(18)

# edit here
OUTFILE="L1TreeGen.root"
#OUTFILE = options.outfile
NEVTS=-1

process.TFileService.fileName=OUTFILE
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(NEVTS) )

# M=1GeV, AV, GEN-SIM-RAW
#readFiles.extend( ['file:///user/ndaci/Data/DarkMonojet/Fall13dr/DarkMatter_Monojet_M-1_AV_Tune4C_13TeV-madgraph/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/08E15E24-04D7-E311-ABEA-002590A3C992.root' ] )

# M=1000GeV, AV, GEN-SIM-RAW
readFiles.extend( ['file:///user/ndaci/Data/DarkMonojet/Fall13dr/DarkMatter_Monojet_M-1000_AV_Tune4C_13TeV-madgraph/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/020E17AC-5ED9-E311-AA77-001E67396581.root' ] )

## processDumpFile = open('config.dump', 'w')
## print >> processDumpFile, process.dumpPython()
