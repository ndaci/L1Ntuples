import FWCore.ParameterSet.Config as cms

l1GenTreeProducer = cms.EDAnalyzer("L1GenTreeProducer",

                                   jetTag                  = cms.untracked.InputTag("ak5GenJets"),
                                   metTagCalo              = cms.untracked.InputTag("genMetCalo"),
                                   metTagTrue              = cms.untracked.InputTag("genMetTrue"),
                                   #  verticesTag             = cms.untracked.InputTag("offlinePrimaryVertices"),
                                   #  tracksTag               = cms.untracked.InputTag("generalTracks"),
                                   
                                   maxJet = cms.uint32(20),
                                   # maxTrk = cms.uint32(100),
                                   # maxVtx = cms.uint32(10),                                   

                                   jetptThreshold = cms.double(10),
                                   excludePdgId   = cms.vuint32()
)
