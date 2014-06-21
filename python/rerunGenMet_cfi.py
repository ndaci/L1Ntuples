# Rerun GenMET production with different input gen particles

# Imports
import FWCore.ParameterSet.Config as cms
from RecoMET.Configuration.GenMETParticles_cff import *
from RecoMET.METProducers.genMetCalo_cfi import *
from RecoMET.METProducers.genMetTrue_cfi import *

# Inputs for genMetCalo
genCandidatesForMET.ignoreParticleIDs = cms.vuint32(
    1000022,
    1000012, 1000014, 1000016,
    2000012, 2000014, 2000016,
    1000039, 5100039,
    4000012, 4000014, 4000016,
    9900012, 9900014, 9900016,
    39, 12, 13, 14, 16, 18
    ) 

# Inputs for genMetTrue
genParticlesForMETAllVisible.ignoreParticleIDs = cms.vuint32(
    1000022,
    1000012, 1000014, 1000016,
    2000012, 2000014, 2000016,
    1000039, 5100039,
    4000012, 4000014, 4000016,
    9900012, 9900014, 9900016,
    39, 12, 14, 16, 18
    )

# Build sequences
genMETParticles = cms.Sequence(genCandidatesForMET+genParticlesForMETAllVisible)
recoGenMET      = cms.Sequence(genMETParticles+genMetCalo+genMetTrue)
