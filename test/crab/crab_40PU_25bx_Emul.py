from WMCore.Configuration import Configuration
config = Configuration()

name = 'Rate_13TeV_40PU_25ns_62X_ReEmul2015_v4_12June2014_Emul'

# GENERAL
config.section_("General")
config.General.requestName = name 
config.General.workArea    = '/user/ndaci/CRABBY/L1Trigger/L1Menu2015/Neutrino13TeV_40PU_25bx/'
config.General.saveLogs    = True
#config.General.transferOutput = True

# JOB TYPE
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../customL1NtupleFromRaw.py'
config.JobType.pyCfgParams = ['reEmulation=True','reEmulMuons=True','reEmulCalos=True','patchNtuple=True','force2012Config=True','customDTTF=True','dttfLutsFile=sqlite:src/L1TriggerDPG/L1Menu/data/dttf_config.db','useStage1Layer2=True','globalTag=POSTLS162_V2::All','runOnMC=True','runOnPostLS1=True','whichPU=40']
#config.JobType.inputFiles = '../../data/dttf_config.db'
config.JobType.allowNonProductionCMSSW = True
config.JobType.maxmemory = 2500

# INPUT DATA
config.section_("Data")
config.Data.inputDataset = '/Neutrino_Pt-2to20_gun/Fall13dr-tsg_PU40bx25_POSTLS162_V2-v1/GEN-SIM-RAW'
config.Data.dbsUrl = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 100
config.Data.publication = True
config.Data.publishDbsUrl = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.publishDataName = name
config.Data.ignoreLocality = False # allows to process inputs on CE != site hosting inputs
#config.Data.lumiMask = 
#config.Data.runRange = 

#A custom string to insert in the output file name inside the CRAB-created directory path to allow organizing groups of tasks.
#config.Data.prefix =  

# USER
config.section_("User")
config.User.email = 'nadir.daci@cern.ch'
#config.User.voRole = 
config.User.voGroup = 'becms'

# GRID
config.section_("Site")
config.Site.storageSite = 'T2_BE_IIHE'
#config.Site.whitelist = 
config.Site.blacklist = ['T1_US_FNAL']
