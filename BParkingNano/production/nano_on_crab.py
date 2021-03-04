
from CRABClient.UserUtilities import config, ClientException
#from input_crab_data import dataset_files
import yaml
import random
import datetime
from fnmatch import fnmatch
from argparse import ArgumentParser

random.seed(2)
production_tag = datetime.date.today().strftime('%Y%b%d')

config = config()
config.section_("General")
config.General.requestName = "privateHNL_nano%s" % production_tag
config.General.workArea = 'crab_privateMCnano'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
#config.JobType.pluginName = 'PrivateMC'
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.eventsPerLumi=100
config.JobType.psetName = '/afs/cern.ch/work/r/ratramon/HNL/NanoProd/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/test/run_HNLnano_cfg.py'
#config.JobType.pyCfgParams = [ 'nThr=1', 'inputFile=BPH-step2_.root', 'outputFile=BPH-step3.root','seedOffset=%i'% random.randint(0,1000)]
config.JobType.disableAutomaticOutputCollection = False
#config.JobType.scriptExe = '/afs/cern.ch/work/r/ratramon/HNL/CMSSW_10_2_15/src/HNLsGen/slurm/V00_v00_n10_njt1/slurm_mass1.0_ctau107.535167851_prod.sh'
#config.JobType.outputFiles = ['MiniAOD.root','GenSimAODSim_step1.log', 'GenSimAODSim_step2.log', 'FrameworkJobReport.xml', 'job.log']
config.JobType.inputFiles = [
#'/afs/cern.ch/work/r/ratramon/HNL/CMSSW_10_2_15/src/HNLsGen/slurm/V00_v00_n10_njt1/BPH_mod_BtoDMuN_NtoE_cfg.py','/afs/cern.ch/work/r/ratramon/HNL/CMSSW_10_2_15/src/HNLsGen/slurm/V00_v00_n10_njt1/step2.py' ,'/afs/cern.ch/work/r/ratramon/HNL/CMSSW_10_2_15/src/HNLsGen/slurm/V00_v00_n10_njt1/step3.py','/afs/cern.ch/work/r/ratramon/HNL/CMSSW_10_2_15/src/HNLsGen/slurm/V00_v00_n10_njt1/step4.py',
'/afs/cern.ch/work/r/ratramon/HNL/NanoProd/CMSSW_10_2_15/src/PhysicsTools/BParkingNano/test/lowPtEleReg_2018_02062020_nv.db' ]
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.outputPrimaryDataset = 'private_BtoDMuN_NtoE_nano'
config.Data.userInputFiles = open("nanoIn.txt").readlines() 
#config.Data.splitting = 'EventBased'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.unitsPerJob = 100000
#config.Data.totalUnits = 50000000
config.Data.publication = False
config.Data.outputDatasetTag = 'privateHNLnano_%s'% production_tag
#config.Data.inputDBS = 'phys03'
## T3 Beijing
#config.Data.ignoreLocality = True
config.Data.outLFNDirBase = '/store/user/ratramon/HNLGen_ntuples/'
## T3 Beijing

config.section_("Data")
config.Site.whitelist = ['T2_IT_Rome']
config.Site.storageSite = 'T2_IT_Rome'


print(config)
