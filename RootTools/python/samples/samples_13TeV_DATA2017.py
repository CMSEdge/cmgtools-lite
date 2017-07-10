import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
kreator = ComponentCreator()

### ----------------------------- Zero Tesla run  ----------------------------------------

dataDir = "$CMSSW_BASE/src/CMGTools/TTHAnalysis/data"  # use environmental variable, useful for instance to run on CRAB
json='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_294297-999999_13TeV_PromptReco_Collisions17_JSON.txt'

run_range = (297047,297678)
label = "_runs%s_%s"%(run_range[0], run_range[1])


DoubleMuon_Run2017B                 = kreator.makeDataComponent("DoubleMuon_Run2017B", "/DoubleMuon/Run2017B-23Jun2017-v1/MINIAOD", "CMS", ".*root", json)
DoubleEG_Run2017B                   = kreator.makeDataComponent("DoubleEG_Run2017B", "/DoubleEG/Run2017B-23Jun2017-v1/MINIAOD", "CMS", ".*root", json)
MuonEG_Run2017B                     = kreator.makeDataComponent("MuonEG_Run2017B", "/MuonEG/Run2017B-23Jun2017-v1/MINIAOD", "CMS", ".*root", json)


### ----------------------------- Run2016B v2 03Feb2017 ----------------------------------------
\

samples = [DoubleMuon_Run2017B]

from CMGTools.TTHAnalysis.setup.Efficiencies import *
dataDir = "$CMSSW_BASE/src/CMGTools/TTHAnalysis/data"

for comp in samples:
    comp.splitFactor = 1000
    comp.isMC = False
    comp.isData = True

if __name__ == "__main__":
    from CMGTools.RootTools.samples.tools import runMain
    runMain(samples)
