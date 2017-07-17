

##########################################################
##       CONFIGURATION FOR SUSY MULTILEPTON TREES       ##
## skim condition: >= 2 loose leptons, no pt cuts or id ##
##########################################################
import PhysicsTools.HeppyCore.framework.config as cfg
import re


#-------- LOAD ALL ANALYZERS -----------

from CMGTools.TTHAnalysis.analyzers.susyCore_modules_cff import *
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption

#-------- SET OPTIONS AND REDEFINE CONFIGURATIONS -----------

is50ns = getHeppyOption("is50ns",False)
analysis = getHeppyOption("analysis","susy")
runData = getHeppyOption("runData",False)
scaleProdToLumi = float(getHeppyOption("scaleProdToLumi",-1)) # produce rough equivalent of X /pb for MC datasets
saveSuperClusterVariables = getHeppyOption("saveSuperClusterVariables",False)
removeJetReCalibration = getHeppyOption("removeJetReCalibration",False)
removeJecUncertainty = getHeppyOption("removeJecUncertainty",True)
doMETpreprocessor = getHeppyOption("doMETpreprocessor",False)
skipT1METCorr = getHeppyOption("skipT1METCorr",False)
#doAK4PFCHSchargedJets = getHeppyOption("doAK4PFCHSchargedJets",False)
forcedSplitFactor = getHeppyOption("splitFactor",-1)
forcedFineSplitFactor = getHeppyOption("fineSplitFactor",-1)
isTest = getHeppyOption("test",None) != None and not re.match("^\d+$",getHeppyOption("test"))
selectedEvents=getHeppyOption("selectEvents","")

sample = "main"
#if runDataQCD or runFRMC: sample="qcd1l"
#sample = "z3l"

if analysis not in ['ttH','susy','SOS']: raise RuntimeError, 'Analysis type unknown'
print 'Using analysis type: %s'%analysis

# Lepton Skimming
ttHLepSkim.minLeptons = 2
ttHLepSkim.maxLeptons = 999

if analysis=='susy':
    susyCoreSequence.insert(susyCoreSequence.index(ttHLepSkim)+1,globalSkim)
    susyCoreSequence.remove(ttHLepSkim)
    globalSkim.selections=["2lep5","1lep5_1tau18", "2tau18","1lep5[maxObj1]"]
#   [ lambda ev: 2<=sum([(lep.miniRelIso<0.4) for lep in ev.selectedLeptons]) ] 
#   ["2lep5[os:!DS_TTW_RA5_sync]_1lep50"]#, "1lep5_1tau18", "2tau18","2lep5_1met50"]

# Run miniIso
lepAna.doMiniIsolation = True
lepAna.packedCandidates = 'packedPFCandidates'
lepAna.miniIsolationPUCorr = 'rhoArea'
lepAna.miniIsolationVetoLeptons = None # use 'inclusive' to veto inclusive leptons and their footprint in all isolation cones
lepAna.doIsolationScan = False

# Lepton Preselection
lepAna.loose_electron_id = "MVA_ID_NonTrig_Spring16_VLooseIdEmu"
isolation = "miniIso"

jetAna.copyJetsByValue = True # do not remove this
metAna.copyMETsByValue = True # do not remove this

if analysis=='susy':
    jetAna.cleanJetsFromLeptons=False
    jetAna.cleanSelectedLeptons=True
    jetAna.storeLowPtJets=True
    jetAna.jetEtaCentral = jetAna.jetEta
    jetAna.mcGT="Spring16_25nsV8_MC"    
    #jetAna.dataGT   = "Spring16_25nsV8BCD_DATA Spring16_25nsV8E_DATA Spring16_25nsV8F_DATA Spring16_25nsV8_DATA"
    jetAna.dataGT   = "Spring16_25nsV8BCD_DATA"
    jetAna.runsDataJEC   = [276811, 277420, 278802]
if not removeJecUncertainty:
    jetAna.addJECShifts = True
    jetAnaScaleDown.copyJetsByValue = True # do not remove this
    metAnaScaleDown.copyMETsByValue = True # do not remove this
    jetAnaScaleUp.copyJetsByValue = True # do not remove this
    metAnaScaleUp.copyMETsByValue = True # do not remove this
    susyCoreSequence.insert(susyCoreSequence.index(jetAna)+1, jetAnaScaleDown)
    susyCoreSequence.insert(susyCoreSequence.index(jetAna)+1, jetAnaScaleUp)
    susyCoreSequence.insert(susyCoreSequence.index(metAna)+1, metAnaScaleDown)
    susyCoreSequence.insert(susyCoreSequence.index(metAna)+1, metAnaScaleUp)
    if analysis=='susy':
        jetAnaScaleDown.cleanJetsFromLeptons=False
        jetAnaScaleDown.cleanSelectedLeptons=True
        jetAnaScaleDown.storeLowPtJets=True
        jetAnaScaleDown.jetEtaCentral = jetAnaScaleDown.jetEta
        jetAnaScaleUp.cleanJetsFromLeptons=False
        jetAnaScaleUp.cleanSelectedLeptons=True
        jetAnaScaleUp.storeLowPtJets=True
        jetAnaScaleUp.jetEtaCentral = jetAnaScaleUp.jetEta
        jetAnaScaleDown.mcGT="Spring16_25nsV8_MC"    
        #jetAnaScaleDown.dataGT   = "Spring16_25nsV8BCD_DATA Spring16_25nsV8E_DATA Spring16_25nsV8F_DATA Spring16_25nsV8_DATA"
        jetAnaScaleDown.dataGT   = "Spring16_25nsV8BCD_DATA"
        jetAnaScaleDown.runsDataJEC   = [276811, 277420, 278802]
        jetAnaScaleUp.mcGT="Spring16_25nsV8_MC"    
        #jetAnaScaleUp.dataGT   = "Spring16_25nsV8BCD_DATA Spring16_25nsV8E_DATA Spring16_25nsV8F_DATA Spring16_25nsV8_DATA"
        jetAnaScaleUp.dataGT   = "Spring16_25nsV8BCD_DATA"
        jetAnaScaleUp.runsDataJEC   = [276811, 277420, 278802]


if analysis in ['SOS']:
## -- SOS preselection settings ---

    lepAna.doDirectionalIsolation = [0.3,0.4]
    lepAna.doFixedConeIsoWithMiniIsoVeto = True

    # Lepton Skimming
    ttHLepSkim.minLeptons = 2
    ttHLepSkim.maxLeptons = 999
    
#    # Jet-Met Skimming
#    ttHJetMETSkim.jetPtCuts = [0,]
    ttHJetMETSkim.metCut    = 50
    susyCoreSequence.append(ttHJetMETSkim)

    # Lepton Preselection
    lepAna.inclusive_muon_pt  = 3
    lepAna.loose_muon_pt  = 3
    lepAna.inclusive_electron_pt  = 5
    lepAna.loose_electron_pt  = 5
    #isolation = "absIso04"
    #isolation = None
    isolation = "Iperbolic"
    lepAna.loose_electron_id = "POG_MVA_ID_Spring15_NonTrig_VLooseIdEmu"

    # Lepton-Jet Cleaning
    jetAna.cleanSelectedLeptons = False
    #jetAna.minLepPt = 20 
    #jetAnaScaleUp.minLepPt = 20 
    #jetAnaScaleDown.minLepPt = 20 
    # otherwise with only absIso cut at 10 GeV and no relIso we risk cleaning away good jets

if isolation == "miniIso": 
    if (analysis=="ttH"):
        lepAna.loose_muon_isoCut     = lambda muon : muon.miniRelIso < 0.4 and muon.sip3D() < 8
        lepAna.loose_electron_isoCut = lambda elec : elec.miniRelIso < 0.4 and elec.sip3D() < 8
    elif analysis=="susy":
        lepAna.loose_muon_isoCut     = lambda muon : muon.miniRelIso < 0.4
        lepAna.loose_electron_isoCut = lambda elec : elec.miniRelIso < 0.4
    else: raise RuntimeError,'analysis field is not correctly configured'
elif isolation == None:
    lepAna.loose_muon_isoCut     = lambda muon : True
    lepAna.loose_electron_isoCut = lambda elec : True
elif isolation == "absIso04":
    lepAna.loose_muon_isoCut     = lambda muon : muon.RelIsoMIV04*muon.pt() < 10 and muon.sip3D() < 8
    lepAna.loose_electron_isoCut = lambda elec : elec.RelIsoMIV04*elec.pt() < 10 and elec.sip3D() < 8
elif isolation == "Iperbolic":
    lepAna.loose_muon_isoCut     = lambda muon : muon.relIso03*muon.pt() < (20+300/muon.pt()) and  abs(muon.ip3D()) < 0.0175 and muon.sip3D() < 2.5
    lepAna.loose_electron_isoCut = lambda elec : elec.relIso03*elec.pt() < (20+300/elec.pt()) and  abs(elec.ip3D()) < 0.0175 and elec.sip3D() < 2.5
else:
    # nothing to do, will use normal relIso03
    pass

# Switch on slow QGL
jetAna.doQG = True
jetAnaScaleUp.doQG = True
jetAnaScaleDown.doQG = True

# Switch off slow photon MC matching
photonAna.do_mc_match = False

# Loose Tau configuration
tauAna.loose_ptMin = 20
tauAna.loose_etaMax = 2.3
tauAna.loose_decayModeID = "decayModeFindingNewDMs"
tauAna.loose_tauID = "decayModeFindingNewDMs"
#if analysis in ["ttH"]: #if cleaning jet-loose tau cleaning
#    jetAna.cleanJetsFromTaus = True
#    jetAnaScaleUp.cleanJetsFromTaus = True
#    jetAnaScaleDown.cleanJetsFromTaus = True


#-------- ADDITIONAL ANALYZERS -----------

if analysis in ['SOS']:
    #Adding LHE Analyzer for saving lheHT
    from PhysicsTools.Heppy.analyzers.gen.LHEAnalyzer import LHEAnalyzer 
    LHEAna = LHEAnalyzer.defaultConfig
    susyCoreSequence.insert(susyCoreSequence.index(ttHCoreEventAna), 
                            LHEAna)

## Event Analyzer for susy multi-lepton (at the moment, it's the TTH one)
from CMGTools.TTHAnalysis.analyzers.ttHLepEventAnalyzer import ttHLepEventAnalyzer
ttHEventAna = cfg.Analyzer(
    ttHLepEventAnalyzer, name="ttHLepEventAnalyzer",
    minJets25 = 0,
    )

## JetTau analyzer, to be called (for the moment) once bjetsMedium are produced
from CMGTools.TTHAnalysis.analyzers.ttHJetTauAnalyzer import ttHJetTauAnalyzer
ttHJetTauAna = cfg.Analyzer(
    ttHJetTauAnalyzer, name="ttHJetTauAnalyzer",
    )

## Insert the FatJet, SV, HeavyFlavour analyzers in the sequence
susyCoreSequence.insert(susyCoreSequence.index(ttHCoreEventAna), 
                        ttHFatJetAna)
susyCoreSequence.insert(susyCoreSequence.index(ttHCoreEventAna), 
                        ttHSVAna)
susyCoreSequence.insert(susyCoreSequence.index(ttHCoreEventAna), 
                        ttHHeavyFlavourHadronAna)

## Insert declustering analyzer
from CMGTools.TTHAnalysis.analyzers.ttHDeclusterJetsAnalyzer import ttHDeclusterJetsAnalyzer
ttHDecluster = cfg.Analyzer(
    ttHDeclusterJetsAnalyzer, name='ttHDecluster',
    lepCut     = lambda lep,ptrel : lep.pt() > 10,
    maxSubjets = 6, # for exclusive reclustering
    ptMinSubjets = 5, # for inclusive reclustering
    drMin      = 0.2, # minimal deltaR(l,subjet) required for a successful subjet match
    ptRatioMax = 1.5, # maximum pt(l)/pt(subjet) required for a successful match
    ptRatioDiff = 0.1,  # cut on abs( pt(l)/pt(subjet) - 1 ) sufficient to call a match successful
    drMatch     = 0.02, # deltaR(l,subjet) sufficient to call a match successful
    ptRelMin    = 5,    # maximum ptRelV1(l,subjet) sufficient to call a match successful
    prune       = True, # also do pruning of the jets 
    pruneZCut       = 0.1, # pruning parameters (usual value in CMS: 0.1)
    pruneRCutFactor = 0.5, # pruning parameters (usual value in CMS: 0.5)
    verbose     = 0,   # print out the first N leptons
    jetCut = lambda jet : jet.pt() > 20,
    mcPartonPtCut = 20,
    mcLeptonPtCut =  5,
    mcTauPtCut    = 15,
    )
susyCoreSequence.insert(susyCoreSequence.index(ttHFatJetAna)+1, ttHDecluster)


from CMGTools.TTHAnalysis.analyzers.treeProducerSusyMultilepton import * 

if analysis=="susy":
    ttHLepSkim.allowLepTauComb = True
    susyCoreSequence.insert(susyCoreSequence.index(ttHCoreEventAna),
                            susyLeptonMatchAna)
    susyCoreSequence.insert(susyCoreSequence.index(ttHCoreEventAna),
                            susyTauMatchAna)
    leptonTypeSusyExtraLight.addVariables([
            NTupleVariable("mcUCSXMatchId", lambda x : x.mcUCSXMatchId if hasattr(x,'mcUCSXMatchId') else -1, mcOnly=True, help="MC truth matching a la UCSX"),
            ])
    tauTypeSusy.addVariables([
            NTupleVariable("mcUCSXMatchId", lambda x : x.mcUCSXMatchId if hasattr(x,'mcUCSXMatchId') else -1, mcOnly=True, help="MC truth matching a la UCSX"),
            ])
    photonAna.do_mc_match = True
    susyMultilepton_collections.update({ # for conversion studies
            "selectedPhotons"    : NTupleCollection("PhoGood", photonTypeSusy, 10, help="Selected photons"),
            }) 
    del susyMultilepton_collections["discardedJets"]
    susyMultilepton_collections.update({"discardedJets"   : NTupleCollection("DiscJet", jetTypeSusySuperLight, 15, help="Jets discarted in the jet-lepton cleaning (JEC)")
                                        })
    #keepJetVars=["pt","eta","phi","mass",
    #             #"etaetaMoment","phiphiMoment",
    #             "mcFlavour","partonFlavour",
    #             "btagCSV"]
    #for var in jetTypeSusyExtraLight.allVars(not runData):
    #    if var.name not in keepJetVars: 
    #        jetTypeSusyExtraLight.removeVariable(var)
    #jetTypeSusyExtraLight.addVariables([
    #        NTupleVariable("etaetaMoment", lambda x : x.etaetaMoment() if hasattr(x,'etaetaMoment') else -1, mcOnly=True, help="eta eta moment"),
    #        NTupleVariable("phiphiMoment", lambda x : x.phiphiMoment() if hasattr(x,'phiphiMoment') else -1, mcOnly=True, help="phi phi moment"),
    #        ])
elif analysis=='SOS':
    # Soft lepton MVA
    ttHCoreEventAna.doLeptonMVASoft = True
    leptonTypeSusyExtraLight.addVariables([
            NTupleVariable("mvaSoftT2tt",    lambda lepton : lepton.mvaValueSoftT2tt, help="Lepton MVA (Soft T2tt version)"),
            NTupleVariable("mvaSoftEWK",    lambda lepton : lepton.mvaValueSoftEWK, help="Lepton MVA (Soft EWK version)"),
            ])

# Spring16 electron MVA - follow instructions on pull request for correct area setup
leptonTypeSusy.addVariables([
        NTupleVariable("mvaIdSpring16HZZ",   lambda lepton : lepton.mvaRun2("Spring16HZZ") if abs(lepton.pdgId()) == 11 else 1, help="EGamma POG MVA ID, Spring16, HZZ; 1 for muons"),
        NTupleVariable("mvaIdSpring16GP",   lambda lepton : lepton.mvaRun2("Spring16GP") if abs(lepton.pdgId()) == 11 else 1, help="EGamma POG MVA ID, Spring16, GeneralPurpose; 1 for muons"),
        ])

if lepAna.doIsolationScan:
    leptonTypeSusyExtraLight.addVariables([
            NTupleVariable("scanAbsIsoCharged005", lambda x : x.ScanAbsIsoCharged005 if hasattr(x,'ScanAbsIsoCharged005') else -999, help="PF abs charged isolation dR=0.05, no pile-up correction"),
            NTupleVariable("scanAbsIsoCharged01", lambda x : x.ScanAbsIsoCharged01 if hasattr(x,'ScanAbsIsoCharged01') else -999, help="PF abs charged isolation dR=0.1, no pile-up correction"),
            NTupleVariable("scanAbsIsoCharged02", lambda x : x.ScanAbsIsoCharged02 if hasattr(x,'ScanAbsIsoCharged02') else -999, help="PF abs charged isolation dR=0.2, no pile-up correction"),
            NTupleVariable("scanAbsIsoCharged03", lambda x : x.ScanAbsIsoCharged03 if hasattr(x,'ScanAbsIsoCharged03') else -999, help="PF abs charged isolation dR=0.3, no pile-up correction"),
            NTupleVariable("scanAbsIsoCharged04", lambda x : x.ScanAbsIsoCharged04 if hasattr(x,'ScanAbsIsoCharged04') else -999, help="PF abs charged isolation dR=0.4, no pile-up correction"),
            NTupleVariable("scanAbsIsoNeutral005", lambda x : x.ScanAbsIsoNeutral005 if hasattr(x,'ScanAbsIsoNeutral005') else -999, help="PF abs neutral+photon isolation dR=0.05, no pile-up correction"),
            NTupleVariable("scanAbsIsoNeutral01", lambda x : x.ScanAbsIsoNeutral01 if hasattr(x,'ScanAbsIsoNeutral01') else -999, help="PF abs neutral+photon isolation dR=0.1, no pile-up correction"),
            NTupleVariable("scanAbsIsoNeutral02", lambda x : x.ScanAbsIsoNeutral02 if hasattr(x,'ScanAbsIsoNeutral02') else -999, help="PF abs neutral+photon isolation dR=0.2, no pile-up correction"),
            NTupleVariable("scanAbsIsoNeutral03", lambda x : x.ScanAbsIsoNeutral03 if hasattr(x,'ScanAbsIsoNeutral03') else -999, help="PF abs neutral+photon isolation dR=0.3, no pile-up correction"),
            NTupleVariable("scanAbsIsoNeutral04", lambda x : x.ScanAbsIsoNeutral04 if hasattr(x,'ScanAbsIsoNeutral04') else -999, help="PF abs neutral+photon isolation dR=0.4, no pile-up correction"),
            NTupleVariable("miniIsoR", lambda x: getattr(x,'miniIsoR',-999), help="miniIso cone size"),
            NTupleVariable("effArea", lambda x: getattr(x,'EffectiveArea03',-999), help="effective area used for PU subtraction"),
            NTupleVariable("rhoForEA", lambda x: getattr(x,'rho',-999), help="rho used for EA PU subtraction")
            ])

# for electron scale and resolution checks
if saveSuperClusterVariables:
    leptonTypeSusyExtraLight.addVariables([
            NTupleVariable("e5x5", lambda x: x.e5x5() if (abs(x.pdgId())==11 and hasattr(x,"e5x5")) else -999, help="Electron e5x5"),
            NTupleVariable("r9", lambda x: x.r9() if (abs(x.pdgId())==11 and hasattr(x,"r9")) else -999, help="Electron r9"),
            NTupleVariable("sigmaIetaIeta", lambda x: x.sigmaIetaIeta() if (abs(x.pdgId())==11 and hasattr(x,"sigmaIetaIeta")) else -999, help="Electron sigmaIetaIeta"),
            NTupleVariable("sigmaIphiIphi", lambda x: x.sigmaIphiIphi() if (abs(x.pdgId())==11 and hasattr(x,"sigmaIphiIphi")) else -999, help="Electron sigmaIphiIphi"),
            NTupleVariable("hcalOverEcal", lambda x: x.hcalOverEcal() if (abs(x.pdgId())==11 and hasattr(x,"hcalOverEcal")) else -999, help="Electron hcalOverEcal"),
            NTupleVariable("full5x5_e5x5", lambda x: x.full5x5_e5x5() if (abs(x.pdgId())==11 and hasattr(x,"full5x5_e5x5")) else -999, help="Electron full5x5_e5x5"),
            NTupleVariable("full5x5_r9", lambda x: x.full5x5_r9() if (abs(x.pdgId())==11 and hasattr(x,"full5x5_r9")) else -999, help="Electron full5x5_r9"),
            NTupleVariable("full5x5_sigmaIetaIeta", lambda x: x.full5x5_sigmaIetaIeta() if (abs(x.pdgId())==11 and hasattr(x,"full5x5_sigmaIetaIeta")) else -999, help="Electron full5x5_sigmaIetaIeta"),
            NTupleVariable("full5x5_sigmaIphiIphi", lambda x: x.full5x5_sigmaIphiIphi() if (abs(x.pdgId())==11 and hasattr(x,"full5x5_sigmaIphiIphi")) else -999, help="Electron full5x5_sigmaIphiIphi"),
            NTupleVariable("full5x5_hcalOverEcal", lambda x: x.full5x5_hcalOverEcal() if (abs(x.pdgId())==11 and hasattr(x,"full5x5_hcalOverEcal")) else -999, help="Electron full5x5_hcalOverEcal"),
            NTupleVariable("correctedEcalEnergy", lambda x: x.correctedEcalEnergy() if (abs(x.pdgId())==11 and hasattr(x,"correctedEcalEnergy")) else -999, help="Electron correctedEcalEnergy"),
            NTupleVariable("eSuperClusterOverP", lambda x: x.eSuperClusterOverP() if (abs(x.pdgId())==11 and hasattr(x,"eSuperClusterOverP")) else -999, help="Electron eSuperClusterOverP"),
            NTupleVariable("ecalEnergy", lambda x: x.ecalEnergy() if (abs(x.pdgId())==11 and hasattr(x,"ecalEnergy")) else -999, help="Electron ecalEnergy"),
            NTupleVariable("superCluster_rawEnergy", lambda x: x.superCluster().rawEnergy() if (abs(x.pdgId())==11 and hasattr(x,"superCluster")) else -999, help="Electron superCluster.rawEnergy"),
            NTupleVariable("superCluster_preshowerEnergy", lambda x: x.superCluster().preshowerEnergy() if (abs(x.pdgId())==11 and hasattr(x,"superCluster")) else -999, help="Electron superCluster.preshowerEnergy"),
            NTupleVariable("superCluster_correctedEnergy", lambda x: x.superCluster().correctedEnergy() if (abs(x.pdgId())==11 and hasattr(x,"superCluster")) else -999, help="Electron superCluster.correctedEnergy"),
            NTupleVariable("superCluster_energy", lambda x: x.superCluster().energy() if (abs(x.pdgId())==11 and hasattr(x,"superCluster")) else -999, help="Electron superCluster.energy"),
            NTupleVariable("superCluster_clustersSize", lambda x: x.superCluster().clustersSize() if (abs(x.pdgId())==11 and hasattr(x,"superCluster")) else -999, help="Electron superCluster.clustersSize"),
            NTupleVariable("superCluster_seed.energy", lambda x: x.superCluster().seed().energy() if (abs(x.pdgId())==11 and hasattr(x,"superCluster")) else -999, help="Electron superCluster.seed.energy"),
])

if not removeJecUncertainty:
    susyMultilepton_globalObjects.update({
            "met_jecUp" : NTupleObject("met_jecUp", metType, help="PF E_{T}^{miss}, after type 1 corrections (JEC plus 1sigma)"),
            "met_jecDown" : NTupleObject("met_jecDown", metType, help="PF E_{T}^{miss}, after type 1 corrections (JEC minus 1sigma)"),
            })
    susyMultilepton_collections.update({
            "cleanJets_jecUp"       : NTupleCollection("Jet_jecUp",     jetTypeSusyExtraLight, 15, help="Cental jets after full selection and cleaning, sorted by pt (JEC plus 1sigma)"),
            "cleanJets_jecDown"     : NTupleCollection("Jet_jecDown",     jetTypeSusyExtraLight, 15, help="Cental jets after full selection and cleaning, sorted by pt (JEC minus 1sigma)"),
            "discardedJets_jecUp"   : NTupleCollection("DiscJet_jecUp", jetTypeSusySuperLight if analysis=='susy' else jetTypeSusyExtraLight, 15, help="Jets discarted in the jet-lepton cleaning (JEC +1sigma)"),
            "discardedJets_jecDown" : NTupleCollection("DiscJet_jecDown", jetTypeSusySuperLight if analysis=='susy' else jetTypeSusyExtraLight, 15, help="Jets discarted in the jet-lepton cleaning (JEC -1sigma)"),
            })

## Tree Producer
treeProducer = cfg.Analyzer(
     AutoFillTreeProducer, name='treeProducerSusyMultilepton',
     vectorTree = True,
     saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
     defaultFloatType = 'F', # use Float_t for floating point
     PDFWeights = PDFWeights,
     globalVariables = susyMultilepton_globalVariables,
     globalObjects = susyMultilepton_globalObjects,
     collections = susyMultilepton_collections,
)
if analysis in ['SOS']:
    del treeProducer.collections["discardedLeptons"]

## histo counter
susyCoreSequence.insert(susyCoreSequence.index(skimAnalyzer), susyCounter)
susyScanAna.doLHE=False # until a proper fix is put in the analyzer

# HBHE new filter
from CMGTools.TTHAnalysis.analyzers.hbheAnalyzer import hbheAnalyzer
hbheAna = cfg.Analyzer(
    hbheAnalyzer, name="hbheAnalyzer", IgnoreTS4TS5ifJetInLowBVRegion=False
    )

susyCoreSequence.insert(susyCoreSequence.index(ttHCoreEventAna),hbheAna)
treeProducer.globalVariables.append(NTupleVariable("hbheFilterNew50ns", lambda ev: ev.hbheFilterNew50ns, int, help="new HBHE filter for 50 ns"))
treeProducer.globalVariables.append(NTupleVariable("hbheFilterNew25ns", lambda ev: ev.hbheFilterNew25ns, int, help="new HBHE filter for 25 ns"))
treeProducer.globalVariables.append(NTupleVariable("hbheFilterIso", lambda ev: ev.hbheFilterIso, int, help="HBHE iso-based noise filter"))
treeProducer.globalVariables.append(NTupleVariable("Flag_badChargedHadronFilter", lambda ev: ev.badChargedHadron, help="bad charged hadron filter decision"))
treeProducer.globalVariables.append(NTupleVariable("Flag_badMuonFilter", lambda ev: ev.badMuon, help="bad muon filter decision"))


#additional MET quantities
metAna.doTkMet = True
treeProducer.globalVariables.append(NTupleVariable("met_trkPt", lambda ev : ev.tkMet.pt() if  hasattr(ev,'tkMet') else  0, help="tkmet p_{T}"))
treeProducer.globalVariables.append(NTupleVariable("met_trkPhi", lambda ev : ev.tkMet.phi() if  hasattr(ev,'tkMet') else  0, help="tkmet phi"))

if not skipT1METCorr:
    if doMETpreprocessor: 
        print "WARNING: you're running the MET preprocessor and also Type1 MET corrections. This is probably not intended."
    jetAna.calculateType1METCorrection = True
    metAna.recalibrate = "type1"
    jetAnaScaleUp.calculateType1METCorrection = True
    metAnaScaleUp.recalibrate = "type1"
    jetAnaScaleDown.calculateType1METCorrection = True
    metAnaScaleDown.recalibrate = "type1"


#ISR jet collection
if analysis=='SOS':
    from CMGTools.TTHAnalysis.analyzers.nIsrAnalyzer import NIsrAnalyzer
    nIsrAnalyzer = cfg.Analyzer(
        NIsrAnalyzer, name='nIsrAnalyzer',
    )
    susyCoreSequence.insert(susyCoreSequence.index(jetAna)+1,nIsrAnalyzer)
    treeProducer.globalVariables.append(NTupleVariable("nISR", lambda ev: ev.nIsr, int, help="number of ISR jets according to SUSY recommendations"))


#-------- SAMPLES AND TRIGGERS -----------
from CMGTools.RootTools.samples.samples_13TeV_DATA2017 import *

selectedComponents = [DoubleMuon_Run2017B, DoubleEG_Run2017B, MuonEG_Run2017B]

from CMGTools.RootTools.samples.triggers_13TeV_DATA2017 import *

triggerFlagsAna.triggerBits = {
    'DoubleMu' : triggers_mumu,
    'DoubleEl' : triggers_ee,
    'MuEG'     : triggers_mue,
    #'MonoJet80MET90' : triggers_Jet80MET90,
    #'MonoJet80MET120' : triggers_Jet80MET120,
    #'METMu5' : triggers_MET120Mu5,
}
triggerFlagsAna.unrollbits = True
triggerFlagsAna.saveIsUnprescaled = True
triggerFlagsAna.checkL1Prescale = True






#from CMGTools.RootTools.samples.samples_13TeV_RunIISpring16MiniAODv1 import *
#from CMGTools.RootTools.samples.samples_13TeV_RunIISpring16MiniAODv2 import *
from CMGTools.HToZZ4L.tools.configTools import printSummary, configureSplittingFromTime, cropToLumi, prescaleComponents, insertEventSelector


if analysis=='susy':
    print 'muy bien' 

elif analysis=='SOS':

    selectedComponents = selectedComponents
    #samples_scans = [SMS_TChiWZ, SMS_T2ttDiLep_mStop_10to80] 
    #samples_privateSig = Higgsino 
    #samples_mainBkg = [TTJets_DiLepton, TBar_tWch_ext, T_tWch_ext] + DYJetsM5to50HT + DYJetsM50HT
    #samples_mainBkgVV = [VVTo2L2Nu, VVTo2L2Nu_ext]
    #samples_fakesBkg = [TTJets_SingleLeptonFromTbar, TTJets_SingleLeptonFromT, WJetsToLNuHT] 
    #samples_rareBkg = [WZTo3LNu, WWToLNuQQ, WZTo1L3Nu, WZTo1L1Nu2Q, ZZTo2L2Q, ZZTo4L, WWW, WZZ, WWZ, ZZZ, T_tch_powheg, TBar_tch_powheg, TToLeptons_sch_amcatnlo, TTHnobb_pow, WWDouble, WpWpJJ, TTWToLNu_ext, TTZToLLNuNu_ext, TTZToLLNuNu_m1to10, TTGJets, WGToLNuG_amcatnlo_ext, ZGTo2LG_ext, TGJets] #WZTo2L2Q,WGToLNuG, #still missing
    #selectedComponents = samples_mainBkg + samples_mainBkgVV
    #configureSplittingFromTime(samples_mainBkg,150,3)

    

elif analysis=="ttH":
    selectedComponents = selectedComponents
#    samples_2l = [ TTWToLNu, TTZToLLNuNu, TTLLJets_m1to10, TTTT_ext, tZq_ll ] + TTHnobb_mWCutfix
#    samples_2l = [WJetsToLNu_LO, WJetsToLNu, DYJetsToLL_M10to50_LO, DYJetsToLL_M10to50, DYJetsToLL_M50, DYJetsToLL_M50_LO, TTJets, TT_pow, TTJets_SingleLeptonFromTbar, TTJets_SingleLeptonFromT, TTJets_DiLepton, TBar_tWch, T_tWch, TToLeptons_tch_amcatnlo, TToLeptons_sch_amcatnlo, TTGJets, WGToLNuG, ZGTo2LG, TGJets, WWDouble, WpWpJJ, TTTT, VHToNonbb, GGHZZ4L,tZq_ll, WZTo3LNu, ZZTo4L, WWTo2L2Nu, WWW, WWZ, WZZ, ZZZ, TTHnobb_pow, TTW_LO, TTZ_LO, TTWToLNu, TTZToLLNuNu, TTLLJets_m1to10] + TTHnobb_mWCutfix
#    samples_1l = [QCD_Mu15] + QCD_Mu5 + [WJetsToLNu_LO,DYJetsToLL_M10to50_LO,DYJetsToLL_M50_LO,TT_pow] + QCDPtEMEnriched + QCDPtbcToE
#    selectedComponents = samples_2l
#    for comp in selectedComponents: comp.splitFactor = 200
#    printSummary(selectedComponents)
#    cropToLumi([TTTT_ext,tZq_ll],200)
#    cropToLumi(TTHnobb_mWCutfix,2000)
#    configureSplittingFromTime(samples_1l,50,3)
#    configureSplittingFromTime(samples_2l,100,3)

if scaleProdToLumi>0: # select only a subset of a sample, corresponding to a given luminosity (assuming ~30k events per MiniAOD file, which is ok for central production)
    target_lumi = scaleProdToLumi # in inverse picobarns
    for c in selectedComponents:
        if not c.isMC: continue
        nfiles = int(min(ceil(target_lumi * c.xSection / 30e3), len(c.files)))
        #if nfiles < 50: nfiles = min(4*nfiles, len(c.files))
        print "For component %s, will want %d/%d files; AAA %s" % (c.name, nfiles, len(c.files), "eoscms" not in c.files[0])
        c.files = c.files[:nfiles]
        c.splitFactor = len(c.files)
        c.fineSplitFactor = 1




if removeJetReCalibration:
    jetAna.recalibrateJets = False
    jetAnaScaleUp.recalibrateJets = False
    jetAnaScaleDown.recalibrateJets = False


#-------- SEQUENCE -----------

sequence = cfg.Sequence(susyCoreSequence+[
        ttHJetTauAna,
        ttHEventAna,
        treeProducer,
    ])
preprocessor = None

#-------- HOW TO RUN -----------
json='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-297723_13TeV_PromptReco_Collisions17_JSON.txt'
#json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
#sequence.remove(jsonAna)
#DoubleMuon = kreator.makeDataComponent("DoubleMuon_Run2017B", "/DoubleMuon/Run2017B-23Jun2017-v1/MINIAOD", "CMS", ".*root", run_range = (294927,297723), triggers = triggers_mumu)
#DoubleEG = kreator.makeDataComponent("DoubleEG_Run2017B", "/DoubleEG/Run2017B-23Jun2017-v1/MINIAOD", "CMS", ".*root", run_range = (294927,297723), triggers = triggers_ee)
#MuonEG = kreator.makeDataComponent("MuonEG_Run2017B", "/MuonEG/Run2017B-23Jun2017-v1/MINIAOD", "CMS", ".*root", run_range = (294927,297723), triggers = triggers_mue)
DoubleMuon_v1 = kreator.makeDataComponent("DoubleMuon_Run2017B", "/DoubleMuon/Run2017B-PromptReco-v1/MINIAOD", "CMS", ".*root", run_range = (294927,297723), triggers = triggers_mumu)
DoubleEG_v1 = kreator.makeDataComponent("DoubleEG_Run2017B", "/DoubleEG/Run2017B-PromptReco-v1/MINIAOD", "CMS", ".*root", run_range = (294927,297723), triggers = triggers_ee)
MuonEG_v1 = kreator.makeDataComponent("MuonEG_Run2017B", "/MuonEG/Run2017B-PromptReco-v1/MINIAOD", "CMS", ".*root", run_range = (294927,297723), triggers = triggers_mue)
#DoubleEG = kreator.makeDataComponent("DoubleEG_Run2016H_run283885", "/DoubleEG/Run2016H-PromptReco-v2/MINIAOD", "CMS", ".*root", run_range = (283885,283885), triggers = triggers_ee)
#DoubleMuon.files = [ 'root://eoscms//eos/cms/store/data/Run2017B/DoubleMuon/MINIAOD/23Jun2017-v1/10000/046A6D49-4859-E711-8CAF-0025904B2C68.root' ]
#DoubleEG.files = [ 'root://eoscms//eos/cms/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/885/00000/743981FC-949D-E611-836E-FA163EC09DF2.root' ]
selectedComponents = [ DoubleMuon_v1, DoubleEG_v1, MuonEG_v1 ]
for comp in selectedComponents:
    comp.json='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-297723_13TeV_PromptReco_Collisions17_JSON.txt'
    comp.splitFactor = 100
    comp.isMC = False
    comp.isData = True
#selectedComponents = [ DoubleMuon, DoubleEG ]




test = getHeppyOption('test')
if test == '1':
    comp = selectedComponents[0]
    comp.files = comp.files[:1]
    comp.splitFactor = 1
    comp.fineSplitFactor = 1
    selectedComponents = [ comp ]
elif test == '92X-Data':
    json='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-297723_13TeV_PromptReco_Collisions17_JSON.txt'
    #json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
    #sequence.remove(jsonAna)
    DoubleMuon = kreator.makeDataComponent("DoubleMuon_Run2017B", "/DoubleMuon/Run2017B-23Jun2017-v1/MINIAOD", "CMS", ".*root", run_range = (294927,297723), triggers = triggers_mumu)
    DoubleEG = kreator.makeDataComponent("DoubleEG_Run2017B", "/DoubleEG/Run2017B-23Jun2017-v1/MINIAOD", "CMS", ".*root", run_range = (294927,297723), triggers = triggers_ee)
    MuonEG = kreator.makeDataComponent("MuonEG_Run2017B", "/MuonEG/Run2017B-23Jun2017-v1/MINIAOD", "CMS", ".*root", run_range = (294927,297723), triggers = triggers_mue)
    #DoubleEG = kreator.makeDataComponent("DoubleEG_Run2016H_run283885", "/DoubleEG/Run2016H-PromptReco-v2/MINIAOD", "CMS", ".*root", run_range = (283885,283885), triggers = triggers_ee)
    #DoubleMuon.files = [ 'root://eoscms//eos/cms/store/data/Run2017B/DoubleMuon/MINIAOD/23Jun2017-v1/10000/046A6D49-4859-E711-8CAF-0025904B2C68.root' ]
    #DoubleEG.files = [ 'root://eoscms//eos/cms/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/283/885/00000/743981FC-949D-E611-836E-FA163EC09DF2.root' ]
    selectedComponents = [ DoubleMuon, DoubleEG, MuonEG ]
    for comp in selectedComponents:
        comp.json='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-297723_13TeV_PromptReco_Collisions17_JSON.txt'
        comp.splitFactor = 100
        comp.isMC = False
        comp.isData = True
#selectedComponents = [ DoubleMuon, DoubleEG ]

## FAST mode: pre-skim using reco leptons, don't do accounting of LHE weights (slow)"
## Useful for large background samples with low skim efficiency
if getHeppyOption("fast"):
    susyCounter.doLHE = False
    from CMGTools.TTHAnalysis.analyzers.ttHFastLepSkimmer import ttHFastLepSkimmer
    fastSkim = cfg.Analyzer(
        ttHFastLepSkimmer, name="ttHFastLepSkimmer2lep",
        muons = 'slimmedMuons', muCut = lambda mu : mu.pt() > 3 and mu.isLooseMuon(),
        electrons = 'slimmedElectrons', eleCut = lambda ele : ele.pt() > 5,
        minLeptons = 2, 
    )
    if jsonAna in sequence:
        sequence.insert(sequence.index(jsonAna)+1, fastSkim)
    else:
        sequence.insert(sequence.index(skimAnalyzer)+1, fastSkim)
if not getHeppyOption("keepLHEweights",False):
    if "LHE_weights" in treeProducer.collections: treeProducer.collections.pop("LHE_weights")
    if lheWeightAna in sequence: sequence.remove(lheWeightAna)
    susyCounter.doLHE = False

## Auto-AAA
from CMGTools.RootTools.samples.autoAAAconfig import *
if not getHeppyOption("isCrab"):
    autoAAA(selectedComponents)

## output histogram
outputService=[]
from PhysicsTools.HeppyCore.framework.services.tfile import TFileService
output_service = cfg.Service(
    TFileService,
    'outputfile',
    name="outputfile",
    fname='treeProducerSusyMultilepton/tree.root',
    option='recreate'
    )    
outputService.append(output_service)


# the following is declared in case this cfg is used in input to the heppy.py script
from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
from CMGTools.TTHAnalysis.tools.EOSEventsWithDownload import EOSEventsWithDownload
event_class = EOSEventsWithDownload if not preprocessor else Events
EOSEventsWithDownload.aggressive = 2 # always fetch if running on Wigner
if getHeppyOption("nofetch") or getHeppyOption("isCrab"):
    event_class = Events
    if preprocessor: preprocessor.prefetch = False
config = cfg.Config( components = selectedComponents,
                     sequence = sequence,
                     services = outputService, 
                     preprocessor = preprocessor, 
                     events_class = event_class)
