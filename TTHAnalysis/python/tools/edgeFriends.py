## if MT2 doesnt work, put this line in the MT2 file: from ROOT.heppy import Davismt2

from CMGTools.TTHAnalysis.treeReAnalyzer import *
from CMGTools.TTHAnalysis.tools.eventVars_MT2 import *
print 'loading stuff for MT2'
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.AutoLibraryLoader.enable()
print 'done loading MT2 stuff.'
import time

#ROOT.gSystem.Load('libCondFormatsBTagObjects') 
ROOT.gSystem.Load("pluginRecoBTagPerformanceDBplugins.so")

#ROOT.gROOT.ProcessLine('.L /afs/cern.ch/work/p/pablom/private/run/pruebaestupida/CMSSW_7_4_12/src/CMGTools/TTHAnalysis/python/tools/BTagCalibrationStandalone.cc+') 

import copy
import math
class edgeFriends:
    def __init__(self,label,tightLeptonSel,cleanJet):
        self.label = "" if (label in ["",None]) else ("_"+label)
        self.tightLeptonSel = tightLeptonSel
        self.cleanJet = cleanJet
        #self.isMC = isMC
        ## with nvtx self.puFile = open("/afs/cern.ch/work/m/mdunser/public/puWeighting/puWeightsVinceLumi1p28.txt","r")
        ##self.puFile = open("/afs/cern.ch/work/m/mdunser/public/puWeighting/puWeightsOfficialPrescription.txt","r")
        ##self.pu_dict = eval(self.puFile.read())
        ##self.puFile.close()
        #self.puFile = ROOT.TFile("/afs/cern.ch/work/m/mdunser/public/puWeighting/2016/pileup_nominalUpDown.root","READ")
        self.puFile = ROOT.TFile("/afs/cern.ch/work/m/mdunser/public/puWeighting/2016/pileup_jul21_nominalUpDown.root","READ")
        self.puHist =copy.deepcopy( self.puFile.Get('weightsNominal') )
        self.puFile.Close()
        ##B-tagging stuff
        self.calib = ROOT.BTagCalibration("csvv2", "/afs/cern.ch/user/p/pablom/public/CSVv2_4invfb_systJuly15.csv")
        self.reader_heavy    = ROOT.BTagCalibrationReader(self.calib, 1, "mujets", "central")
        self.reader_heavy_UP = ROOT.BTagCalibrationReader(self.calib, 1, "mujets", "up")
        self.reader_heavy_DN = ROOT.BTagCalibrationReader(self.calib, 1, "mujets", "down")
        self.reader_light    = ROOT.BTagCalibrationReader(self.calib, 1, "incl"  , "central")
        self.reader_light_UP = ROOT.BTagCalibrationReader(self.calib, 1, "incl"  , "up")
        self.reader_light_DN = ROOT.BTagCalibrationReader(self.calib, 1, "incl"  , "down")

        self.calibFASTSIM = ROOT.BTagCalibration("csvv2", "/afs/cern.ch/user/p/pablom/public/CSV_13TEV_Combined_20_11_2015.csv")
        self.reader_heavyFASTSIM    = ROOT.BTagCalibrationReader(self.calibFASTSIM, 1, "fastsim", "central")
        self.reader_heavy_UPFASTSIM = ROOT.BTagCalibrationReader(self.calibFASTSIM, 1, "fastsim", "up")
        self.reader_heavy_DNFASTSIM = ROOT.BTagCalibrationReader(self.calibFASTSIM, 1, "fastsim", "down")
        self.reader_lightFASTSIM    = ROOT.BTagCalibrationReader(self.calibFASTSIM, 1, "fastsim"  , "central")
        self.reader_light_UPFASTSIM = ROOT.BTagCalibrationReader(self.calibFASTSIM, 1, "fastsim"  , "up")
        self.reader_light_DNFASTSIM = ROOT.BTagCalibrationReader(self.calibFASTSIM, 1, "fastsim"  , "down")

        self.f_btag_eff      = ROOT.TFile("/afs/cern.ch/user/p/pablom/public/btageff__ttbar_powheg_pythia8_25ns_16.root")
        self.h_btag_eff_b    = copy.deepcopy(self.f_btag_eff.Get("h2_BTaggingEff_csv_med_Eff_b"   ))
        self.h_btag_eff_c    = copy.deepcopy(self.f_btag_eff.Get("h2_BTaggingEff_csv_med_Eff_c"   ))
        self.h_btag_eff_udsg = copy.deepcopy(self.f_btag_eff.Get("h2_BTaggingEff_csv_med_Eff_udsg"))
        self.f_btag_eff.Close()

        # el sergio #######
        fSFMuon_FullFast_ID  = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/SFs/forICHEP/muons/FullSimFastSim/sf_mu_medium.root','read')
        fSFMuon_FullFast_ISO = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/SFs/forICHEP/muons/FullSimFastSim/sf_mu_mediumID_mini02.root','read')
        fSFMuon_FullFast_IP  = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/SFs/forICHEP/muons/FullSimFastSim/sf_mu_tightIP2D.root','read') 
        fSFMuon_FullData_ID  = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/SFs/forICHEP/muons/dataFullSim/TnP_MuonID_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root','read') 
        fSFMuon_FullData_ISO = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/SFs/forICHEP/muons/dataFullSim/TnP_MuonID_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root','read')
        fSFMuon_FullData_IP  = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/SFs/forICHEP/muons/dataFullSim/TnP_MuonID_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root','read')
        fSFMuon_Tracking     = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/SFs/forICHEP/muons/ratios_early.root')

        self.hMuonDataFull_ID  = copy.deepcopy(fSFMuon_FullData_ID .Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0"))
        self.hMuonDataFull_ISO = copy.deepcopy(fSFMuon_FullData_ISO.Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass"))
        self.hMuonDataFull_IP  = copy.deepcopy(fSFMuon_FullData_IP .Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass"))
        self.hMuonFullFast_ID  = copy.deepcopy(fSFMuon_FullFast_ID .Get("histo2D"))
        self.hMuonFullFast_ISO = copy.deepcopy(fSFMuon_FullFast_ISO.Get("histo2D"))
        self.hMuonFullFast_IP  = copy.deepcopy(fSFMuon_FullFast_IP .Get("histo2D"))
        self.hMuonTracking     = copy.deepcopy(fSFMuon_Tracking    .Get("ratio_eta"))

        fSFMuon_FullFast_ID .Close()
        fSFMuon_FullFast_ISO.Close()
        fSFMuon_FullFast_IP .Close()
        fSFMuon_FullData_ID .Close()
        fSFMuon_FullData_ISO.Close()
        fSFMuon_FullData_IP .Close()
        fSFMuon_Tracking    .Close()

        print 'muon', self.hMuonDataFull_ID, self.hMuonDataFull_ISO, self.hMuonDataFull_IP
        print 'muon', self.hMuonFullFast_ID, self.hMuonFullFast_ISO, self.hMuonFullFast_IP
        fSFElec_FullFast_ID  = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/SFs/forICHEP/electrons/FullFast/sf_el_tight2d3d.root') 
        fSFElec_FullFast_ISO = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/SFs/forICHEP/electrons/FullFast/sf_el_mini01.root')
        fSFElec_FullFast_IP  = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/SFs/forICHEP/electrons/FullFast/sf_el_inhit_eq0.root')  #IP in electrons is not IP but conversion veto (IP included in id scale factors)
        fSFElec_FullData     = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/SFs/forICHEP/electrons/scaleFactors.root')  
        fSFElec_Tracking     = ROOT.TFile('/afs/cern.ch/work/s/sesanche/public/SFs/forICHEP/electrons/egammaEffi.txt_SF2D.root')  

        self.hElecDataFull_ID  = copy.deepcopy( fSFElec_FullData.Get("GsfElectronToTight2D3D"))
        self.hElecDataFull_ISO = copy.deepcopy( fSFElec_FullData.Get("MVAVLooseElectronToMini"))
        self.hElecDataFull_IP  = copy.deepcopy( fSFElec_FullData.Get("MVATightElectronToConvIHit0"))
        self.hElecFullFast_ID  = copy.deepcopy( fSFElec_FullFast_ID .Get("histo2D"))
        self.hElecFullFast_ISO = copy.deepcopy( fSFElec_FullFast_ISO.Get("histo2D"))
        self.hElecFullFast_IP  = copy.deepcopy( fSFElec_FullFast_IP .Get("histo2D"))
        self.hElecTracking     = copy.deepcopy( fSFElec_Tracking    .Get("EGamma_SF2D"))


        fSFElec_FullFast_ID .Close()
        fSFElec_FullFast_ISO.Close()
        fSFElec_FullFast_IP .Close()
        fSFElec_FullData    .Close()

        print 'elec', self.hElecDataFull_ID, self.hElecDataFull_ISO, self.hElecDataFull_IP
        print 'elec', self.hElecFullFast_ID, self.hElecFullFast_ISO, self.hElecFullFast_IP


        ## =================
        ## pdf things
        ## =================
        self.an_file = ROOT.TFile("/afs/cern.ch/work/m/mdunser/public/pdfsForLikelihood/pdfs_version5_80X_2016Data_savingTheWorkspace_withSFPDFs_12p9invfb.root")
        ## file for 7.65 self.an_file = ROOT.TFile("/afs/cern.ch/work/m/mdunser/public/pdfsForLikelihood/pdfs_version4_80X_2016Data_savingTheWorkspace_withSFPDFs_7p65invfb.root")
        ## file used before topup to 7.65 self.an_file = ROOT.TFile("/afs/cern.ch/work/m/mdunser/public/pdfsForLikelihood/pdfs_version1_80X_2016Data_savingTheWorkspace.root")
        ## file with SF PDFs self.an_file = ROOT.TFile("/afs/cern.ch/work/m/mdunser/public/pdfsForLikelihood/pdfs_version3_80X_2016Data_savingTheWorkspace_withSFPDFs.root")
        self.wspace = copy.deepcopy( self.an_file.Get('w') )
        # data
        for t in ['DA','MC']:#,'MC_SF']:
            for var in [['mlb','sum_mlb_Edge'],['met','met_Edge'],['zpt','lepsZPt_Edge'],['ldr','lepsDR_Edge'],['a3d','d3D_Edge'],['ldp','lepsDPhi_Edge']]:
                print 'loading likelihoods for variable %s in %s'%(var[0],t)
                setattr(self,'h_lh_ana_%s_%s' %(var[0],t), self.wspace.pdf('%s_analyticalPDF_%s'%(var[0],t)))
                setattr(self,'var_ana_%s_%s'  %(var[0],t), self.wspace.var(var[1]))
                setattr(self,'frame_ana_%s_%s'%(var[0],t),getattr(self,'var_ana_%s_%s'%(var[0],t)).frame())
                getattr(self,'h_lh_ana_%s_%s' %(var[0],t)).plotOn(getattr(self,'frame_ana_%s_%s'%(var[0],t)))
                setattr(self,'obs_ana%s_%s'   %(var[0],t), ROOT.RooArgSet(self.wspace.var(var[1])))
        ## =================

        self.susymasslist = ['GenSusyMScan1'     , 'GenSusyMScan2'      , 'GenSusyMScan3'      , 'GenSusyMScan4'      ,
                             'GenSusyMGluino'    , 'GenSusyMGravitino'  , 'GenSusyMStop'       , 'GenSusyMSbottom'    ,
                             'GenSusyMStop2'     , 'GenSusyMSbottom2'   , 'GenSusyMSquark'     ,
                             'GenSusyMNeutralino', 'GenSusyMNeutralino2', 'GenSusyMNeutralino3', 'GenSusyMNeutralino4',
                             'GenSusyMChargino'  , 'GenSusyMChargino2']

        self.triggerlist = ['HLT_mu17mu8',
                            'HLT_mu17mu8_dz',
                            'HLT_mu17tkmu8',
                            'HLT_mu17tkmu8_dz',
                            'HLT_mu17el12',
                            'HLT_mu23el12',
                            'HLT_mu23el8',
                            'HLT_mu27tkmu8',
                            'HLT_el17el12_dz',
                            'HLT_el23el12_dz',
                            'HLT_mu8el17',
                            'HLT_mu8el23',
                            'HLT_mu30tkmu11_noniso',
                            'HLT_mu30ele30_noniso',
                            'HLT_doubleele33_noniso',
                            'HLT_doubleele33_MW_noniso',
                            'HLT_htall',
                            'HLT_atall',
                            'HLT_htmet' ]



        self.btagMediumCut = 0.800
        self.btagLooseCut  = 0.460

    def listBranches(self):
        label = self.label
        biglist = [ ("evt"+label, "D"),
                    ("run"+label, "I"),
                    ("lumi"+label, "I"),
                    ("nVert"+label, "I"),
                    ("nLepTight"+label, "I"),
                    ("nLepLoose"+label, "I"),
                    ("nJetSel"+label, "I"), ("nJetSel_jecUp"+label, "I"), ("nJetSel_jecDn"+label, "I"),
                    ("bestMjj"+label, "F"),
                    ("minMjj"+label, "F"),
                    ("maxMjj"+label, "F"),
                    ("hardMjj"+label, "F"),
                    ("hardJJDphi"+label, "F"),
                    ("hardJJDR"+label, "F"),
                    ("j1MetDPhi"+label, "F"),
                    ("j2MetDPhi"+label, "F"),
                    ## not really needed for much ("isLefthanded"+label, "I"),
                    ## not really needed for much ("isRighthanded"+label, "I"),
                    ("nPairLep"+label, "I"),
                    ("iLT"+label,"I",20,"nLepTight"+label), 
                    ("iJ"+label,"I",20,"nJetSel"+label), # index >= 0 if in Jet; -1-index (<0) if in DiscJet
                    ("nLepGood20"+label, "I"),
                    ("nLepGood20T"+label, "I"),
                    ("nJet35"+label        , "I") , ("nJet35_jecUp"+label        , "I") , ("nJet35_jecDn"+label        , "I") ,
                    ("htJet35j"+label      , "F") , ("htJet35j_jecUp"+label      , "F") , ("htJet35j_jecDn"+label      , "F") ,
                    ("nBJetMedium25"+label , "I") , ("nBJetMedium25_jecUp"+label , "I") , ("nBJetMedium25_jecDn"+label , "I") ,
                    ("nBJetLoose35"+label  , "I") , ("nBJetLoose35_jecUp"+label  , "I") , ("nBJetLoose35_jecDn"+label  , "I") ,
                    ("nBJetMedium35"+label , "I") , ("nBJetMedium35_jecUp"+label , "I") , ("nBJetMedium35_jecDn"+label , "I") ,
                    ("iL1T"+label, "I"),
                    ("iL2T"+label, "I"), 
                    ("lepsMll"+label, "F"),
                    ("lepsJZB"+label, "F"), 
                    ("lepsJZB_raw"+label, "F"),
                    ("lepsJZB_recoil"+label, "F"),
                    ("lepsDR"+label, "F"),
                    ("lepsMETRec"+label, "F"),
                    ("lepsZPt"+label, "F"),
                    ("metl1DPhi"+label, "F"),
                    ("metl2DPhi"+label, "F"),
                    ("met"+label, "F"), ("met_phi"+label, "F"), ("met_jecUp"+label, "F"), ("met_jecDn"+label, "F"), ("met_raw"+label, "F"),
                    ("lepsDPhi"+label, "F"),
                    ("Lep1_pt"+label, "F"), 
                    ("Lep1_eta"+label, "F"), 
                    ("Lep1_phi"+label, "F"),
                    ("Lep1_miniRelIso"+label, "F"),
                    ("Lep1_relIso03"+label, "F"),
                    ("Lep1_relIso04"+label, "F"),
                    ("Lep1_dxy"+label, "F"),
                    ("Lep1_dz"+label, "F"),
                    ("Lep1_sip3d"+label, "F"),
                    ("Lep1_pdgId"+label, "I"), 
                    ("Lep1_tightCharge"+label, "F"), 
                    ("Lep1_mvaIdSpring15"+label, "F"),
                    ("Lep1_mcMatchId"+label, "F"),
                    ("Lep1_minTauDR"+label, "F"),
                    ("Lep2_pt"+label, "F"), 
                    ("Lep2_eta"+label, "F"),
                    ("Lep2_phi"+label, "F"),
                    ("Lep2_miniRelIso"+label, "F"),
                    ("Lep2_relIso03"+label, "F"),
                    ("Lep2_relIso04"+label, "F"),
                    ("Lep2_dxy"+label, "F"),
                    ("Lep2_dz"+label, "F"),
                    ("Lep2_sip3d"+label, "F"),
                    ("Lep2_pdgId"+label, "I"),
                    ("Lep2_tightCharge"+label, "F"),
                    ("Lep2_mvaIdSpring15"+label, "F"),
                    ("Lep2_mcMatchId"+label, "F"),
                    ("Lep2_minTauDR"+label, "F"),
                    ("PileupW"+label, "F"), 
                    ("min_mlb1"+label, "F"),
                    ("min_mlb2"+label, "F"),
                    ("sum_mlb"+label, "F"), 
                    ("sum_mlb_jecUp"+label, "F"), 
                    ("sum_mlb_jecDn"+label, "F"), 
                    ("st"+label,"F"), 
                    ("srID"+label, "I"), 
                    ("mt2"+label, "F"), ("mt2_jecUp"+label, "F"), ("mt2_jecDn"+label, "F"),
                    ("lh_ana_zpt_data"+label, "F") ,
                    ("lh_ana_a3d_data"+label, "F") ,
                    ("lh_ana_met_data"+label, "F") ,
                    ("lh_ana_mlb_data"+label, "F") ,
                    ("lh_ana_mlb_jecUp_data"+label, "F") , 
                    ("lh_ana_mlb_jecDn_data"+label, "F") , 
                    ("lh_ana_ldr_data"+label, "F") ,
                    ("lh_ana_ldp_data"+label, "F") ,
                    ("lh_ana_zpt_mc"+label  , "F") , 
                    ("lh_ana_met_mc"+label  , "F") , 
                    ("lh_ana_mlb_mc"+label  , "F") , 
                    ("lh_ana_mlb_jecUp_mc"+label  , "F") , 
                    ("lh_ana_mlb_jecDn_mc"+label  , "F") , 
                    ("lh_ana_ldr_mc"+label  , "F") ,
                    ("lh_ana_a3d_mc"+label  , "F") ,
                    ("lh_ana_ldp_mc"+label  , "F") ,
                    ("lh_ana_zpt_mc_sf"+label  , "F") , 
                    ("lh_ana_met_mc_sf"+label  , "F") , 
                    ("lh_ana_mlb_mc_sf"+label  , "F") , 
                    ("lh_ana_mlb_jecUp_mc_sf"+label  , "F") , 
                    ("lh_ana_mlb_jecDn_mc_sf"+label  , "F") , 
                    ("lh_ana_ldr_mc_sf"+label  , "F") ,
                    ("lh_ana_a3d_mc_sf"+label  , "F") ,
                    ("lh_ana_ldp_mc_sf"+label  , "F") ,
                    ("nll"+label, "F"), ("nll_jecUp"+label, "F"), ("nll_jecDn"+label, "F"),
                    ("nll_mc"+label, "F"), ("nll_mc_jecUp"+label, "F"), ("nll_mc_jecDn"+label, "F"),
                    ("nll_mc_sf"+label, "F"),
                    ("weight_trigger"+label  , "F") ,
                    ("weight_btagsf"+label  , "F") ,
                    ("weight_btagsf_heavy_UP"+label, "F") ,
                    ("weight_btagsf_heavy_DN"+label, "F") ,
                    ("weight_btagsf_light_UP"+label, "F") ,
                    ("weight_btagsf_light_DN"+label, "F") ,
                    ("d3D" + label, "F"),
                    ("parPt" + label, "F"),
                    ("ortPt" + label, "F"),
                    ("dTheta" + label, "F"),
                    ('passesFilters' +label, 'I'),
                    ('genWeight' +label, 'F'),
                    ('weight_LepSF'+label,'F'),
                    ('weight_LepSF_MuUp'+label,'F'),
                    ('weight_LepSF_MuDn'+label,'F'),
                    ('weight_LepSF_ElUp'+label,'F'),
                    ('weight_LepSF_ElDn'+label,'F'),
                 ]
        for trig in self.triggerlist:
            biglist.append( ( '{tn}{lab}'.format(lab=label, tn=trig)) )
        for mass in self.susymasslist:
            biglist.append( ( '{tn}{lab}'.format(lab=label, tn=mass)) )
        ## for lfloat in 'pt eta phi miniRelIso pdgId'.split():
        ##     if lfloat == 'pdgId':
        ##         biglist.append( ("Lep"+label+"_"+lfloat,"I", 10, "nPairLep"+label) )
        ##     else:
        ##         biglist.append( ("Lep"+label+"_"+lfloat,"F", 10, "nPairLep"+label) )
        for jfloat in "pt eta phi mass btagCSV rawPt".split():
            biglist.append( ("JetSel"+label+"_"+jfloat,"F",20,"nJetSel"+label) )
        #if self.isMC:
        biglist.append( ("JetSel"+label+"_mcPt",     "F",20,"nJetSel"+label) )
        biglist.append( ("JetSel"+label+"_mcFlavour","I",20,"nJetSel"+label) )
        biglist.append( ("JetSel"+label+"_mcMatchId","I",20,"nJetSel"+label) )
        return biglist
    def __call__(self,event):
        t0 = time.time()
        leps  = [l for l in Collection(event,"LepGood","nLepGood")]
        lepso = [l for l in Collection(event,"LepOther","nLepOther")]
        jetsc = [j for j in Collection(event,"Jet","nJet")]
        jetsd = [j for j in Collection(event,"DiscJet","nDiscJet")]
        jetsc_jecUp = [j for j in Collection(event,"Jet_jecUp","nJet_jecUp")]
        jetsd_jecUp = [j for j in Collection(event,"DiscJet_jecUp","nDiscJet_jecUp")]
        jetsc_jecDn = [j for j in Collection(event,"Jet_jecDown","nJet_jecDown")]
        jetsd_jecDn = [j for j in Collection(event,"DiscJet_jecDown","nDiscJet_jecDown")]
        #metco = [m for m in Collection(event,"metcJet","nDiscJet")]
        (met, metphi)  = event.met_pt, event.met_phi
        (met_raw, metphi_raw)  = event.met_rawPt, event.met_rawPhi
        isData = event.isData
        if not isData:
            gentaus  = [t for t in Collection(event,"genTau","ngenTau")]
            ntrue = event.nTrueInt
        ## nvtx = event.nVert
        metp4 = ROOT.TLorentzVector()
        metp4.SetPtEtaPhiM(met,0,metphi,0)
        metp4_raw = ROOT.TLorentzVector()
        metp4_raw.SetPtEtaPhiM(met_raw,0,metphi_raw,0)
        ret = {}; jetret = {}; 
        lepret  = {}
        trigret = {}
        ret['met'] = met; ret['met_phi'] = metphi; ret['met_raw'] = met_raw;
        ret['met_jecUp'] = event.met_jecUp_pt; ret['met_jecDn'] = event.met_jecDown_pt 
        ret['run'] = event.run
        ret['lumi'] = event.lumi
        ret['evt'] = long(event.evt)
        ret['nVert'] = event.nVert
        
        t01 = time.time()
        ## copy the triggers, susy masses and filters!!
        for mass in self.susymasslist:
            ret[mass] = (-1 if not hasattr(event, mass) else getattr(event, mass) )
        self.isSMS =  (ret['GenSusyMScan1'] > 0)

        for trig in self.triggerlist:
            ##if not isData:
            ##    trigret[trig] = -1
            ##else:
            trigret[trig] = (-1 if not hasattr(event, trig) else getattr(event, trig) )
        t1 = time.time()

        ret['genWeight']          = ( 1. if not hasattr(event, 'genWeight'         ) else getattr(event, 'genWeight') )
        ## twiki:
        ##Flag_HBHENoiseFilter:Flag_HBHENoiseIsoFilter:Flag_EcalDeadCellTriggerPrimitiveFilter:Flag_goodVertices:Flag_eeBadScFilter:Flag_globalTightHalo2016Filter:badMuon:badChargedhadron
        if isData:
            ret['passesFilters'] = int(event.Flag_HBHENoiseFilter > 0 and 
                                    event.Flag_HBHENoiseIsoFilter > 0 and 
                                    event.Flag_EcalDeadCellTriggerPrimitiveFilter > 0 and 
                                    event.Flag_goodVertices > 0 and 
                                    event.Flag_eeBadScFilter > 0 and 
                                    event.Flag_globalTightHalo2016Filter > 0 and
                                    event.Flag_badChargedHadronFilter > 0 and
                                    event.Flag_badMuonFilter)# and 
                                    #event.Flag_METFilters)
        else:
            ret['passesFilters'] = 1

        ## this will be slow
        ## ret['isLefthanded' ] = 0
        ## ret['isRighthanded'] = 0
        ## if not isData:
        ##     genparts = [p for p in Collection(event,"GenPart","nGenPart")]
        ##     for p in genparts:
        ##         if p.status == 62:
        ##             if   abs(p.pdgId)-1000000 in [11, 13]: ret['isLefthanded' ] = 1
        ##             elif abs(p.pdgId)-2000000 in [11, 13]: ret['isRighthanded'] = 1
        ##             break

        t2 = time.time()

        #
        ### Define tight leptons
        ret["iLT"] = []; ret["nLepGood20T"] = 0

        # ====================
        # do pileupReweighting
        # ====================
        #puWt = self.pu_dict[int(ntrue)] if not isData else 1.
        puWt = self.puHist.GetBinContent(self.puHist.FindBin(ntrue)) if not isData else 1.
        #if puWt > 10: puWt = 10.
        ret["PileupW"] = puWt
        t21 = time.time()

        # ===============================
        # new, simpler sorting of leptons
        # ===============================
        nLepLoose = 0
        for il,lep in enumerate(leps):
            if not self._susyEdgeLoose(lep): continue
            nLepLoose+= 1
            if not self.tightLeptonSel(lep): continue
            ret["iLT"].append(il)
            ret["nLepGood20T"] += 1
        # other leptons, negative indices
        for il,lep in enumerate(lepso):
            if not self._susyEdgeLoose(lep): continue
            nLepLoose+= 1
            if not self.tightLeptonSel(lep): continue
            ret["iLT"].append(-1-il)
            ret["nLepGood20T"] += 1
        ret["nLepLoose"] = nLepLoose
        ret["nLepTight"] = len(ret["iLT"])

        t22 = time.time()
        #
        # sort the leptons by pT:
        ret["iLT"].sort(key = lambda idx : leps[idx].pt if idx >= 0 else lepso[-1-idx].pt, reverse = True)
        t23 = time.time()

        ## search for the lepton pair
        #lepst  = [ leps [il] for il in ret["iLT"] ]

        lepst = []
        for il in ret['iLT']:
            if il >=0: 
                lepst.append(leps[il])
            else: 
                lepst.append(lepso[-1-il])
        #
        iL1iL2 = self.getPairVariables(lepst, metp4, metp4_raw)
        t24 = time.time()
        ret['iL1T'] = ret["iLT"][ iL1iL2[0] ] if (len(ret["iLT"]) >=1 and iL1iL2[0] != -999) else -999
        ret['iL2T'] = ret["iLT"][ iL1iL2[1] ] if (len(ret["iLT"]) >=2 and iL1iL2[1] != -999) else -999
        ret['lepsMll'] = iL1iL2[2] 
        ret['lepsJZB'] = iL1iL2[3] 
        ret['lepsJZB_raw'] = iL1iL2[4] 
        ret['lepsDR'] = iL1iL2[5] 
        ret['lepsMETRec'] = iL1iL2[6] 
        ret['lepsZPt'] = iL1iL2[7] 
        ret['lepsDPhi'] = iL1iL2[8]
        ret['d3D']      = iL1iL2[9]
        ret['parPt']    = iL1iL2[10]
        ret['ortPt']    = iL1iL2[11]
        ret['dTheta']    = iL1iL2[12]
        t3 = time.time()

        #print 'new event =================================================='
        l1 = ROOT.TLorentzVector()
        l2 = ROOT.TLorentzVector()
        ltlvs = [l1, l2]
        lepvectors = []

        for lfloat in 'pt eta phi miniRelIso pdgId mvaIdSpring15 dxy dz sip3d relIso03 relIso04 tightCharge mcMatchId'.split():
            if lfloat == 'pdgId':
                lepret["Lep1_"+lfloat+self.label] = -99
                lepret["Lep2_"+lfloat+self.label] = -99
            else:
                lepret["Lep1_"+lfloat+self.label] = -42.
                lepret["Lep2_"+lfloat+self.label] = -42.
        if ret['iL1T'] != -999 and ret['iL2T'] != -999:
            ret['nPairLep'] = 2
            # compute the variables for the two leptons in the pair
            lcount = 1
            for idx in [ret['iL1T'], ret['iL2T']]:
                lep = leps[idx] if idx >= 0 else lepso[-1-idx]
                minDRTau = 99.
                if not isData:
                    for tau in gentaus:
                        tmp_dr = deltaR(lep, tau)
                        if tmp_dr < minDRTau:
                            minDRTau = tmp_dr
                for lfloat in 'pt eta phi miniRelIso pdgId mvaIdSpring15 dxy dz sip3d relIso03 relIso04 tightCharge mcMatchId'.split():
                    if lfloat == 'mcMatchId' and isData:
                        lepret["Lep"+str(lcount)+"_"+lfloat+self.label] = 1
                    else:
                        lepret["Lep"+str(lcount)+"_"+lfloat+self.label] = getattr(lep,lfloat)
                lepvectors.append(lep)
                lepret['metl'+str(lcount)+'DPhi'+self.label] = abs( deltaPhi( getattr(lep, 'phi'), metphi ))
                lepret["Lep"+str(lcount)+"_"+"minTauDR"+self.label] = minDRTau
                ltlvs[lcount-1].SetPtEtaPhiM(lep.pt, lep.eta, lep.phi, 0.0005 if lep.pdgId == 11 else 0.106)
                lcount += 1
                #print 'good lepton', getattr(lep,'pt'), getattr(lep,'eta'), getattr(lep,'phi'), getattr(lep,'pdgId') 
        else:
            ret['nPairLep'] = 0
        t4 = time.time()

        mt2 = -1.; mt2_jecUp = -1.; mt2_jecDn = -1.
        if ret['nPairLep'] == 2:
            l1mt2 = ROOT.reco.Particle.LorentzVector(lepvectors[0].p4().Px(), lepvectors[0].p4().Py(),lepvectors[0].p4().Pz(),lepvectors[0].p4().Energy())
            l2mt2 = ROOT.reco.Particle.LorentzVector(lepvectors[1].p4().Px(), lepvectors[1].p4().Py(),lepvectors[1].p4().Pz(),lepvectors[1].p4().Energy())
            metp4obj = ROOT.reco.Particle.LorentzVector(met*cos(metphi),met*sin(metphi),0,met)
            metp4obj_jecUp = ROOT.reco.Particle.LorentzVector(ret['met_jecUp']*cos(metphi),ret['met_jecUp']*sin(metphi),0,ret['met_jecUp'])
            metp4obj_jecDn = ROOT.reco.Particle.LorentzVector(ret['met_jecDn']*cos(metphi),ret['met_jecDn']*sin(metphi),0,ret['met_jecDn'])
            mt2       = computeMT2(l1mt2, l2mt2, metp4obj)
            mt2_jecUp = computeMT2(l1mt2, l2mt2, metp4obj_jecUp)
            mt2_jecDn = computeMT2(l1mt2, l2mt2, metp4obj_jecDn)
            del metp4obj, metp4obj_jecUp, metp4obj_jecDn
        ret['mt2'] = mt2
        ret['mt2_jecUp'] = mt2_jecUp
        ret['mt2_jecDn'] = mt2_jecDn
        t5 = time.time()
            
        ### Define jets
        ret["iJ"] = []
        jetsc       = self.setJetCollection(jetsc, lepst)      ; jetsd       = self.setJetCollection(jetsd, lepst);
        jetsc_jecUp = self.setJetCollection(jetsc_jecUp, lepst); jetsd_jecUp = self.setJetCollection(jetsd_jecUp, lepst);
        jetsc_jecDn = self.setJetCollection(jetsc_jecDn, lepst); jetsd_jecDn = self.setJetCollection(jetsd_jecDn, lepst);

        (ret["iJ"]      , nb25      , nb35      , nl35      , n35      , ht35      , theJets      , theBJets      ) = self.countJets(jetsc      , jetsd      )
        (ijlist_jecup   , nb25_jecUp, nb35_jecUp, nl35_jecUp, n35_jecUp, ht35_jecUp, theJets_jecUp, theBJets_jecUp) = self.countJets(jetsc_jecUp, jetsd_jecUp)
        (ijlist_jecdn   , nb25_jecDn, nb35_jecDn, nl35_jecDn, n35_jecDn, ht35_jecDn, theJets_jecDn, theBJets_jecDn) = self.countJets(jetsc_jecDn, jetsd_jecDn)

        ret['nJet35']        = n35 ; ret['nJet35_jecUp']        = n35_jecUp ; ret['nJet35_jecDn']        = n35_jecDn
        ret['nBJetMedium25'] = nb25; ret['nBJetMedium25_jecUp'] = nb25_jecUp; ret['nBJetMedium25_jecDn'] = nb25_jecDn
        ret['nBJetMedium35'] = nb35; ret['nBJetMedium35_jecUp'] = nb35_jecUp; ret['nBJetMedium35_jecDn'] = nb35_jecDn
        ret['nBJetLoose35']  = nl35; ret['nBJetLoose35_jecUp']  = nl35_jecUp; ret['nBJetLoose35_jecDn']  = nl35_jecDn
        ret["htJet35j"]      = ht35; ret["htJet35j_jecUp"]      = ht35_jecUp; ret["htJet35j_jecDn"]      = ht35_jecDn

        # 2. compute the jet list
        ret['nJetSel']       = len(ret["iJ"])
        ret['nJetSel_jecUp'] = len(ijlist_jecup)
        ret['nJetSel_jecDn'] = len(ijlist_jecdn)
        

        # 3. sort the jets by pt
        
        ret["iJ"].sort(key = lambda idx : jetsc[idx].pt if idx >= 0 else jetsd[-1-idx].pt, reverse = True)

        # 4. compute the variables
        
        for jfloat in "pt eta phi mass btagCSV rawPt".split():
            jetret[jfloat] = []
        if not isData:
            for jmc in "mcPt mcFlavour mcMatchId".split():
                jetret[jmc] = []
        for idx in ret["iJ"]:
            jet = jetsc[idx] if idx >= 0 else jetsd[-1-idx]
            for jfloat in "pt eta phi mass btagCSV rawPt".split():
                jetret[jfloat].append( getattr(jet,jfloat) )
            if not isData:
                for jmc in "mcPt mcFlavour mcMatchId".split():
                    jetret[jmc].append( getattr(jet,jmc) if not isData else -1.)
        t6 = time.time()
        
        totalRecoil = ROOT.TLorentzVector()
        for j in theJets:
            jet = ROOT.TLorentzVector()
            jet.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)
            totalRecoil = totalRecoil + jet
            

        ## compute mlb for the two lepton  

        theJets  = sorted(theJets , key = lambda j : j.pt, reverse = True)
        theBJets = sorted(theBJets, key = lambda j : j.pt, reverse = True)

        theJets_jecUp  = sorted(theJets_jecUp , key = lambda j : j.pt, reverse = True)
        theBJets_jecUp = sorted(theBJets_jecUp, key = lambda j : j.pt, reverse = True)
        theJets_jecDn  = sorted(theJets_jecDn , key = lambda j : j.pt, reverse = True)
        theBJets_jecDn = sorted(theBJets_jecDn, key = lambda j : j.pt, reverse = True)
        ret['lepsJZB_recoil'] = totalRecoil.Pt() - ret['lepsZPt']
        ret['bestMjj'] = self.getBestMjj(theJets)
        ret['minMjj']  = self.getMinMjj (theJets)
        ret['maxMjj']  = self.getMaxMjj (theJets)
        ret['hardMjj'] = self.getHardMjj(theJets)
        ret['hardJJDphi'] = self.getHardMjj(theJets, True)
        ret['hardJJDR'] = self.getHardMjj(theJets, True, True)
        ret['j1MetDPhi'] = deltaPhi(metphi, theJets[0].phi) if len(theJets) > 0 else -99.
        ret['j2MetDPhi'] = deltaPhi(metphi, theJets[1].phi) if len(theJets) > 1 else -99.
         
        t7 = time.time()
        [wtbtag, wtbtagUp_heavy, wtbtagUp_light, wtbtagDown_heavy, wtbtagDown_light] = (self.getWeightBtag(theJets) if not isData else [1., 1., 1., 1., 1.])

        ret['weight_trigger'] = 1.
        if not isData:
            if abs(lepret["Lep1_pdgId"+self.label] * lepret["Lep2_pdgId"+self.label]) == 169: ret['weight_trigger'] = 0.94
            if abs(lepret["Lep1_pdgId"+self.label] * lepret["Lep2_pdgId"+self.label]) == 143: ret['weight_trigger'] = 0.89
            if abs(lepret["Lep1_pdgId"+self.label] * lepret["Lep2_pdgId"+self.label]) == 121: ret['weight_trigger'] = 0.97

        ret['weight_btagsf'] = wtbtag
        ret['weight_btagsf_heavy_UP'] = wtbtagUp_heavy
        ret['weight_btagsf_heavy_DN'] = wtbtagDown_heavy
        ret['weight_btagsf_light_UP'] = wtbtagUp_light
        ret['weight_btagsf_light_DN'] = wtbtagDown_light
                                               
        [lep1SF, lep1SF_MuUp, lep1SF_MuDn, lep1SF_ElUp, lep1SF_ElDn] = (self.getLepSF(lepret["Lep1_eta"+self.label],lepret["Lep1_pt"+self.label],lepret["Lep1_pdgId"+self.label]) if not isData else [1., 1.,1.,1.,1.])
        [lep2SF, lep2SF_MuUp, lep2SF_MuDn, lep2SF_ElUp, lep2SF_ElDn] = (self.getLepSF(lepret["Lep2_eta"+self.label],lepret["Lep2_pt"+self.label],lepret["Lep2_pdgId"+self.label]) if not isData else [1., 1.,1.,1.,1.])
        ret['weight_LepSF']      = lep1SF      * lep2SF
        ret['weight_LepSF_MuUp'] = lep1SF_MuUp * lep2SF_MuUp
        ret['weight_LepSF_MuDn'] = lep1SF_MuDn * lep2SF_MuDn
        ret['weight_LepSF_ElUp'] = lep1SF_ElUp * lep2SF_ElUp
        ret['weight_LepSF_ElDn'] = lep1SF_ElDn * lep2SF_ElDn
        t8 = time.time()
        ##print 'njets: %.0d nbjets35medium: %.0d / %.0d'%(ret["nJet35"], len(theBJets), ret["nBJetMedium35"])
	
        jet = ROOT.TLorentzVector()
        min_mlb = 1e6; min_mlb_jecUp = 1e6; min_mlb_jecDn = 1e6
        max_mlb = 1e6; max_mlb_jecUp = 1e6; max_mlb_jecDn = 1e6
        _lmin, _jmin, _lmin_jecUp, _lmin_jecDn, _jmin_jecUp, _jmin_jecDn = -1, -1, -1, -1, -1, -1
        _lmax, _jmax, _lmax_jecUp, _lmax_jecDn, _jmax_jecUp, _jmax_jecDn = -1, -1, -1, -1, -1, -1
        leplist = [l1, l2]
        ## DO MLB CALCULATION HERE
        # find the global minimum mlb (or mlj)
        # new mlb calculations
        jet1coll = (theBJets if len(theBJets) >= 1 else theJets)
        jet2coll = (theBJets if len(theBJets) >= 2 else theJets)

        jet1coll_jecUp = (theBJets_jecUp if len(theBJets_jecUp) >= 1 else theJets_jecUp)
        jet2coll_jecUp = (theBJets_jecUp if len(theBJets_jecUp) >= 2 else theJets_jecUp)
        jet1coll_jecDn = (theBJets_jecDn if len(theBJets_jecDn) >= 1 else theJets_jecDn)
        jet2coll_jecDn = (theBJets_jecDn if len(theBJets_jecDn) >= 2 else theJets_jecDn)
        if ret['nPairLep'] > 1:
            for _il,lep in enumerate(leplist):
                for _ij,j in enumerate(jet1coll):
                    jet.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)           
                    tmp_mlb = (lep+jet).M()
                    if tmp_mlb < min_mlb:
                        min_mlb = tmp_mlb
                        _lmin = _il
                        _jmin = _ij

                for _ij, j in enumerate(jet1coll_jecUp):
                    jet.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)           
                    tmp_mlb = (lep+jet).M()
                    if tmp_mlb < min_mlb_jecUp:
                        min_mlb_jecUp = tmp_mlb
                        _lmin_jecUp = _il
                        _jmin_jecUp = _ij

                for _ij, j in enumerate(jet1coll_jecDn):
                    jet.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)           
                    tmp_mlb = (lep+jet).M()
                    if tmp_mlb < min_mlb_jecDn:
                        min_mlb_jecDn = tmp_mlb
                        _lmin_jecDn = _il
                        _jmin_jecDn = _ij

            for _il,lep in enumerate(leplist):
                if _il != _lmin:
                    for _ij,j in enumerate(jet2coll):
                        if len(theBJets) == 1 and j.btagCSV >= self.btagMediumCut:
                            continue
                        if (len(theBJets) == 0 or len(theBJets) >= 2) and _ij == _jmin: continue
                        jet.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)           
                        tmp_mlb = (lep+jet).M()
                        if tmp_mlb < max_mlb:
                            max_mlb = tmp_mlb
                            _lmax = _il
                            _jmax = _ij
                if _il != _lmin_jecUp:
                    for _ij,j in enumerate(jet2coll_jecUp):
                        if len(theBJets_jecUp) == 1 and j.btagCSV >= self.btagMediumCut: continue
                        if (len(theBJets_jecUp) == 0 or len(theBJets_jecUp) >= 2) and _ij == _jmin_jecUp: continue
                        jet.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)           
                        tmp_mlb = (lep+jet).M()
                        if tmp_mlb < max_mlb_jecUp:
                            max_mlb_jecUp = tmp_mlb
                            _lmax_jecUp = _il
                            _jmax_jecUp = _ij

                if _il != _lmin_jecDn:
                    for _ij,j in enumerate(jet2coll_jecDn):
                        if len(theBJets_jecDn) == 1 and j.btagCSV >= self.btagMediumCut: continue
                        if (len(theBJets_jecDn) == 0 or len(theBJets_jecDn) >= 2) and _ij == _jmin_jecDn: continue
                        jet.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)           
                        tmp_mlb = (lep+jet).M()
                        if tmp_mlb < max_mlb_jecDn:
                            max_mlb_jecDn = tmp_mlb
                            _lmax_jecDn = _il
                            _jmax_jecDn = _ij            
            
        ret["min_mlb1"] = min_mlb if min_mlb < 1e6  else -1.
        ret["min_mlb2"] = max_mlb if max_mlb < 1e6  else -1.
        min_mlb1_jecUp = min_mlb_jecUp if min_mlb_jecUp < 1e6  else -1.
        min_mlb2_jecUp = max_mlb_jecUp if max_mlb_jecUp < 1e6  else -1.
        min_mlb1_jecDn = min_mlb_jecDn if min_mlb_jecDn < 1e6  else -1.
        min_mlb2_jecDn = max_mlb_jecDn if max_mlb_jecDn < 1e6  else -1.

        ret["sum_mlb"] = (ret["min_mlb1"] + ret["min_mlb2"]) if ret["min_mlb1"] > 0. and ret["min_mlb2"] > 0. else -1.
        ret["sum_mlb_jecUp"] = (min_mlb1_jecUp  + min_mlb2_jecUp) if min_mlb1_jecUp > 0. and min_mlb2_jecUp > 0. else -1.
        ret["sum_mlb_jecDn"] = (min_mlb1_jecDn  + min_mlb2_jecDn) if min_mlb1_jecDn > 0. and min_mlb2_jecDn > 0. else -1.

        
        ret["st"] = met+lepret["Lep1_pt"+self.label]+lepret["Lep2_pt"+self.label]
        t9 = time.time()

        ## get the SR id which is 1xx for central and 2xx for forward. the 10 digit is the number of 
        ## b-tags and the signle digit is the mll region going from 1-5
        isBasicSREvent = (ret['nPairLep'] > 0 and ret["lepsDR"] > 0.1 and lepret["Lep1_pt"+self.label] > 20. and lepret["Lep2_pt"+self.label] > 20. and ret['lepsMll'] > 20.)
        isBasicSREvent = isBasicSREvent * (abs(lepret["Lep1_eta"+self.label] - 1.5) > 0.1 and abs(lepret["Lep2_eta"+self.label] - 1.5) > 0.1)
        isBasicSREvent = isBasicSREvent * (met > 150 and ret['nJetSel'] >= 2 ) ## adapting to 2016 selection
        
        if isBasicSREvent:
            srID = self.getSRID(ret['lepsMll'], lepret["Lep1_eta"+self.label], lepret["Lep2_eta"+self.label], ret["nBJetMedium35"])
            ret["srID"] = srID
            for t in  ['data','mc']:#, 'mc_sf']:
                if t == 'data'  : nam = 'DA'
                if t == 'mc'    : nam = 'MC'
                if t == 'mc_sf' : nam = 'MC_SF'
                for u in ['_ana']:
                    
                    for var in [['mlb',ret['sum_mlb'],'sum_mlb_Edge'],['met',met,'met_Edge'],
                                ['zpt',ret['lepsZPt'],'lepsZPt_Edge'],['ldr',ret['lepsDR'],'lepsDR_Edge'],
                                ['a3d',ret['d3D'],'d3D_Edge'],['ldp',ret['lepsDPhi'],'lepsDPhi_Edge']]:
                                
                        
                        self.wspace.var(var[2]).setVal(var[1])
                        ret["lh%s_%s_%s"%(u,var[0],t)] = getattr(self,'h_lh_ana_%s_%s'%(var[0],nam)).getVal(getattr(self,'obs_ana%s_%s'%(var[0],nam)))

                    self.wspace.var('sum_mlb_Edge').setVal(ret['sum_mlb_jecUp'])
                    ret["lh%s_mlb_jecUp_%s"%(u,t)] = getattr(self, 'h_lh_ana_mlb_%s'%nam).getVal(getattr(self,'obs_anamlb_%s'%nam))
                    self.wspace.var('sum_mlb_Edge').setVal(ret['sum_mlb_jecDn'])
                    ret["lh%s_mlb_jecDn_%s"%(u,t)] = getattr(self, 'h_lh_ana_mlb_%s'%nam).getVal(getattr(self,'obs_anamlb_%s'%nam))

                    if not ret["lh%s_mlb_%s"%(u,t)]: ret["lh%s_mlb_%s"%(u,t)] = 1e-50
                    if not ret["lh%s_mlb_jecUp_%s"%(u,t)]: ret["lh%s_mlb_jecUp_%s"%(u,t)] = 1e-50
                    if not ret["lh%s_mlb_jecDn_%s"%(u,t)]: ret["lh%s_mlb_jecDn_%s"%(u,t)] = 1e-50
                    if not ret["lh%s_ldr_%s"%(u,t)]: ret["lh%s_ldr_%s"%(u,t)] = 1e-50
                    if not ret["lh%s_met_%s"%(u,t)]: ret["lh%s_met_%s"%(u,t)] = 1e-50
                    if not ret["lh%s_zpt_%s"%(u,t)]: ret["lh%s_zpt_%s"%(u,t)] = 1e-50
                    if not ret["lh%s_a3d_%s"%(u,t)]: ret["lh%s_a3d_%s"%(u,t)] = 1e-50
                    if not ret["lh%s_ldp_%s"%(u,t)]: ret["lh%s_ldp_%s"%(u,t)] = 1e-50
            
            ret['nll']       = -1.*math.log(ret["lh_ana_mlb_data"] *ret["lh_ana_met_data"] *ret["lh_ana_zpt_data"] *ret["lh_ana_ldp_data"] )
            ret['nll_jecUp'] = -1.*math.log(ret["lh_ana_mlb_jecUp_data"] *ret["lh_ana_met_data"] *ret["lh_ana_zpt_data"] *ret["lh_ana_ldp_data"] )
            ret['nll_jecDn'] = -1.*math.log(ret["lh_ana_mlb_jecDn_data"] *ret["lh_ana_met_data"] *ret["lh_ana_zpt_data"] *ret["lh_ana_ldp_data"] )
            ret['nll_mc']    = -1.*math.log(ret["lh_ana_mlb_mc"]   *ret["lh_ana_met_mc"]   *ret["lh_ana_zpt_mc"]   *ret["lh_ana_ldp_mc"]   )
            #ret['nll_mc_sf'] = -1.*math.log(ret["lh_ana_mlb_mc_sf"]*ret["lh_ana_met_mc_sf"]*ret["lh_ana_zpt_mc_sf"]*ret["lh_ana_ldp_mc_sf"])
            # print 'min mlb', ret["sum_mlb"], ret["sum_mlb_jecUp"], ret["sum_mlb_jecDn"]
            # print 'nll    ', ret['nll'], ret['nll_jecUp'], ret['nll_jecDn']
            # print 'the bje', len(theBJets_jecUp), len(theBJets_jecUp), len(theBJets_jecDn)
            # print 'jet min', _jmin, _jmin_jecUp, _jmin_jecDn
            # print 'jet max', _jmax, _jmax_jecUp, _jmax_jecDn
        else:
            ret["srID"]      = -99
            for t in ['data', 'mc']:#, 'mc_sf']:
                for u in ['_ana']:
                    ret["lh%s_mlb_%s"%(u,t)] = -999.
                    ret["lh%s_ldr_%s"%(u,t)] = -999.
                    ret["lh%s_met_%s"%(u,t)] = -999.
                    ret["lh%s_zpt_%s"%(u,t)] = -999.
                    ret["lh%s_a3d_%s"%(u,t)] = -999.
                    ret["lh%s_ldp_%s"%(u,t)] = -999.
            ret['nll']       = 0.
            ret['nll_mc']    = 0.
            ret['nll_mc_sf'] = 0.
        t10 = time.time()

        ## print 'time from start to pre trigloaded: %.3f mus'%( (t01-t0)*1000000. )
        ## print 'time from pretrig to posttrigload: %.3f mus'%( (t1 -t01)*1000000. )
        ## print 'time from trigger loaded to l-r-h: %.3f mus'%( (t2 -t1)*1000000. )
        ## print '  time for puW: %.3f'%( (t21-t2)*1000000.  )
        ## print '  time for lepstuff: %.3f'%( (t22-t21)*1000000.  )
        ## print '  time for lepsort : %.3f'%( (t23-t22)*1000000.  )
        ## print '  time for npairlep: %.3f'%( (t24-t23)*1000000.  )
        ## print 'time from l-r-h to done with leps: %.3f mus'%( (t3 -t2)*1000000. )
        ## print 'time from done with lep to premt2: %.3f mus'%( (t4 -t3)*1000000. )
        ## print 'time from premt2 to post mt2     : %.3f mus'%( (t5 -t4)*1000000. )
        ## print 'time from post mt2 to jet filled : %.3f mus'%( (t6 -t5)*1000000. )
        ## print 'time from jet filled to pre btag : %.3f mus'%( (t7 -t6)*1000000. )
        ## print 'time from prebtag to post btag   : %.3f mus'%( (t8 -t7)*1000000. )
        ## print 'time from post btag to poost mlbe: %.3f mus'%( (t9 -t8)*1000000. )
        ## print 'time from post mlb etc to post LH: %.3f mus'%( (t10-t9)*1000000. )
        
        fullret = {}
        for k,v in ret.iteritems(): 
            fullret[k+self.label] = v
        for k,v in jetret.iteritems(): 
            fullret["JetSel%s_%s" % (self.label,k)] = v
        #for k,v in lepret.iteritems(): 
        #    fullret["Lep%s_%s" % (self.label,k)] = v
        for k,v in lepret.iteritems(): 
            fullret[k] = v
        for k,v in trigret.iteritems(): 
            fullret[k+self.label] = v
        return fullret

    def setJetCollection(self, jetcoll, lepst):
        for j in jetcoll:
            j._clean = True
            if abs(j.eta) > 2.4 or j.pt < 25.:
                j._clean = False
                continue
            if j.pt < 35 and j.btagCSV < self.btagMediumCut: 
                j._clean = False
                continue
            for l in lepst:
                #lep = leps[l]
                if deltaR(l,j) < 0.4:
                    j._clean = False
        return jetcoll

    def countJets(self, coll1, coll2):
        nb25 = 0; nb25 = 0; nb35 = 0; ht35 = 0.; nl35 = 0; n35 = 0
        thejets = []; thebjets = []
        retlist = []
        for ijc,j in enumerate(coll1):
            if not j._clean: continue
            bt = j.btagCSV
            pt = j.pt
            if pt > 25 and bt > self.btagMediumCut: nb25 += 1
            #if pt > 25 and bt > self.btagLooseCut : nl25 += 1
            #if pt > 35 and bt > self.btagMediumCut: nb35 += 1
            if pt > 35:
                thejets.append(j)
                n35 += 1; ht35 += pt
                retlist.append(ijc)
                if bt > self.btagMediumCut:
                    nb35 += 1
                    thebjets.append(j)
                if bt > self.btagLooseCut:
                    nl35 += 1
        for ijd,j in enumerate(coll2):
            if not j._clean: continue
            bt = j.btagCSV
            pt = j.pt
            if pt > 25 and bt > self.btagMediumCut: nb25 += 1
            #if pt > 25 and bt > self.btagLooseCut : nl25 += 1
            #if pt > 35 and bt > self.btagMediumCut: nb35 += 1
            if pt > 35:
                thejets.append(j)
                n35 += 1; ht35 += pt
                retlist.append(-1-ijd)
                if bt > self.btagMediumCut:
                    nb35 += 1
                    thebjets.append(j)
                if bt > self.btagLooseCut:
                    nl35 += 1
        return retlist, nb25, nb35, nl35, n35, ht35, thejets, thebjets

    def getMll_JZB(self, l1, l2, met, met_raw):
        metrecoil = (met + l1 + l2).Pt()
        metrawrecoil = (met_raw + l1 + l2).Pt() 
        zpt = (l1 + l2).Pt()
        jzb = metrecoil - zpt
        jzb_raw = metrawrecoil - zpt
        v1 = l1.Vect()
        v2 = l2.Vect()
        return ((l1+l2).M(), jzb, jzb_raw, l1.DeltaR(l2), metrecoil, zpt, abs( deltaPhi( l1.Phi(), l2.Phi() )) , v1.Angle(v2))
    def getParOrtPt(self, l1, l2):
        if l1.Pt() > l2.Pt():
            v1 = l1.Vect()
            v2 = l2.Vect()
        else:
            v1 = l2.Vect()
            v2 = l1.Vect()
        u1 = v1.Unit()                              # direction of the harder lepton
        p1 = math.cos(v1.Angle(v2)) * v2.Mag() * u1 # projection of the softer lepton onto the harder
        o1 = v1 - p1                                # orthogonal to the projection of the softer onto the harder
        return  (p1.Perp(), o1.Perp())

    def getPairVariables(self,lepst, metp4, metp4_raw):
        ret = (-999,-999,-99., -9000., -9000, -99., -99., -99., -99.,-99.,-99.,-99.,-99.,-99.)
        if len(lepst) >= 2:
            [mll, jzb, jzb_raw, dr, metrec, zpt, dphi, d3D] = self.getMll_JZB(lepst[0].p4(), lepst[1].p4(), metp4, metp4_raw)
            [parPt, ortPt] = self.getParOrtPt(lepst[0].p4(),lepst[1].p4())
            ret = (0, 1, mll, jzb, jzb_raw, dr, metrec, zpt, dphi, d3D, parPt, ortPt, lepst[0].p4().Theta() - lepst[1].p4().Theta())
        return ret

    def getSRID(self, mll, eta1, eta2, nb):
        mllid, bid, etaid = -1, -1, -1
        if    20. < mll <  70.:
            mllid = 1
        elif  70. < mll <  81.:
            mllid = 2
        elif  81. < mll < 101.:
            mllid = 3
        elif 101. < mll < 120.:
            mllid = 4
        elif 120. < mll:
            mllid = 5
            
        if abs(eta1) < 1.4 and abs(eta2) < 1.4:
            etaid = 1
        else:
            etaid = 2

        return (100*etaid + 10*nb + mllid)

    def getBestMjj(self, jetsel):
        if len(jetsel) < 2: return -99.
        bestmjj = 1e6
        for jeti in jetsel:
            for jetj in jetsel:
                if jeti == jetj: continue
                jet1 = ROOT.TLorentzVector()
                jet1.SetPtEtaPhiM(jeti.pt, jeti.eta, jeti.phi, jeti.mass)
                jet2 = ROOT.TLorentzVector()
                jet2.SetPtEtaPhiM(jetj.pt, jetj.eta, jetj.phi, jetj.mass)
                dijetmass = (jet1+jet2).M()
                if abs(dijetmass - 80.385) < abs(bestmjj - 80.385):
                    bestmjj = dijetmass
        return bestmjj
    def getMinMjj(self, jetsel):
        if len(jetsel) < 2: return -99.
        minmjj = 1e6
        for jeti in jetsel:
            for jetj in jetsel:
                if jeti == jetj: continue
                jet1 = ROOT.TLorentzVector()
                jet1.SetPtEtaPhiM(jeti.pt, jeti.eta, jeti.phi, jeti.mass)
                jet2 = ROOT.TLorentzVector()
                jet2.SetPtEtaPhiM(jetj.pt, jetj.eta, jetj.phi, jetj.mass)
                dijetmass = (jet1+jet2).M()
                if dijetmass < minmjj:
                    minmjj = dijetmass
        return minmjj
    def getMaxMjj(self, jetsel):
        if len(jetsel) < 2: return -99.
        maxmjj = -99.
        for jeti in jetsel:
            for jetj in jetsel:
                if jeti == jetj: continue
                jet1 = ROOT.TLorentzVector()
                jet1.SetPtEtaPhiM(jeti.pt, jeti.eta, jeti.phi, jeti.mass)
                jet2 = ROOT.TLorentzVector()
                jet2.SetPtEtaPhiM(jetj.pt, jetj.eta, jetj.phi, jetj.mass)
                dijetmass = (jet1+jet2).M()
                if dijetmass > maxmjj:
                    maxmjj = dijetmass
        return maxmjj
    def getHardMjj(self, jetsel, _dphi = False, _dr = False):
        if len(jetsel) < 2: return -99.
        if not _dphi:
            jet1 = ROOT.TLorentzVector()
            jet2 = ROOT.TLorentzVector()
            jet1.SetPtEtaPhiM(jetsel[0].pt, jetsel[0].eta, jetsel[0].phi, jetsel[0].mass)
            jet2.SetPtEtaPhiM(jetsel[1].pt, jetsel[1].eta, jetsel[1].phi, jetsel[1].mass)
            retval = (jet1+jet2).M()
        else:         
            if not _dr: retval = deltaPhi( jetsel[0].phi, jetsel[1].phi)
            else:       retval = deltaR(jetsel[0], jetsel[1])
        return retval

   
    #############Pablin
    def get_SF_btag(self, pt, eta, mcFlavour):

       flavour = 2
       if abs(mcFlavour) == 5: flavour = 0
       elif abs(mcFlavour)==4: flavour = 1
  
       pt_cutoff  = max(30. , min(669., pt))
       eta_cutoff = min(2.39, abs(eta))

       if flavour == 2:
          SF = self.reader_light.eval(flavour,eta_cutoff, pt_cutoff);
          SFup = self.reader_light_UP.eval(flavour,eta_cutoff, pt_cutoff);
          SFdown = self.reader_light_DN.eval(flavour,eta_cutoff, pt_cutoff);
          SFcorr = self.reader_lightFASTSIM.eval(flavour,eta_cutoff, pt_cutoff);
          SFupcorr = self.reader_light_UPFASTSIM.eval(flavour,eta_cutoff, pt_cutoff);
          SFdowncorr = self.reader_light_DNFASTSIM.eval(flavour,eta_cutoff, pt_cutoff);
       else:
          SF = self.reader_heavy.eval(flavour,eta_cutoff, pt_cutoff)
          SFup  = self.reader_heavy_UP.eval(flavour,eta_cutoff, pt_cutoff)
          SFdown = self.reader_heavy_DN.eval(flavour,eta_cutoff, pt_cutoff)
          SFcorr = self.reader_heavyFASTSIM.eval(flavour,eta_cutoff, pt_cutoff)
          SFupcorr  = self.reader_heavy_UPFASTSIM.eval(flavour,eta_cutoff, pt_cutoff)
          SFdowncorr = self.reader_heavy_DNFASTSIM.eval(flavour,eta_cutoff, pt_cutoff)

       return [SF*SFcorr, SFup*SFupcorr, SFdown*SFdowncorr]


    def getBtagEffFromFile(self, pt, eta, mcFlavour):
    
       pt_cutoff = max(20.,min(399., pt))
       if (abs(mcFlavour) == 5): 
           h = self.h_btag_eff_b
           #use pt bins up to 600 GeV for b
           pt_cutoff = max(20.,min(599., pt))
       elif (abs(mcFlavour) == 4):
           h = self.h_btag_eff_c
       else:
           h = self.h_btag_eff_udsg
    
       binx = h.GetXaxis().FindBin(pt_cutoff)
       biny = h.GetYaxis().FindBin(fabs(eta))

       return h.GetBinContent(binx,biny)


    def getWeightBtag(self, jets):

        mcTag = 1.
        mcNoTag = 1.
        dataTag = 1.
        dataNoTag = 1.
        errHup   = 0
        errHdown = 0
        errLup   = 0
        errLdown = 0
    
        for jet in jets:
        
            csv = jet.btagCSV
            mcFlavor = (jet.hadronFlavour if hasattr(jet, 'hadronFlavour') else jet.mcFlavour)
            eta = jet.eta
            pt = jet.pt

            if(eta > 2.5): continue
            if(pt < 20): continue
            eff = self.getBtagEffFromFile(pt, eta, mcFlavor)

            istag = csv > self.btagMediumCut and eta < 2.5 and pt > 20
            SF = self.get_SF_btag(pt, eta, mcFlavor)
            if(istag):
                 mcTag = mcTag * eff
                 dataTag = dataTag * eff * SF[0]
                 if(mcFlavor == 5 or mcFlavor ==4):
	             errHup  = errHup + (SF[1] - SF[0]  )/SF[0]
	             errHdown = errHdown = (SF[0] - SF[2])/SF[0]
                 else: 
	             errLup = errLup + (SF[1] - SF[0])/SF[0]
                     errLdown = errLdown + (SF[0] - SF[2])/SF[0]
            else: 
                 mcNoTag = mcNoTag * (1 - eff)
                 dataNoTag = dataNoTag * (1 - eff*SF[0])
                 if mcFlavor==5 or mcFlavor==4:
	             errHup = errHup * -eff*(SF[1] - SF[0]  )/(1-eff*SF[0])
                     errHdown = errHdown * -eff*(SF[0] - SF[2])/(1-eff*SF[0])	
                 else:
	             errLup = errLup * -eff*(SF[1] - SF[0]  )/(1-eff*SF[0])
	             errLdown = errLdown * -eff*(SF[0] - SF[2])/(1-eff*SF[0]);	


        wtbtag = (dataNoTag * dataTag ) / ( mcNoTag * mcTag )
        wtbtagUp_heavy   = wtbtag*( 1 + errHup   )
        wtbtagUp_light   = wtbtag*( 1 + errLup   )
        wtbtagDown_heavy = wtbtag*( 1 - errHdown )
        wtbtagDown_light = wtbtag*( 1 - errLdown )

        return [wtbtag, wtbtagUp_heavy, wtbtagUp_light, wtbtagDown_heavy, wtbtagDown_light]

    def selfNewMediumMuonId(self, muon):
        if not hasattr(muon, 'isGlobalMuon'):
            return (muon.mediumMuonId == 1)
        goodGlob = (muon.isGlobalMuon and 
                    muon.globalTrackChi2 < 3 and
                    muon.chi2LocalPosition < 12 and
                    muon.trkKink < 20)
        isMedium = (muon.innerTrackValidHitFraction > 0.49 and
                    muon.segmentCompatibility > (0.303 if goodGlob else  0.451) )
        return isMedium
        #muon.segmentCompatibility < 0.49: return False

    def _susyEdgeLoose(self, lep):
            if lep.pt <= 10.: return False
            if abs(lep.dxy) > 0.05: return False
            if abs(lep.dz ) > 0.10: return False
            if lep.sip3d > 8: return False
            lepeta = abs(lep.eta)
            if lep.miniRelIso > 0.4: return False
            ## muons
            if abs(lep.pdgId) == 13:
              if lepeta > 2.4: return False
              #if lep.mediumMuonId != 1: return False
              if not self.selfNewMediumMuonId(lep): return False
            ## electrons
            if abs(lep.pdgId) == 11:
              if lepeta > 2.5: return False
              if (lep.convVeto == 0) or (lep.lostHits > 0) : return False
              if (lepeta < 0.8   and lep.mvaIdSpring15 < -0.70) : return False
              if (lepeta > 0.8   and lepeta < 1.479 and lep.mvaIdSpring15 < -0.83) : return False
              if (lepeta > 1.479 and lep.mvaIdSpring15 < -0.92) : return False
              if hasattr(lep, 'idEmuTTH'):
                if lep.idEmuTTH == 0: return False
            return True

    def getLepSF(self, eta, pt, id):
        result = [1.,1.,1.,1.,1.]
        if abs(id) == 11:
            #electrones
            pt  = min(199,pt)
            sf1 = self.hElecDataFull_ID .GetBinContent(self.hElecDataFull_ID .FindBin(pt, abs(eta)))
            sf2 = self.hElecDataFull_ISO.GetBinContent(self.hElecDataFull_ISO.FindBin(pt, abs(eta)))
            sf3 = self.hElecDataFull_IP .GetBinContent(self.hElecDataFull_IP .FindBin(pt, abs(eta)))
            sf7 = self.hElecTracking    .GetBinContent(self.hElecTracking    .FindBin(eta, pt)     ) 

            sf1_e = self.hElecDataFull_ID .GetBinError(self.hElecDataFull_ID .FindBin(pt, abs(eta)))
            sf2_e = self.hElecDataFull_ISO.GetBinError(self.hElecDataFull_ISO.FindBin(pt, abs(eta)))
            sf3_e = self.hElecDataFull_IP .GetBinError(self.hElecDataFull_IP .FindBin(pt, abs(eta)))
            sf7_e = self.hElecTracking    .GetBinError(self.hElecTracking    .FindBin(eta,pt)      )

            # fullsim by defult
            sf4 = 1.
            sf5 = 1.
            sf6 = 1.
            sf4_e = 0.
            sf5_e = 0.
            sf6_e = 0.

            # if fastsim
            if self.isSMS:
                sf4 = self.hElecFullFast_ID .GetBinContent(self.hElecFullFast_ID .FindBin(pt, abs(eta)))
                sf5 = self.hElecFullFast_ISO.GetBinContent(self.hElecFullFast_ISO.FindBin(pt, abs(eta)))
                sf6 = self.hElecFullFast_IP .GetBinContent(self.hElecFullFast_IP .FindBin(pt, abs(eta)))
                sf4_e = 0.02 * sf4 # as prescribed in the twiki
                sf5_e = 0.02 * sf5
                sf6_e = 0.02 * sf6

            elVar = sqrt( (sf1_e*sf2*sf3*sf4*sf5*sf6*sf7)**2 + (sf1*sf2_e*sf3*sf4*sf5*sf6*sf7)**2
                          +(sf1*sf2*sf3_e*sf4*sf5*sf6*sf7)**2 +(sf1*sf2*sf3*sf4_e*sf5*sf6*sf7)**2
                          +(sf1*sf2*sf3*sf4*sf5_e*sf6*sf7)**2 +(sf1*sf2*sf3*sf4*sf5*sf6_e*sf7)**2
                          +(sf1*sf2*sf3*sf4*sf5*sf6*sf7_e)**2 )
            elSF  = sf1*sf2*sf3*sf4*sf5*sf6*sf7

            return [ elSF, elSF, elSF, elSF+elVar, elSF-elVar]

        else:
            pt  = min(119, pt)
            sf1 = self.hMuonDataFull_ID .GetBinContent(self.hMuonDataFull_ID .FindBin(pt, abs(eta)))
            sf2 = self.hMuonDataFull_ISO.GetBinContent(self.hMuonDataFull_ISO.FindBin(pt, abs(eta)))
            sf3 = self.hMuonDataFull_IP .GetBinContent(self.hMuonDataFull_IP .FindBin(pt, abs(eta)))
#            sf7 = self.hMuonTracking    .GetBinContent(self.hMuonTracking    .FindBin(eta, pt)     )
            sf7 = self.hMuonTracking    .Eval(eta)

            sf1_e = self.hMuonDataFull_ID .GetBinError(self.hMuonDataFull_ID .FindBin(pt, abs(eta)))
            sf2_e = self.hMuonDataFull_ISO.GetBinError(self.hMuonDataFull_ISO.FindBin(pt, abs(eta)))
            sf3_e = self.hMuonDataFull_IP .GetBinError(self.hMuonDataFull_IP .FindBin(pt, abs(eta)))
#            sf7_e = self.hMuonTracking    .Eval(eta)
            sf7_e = 0
            # fullsim by defult
            sf4 = 1.
            sf5 = 1.
            sf6 = 1.
            sf4_e = 0.
            sf5_e = 0.
            sf6_e = 0.
            if self.isSMS:
                sf4 = self.hMuonFullFast_ID .GetBinContent(self.hMuonFullFast_ID .FindBin(pt, abs(eta)))
                sf5 = self.hMuonFullFast_ISO.GetBinContent(self.hMuonFullFast_ISO.FindBin(pt, abs(eta)))
                sf6 = self.hMuonFullFast_IP .GetBinContent(self.hMuonFullFast_IP .FindBin(pt, abs(eta)))
                sf4_e = 0.02 * sf4 # as prescribed in the twiki
                sf5_e = 0.02 * sf5
                sf6_e = 0.02 * sf6

            muVar = sqrt( (sf1_e*sf2*sf3*sf4*sf5*sf6*sf7)**2 + (sf1*sf2_e*sf3*sf4*sf5*sf6*sf7)**2
                          +(sf1*sf2*sf3_e*sf4*sf5*sf6*sf7)**2 +(sf1*sf2*sf3*sf4_e*sf5*sf6*sf7)**2
                          +(sf1*sf2*sf3*sf4*sf5_e*sf6*sf7)**2 +(sf1*sf2*sf3*sf4*sf5*sf6_e*sf7)**2
                          +(sf1*sf2*sf3*sf4*sf5*sf6*sf7_e)**2 )
            muSF  = sf1*sf2*sf3*sf4*sf5*sf6*sf7
            return [ muSF, muSF+muVar, muSF-muVar, muSF, muSF]

            
def newMediumMuonId(muon):
    if not hasattr(muon, 'isGlobalMuon'):
        return (muon.mediumMuonId == 1)
    goodGlob = (muon.isGlobalMuon and 
                muon.globalTrackChi2 < 3 and
                muon.chi2LocalPosition < 12 and
                muon.trkKink < 20)
    isMedium = (muon.innerTrackValidHitFraction > 0.49 and
                muon.segmentCompatibility > (0.303 if goodGlob else  0.451) )
    return isMedium
    #muon.segmentCompatibility < 0.49: return False
 
def _susyEdgeTight(lep):
        if lep.pt <= 20.: return False
        eta = abs(lep.eta)
        if eta          > 2.4: return False
        if abs(lep.dxy) > 0.05: return False
        if abs(lep.dz ) > 0.10: return False
        if eta > 1.4 and eta < 1.6: return False
        if abs(lep.pdgId) == 13:
          ## old medium ID if lep.mediumMuonId != 1: return False
          if not newMediumMuonId(lep): return False
          if lep.miniRelIso > 0.2: return False
        #if abs(lep.pdgId) == 11 and (lep.tightId < 1 or (abs(lep.etaSc) > 1.4442 and abs(lep.etaSc) < 1.566)) : return False
        if abs(lep.pdgId) == 11:
          etatest = (abs(lep.etaSc) if hasattr(lep, 'etaSc') else abs(lep.eta))
          if (etatest > 1.4442 and etatest < 1.566) : return False
          if (lep.convVeto == 0) or (lep.lostHits > 0) : return False
          if (eta < 0.8 and lep.mvaIdSpring15 < 0.87) : return False
          if (eta > 0.8 and eta < 1.479 and lep.mvaIdSpring15 < 0.60) : return False
          if (eta > 1.479 and lep.mvaIdSpring15 < 0.17) : return False
          if lep.miniRelIso > 0.1: return False
        return True

if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("tree")
    tree.vectorTree = True
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            self.sf1 = edgeFriends("Edge", 
                lambda lep : _susyEdgeTight(lep),
                cleanJet = lambda lep,jet,dr : (jet.pt < 35 and dr < 0.4 and abs(jet.eta) > 2.4))
        def analyze(self,ev):
            print "\nrun %6d lumi %4d event %d: leps %d" % (ev.run, ev.lumi, ev.evt, ev.nLepGood)
            print self.sf1(ev)
            print self.sf2(ev)
            print self.sf3(ev)
            print self.sf4(ev)
            print self.sf5(ev)
    el = EventLoop([ Tester("tester") ])
    el.loop([tree], maxEvents = 50)
