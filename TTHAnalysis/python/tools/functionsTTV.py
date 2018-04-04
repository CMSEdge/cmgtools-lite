import itertools

def elMvaIDRun2(lep, name, wp):
     eta = abs(lep.etaSc)
     if name == "Fall17noIso":
         if wp == 'Loose':
             if lep.pt <= 10:
                 if   eta < 0.8  : return lep.mvaIdFall17noIso > -0.13285867293779202
                 elif eta < 1.479: return lep.mvaIdFall17noIso > -0.31765300958836074
                 else            : return lep.mvaIdFall17noIso > -0.0799205914718861
             else:
                 if   eta < 0.8  : return lep.mvaIdFall17noIso > -0.856871961305474
                 elif eta < 1.479: return lep.mvaIdFall17noIso > -0.8107642141584835
                 else            : return lep.mvaIdFall17noIso > -0.7179265933023059
         elif wp == 'wp90':
             if lep.pt<=10 and eta<0.8:
                 c = 0.9165112826974601
                 tau = 2.7381703555094217
                 A = 1.03549199648109
             elif lep.pt>10 and eta<0.8:
                 c = 0.9616542816132922
                 tau = 8.757943837889817
                 A = 3.1390200321591206
             elif lep.pt<=10 and eta<1.479:
                 c = 0.8655738322220173
                 tau = 2.4027944652597073
                 A = 0.7975615613282494
             elif lep.pt>10 and eta<1.479:
                 c = 0.9319258011430132
                 tau = 8.846057432565809
                 A = 3.5985063793347787
             elif lep.pt<=10:
                 c = -3016.035055227131
                 tau = -52140.61856333602
                 A = -3016.3029387236506
             elif lep.pt>10:
                 c = 0.8899260780999244
                 tau = 10.124234115859881
                 A = 4.352791250718547
             return lep.mvaIdFall17noIso > c-A*exp(-lep.pt/tau)
         elif wp == 'wp80':
             if lep.pt<=10 and eta<0.8:
                 c = 0.9530240956555949
                 tau = 2.7591425841003647
                 A = 0.4669644718545271
             elif lep.pt>10 and eta<0.8:
                 c = 0.9825268564943458
                 tau = 8.702601455860762
                 A = 1.1974861596609097
             elif lep.pt<=10 and eta<1.479:
                 c = 0.9336564763961019
                 tau = 2.709276284272272
                 A = 0.33512286599215946
             elif lep.pt>10 and eta<1.479:
                 c = 0.9727509457929913
                 tau = 8.179525631018565
                 A = 1.7111755094657688
             elif lep.pt<=10:
                 c = 0.9313133688365339
                 tau = 1.5821934800715558
                 A = 3.8889462619659265
             elif lep.pt>10:
                 c = 0.9562619539540145
                 tau = 8.109845366281608
                 A = 3.013927699126942
             return lep.mvaIdFall17noIso > c-A*exp(-lep.pt/tau)

     elif name == "Fall17Iso":
         if wp == 'Loose':
             if lep.pt <= 10:
                 if   eta < 0.8  : return lep.mvaIdFall17Iso > -0.09564086146419018
                 elif eta < 1.479: return lep.mvaIdFall17Iso > -0.28229916981926795
                 else            : return lep.mvaIdFall17Iso > -0.05466682296962322
             else:
                 if   eta < 0.8  : return lep.mvaIdFall17Iso > -0.833466688584422
                 elif eta < 1.479: return lep.mvaIdFall17Iso > -0.7677000247570116
                 else            : return lep.mvaIdFall17Iso > -0.6917305995653829
         elif wp == 'wp90':
             if lep.pt<=10 and eta<0.8:
                 c = 0.9387070396095831
                 tau = 2.6525585228167636
                 A = 0.8222647164151365
             elif lep.pt>10 and eta<0.8:
                 c = 0.9717674837607253
                 tau = 8.912850985100356
                 A = 1.9712414940437244
             elif lep.pt<=10 and eta<1.479:
                 c = 0.8948802925677235
                 tau = 2.7645670358783523
                 A = 0.4123381218697539
             elif lep.pt>10 and eta<1.479:
                 c = 0.9458745023265976
                 tau = 8.83104420392795
                 A = 2.40849932040698
             elif lep.pt<=10:
                 c = -1830.8583661119892
                 tau = -36578.11055382301
                 A = -1831.2083578116517
             elif lep.pt>10:
                 c = 0.8979112012086751
                 tau = 9.814082144168015
                 A = 4.171581694893849
             return lep.mvaIdFall17Iso > c-A*exp(-lep.pt/tau)
         elif wp == 'wp80':
             if lep.pt<=10 and eta<0.8:
                 c = 0.9725509559754997
                 tau = 2.976593261509491
                 A = 0.2653858736397496
             elif lep.pt>10 and eta<0.8:
                 c = 0.9896562087723659
                 tau = 10.342490511998674
                 A = 0.40204156417414094
             elif lep.pt<=10 and eta<1.479:
                 c = 0.9508038141601247
                 tau = 2.6633500558725713
                 A = 0.2355820499260076
             elif lep.pt>10 and eta<1.479:
                 c = 0.9819232656533827
                 tau = 9.05548836482051
                 A = 0.772674931169389
             elif lep.pt<=10:
                 c = 0.9365037167596238
                 tau = 1.5765442323949856
                 A = 3.067015289215309
             elif lep.pt>10:
                 c = 0.9625098201744635
                 tau = 8.42589315557279
                 A = 2.2916152615134173
             return lep.mvaIdFall17Iso > c-A*exp(-lep.pt/tau)
             

def _ttV_idEmu_cuts_E2(lep):
    if (abs(lep.pdgId)!=11): return True
    if (lep.hadronicOverEm>=(0.10-0.03*(abs(lep.etaSc)>1.479))): return False
    if (abs(lep.dEtaScTrkIn)>=(0.01-0.002*(abs(lep.etaSc)>1.479))): return False
    if (abs(lep.dPhiScTrkIn)>=(0.04+0.03*(abs(lep.etaSc)>1.479))): return False
    if (lep.eInvMinusPInv<=-0.05): return False
    if (lep.eInvMinusPInv>=(0.01-0.005*(abs(lep.etaSc)>1.479))): return False
    if (lep.sigmaIEtaIEta>=(0.011+0.019*(abs(lep.etaSc)>1.479))): return False
    return True
def _ttV_idEmu_cuts_E2_obj(lep):
    if (abs(lep.pdgId())!=11): return True
    etasc = lep.etaSc
    if (lep.hadronicOverEm()>=(0.10-0.03*(abs(etasc)>1.479))): return False
    if (abs(lep.deltaEtaSuperClusterTrackAtVtx())>=(0.01-0.002*(abs(etasc)>1.479))): return False
    if (abs(lep.deltaPhiSuperClusterTrackAtVtx())>=(0.04+0.03*(abs(etasc)>1.479))): return False
    eInvMinusPInv = (1.0/lep.ecalEnergy() - lep.eSuperClusterOverP()/lep.ecalEnergy()) if lep.ecalEnergy()>0. else 9e9
    if (eInvMinusPInv<=-0.05): return False
    if (eInvMinusPInv>=(0.01-0.005*(abs(etasc)>1.479))): return False
    if (lep.full5x5_sigmaIetaIeta()>=(0.011+0.019*(abs(etasc)>1.479))): return False
    return True

def _soft_MuonId_2016ICHEP(lep):
    if (abs(lep.pdgId())!=13): return False
    if not lep.muonID("TMOneStationTight"): return False #TMOneStationTightMuonId
    if not lep.track().hitPattern().trackerLayersWithMeasurement() > 5: return False
    if not lep.track().hitPattern().pixelLayersWithMeasurement() > 0: return False
    if not (abs(lep.dxy())<0.3 and abs(lep.dz())<20): return False
    return True

def _medium_MuonId_2016ICHEP(lep):
    if (abs(lep.pdgId())!=13): return False
    if not (lep.physObj.isGlobalMuon() or lep.physObj.isTrackerMuon()): return False
    if not (lep.innerTrack().validFraction()>0.49): return False
    if lep.segmentCompatibility()>0.451: return True
    else:
        if not lep.globalTrack().isNonnull(): return False
        if not lep.isGlobalMuon: return False
        if not lep.globalTrack().normalizedChi2()<3: return False
        if not lep.combinedQuality().chi2LocalPosition<12: return False
        if not lep.combinedQuality().trkKink<20: return False 
        if not lep.segmentCompatibility()>0.303: return False

    return True


def _ttV_looseTTVLepId(lep):
    if abs(lep.pdgId)==13:
        if abs(lep.eta) > 2.4: return False
        if lep.pt < 5: return False
        if abs(lep.dz)>=0.1: return False
        if abs(lep.dxy)>=0.05: return False
        if lep.sip3d>=8: return False
        if lep.miniRelIso>=0.4: return False
        if lep.looseIdOnly==0: return False
        return True

    elif abs(lep.pdgId)==11:
        etasc = lep.etaSc
        if abs(etasc) > 2.5: return False 
        if lep.pt < 7: return False
        if abs(lep.dz)>=0.1: return False
        if abs(lep.dxy)>=0.05: return False
        if lep.sip3d>=8: return False
        if lep.miniRelIso>=0.4: return False
        if lep.looseIdOnly==0: return False
        if not elMvaIDRun2(lep,"Fall17noIso","Loose"): return False
        #if abs(etasc)<0.8:
        #    if lep.mvaIdSpring16HZZ<0: return False
        #elif abs(etasc)<1.479:
        #    if lep.mvaIdSpring16HZZ<0: return False
        #else:
        #    if lep.mvaIdSpring16HZZ<0.7: return False
        if lep.lostHits>=2: return False
        return True
    else:
        return False


def _ttV_fakableTTVLepId(lep):
    if not _ttV_looseTTVLepId(lep): return False
    if abs(lep.pdgId)==13:
        if lep.jetBTagCSV > 0.8484: return False
        id1=(lep.mvaTTH>0.75)
        id2=(lep.jetPtRatiov2>0.5 and lep.jetBTagCSV<0.3 and lep.segmentCompatibility>0.3)
        if not (id1 or id2): return False
        return True
    elif abs(lep.pdgId)==11:
        if lep.jetBTagCSV > 0.8484: return False
        id1=(lep.mvaTTH>0.30)
        id2=(lep.jetPtRatiov2>0.5 and lep.jetBTagCSV<0.3)
        if not (id1 or id2): return False
        if not _ttV_idEmu_cuts_E2(lep): return False
        if lep.lostHits!=0: return False 
        return True
    else:
        return False


def _ttV_tightTTVLepId(lep):
    if not _ttV_fakableTTVLepId(lep): return False
    if abs(lep.pdgId)==13:
        if abs(lep.eta)> 2.4: return False
        if lep.mediumMuonID2016==0: return False
        if lep.mvaTTH<0.75: return False
        if lep.tightCharge<2: return False
        return True
    elif abs(lep.pdgId)==11:
        if abs(lep.eta)> 2.5: return False
        if lep.convVeto!=1: return False
        if lep.mvaTTH<0.30: return False
        return True
    else:
        return False


def _ttV_conePt(lep):
    if abs(lep.pdgId)!=11 and abs(lep.pdgId)!=13: return lep.pt
    if (abs(lep.pdgId)!=13 or lep.mediumMuonID2016>0) and \
            lep.mvaTTH>0.75: return lep.pt
    else: return 0.750*lep.pt/lep.jetPtRatiov2


def getLepZIndices(leps, lepIdxs):
    
    #build Zcandidate first
    if len(leps)<2: return [-1,-1]

    lepPair=[-1,-1]
    tmp=10000
    
    for x in itertools.combinations(lepIdxs,2):
        mz = (leps[x[0]].p4() + leps[x[1]].p4()).M()
        diff=abs(mz-91)
        if diff<tmp:
            tmp=abs(mz-91)
            lepPair[0]=lepIdxs[x[0]]
            lepPair[1]=lepIdxs[x[1]]
    return lepPair
