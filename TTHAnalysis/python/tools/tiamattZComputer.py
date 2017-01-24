from CMGTools.TTHAnalysis.treeReAnalyzer import *
from PhysicsTools.HeppyCore.utils.deltar import matchObjectCollection3
import ROOT

#sys.path.append("/afs/cern.ch/work/m/mmarionn/private/Keras/")
#from NNFunctions import evaluateModel, loadModel

sys.path.append("/afs/cern.ch/work/m/mmarionn/private/NeuralHEP/utils/")
from KerasSuperLayers import *
from KerasUtils import evaluateModel, loadModel, getLayerCollection

from CMGTools.TTHAnalysis.tools.functionsTTV import getLepZIndices

class MyVarProxy:
    def __init__(self,lep):
        self._ob = lep
    def __getitem__(self,name):
        return self.__getattr__(name)
    def __getattr__(self,name):
        if name in self.__dict__: return self.__dict__[name]
        else: return getattr(self._ob,name)
    def eta(self): return self._ob.eta
    def phi(self): return self._ob.phi
    def pt(self): return self._ob.pt
    def pdgId(self): return self._ob.pdgId


def buildZCandidate(leps):
    lep1=ROOT.TLorentzVector()
    lep2=ROOT.TLorentzVector()
    lep1.SetPtEtaPhiM(leps[0][0],leps[0][1],leps[0][2],0.0005 if abs(leps[0][3])==11 else 0.105)
    lep2.SetPtEtaPhiM(leps[1][0],leps[1][1],leps[1][2],0.0005 if abs(leps[1][3])==11 else 0.105)

    Z=lep1+lep2
    output=[Z.Pt(),Z.Eta(),Z.Phi(),Z.M(),0,0]  
    return output

def buildWCandidate(lep,met):
    lep4v=ROOT.TLorentzVector()
    met4v=ROOT.TLorentzVector()
    lep4v.SetPtEtaPhiM(lep[0],lep[1],lep[2],0.0005 if abs(lep[3])==11 else 0.105)
    met4v.SetPtEtaPhiM(met[0],0,met[1],0)

    WPt=(lep4v+met4v).Pt()
    WPhi=(lep4v+met4v).Phi()
    dPhi=lep4v.DeltaPhi(met4v)
    WMT=np.sqrt(2*lep[0]*met[0]*(1-np.cos(dPhi)))

    output=[WPt,0,WPhi,WMT,0,0]
    return output 


class TiamattZComputer:

    def __init__(self,label,training,preprocessing,layers=None, inputLabel="", indicesFromMPAF=True):
        self.label = "" if (label in ["",None]) else ("_"+label)
        self.training = training
        self.systsJEC = ["","_Up","_Do"]
        self.debugprinted = False
        self.layers=layers
        self.inputlabel=inputLabel
        self.indicesFromMPAF=indicesFromMPAF
        self.checkLayers=(layers!=None)
        #print " ====== begin "
        self.model=loadModel(weightFile=training,customLayers=superLayers)
        if self.checkLayers:
            self.models=getLayerCollection(self.model, layers)
        self.preprocessModel=loadModel(weightFile=preprocessing,customLayers=superLayers)
        #print " ====== end "
    def listBranches(self):
        label = self.label
        
        #biglist = [
        #    ("nNNScores"+label,"I"),
        #    ("NNRes"+label, "F",5,"nNNScores"+label)
        #    ]
        biglist=[]
        if not self.checkLayers:
            for key in self.systsJEC:
                biglist.extend([
                        ("nNNScores"+label+key,"I"),
                        ("NNRes"+label+key, "F",5,"nNNScores"+label+key)
                        ])
        
        if self.checkLayers:
            for layer in self.layers:
                biglist.extend([
                        ("nOutput"+layer+label,"I"),
                        ("outputLayer"+layer+label, "F",256,"nOutput"+layer+label)
                        ])


        return biglist

    def getTiamattZPrediction(self, ilepsZ, ilepW,ijets,imet,imetphi, varJEC=""):
        lepsZ=[]
        lepW=[]
        met=[]
        jets=[]
        nJMax=10
        #jets
        
        if len(ijets)<2 or len(ilepsZ)!=2 or ilepW==None:
            if not self.checkLayers:
                return [], 0
            else:
                res={}
                for layer in self.models.keys(): 
                    res[layer]=[]
                return res
        
        for k in range(0,min(len(ijets),nJMax)):
            j=ijets[k]
            pt=-100
            if varJEC=="":
                pt=j.pt
            if varJEC=="_Up":
                pt=j.pt*getattr(j,"corr_JECUp")/getattr(j,"CorrFactor_L1L2L3Res")
            if varJEC=="_Do":
                pt=j.pt*getattr(j,"corr_JECDown")/getattr(j,"CorrFactor_L1L2L3Res")
                
            jet=[pt,j.eta,j.phi,getattr(j,"mass"),getattr(j,"btagCSV"),getattr(j,"qgl") ]
            jets.append(jet)
            
        if len(jets)<nJMax:
            for j in range(len(jets),nJMax):
                jets.append([0,0,0,0,0,0]) 
            
        for lep in ilepsZ:
            lepsZ.append([lep.pt, lep.eta, lep.phi, 0.0005 if abs(lep.pdgId)==11 else 0.105, 0,0])
            
        lepW.extend([ilepW.pt, ilepW.eta, ilepW.phi, 0.0005 if abs(ilepW.pdgId)==11 else 0.105,0,0])
        if varJEC=="":
            met.extend([ imet[0], 0, imetphi[0], 0, 0, 0 ])
        if varJEC=="_Up":
            met.extend([ imet[1], 0, imetphi[1], 0, 0, 0 ])
        if varJEC=="_Do":
            met.extend([ imet[-1], 0, imetphi[-1], 0, 0, 0 ])
            
        ZCand=buildZCandidate(lepsZ)
        WCand=buildWCandidate(lepW,met)
        WCand2=buildWCandidate(lepsZ[0],met)
        WCand3=buildWCandidate(lepsZ[1],met)
            
        evt=[lepsZ[0], lepsZ[1], lepW, met]
        evt.extend(jets)
        evt.extend([ZCand, WCand, WCand2, WCand3 ])
        
        #print "---------------------------- ",self.checkLayers
        #print " ============== ", np.array(evt).shape
        if not self.checkLayers:
            preprocEvt=evaluateModel(np.array(evt),self.preprocessModel)
            res=evaluateModel(preprocEvt,self.model)
            NNresult=[]
        
            for r in res[0]:
                NNresult.append(r)
            return NNresult, len(NNresult)
        else:
            res={}
            for layer in self.models.keys():
                preprocEvt=evaluateModel(np.array(evt),self.preprocessModel)
                tmpRes=evaluateModel(preprocEvt, self.models[layer] )
                res[layer]=tmpRes[0].tolist()
                #print "==>> ", res
            #print "============================================ ",res
            return res


    def __call__(self,event):
        self.ev = event
        fullret = {}
        leps = [l for l in Collection(event,"LepGood","nLepGood")]

        lepsZ=[]
        lepW=None
        jets=[]
        if self.indicesFromMPAF:
            lepsZ=[]
            lepW=None
            jets=[]
            lepsZ=[leps[il] for il in getattr(event,"lepZIdx") if il!=-1]#[0:2]
            lepW =leps[getattr(event,"lepWIdx") ] if getattr(event,"lepWIdx") !=-1 else None
            allJets=[j for j in Collection(event,"Jet","nJet")]
            jets=[allJets[ij] for ij in getattr(event,"jetIdx") if ij!=-1]#[0:len(getattr(event,"jetIdx"))]
        else:
            lepsZ=[]
            lepW=None
            jets=[]
            lepsZIdxs=[]
            lepsIdx=[il for il in getattr(event,"iT"+self.inputlabel)]
            if getattr(event,"nLepTight"+self.inputlabel)>=3:
                lepst=[leps[il] for il in getattr(event,"iT"+self.inputlabel)][0:getattr(event,"nLepTight"+self.inputlabel)]
                lepsZIdxs=getLepZIndices(lepst,lepsIdx)
                lepsZ=[leps[il] for il in lepsZIdxs]
                lepW=None
                for idx in lepsIdx:
                    if idx not in lepsZIdxs:
                        lepW = leps[idx]
                        break
            
            allJets=[j for j in Collection(event,"Jet","nJet")]
            if getattr(event,"nJetSel"+self.inputlabel)>=2:
                jets=[allJets[ij] for ij in getattr(event,"iJSel"+self.inputlabel)][0:getattr(event,"nJetSel"+self.inputlabel)]
 
        met={}
        metphi={}
        met[0]=event.met_pt
        met[1]=getattr(event,"met_jecUp_pt",event.met_pt) if hasattr(event,"met_jecUp_pt") else event.met_pt #fixme!!!!
        met[-1]=getattr(event,"met_jecDown_pt",event.met_pt) if hasattr(event,"met_jecDown_pt") else event.met_pt
        metphi[0]= event.met_phi
        metphi[1]= getattr(event,"met_jecUp_phi",event.met_phi) if hasattr(event,"met_jecUp_phi") else event.met_phi
        metphi[-1]= getattr(event,"met_jecDown_phi",event.met_phi) if  hasattr(event,"met_jecDown_phi") else event.met_phi
            
        if not self.checkLayers:
            for var in self.systsJEC:
                fullret["NNRes"+self.label+var], fullret["nNNScores"+self.label+var]=self.getTiamattZPrediction(lepsZ, lepW, jets, met, metphi, varJEC=var)
        else:
            ret=self.getTiamattZPrediction(lepsZ, lepW, jets, met, metphi,varJEC="")
            #print ret
            for layer in ret.keys():
                fullret["nOutput"+layer+self.label]=len(ret[layer])
                fullret["outputLayer"+layer+self.label]=ret[layer]

        return fullret



if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("tree")
    tree.vectorTree = True
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
        def analyze(self,ev):
            print "\nrun %6d lumi %4d event %d: leps %d" % (ev.run, ev.lumi, ev.evt, ev.nLepGood)
            print self.sf1(ev)
            print self.sf2(ev)
            print self.sf3(ev)
            print self.sf4(ev)
            print self.sf5(ev)
    el = EventLoop([ Tester("tester") ])
    el.loop([tree], maxEvents = 50)

        
