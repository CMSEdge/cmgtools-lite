import FWCore.ParameterSet.Config as cms

from CMGTools.H2TauTau.objects.eleCuts_cff import getEleCuts
from CMGTools.H2TauTau.objects.muCuts_cff import getMuCuts
# from CMGTools.H2TauTau.objects.eleSkimCuts_cff import * 
# from CMGTools.H2TauTau.objects.tauSkimCuts_cff import * 

muEleCuts = cms.PSet(
    baseline = cms.PSet(
         #mass = cms.string('mass()>10'),
         muLeg = getMuCuts('leg1','muEle').clone(),
         eleLeg = getEleCuts('leg2', 'muEle').clone()
         ),
    )

if __name__ == '__main__':

    print muEleCuts
