MODULES = []

 
from CMGTools.TTHAnalysis.tools.edgeFriends import edgeFriends, _susyEdgeTight


MODULES.append( ('edgeFriends', edgeFriends("Edge",  
                                lambda lep : _susyEdgeTight(lep),
                                cleanJet = lambda lep,jet,dr : (jet.pt < 35 and dr < 0.4)) ) )
