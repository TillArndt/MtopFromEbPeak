#!/usr/bin/env python

import optparse
import os,sys
import json
import pickle
import ROOT
from subprocess import Popen, PIPE

"""
Perform the analysis on a single file
"""
def runBJetEnergyPeak(inFileURL, outFileURL, xsec=None):

    print '...analysing %s' % inFileURL

    #load btag effs
    cachefile = open('%s/src/UserCode/MtopFromEbPeak/analyzeNplot/data/btagefficiencies.pck' % os.environ['CMSSW_BASE'], 'r')
    btagEffs = pickle.load(cachefile)
    cachefile.close()


    ROOT.gSystem.Load('libCondFormatsBTauObjects') 
    calib = ROOT.BTagCalibration("CSVv2", "CSVv2.csv")
    reader = ROOT.BTagCalibrationReader(calib,0,"incl","central")  # 0 is for loose op , central is for the systematic
    #book some histograms
    histos={ 
        'nvtx'  :ROOT.TH1F('nvtx',';Vertex multiplicity; Events',30,0,30),
        'nbtags':ROOT.TH1F('nbtags',';b-tag multiplicity; Events',5,0,5),
        'bjeten':ROOT.TH1F('bjeten',';Energy [GeV]; Jets',30,0,300),
        'bjetenls':ROOT.TH1F('bjetenls',';log(E);  1/E dN_{b jets}/dlog(E)',20,3.,7.,)
        }
    for key in histos:
        histos[key].Sumw2()
        histos[key].SetDirectory(0)

    #open file and loop over events tree
    fIn=ROOT.TFile.Open(inFileURL)
    tree=fIn.Get('data')
    totalEntries=tree.GetEntriesFast()
    for i in xrange(0,totalEntries):

        tree.GetEntry(i)
        if i%100==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
        #require at least two jets
        nJets, nBtags = 0, 0
        taggedJetsP4=[]
        
        bwgt_mc = 1.
        bwgt_data= 1.

        #BTagEffs
        for ij in xrange(0,tree.nJet):
            if tree.nGenWeight>0 :
                b_sf = reader.eval(0, tree.Jet_eta[ij], tree.Jet_pt[ij])
                c_sf = reader.eval(1, tree.Jet_eta[ij], tree.Jet_pt[ij])
                l_sf = reader.eval(2, tree.Jet_eta[ij], tree.Jet_pt[ij])
                b_eff = btagEffs['b']['loose'].Eval(tree.Jet_pt[ij])
                c_eff = btagEffs['c']['loose'].Eval(tree.Jet_pt[ij])
                l_eff = btagEffs['udsg']['loose'].Eval(tree.Jet_pt[ij])

                if  tree.Jet_CombIVF[ij]>0.460:
                    if abs(tree.Jet_flavour[ij]) ==5:
                        bwgt_data *=b_sf
                        bwgt_data*= b_eff
                        bwgt_mc *= b_eff
                    elif abs(tree.Jet_flavour[ij]) ==4:
                        bwgt_data *= c_sf
                        bwgt_data*= c_eff
                        bwgt_mc *= c_eff
                    else:
                        bwgt_data *= l_sf
                        bwgt_data*= l_eff
                        bwgt_mc *= l_eff
                else :
                    if abs(tree.Jet_flavour[ij]) ==5:
                        bwgt_data *=(1- (b_sf*b_eff))
                        bwgt_mc *=(1-b_eff)
                    elif abs(tree.Jet_flavour[ij]) ==4:
                        bwgt_data *=(1-(c_sf*c_eff))
                        bwgt_mc *=(1-c_eff)
                    else:
                        bwgt_data *=(1-(l_sf*l_eff))
                        bwgt_mc *=(1- l_eff)

        

        for ij in xrange(0,tree.nJet):

            #get the kinematics and select the jet
            jp4=ROOT.TLorentzVector()
            jp4.SetPtEtaPhiM(tree.Jet_pt[ij],tree.Jet_eta[ij],tree.Jet_phi[ij],tree.Jet_mass[ij])
            if jp4.Pt()<30 or ROOT.TMath.Abs(jp4.Eta())>2.4 : continue

            #count selected jet
            nJets +=1

            #save P4 for b-tagged jet
            if tree.Jet_CombIVF[ij]>0.460:
                nBtags+=1
                taggedJetsP4.append(jp4)
        
        if nJets<2 : continue

        #generator level weight only for MC
        evWgt=1.0
        if xsec              : evWgt  = xsec*tree.LepSelEffWeights[0]*tree.PUWeights[0]
        if tree.nGenWeight>0 : evWgt *= tree.GenWeights[0]

        #ready to fill the histograms
        histos['nvtx'].Fill(tree.nPV,evWgt)
        histos['nbtags'].Fill(nBtags,evWgt)


        #use up to two leading b-tagged jets
        for ij in xrange(0,len(taggedJetsP4)):
            if ij>1 : break
            bevtWgt = bwgt_data/bwgt_mc
            bevtWgt*=evWgt
            histos['bjeten'].Fill(taggedJetsP4[ij].E(),bevtWgt)
            histos['bjetenls'].Fill(ROOT.TMath.Log(taggedJetsP4[ij].E()),bevtWgt/taggedJetsP4[ij].E())
        
    #all done with this file
    fIn.Close()

    #save histograms to file
    fOut=ROOT.TFile.Open(outFileURL,'RECREATE')
    for key in histos: histos[key].Write()
    fOut.Close()


"""
Wrapper to be used when run in parallel
"""
def runBJetEnergyPeakPacked(args):
    
    try:
        return runBJetEnergyPeak(inFileURL=args[0],
                                 outFileURL=args[1],
                                 xsec=args[2])
    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],args[0])
        print 50*'<'
        return False


"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',      default=None,        type='string')
    parser.add_option('-i', '--inDir',       dest='inDir',       help='input directory with files',   default=None,        type='string')
    parser.add_option('-o', '--outDir',      dest='outDir',      help='output directory',             default='analysis',  type='string')
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel',    default=0,           type='int')
    (opt, args) = parser.parse_args()

    #read list of samples
    jsonFile = open(opt.json,'r')
    samplesList=json.load(jsonFile,encoding='utf-8').items()
    jsonFile.close()

    #prepare output
    if len(opt.outDir)==0    : opt.outDir='./'
    os.system('mkdir -p %s' % opt.outDir)
     
    #create the analysis jobs
    taskList = []
    for sample, sampleInfo in samplesList: 
        inFileURL  = 'file://%s/%s.root' % (opt.inDir,sample)
        #if not os.path.isfile(inFileURL): continue
        xsec=sampleInfo[0] if sampleInfo[1]==0 else None        
        outFileURL = '%s/%s.root' % (opt.outDir,sample)
        taskList.append( (inFileURL,outFileURL,xsec) )

    #run the analysis jobs
    if opt.njobs == 0:
        for inFileURL, outFileURL, xsec in taskList:
            runBJetEnergyPeak(inFileURL=inFileURL, outFileURL=outFileURL, xsec=xsec)
    else:
        from multiprocessing import Pool
        pool = Pool(opt.njobs)
        pool.map(runBJetEnergyPeakPacked,taskList)

    #all done here
    print 'Analysis results are available in %s' % opt.outDir
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
