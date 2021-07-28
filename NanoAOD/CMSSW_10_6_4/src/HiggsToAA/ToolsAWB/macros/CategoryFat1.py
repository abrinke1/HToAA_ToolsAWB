#! /usr/bin/env python

#############################################################################################################
###  CategoryFat1.py : a script to study signal and background yields in the 3- or 4-merged jet category  ###
#############################################################################################################

import ROOT as R
R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

import os
import sys
import math
import copy
import numpy
import pickle
import subprocess

sys.path.insert(0, '%s/config/CategoryFat1' % os.getcwd())
# import NoAK4_SelV1
# CUTS = NoAK4_SelV1.Cuts()
import NoAK4_SelV1_JetHT_AB
CUTS = NoAK4_SelV1_JetHT_AB.Cuts()

sys.path.insert(0, '%s/python' % os.getcwd())
import EventWeights
EVT_WGTS = EventWeights.GetWeights()
import EventSelection
EVT_SEL = EventSelection.GetSelection()


## User-defined constants
MAX_EVT =     -1  ## Maximum number of events to process
PRT_EVT =  10000  ## Print to screen every Nth event
VERBOSE = False   ## Print extra info about each event

SIG_MODEL  = ''
BOSON_ID   = 36       ## PDG ID of scalar 's' or pseudo-scalar 'a' (can be 36, XXX, or 9000006)
BOSON_MASS = 20       ## "a" or "s" boson mass
# SIG_MODEL  = 'VBF_MS-40_ctauS-0'
# BOSON_ID   = 9000006  ## PDG ID of scalar 's' or pseudo-scalar 'a' (can be 36, XXX, or 9000006)
# BOSON_MASS = 40       ## "a" or "s" boson mass
MASS_B     = 4.18     ## Bottom quark mass in GeV
MASS_MU    = 0.10566  ## Muon mass in GeV

MAX_DR_JET   = 0.4  ## Maximum dR for matching between GEN and RECO AK4 jets
MAX_DR_FAT   = 0.8  ## Maximum dR for matching between GEN and RECO AK8 jets
MIN_PT_GEN   = 5    ## Minimum pT for GEN b-quarks
MAX_ETA_GEN  = 2.4  ## Maximum |eta| value for GEN b-quarks

SEL_MU  = False ## Select only events with a muon matched to one of the selected b-jets (and scale by expected luminosity)
SEL_4B  = False ## Select only events with all 4 GEN b-quarks passing the pT and eta cuts

XSEC = 43920.  ## ggH cross section in fb: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV2014#s_13_0_TeV
N_MC = 999000  ## ggH MC generated with 999k events

DO_NN  = False  ## Get NN output scores from data/NN_outputs
MIN_NN = -99    ## Minimum NN score for selected events



def main():

###################
## Initialize files
###################

    ## Location of input files
    in_file_names = []
    # in_dir = os.getcwd()+'/skims/'
    # in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/MC/VBFHToAATo4B/'
    # in_dir = '/cms/data/store/user/mcmaster/CrabTestGGHK1M_2/CrabTestGGHK1M_5/200303_211319/0000/'
    # in_dir = '/cms/data/store/user/mcmaster/ParkingBPH4/CRAB3_2018D_Parked_promptD-v1/200218_214714/0001/'

    # in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/data/ParkingBPH4/Nano25Oct2019/Run2018D-05May2019promptD-v1/SkimsFatJet/'
    # in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/data/JetHT/Nano25Oct2019/2018/SkimsFatJet/'

    # in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/MC/GGHToAATo4B/Mass20/CrabTestGGHK1M_2_CrabTestGGHK1M_5/200303_211319/0000/SkimsFatJet/'
    in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/MC/GGHToAATo4B/Mass20_HT100/CrabTestGGHK1M_HPT_5/200705_000519/0000/SkimsFatJet/'
    # in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/MC/VBFH_HToSSTo4b/VBFH_HToSSTo4b_MH-125_TuneCP5_13TeV-powheg-pythia8/Nano25Oct2019/SkimsFatJet/'
    
    # in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/MC/QCD/QCD_BGenFilter/Nano25Oct2019/SkimsFatJet/'
    # in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/MC/QCD/QCD_bEnriched/Nano25Oct2019/SkimsFatJet/'
    # in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/MC/QCD/QCD_Inclusive/Nano25Oct2019/SkimsFatJet/'
    # in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/MC/TTJets/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/Nano25Oct2019/SkimsFatJet/'
    # in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/MC/ZJets/ZJetsToQQ_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8/Nano25Oct2019/SkimsFatJet/'
    # in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/MC/WJets/WJetsToQQ_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/Nano25Oct2019/SkimsFatJet/'

    # in_dir = '/cms/data/store/user/abrinke1/NanoAOD/2018/MC/TTJets/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/Nano25Oct2019/SkimsTTBar/'

    for file_name in subprocess.check_output(['ls', in_dir]).splitlines():
        # if not 'VBHXLK4_198k.root' in file_name: continue
        if not '.root' in file_name: continue
        # if not 'nFat1_999k.root' in file_name: continue
        # if not '5_18p6M.root' in file_name: continue
        # if '214714_000' in file_name: continue
        # if not 'massH_70_msoft_70' in file_name: continue
        # if not 'nFat1_doubB_0p8' in file_name: continue
        if not 'msoft_90' in file_name: continue
        if 'Mu_pT_6' in file_name: continue
        # if not 'GenBVeto' in file_name: continue
        if 'Veto_HT' in file_name: continue
        # if not '2018_JetHT' in file_name: continue
        # if not '400toInf' in file_name: continue
        # if not 'oneAK8_bloose' in file_name: continue
        in_file_names.append(in_dir+file_name)
        print 'Opening file: '+in_file_names[-1]

    ## Chain together trees from input files
    ch = R.TChain('Events')
    for i in range(len(in_file_names)):
        ch.Add( in_file_names[i] )

    BASE_STR = 'CategoryFat1'
    BASE_STR += '_%d' % BOSON_MASS
    if   'Parking'        in in_file_names[-1]: BASE_STR += '_Data'
    elif 'JetHT'          in in_file_names[-1]: BASE_STR += '_Data'
    elif 'QCD_Inclusive'  in in_file_names[-1]: BASE_STR += '_QCD_Inclusive'
    elif 'QCD_BGenFilter' in in_file_names[-1]: BASE_STR += '_QCD_BGenFilter'
    elif 'QCD_bEnriched'  in in_file_names[-1]: BASE_STR += '_QCD_bEnriched'
    elif 'Mass20_HT100'   in in_file_names[-1]: BASE_STR += '_Sig_MC_HT100'
    elif 'GGH'            in in_file_names[-1]: BASE_STR += '_Sig_MC'
    elif 'VBF'            in in_file_names[-1]: BASE_STR += '_Sig_MC_%s' % SIG_MODEL
    elif 'TTJets'         in in_file_names[-1]: BASE_STR += '_TTJets'
    elif 'ZJets'          in in_file_names[-1]: BASE_STR += '_ZJets'
    elif 'WJets'          in in_file_names[-1]: BASE_STR += '_WJets'
    else:
        print '\n\n*** Sample for %s not recognized! Quitting. ***'
        sys.exit()

    ## For ParkingBPH dataset, triggers are only appied to MC (if at all)
    if CUTS.TRG_BPH and not 'Data' in BASE_STR: BASE_STR += '_TrgBPH'
    ## For JetHT dataset, apply triggers to both data and MC
    if CUTS.TRG_JET_HT == 'All': BASE_STR += '_TrgJetHT'
    if CUTS.TRG_JET_HT == 'ABC': BASE_STR += '_TrgJetHT_ABC'
    if CUTS.TRG_JET_HT == 'AB' : BASE_STR += '_TrgJetHT_AB'
    if CUTS.TRG_JET_HT == 'CX' : BASE_STR += '_TrgJetHT_CX'
    BASE_STR += '_%s' % CUTS.SEL_NAME
    if not SEL_MU: BASE_STR += '_noMu'

    if MAX_EVT > 0:
        BASE_STR += '_%dk' % int(MAX_EVT / 1000)

    if MIN_NN > 0:
        BASE_STR += '_NN_0p%d' % int(MIN_NN*10)

    ## Bools to help switches in the processing
    isSig      = ('Sig_MC'         in BASE_STR)
    isQCD_Incl = ('QCD_Inclusive'  in BASE_STR)
    isQCD_BGen = ('QCD_BGenFilter' in BASE_STR)
    isQCD_bEnr = ('QCD_bEnriched'  in BASE_STR)
    isTTBar    = ('TTJets'         in BASE_STR)
    isZJets    = ('ZJets'          in BASE_STR)
    isWJets    = ('WJets'          in BASE_STR)
    isData     = ('Data'           in BASE_STR)
    isQCD = (isQCD_Incl or isQCD_BGen or isQCD_bEnr)
    isMC  = (isSig or isQCD or isTTBar or isZJets or isWJets)

    print '\nisSig = %s, isQCD = %s, isData = %s, isMC = %s' % (isSig, isQCD, isData, isMC)

    if not (isMC or isData):
        print '\n\n***** NEITHER DATA NOR MC IN %s??? Bizzare! Quitting. *****\n\n' % BASE_STR

    ## Set output directories (create if they do not exist)
    if not os.path.exists('plots/png/%s/' % BASE_STR):
        os.makedirs('plots/png/%s/' % BASE_STR)
    if not os.path.exists('data/NN_inputs/'):
        os.makedirs('data/NN_inputs/')

    ## Create output ROOT file for histograms
    out_file = R.TFile('plots/%s.root' % BASE_STR, 'recreate')
    ## Create output directory for png images of plots
    png_dir  = 'plots/png/%s/' % BASE_STR


    ## Loop over neural networks
    NN_list = ['BMM_2020_09_29', 'BMM_2020_09_29_btag', 'BMM_2020_09_29_phys']  ## 'BMM_2020_06_15'
    ## Dictionaries for NN input and output values by event
    NN_input_dict  = {}
    NN_output_dict = {}
    for NN in NN_list:
        ## Create dictionary for NN input values by event
        NN_input_dict[NN] = {}

        if not DO_NN: continue

        ## Access dictionaries of pre-computed NN output values by event, if they exist
        print '\nLoading NN weights from data/NN_outputs/%s_%s.pkl' % (BASE_STR, NN)
        with open(('data/NN_outputs/%s_%s.pkl' % (BASE_STR.replace('_1k',''), NN)), 'rb') as NN_out_file:
            NN_output_dict[NN] = pickle.load(NN_out_file)
    ## End loop over NN

    ## Signal 10% quantile boundaries
    NN_quant = {}
    NN_quant['BMM_2020_06_15']      = [0.0, 0.21008, 0.32144, 0.41611, 0.49704, 0.56477, 0.61700, 0.67510, 0.72722, 0.75547, 1.0]
    NN_quant['BMM_2020_09_29']      = [0.0, 0.62242, 0.72669, 0.79669, 0.84072, 0.87201, 0.89924, 0.91971, 0.93749, 0.95762, 1.0]
    NN_quant['BMM_2020_09_29_btag'] = [0.0, 0.86342, 0.89955, 0.92326, 0.93792, 0.94842, 0.95667, 0.96306, 0.96952, 0.97517, 1.0]
    NN_quant['BMM_2020_09_29_phys'] = [0.0, 0.88571, 0.89499, 0.90197, 0.90867, 0.91455, 0.92090, 0.92750, 0.93561, 0.94921, 1.0]

    # ## Load ggH pT re-weighting file
    # if isSig and not 'VBF' in SIG_MODEL:
    #     ggH_pt_file = R.TFile('data/GGH_pt/GGH_pt_MadGraph_vs_Powheg_MINLO_weights.root', 'read')
    #     print '\nLoading Higgs pT re-weighting from %s' % ggH_pt_file.GetName()
    #     h_ggH_pt_wgt = ggH_pt_file.Get('pt_diff_log2_ratio')
    #     h_ggH_pt_wgt.SetDirectory(0)  ## Save to local memory for this job
    #     ggH_pt_file.Close()
    print '\n\n*** For now, just implementing functional form of ggH pT re-weighting. ***\n\n'

    ## 2018 DeepCSV and DeepFlavB cuts: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    DeepB = {'0': -999, 'L': 0.1241, 'M': 0.4184, 'T': 0.7527}
    FlavB = {'0': -999, 'L': 0.0494, 'M': 0.2770, 'T': 0.7264}


#############
## Histograms
#############

    ## Histogram bins: [# of bins, minimum x, maximum x]
    ptS_bins  = [ 30,    0,  30]  ## Bins for pT (small)
    ptL_bins  = [100,    0, 1000] ## Bins for pT (large)
    log_bins  = [ 45,    2,  11]  ## Bins for log2(pT)
    eta_bins  = [ 30, -3.0, 3.0]  ## Bins for eta distributions
    phi_bins  = [ 32, -3.2, 3.2]  ## Bins for phi distributions
    IP_bins   = [ 20,    0,  10]  ## Bins for log2(sig_dXY)
    dRs_bins  = [ 30,    0, 0.6]  ## Bins for dR separation (small)
    dRm_bins  = [ 24,    0, 1.2]  ## Bins for dR separation (medium)
    dRl_bins  = [ 60,    0,   6]  ## Bins for dR separation (large)
    dR_bins   = [ 40,    0, 4.0]  ## Bins for dR separation (default)
    dEta_bins = [ 26, -5.1, 5.1]  ## Bins for dEta separation
    dPhi_bins = [ 32,    0, 3.2]  ## Bins for dPhi separation
    mass_bins = [ 40,    0, 200]  ## Bins for invariant mass
    # mass_bins = [ 80,    0, 200]  ## Bins for invariant mass
    # massH_bins = [40,    0, 200]  ## Bins for invariant mass of Higgs
    massA_bins = [50,    0,  50]  ## Bins for invariant mass of "a" boson
    massV_bins = [30,    0,  15]  ## Bins for log2 mass of VBF di-jet system
    tag_bins  = [  5, -0.5, 4.5]  ## Bins for number of b-tags
    PV_bins   = [ 50, -1.5, 98.5] ## Bins for number of primary vertices
    SV_bins   = [ 10, -1.5, 8.5]  ## Bins for number of secondary vertices
    NN_bins   = [ 12, -0.1, 1.1]  ## Bins for Neural Network output score


    ## Book 1D histograms
    ## Important to use '1D' instead of '1F' when dealing with large numbers of entries, and weighted events (higher precision)

    hst = {}  ## Dictionary of histograms

    ## Event weights
    hst['wgt_evt']         = R.TH1D('h_wgt_evt',        'log_{10}(Event weight), unweighted', 1000, -5, 5)
    hst['wgt_evt_wgt']     = R.TH1D('h_wgt_evt_wgt',    'log_{10}(Event weight), weighted',   1000, -5, 5)
    hst['wgt_cand']        = R.TH1D('h_wgt_cand',       'log_{10}(Cand. weight), unweighted', 1000, -5, 5)
    hst['wgt_cand_wgt']    = R.TH1D('h_wgt_cand_wgt',   'log_{10}(Cand. weight), weighted',   1000, -5, 5)
    hst['wgt_lumi']        = R.TH1D('h_wgt_lumi',       'Effective luminosity, unweighted',    650, 0, 65)
    hst['wgt_lumi_wgt']    = R.TH1D('h_wgt_lumi_wgt',   'Effectuve luminosity, weighted',      650, 0, 65)
    hst['wgt_PU']          = R.TH1D('h_wgt_PU',         'Pileup wieght, unweighted',          1000, 0, 10)
    hst['wgt_PU_wgt']      = R.TH1D('h_wgt_PU_wgt',     'Pileup weight, weighted',            1000, 0, 10)
    hst['wgt_trg']         = R.TH1D('h_wgt_trg',        'Trigger efficiency SF, unweighted',   200, 0,  2)
    hst['wgt_trg_wgt']     = R.TH1D('h_wgt_trg_wgt',    'Trigger efficiency SF, weighted',     200, 0,  2)
    hst['wgt_ggH_pt']      = R.TH1D('h_wgt_ggH_pt',     'ggH p_{T} weight, unweighted',        150, 0, 1.5)
    hst['wgt_ggH_pt_wgt']  = R.TH1D('h_wgt_ggH_pt_wgt', 'ggH p_{T} weight, weighted',          150, 0, 1.5)
    hst['wgt_Xsec']        = R.TH1D('h_wgt_Xsec',       'log_{10}(Xsec weight), unweighted',   600, -4, 2)
    hst['wgt_Xsec_wgt']    = R.TH1D('h_wgt_Xsec_wgt',   'log_{10}(Xsec weight), weighted',     600, -4, 2)
    hst['wgt_HEM']         = R.TH1D('h_wgt_HEM',        'HEM 15/16 veto, unweighted',          110, -0.05, 1.05)
    hst['wgt_HEM_wgt']     = R.TH1D('h_wgt_HEM_wgt',    'HEM 15/16 veto, weighted',            110, -0.05, 1.05)


    ## *** Trigger efficiency (unique w.r.t. HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4, which collected 54.54/59.83 fb-1 of data in 2018) *** ##
    ## *** See brilcalc procedures and instructions at abrinke1@lxplus.cern.ch:~/BrilCalc/info.txt *** ##
    ## 64.8% ( 0.0%), 54.54 fb-1, 100% for pT > 350: JetHT  : [L1_SingleJet180                                                and HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4]
    ## 45.9% (46.2%), 54.54 fb-1, DROPS at high pT : JetHT  : [L1_DoubleJet112er2p3_dEta_Max1p6 or L1_DoubleJet150er2p5       and HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71]
    ## 12.5% (11.4%), 59.83 fb-1, DROPS at high pT : JetHT  : [L1_HTT360er                                                    and HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5]
    ## 30.7% (31.4%),  7 - 35 fb-1, flat in pT     : ParkBP : [L1_SingleMu6er1p5 or L1_SingleMu22                             and HLT_Mu7_IP4 or HLT_Mu8_IP3]
    ## 41.4% (34.5%), 16 - 20 fb-1, flat in pT     : BTagMu : [L1_Mu3_Jet120er2p5_dR_Max0p8 or L1_Mu3_Jet120er2p5_dR_Max0p4   and HLT_BTagMu_AK8DiJet170_Mu5_noalgo
    ##                                                                                                                         or HLT_BTagMu_AK4DiJet170_Mu5_noalgo]
    ## 27.9% ( 4.4%), 31.84 fb-1,  50% for pT > 350: BTagMu : [L1_SingleJet200                                                and HLT_BTagMu_AK8Jet300_Mu5_noalgo]
    ##  7.2% ( 8.1%), 54.55 fb-1, prefers lower pT : JetHT  : [L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6 and HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71]
    ##  6.0% ( 2.5%), 30.08 fb-1,  higher in pT    : JetHT  : [L1_Mu3_Jet120er2p5_dR_Max0p4                                   and HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71]
    ## 21.2% ( 1.4%), 59.83 fb-1, 100% for pT > 550: JetHT  : [L1_SingleJet180                                                and HLT_AK8PFJet500]
    ## https://cmsoms.cern.ch/cms/triggers/hlt_report?cms_run=324970
    ## https://cmsoms.cern.ch/cms/triggers/prescale?cms_run=324970&cms_run_sequence=GLOBAL-RUN
    hst['trig_eff'] = R.TH1D('h_trig_eff', 'h_trig_eff', 9, -0.5, 8.5)

    ## Number of Higgs candidates and jets in each event
    hst['nCandsSel']    = R.TH1D('h_nCandsSel',    'h_nCandsSel',    11, -0.5, 10.5)
    # hst['nCandsSV']     = R.TH1D('h_nCandsSV',     'h_nCandsSV',     11, -0.5, 10.5)
    # hst['nCandsHT']     = R.TH1D('h_nCandsHT',     'h_nCandsHT',     11, -0.5, 10.5)
    hst['nJet']         = R.TH1D('h_nJet',         'h_nJet',         21, -0.5, 20.5)

    ## Number of GEN b quarks matched to the Higgs candidate
    hst['nGenB_Higg'] = R.TH1D('h_nGenB_Higg', 'h_nGenB_Higg', 7, -0.5, 6.5)
    hst['nGenB_fatH'] = R.TH1D('h_nGenB_fatH', 'h_nGenB_fatH', 7, -0.5, 6.5)
    hst['nGenB_jetH'] = R.TH1D('h_nGenB_jetH', 'h_nGenB_jetH', 7, -0.5, 6.5)

    ## Muon and pileup properties
    hst['mu_pt_log'] = R.TH1D('h_mu_pt_log', 'h_mu_pt_log', log_bins[0], log_bins[1], log_bins[2])
    hst['mu_pt']     = R.TH1D('h_mu_pt',     'h_mu_pt',     ptS_bins[0], ptS_bins[1], ptS_bins[2])
    hst['mu_eta']    = R.TH1D('h_mu_eta',    'h_mu_eta',    eta_bins[0], eta_bins[1], eta_bins[2])
    hst['mu_IP_log'] = R.TH1D('h_mu_IP_log', 'h_mu_IP_log',  IP_bins[0],  IP_bins[1],  IP_bins[2])
    hst['mu_IP']     = R.TH1D('h_mu_IP',     'h_mu_IP',              20,           0,          20)

    hst['nMuon']         = R.TH1D('h_nMuon',         'h_nMuon',         9, -0.5, 8.5)
    hst['mu_charge']     = R.TH1D('h_mu_charge',     'h_mu_charge',     3, -1.5, 1.5)
    hst['mu_set_charge'] = R.TH1D('h_mu_set_charge', 'h_mu_set_charge', 7, -3.5, 3.5)

    hst['nPV_all']  = R.TH1D('h_nPV_all',  'h_nPV_all',  PV_bins[0], PV_bins[1], PV_bins[2])
    hst['nPV_good'] = R.TH1D('h_nPV_good', 'h_nPV_good', PV_bins[0], PV_bins[1], PV_bins[2])

    ## 'a' boson and Higgs pT, mass, and eta
    hst['pt_log_fatH'] = R.TH1D('h_pt_log_fatH', 'h_pt_log_fatH', log_bins[0], log_bins[1], log_bins[2])
    hst['pt_log_jetH'] = R.TH1D('h_pt_log_jetH', 'h_pt_log_jetH', log_bins[0], log_bins[1], log_bins[2])
    hst['pt_log_Higg'] = R.TH1D('h_pt_log_Higg', 'h_pt_log_Higg', log_bins[0], log_bins[1], log_bins[2])
    hst['pt_log_ISR']  = R.TH1D('h_pt_log_ISR',  'h_pt_log_ISR',  log_bins[0], log_bins[1], log_bins[2])
    
    hst['pt_fatH']   = R.TH1D('h_pt_fatH',   'h_pt_fatH',   ptL_bins[0], ptL_bins[1], ptL_bins[2])
    hst['pt_jetH']   = R.TH1D('h_pt_jetH',   'h_pt_jetH',   ptL_bins[0], ptL_bins[1], ptL_bins[2])
    hst['pt_jetInA'] = R.TH1D('h_pt_jetInA', 'h_pt_jetInA', ptL_bins[0], ptL_bins[1], ptL_bins[2])
    hst['pt_jetInB'] = R.TH1D('h_pt_jetInB', 'h_pt_jetInB', ptL_bins[0], ptL_bins[1], ptL_bins[2])
    hst['pt_Higg']   = R.TH1D('h_pt_Higg',   'h_pt_Higg',   ptL_bins[0], ptL_bins[1], ptL_bins[2])
    hst['pt_ISR']    = R.TH1D('h_pt_ISR',    'h_pt_ISR',    ptL_bins[0], ptL_bins[1], ptL_bins[2])
    
    hst['mass_fatH']  = R.TH1D('h_mass_fatH',  'h_mass_fatH',  mass_bins[0], mass_bins[1], mass_bins[2])
    hst['msoft_fatH'] = R.TH1D('h_msoft_fatH', 'h_msoft_fatH', mass_bins[0], mass_bins[1], mass_bins[2])
    hst['mass_fatA']  = R.TH1D('h_mass_fatA',  'h_mass_fatA',  mass_bins[0], mass_bins[1], mass_bins[2])
    hst['msoft_fatA'] = R.TH1D('h_msoft_fatA', 'h_msoft_fatA', mass_bins[0], mass_bins[1], mass_bins[2])
    hst['mass_jetH']  = R.TH1D('h_mass_jetH',  'h_mass_jetH',  mass_bins[0], mass_bins[1], mass_bins[2])
    hst['mass_Higg']  = R.TH1D('h_mass_Higg',  'h_mass_Higg',  mass_bins[0], mass_bins[1], mass_bins[2])
    hst['msoft_Higg'] = R.TH1D('h_msoft_Higg', 'h_msoft_Higg', mass_bins[0], mass_bins[1], mass_bins[2])
    hst['mass_VBF']   = R.TH1D('h_mass_VBF',   'h_mass_VBF',   massV_bins[0], massV_bins[1], massV_bins[2])
    
    hst['eta_fatH'] = R.TH1D('h_eta_fatH', 'h_eta_fatH', eta_bins[0], eta_bins[1], eta_bins[2])
    hst['eta_jetH'] = R.TH1D('h_eta_jetH', 'h_eta_jetH', eta_bins[0], eta_bins[1], eta_bins[2])
    hst['eta_Higg'] = R.TH1D('h_eta_Higg', 'h_eta_Higg', eta_bins[0], eta_bins[1], eta_bins[2])
    hst['eta_ISR']  = R.TH1D('h_eta_ISR',  'h_eta_ISR',  eta_bins[0], eta_bins[1], eta_bins[2])

    hst['phi_fatH_Pos'] = R.TH1D('h_phi_fatH_Pos', 'h_phi_fatH_Pos', phi_bins[0], phi_bins[1], phi_bins[2])
    hst['phi_ISR_Pos']  = R.TH1D('h_phi_ISR_Pos',  'h_phi_ISR_Pos',  phi_bins[0], phi_bins[1], phi_bins[2])
    hst['phi_fatH_Neg'] = R.TH1D('h_phi_fatH_Neg', 'h_phi_fatH_Neg', phi_bins[0], phi_bins[1], phi_bins[2])
    hst['phi_ISR_Neg']  = R.TH1D('h_phi_ISR_Neg',  'h_phi_ISR_Neg',  phi_bins[0], phi_bins[1], phi_bins[2])

    ## Separation betweeen AK4 and AK8 jets in Higgs, or ISR jet and Higgs
    hst['dR_Higg']   = R.TH1D('h_dR_Higg',   'h_dR_Higg',     dR_bins[0],   dR_bins[1],   dR_bins[2])
    hst['dEta_Higg'] = R.TH1D('h_dEta_Higg', 'h_dEta_Higg', dEta_bins[0], dEta_bins[1], dEta_bins[2])
    hst['dPhi_Higg'] = R.TH1D('h_dPhi_Higg', 'h_dPhi_Higg', dPhi_bins[0], dPhi_bins[1], dPhi_bins[2])
    hst['dR_ISR']    = R.TH1D('h_dR_ISR',    'h_dR_ISR',      dR_bins[0],   dR_bins[1],   dR_bins[2])
    hst['dEta_ISR']  = R.TH1D('h_dEta_ISR',  'h_dEta_ISR',  dEta_bins[0], dEta_bins[1], dEta_bins[2])
    hst['dPhi_ISR']  = R.TH1D('h_dPhi_ISR',  'h_dPhi_ISR',  dPhi_bins[0], dPhi_bins[1], dPhi_bins[2])

    ## Minimum dR from triggering muon to one of the b-quark or ISR jets
    hst['dR_mu_Higg'] = R.TH1D('h_dR_mu_Higg', 'h_dR_mu_Higg', dRm_bins[0], dRm_bins[1], dRm_bins[2])
    hst['dR_mu_ISR']  = R.TH1D('h_dR_mu_ISR',  'h_dR_mu_ISR',  dRm_bins[0], dRm_bins[1], dRm_bins[2])

    # ## Full event info (Higgs + up to 2 extra jets)
    # hst['evt_mass']   = R.TH1D('h_evt_mass', 'h_evt_mass', log_bins[0], log_bins[1], log_bins[2])
    # hst['evt_MT']     = R.TH1D('h_evt_MT',   'h_evt_MT',   log_bins[0], log_bins[1], log_bins[2])
    # hst['evt_HT']     = R.TH1D('h_evt_HT',   'h_evt_HT',   log_bins[0], log_bins[1], log_bins[2])
    hst['LHE_HT']     = R.TH1D('h_LHE_HT',     'LHE HT',        250,   0, 2500)
    hst['LHE_HT_log'] = R.TH1D('h_LHE_HT_log', 'log_{2} LHE HT', 80, 3.5, 11.5)

    ## Properties of signal candidate jets
    hst['QGL_jetH']     = R.TH1D('h_QGL_jetH',     'h_QGL_jetH',     40, -1.0, 1.0)
    hst['QGL_ISR']      = R.TH1D('h_QGL_ISR',      'h_QGL_ISR',      40, -1.0, 1.0)
    hst['jetID_jetH']   = R.TH1D('h_jetID_jetH',   'h_jetID_jetH',    7, -0.5, 6.5)
    hst['jetID_ISR']    = R.TH1D('h_jetID_ISR',    'h_jetID_ISR',     7, -0.5, 6.5)
    hst['jetPU_jetH']   = R.TH1D('h_jetPU_jetH',   'h_jetPU_jetH',    8, -0.5, 7.5)
    hst['jetPU_ISR']    = R.TH1D('h_jetPU_ISR',    'h_jetPU_ISR',     8, -0.5, 7.5)
    hst['deepB_fatH']   = R.TH1D('h_deepB_fatH',   'h_deepB_fatH',   22, -0.1, 1.0)
    hst['deepB_fatA']   = R.TH1D('h_deepB_fatA',   'h_deepB_fatA',   22, -0.1, 1.0)
    hst['deepB_jetH']   = R.TH1D('h_deepB_jetH',   'h_deepB_jetH',   22, -0.1, 1.0)
    hst['deepB_jetInA'] = R.TH1D('h_deepB_jetInA', 'h_deepB_jetInA', 22, -0.1, 1.0)
    hst['deepB_jetInB'] = R.TH1D('h_deepB_jetInB', 'h_deepB_jetInB', 22, -0.1, 1.0)
    hst['deepB_ISR']    = R.TH1D('h_deepB_ISR',    'h_deepB_ISR',    22, -0.1, 1.0)
    hst['flavB_jetH']   = R.TH1D('h_flavB_jetH',   'h_flavB_jetH',   20,  0.0, 1.0)
    hst['flavB_ISR']    = R.TH1D('h_flavB_ISR',    'h_flavB_ISR',    20,  0.0, 1.0)
    hst['doubB_fatH']   = R.TH1D('h_doubB_fatH',   'h_doubB_fatH',   60, -0.1, 1.0)
    hst['doubB_fatA']   = R.TH1D('h_doubB_fatA',   'h_doubB_fatA',   60, -0.1, 1.0)
    hst['bbVsLF_fatH']  = R.TH1D('h_bbVsLF_fatH',  'h_bbVsLF_fatH', 120, -0.1, 1.0)
    hst['ZvsQCD_fatH']  = R.TH1D('h_ZvsQCD_fatH',  'h_ZvsQCD_fatH',  60, -0.1, 1.0)
    hst['HvsQCD_fatH']  = R.TH1D('h_HvsQCD_fatH',  'h_HvsQCD_fatH',  60, -0.1, 1.0)
    hst['XvsQCD_fatH']  = R.TH1D('h_XvsQCD_fatH',  'h_XvsQCD_fatH',  60, -0.1, 1.0)
    hst['fourQ_fatH']   = R.TH1D('h_fourQ_fatH',   'h_fourQ_fatH',   22, -0.1, 1.0)
    hst['fourQ_fatA']   = R.TH1D('h_fourQ_fatA',   'h_fourQ_fatA',   22, -0.1, 1.0)
    hst['nConst_jetH']  = R.TH1D('h_nConst_jetH',  'h_nConst_jetH',  50,    0,  100)
    hst['nConst_ISR']   = R.TH1D('h_nConst_ISR',   'h_nConst_ISR',   50,    0,  100)
    hst['nSV_fatH']     = R.TH1D('h_nSV_fatH',     'h_nSV_fatH',     SV_bins[0], SV_bins[1], SV_bins[2])
    hst['nSV_fatA']     = R.TH1D('h_nSV_fatA',     'h_nSV_fatA',     SV_bins[0], SV_bins[1], SV_bins[2])
    hst['nSV_jetH']     = R.TH1D('h_nSV_jetH',     'h_nSV_jetH',     SV_bins[0], SV_bins[1], SV_bins[2])
    hst['nSV_Higg']     = R.TH1D('h_nSV_Higg',     'h_nSV_Higg',     SV_bins[0], SV_bins[1], SV_bins[2])
    hst['nSV_ISR']      = R.TH1D('h_nSV_ISR',      'h_nSV_ISR',      SV_bins[0], SV_bins[1], SV_bins[2])
    hst['nSV_extra']    = R.TH1D('h_nSV_extra',    'h_nSV_extra',    SV_bins[0], SV_bins[1], SV_bins[2])

    hst['tau1_fatH'] = R.TH1D('h_tau1_fatH',  'h_tau1_fatH', 28, -0.05, 0.65)
    hst['tau2_fatH'] = R.TH1D('h_tau2_fatH',  'h_tau2_fatH', 20, -0.05, 0.35)
    hst['tau3_fatH'] = R.TH1D('h_tau3_fatH',  'h_tau3_fatH', 30, -0.05, 0.25)
    hst['tau4_fatH'] = R.TH1D('h_tau4_fatH',  'h_tau4_fatH', 20, -0.05, 0.15)

    hst['n2b1_fatH'] = R.TH1D('h_n2b1_fatH', 'h_n2b1_fatH', 30, -0.1, 0.5)
    hst['n3b1_fatH'] = R.TH1D('h_n3b1_fatH', 'h_n3b1_fatH', 42, -0.1, 4.1)

    hst['sub1_pt_frac_fatH'] = R.TH1D('h_sub1_pt_frac_fatH', 'h_sub1_pt_frac_fatH', 20, 0, 1)
    hst['sub1_mass_fatH']    = R.TH1D('h_sub1_mass_fatH',    'h_sub1_mass_fatH',    50, 0, 100)
    hst['sub1_deepB_fatH']   = R.TH1D('h_sub1_deepB_fatH',   'h_sub1_deepB_fatH',   24, -0.1, 1.1)
    hst['sub1_tau1_fatH']    = R.TH1D('h_sub1_tau1_fatH',    'h_sub1_tau1_fatH',    28, -0.05, 0.65)
    hst['sub1_tau2_fatH']    = R.TH1D('h_sub1_tau2_fatH',    'h_sub1_tau2_fatH',    20, -0.05, 0.35)
    hst['sub1_tau3_fatH']    = R.TH1D('h_sub1_tau3_fatH',    'h_sub1_tau3_fatH',    30, -0.05, 0.25)
    hst['sub1_tau4_fatH']    = R.TH1D('h_sub1_tau4_fatH',    'h_sub1_tau4_fatH',    20, -0.05, 0.15)
    hst['sub1_n2b1_fatH']    = R.TH1D('h_sub1_n2b1_fatH',    'h_sub1_n2b1_fatH',    20,  0.15, 0.55)
    hst['sub1_n3b1_fatH']    = R.TH1D('h_sub1_n3b1_fatH',    'h_sub1_n3b1_fatH',    42, -0.1, 4.1)

    hst['sub2_pt_frac_fatH'] = R.TH1D('h_sub2_pt_frac_fatH', 'h_sub2_pt_frac_fatH', 20, 0, 1)
    hst['sub2_mass_fatH']    = R.TH1D('h_sub2_mass_fatH',    'h_sub2_mass_fatH',    50, 0, 100)
    hst['sub2_deepB_fatH']   = R.TH1D('h_sub2_deepB_fatH',   'h_sub2_deepB_fatH',   24, -0.1, 1.1)
    hst['sub2_tau1_fatH']    = R.TH1D('h_sub2_tau1_fatH',    'h_sub2_tau1_fatH',    28, -0.05, 0.65)
    hst['sub2_tau2_fatH']    = R.TH1D('h_sub2_tau2_fatH',    'h_sub2_tau2_fatH',    20, -0.05, 0.35)
    hst['sub2_tau3_fatH']    = R.TH1D('h_sub2_tau3_fatH',    'h_sub2_tau3_fatH',    30, -0.05, 0.25)
    hst['sub2_tau4_fatH']    = R.TH1D('h_sub2_tau4_fatH',    'h_sub2_tau4_fatH',    20, -0.05, 0.15)
    hst['sub2_n2b1_fatH']    = R.TH1D('h_sub2_n2b1_fatH',    'h_sub2_n2b1_fatH',    30,  0.15, 0.55)
    hst['sub2_n3b1_fatH']    = R.TH1D('h_sub2_n3b1_fatH',    'h_sub2_n3b1_fatH',    42, -0.1, 4.1)

    hst['sub12_pt_frac_fatH']   = R.TH1D('h_sub12_pt_frac_fatH',   'h_sub12_pt_frac_fatH',   40, 0.8, 1.2)
    hst['sub12_mass_frac_fatH'] = R.TH1D('h_sub12_mass_frac_fatH', 'h_sub12_mass_frac_fatH', 80, 0.4, 1.2)

    hst['sub12_dR']       = R.TH1D('h_sub12_dR',      'h_sub12_dR',      dRm_bins[0], dRm_bins[1], dRm_bins[2])
    hst['sub12_dR_fatH']  = R.TH1D('h_sub12_dR_fatH', 'h_sub12_dR_fatH', dRs_bins[0], dRs_bins[1], dRs_bins[2])

    ## Neural network output scores
    for NN in NN_list:
        hst[NN] = R.TH1D('h_NN_%s' % NN, 'h_NN_%s' % NN, NN_bins[0], NN_bins[1], NN_bins[2])

    ## Make sure proper uncertainties are stored and propagated
    for key in hst.keys():
        hst[key].Sumw2()

    ## Get "Model" branch(es)
    if isSig and SIG_MODEL != '':
        VBF_model = ch.GetBranch('GenModel_VBFH_HToSSTobbbb_MH-125_%s_TuneCP5_13TeV-powheg-pythia8' % SIG_MODEL.replace('VBF_',''))
        VBF_model.SetName('VBF_model')

    #############
    ## Event loop
    #############

    ## Keep track of all passing events, to make sure there are no duplicates
    event_list = []

    nEvtFatJet  = 0
    nEvtMuons   = 0
    nEvtTrig    = 0
    nEvtGen4B   = 0
    nEvtCatFat1 = 0
    nEvtSel     = 0

    nLostNN = 0

    ## Compute prescale to only access every Nth event
    PS = ( 1 if (MAX_EVT >= ch.GetEntries() or MAX_EVT < 1) else int(math.floor(ch.GetEntries() / MAX_EVT)) )
    # PS = 1

    print '\nAbout to run over %d events in chain, up to a maximum of %d (prescale = %d)\n' % (ch.GetEntries(), MAX_EVT, PS)

    iPrintSig = 0  ## Count passed signal events for printing
    for iEvt in range(ch.GetEntries() / PS):

        if iEvt >= MAX_EVT and MAX_EVT > 0: break
        if iEvt % PRT_EVT is 0: print 'Event #', iEvt

        ch.GetEntry(iEvt*PS)

        # if not (ch.run == 320010 and ch.luminosityBlock == 159 and ch.event == 246071540):
        #     continue
        # print 'Run = %d, LS = %d, Event = %d' % (ch.run, ch.luminosityBlock, ch.event)

        ####################################################
        ## Pre-select events based on simple RECO quantities
        ####################################################

        nMuon = ch.nMuon    ## Number of RECO muons in the event
        nJet  = ch.nJet     ## Number of RECO AK4 jets (dR = 0.4) in the event
        nFat  = ch.nFatJet  ## Number of RECO AK8 jets (dR = 0.8) in the event

        ## Quickly eliminate events not passing basic multiplicity cuts
        if nFat  < CUTS.NUM_FATS:  continue
        if nMuon < CUTS.NUM_MUONS: continue

        ## Remove data events not in Golden JSON
        if isData:
            pass_JSON = EVT_SEL.PassJSON(ch.run, ch.luminosityBlock)
            if not pass_JSON:
                # print 'Run %d, LS %d, event %d not in Golden JSON' % (ch.run, ch.luminosityBlock, ch.event)
                continue

        ## Remove events with no good collision vertices
        if ch.PV_npvsGood < 1:
            # print 'Run %d, LS %d, event %d has %d good primary vertices' % (ch.run, ch.luminosityBlock,
            #                                                                 ch.event, ch.PV_npvsGood)
            continue

        # print 'Here1'

        ## Only consider events containing at least 1 fat jet with pT > X and |eta| < Y
        nGoodFat = 0
        for iFat in range(nFat):
            if nGoodFat == CUTS.MAX_COMB_FAT: break
            if ch.FatJet_pt[iFat] <= CUTS.MIN_PT_FAT: break
            if ch.FatJet_pt[iFat] >  CUTS.MAX_PT_FAT: continue
            if abs(ch.FatJet_eta[iFat]) >= CUTS.MAX_ETA_FAT: continue
            nGoodFat += 1

        if nGoodFat < CUTS.NUM_FATS: continue
        nEvtFatJet += 1

        # print 'Here1a'

        ## Only consider events containing at least N muons passing softID, pT > 7 GeV, and |eta| < 2.4
        nGoodMu = 0
        for iMu in range(nMuon):
            if  ch.Muon_softId[iMu] != 1:   continue
            if      ch.Muon_pt[iMu]  < 7:   continue
            if abs(ch.Muon_eta[iMu]) > 2.4: continue
            nGoodMu += 1

        if nGoodMu < CUTS.NUM_MUONS: continue
        nEvtMuons += 1

        # print 'Here1b'

        ## Select only events that would pass HLT_Mu7_IP4 or HLT_Mu8_IP3
        if CUTS.TRG_BPH and isMC:
            if ch.HLT_Mu7_IP4_part0 == 0 and ch.HLT_Mu8_IP3_part0 == 0: continue
        ## Select events passing any of the relevant JetHT paths
        if CUTS.TRG_JET_HT == 'All' or CUTS.TRG_JET_HT == 'ABC':
            if not ( ( ch.L1_SingleJet180                      and ch.HLT_AK8PFJet500 ) or
                     ( ch.L1_SingleJet180                      and ch.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 ) or
                     ( ( ch.L1_DoubleJet150er2p5 or
                         ch.L1_DoubleJet112er2p3_dEta_Max1p6 ) and ch.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 ) ): continue
        ## Select only events passing the AK8 jet triggers
        if CUTS.TRG_JET_HT == 'AB':
            if not ( ( ch.L1_SingleJet180 and ch.HLT_AK8PFJet500 ) or
                     ( ch.L1_SingleJet180 and ch.HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4 ) ): continue
        ## Select only events passing the di-jet trigger
        if CUTS.TRG_JET_HT == 'CX':
            if not ( ( ch.L1_DoubleJet150er2p5 or
                       ch.L1_DoubleJet112er2p3_dEta_Max1p6 ) and ch.HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 ): continue

        nEvtTrig += 1

        # print 'Here2'

        ##############################################
        ## Retrieve GEN-level quantities for signal MC
        ##############################################

        nGen = (0 if isData else len(ch.GenPart_pdgId))  ## Number of GEN particles in the event

        if VERBOSE: print '\nIn event %d we find %d GEN particles, %d muons, %d jets, %d fat jets' % (iEvt, nGen, nMuon, nJet, nFat)

        genAIdx   = []  ## Indices of all GEN 'a' or 's' bosons
        genBIdx   = []  ## Indices of all GEN b-quarks not radiated from other b-quarks

        if isSig and SIG_MODEL != '':
            if not ch.VBF_model == 1: continue

        hasGenB_15 = False  ## Track whether event contains at least one GEN-level b-quark with pT > 15 GeV
        for iGen in range(nGen):
            genId = ch.GenPart_pdgId[iGen]  ## GEN particle PDG ID

            if abs(genId) != 5 and abs(genId) != BOSON_ID: continue  ## Only care about b-quarks and 'a'/'s' bosons
            if abs(genId) == BOSON_ID: genAIdx.append(iGen)          ## Store indices of 'a'/'s' bosons from Higgs decay

            if abs(genId) == 5 and ch.GenPart_pt[iGen] > 15: hasGenB_15 = True

            iGenMom  = ch.GenPart_genPartIdxMother[iGen]
            genMomId = -9999
            if iGenMom != -1:
                genMomId = ch.GenPart_pdgId[iGenMom]
            if VERBOSE: print '  * iGen = %d, genId = %d, iGenMom = %d, genMomId = %d' % (iGen, genId, iGenMom, genMomId)

            # if abs(genId) == 5 and abs(genMomId) == BOSON_ID:  ## Must come from 'a'/'s' boson decay
            #     genBIdx.append(iGen)
            ## Must not be radiated from another b-quark
            if abs(genId) == 5 and genId != genMomId and ch.GenPart_pt[iGen] > 0:
                ## If signal, must not be radiated from a light flavor quark, photon, or gluon
                if isSig and (abs(genMomId) <= 5 or genMomId == 21 or genMomId == 22): continue
                ## In signal, also exclude extra b-quarks with status = 33,41,43-44,51-52, or 71-73 and no mother
                if isSig and ch.GenPart_status[iGen] in [33,41,43,44,51,52,71,72,73] and iGenMom < 0: continue
                ## In signal, make a special note of extra b-quarks with status = 62 and no mother
                if isSig and ch.GenPart_status[iGen] == 62 and iGenMom < 0:
                    print '\n\n*** Interesting case in event %d!  b-quark with status = 62, no mother, eta = %.2f, pT = %.2f ***\n\n' % (ch.event, ch.GenPart_pt[iGen], ch.GenPart_eta[iGen])
                    continue
                ## Add this b-quark to the list
                genBIdx.append(iGen)

        if (isSig and (len(genAIdx) != 2 or len(genBIdx) != 4)) or (isTTBar and len(genBIdx) < 2):
            print '\n\nEvent %d bizzare error!!! %d genA, %d genB. Quitting!' % (ch.event, len(genAIdx), len(genBIdx))
            sys.exit()

        ## Select only events that have the chosen boson mass
        if isSig and ch.GenPart_mass[genAIdx[0]] != BOSON_MASS: continue

        ## Inclusive QCD sample must not contain any high-pT b-hadrons
        if isQCD_Incl and hasGenB_15: continue
        ## bEnriched and BGenFilter QCD samples must contain at least one high-pT b-hadron
        if (isQCD_bEnr or isQCD_BGen) and not hasGenB_15: continue

        ## Sort 'a' bosons by pT, b-quarks by pT after sorting by parent pT (for signal)
        iGenA = []
        if isSig:
            if ch.GenPart_pt[genAIdx[0]] >= ch.GenPart_pt[genAIdx[1]]: iGenA = [genAIdx[0], genAIdx[1]]
            if ch.GenPart_pt[genAIdx[1]] >  ch.GenPart_pt[genAIdx[0]]: iGenA = [genAIdx[1], genAIdx[0]]
        iGenB = []
        for iGen in genBIdx:  ## Double-loop over 2 b-quarks at the same time
            if not isSig:  ## For QCD background, just append all b-quarks
                iGenB.append(iGen)
                continue
            for jGen in genBIdx:
                if iGen == jGen: continue  ## Can't be 2 identical b-quarks
                if iGen in iGenB or jGen in iGenB: continue  ## Can't add in twice
                if ch.GenPart_genPartIdxMother[iGen] == iGenA[0] and \
                   ch.GenPart_genPartIdxMother[jGen] == iGenA[0]:  ## First consider b's from higher-pT 'a'
                    if ch.GenPart_pt[iGen] > ch.GenPart_pt[jGen]:  ## If first b is higher-pT, add to front of list
                        iGenB.insert(0, jGen)
                        iGenB.insert(0, iGen)
                    else:  ## If second b is higher-pT, add that to front of list
                        iGenB.insert(0, iGen)
                        iGenB.insert(0, jGen)
                elif ch.GenPart_genPartIdxMother[iGen] == iGenA[1] and \
                     ch.GenPart_genPartIdxMother[jGen] == iGenA[1]:  ## Next consider b's from lower-pT 'a'
                    if ch.GenPart_pt[iGen] < ch.GenPart_pt[jGen]:    ## If first b is lower-pT, add to back of list
                        iGenB.append(jGen)
                        iGenB.append(iGen)
                    else:  ## If second b is lower-pT, add that to back of list
                        iGenB.append(iGen)
                        iGenB.append(jGen)
            ## End loop: for jGen in genBIdx
        ## End loop: for iGen in genBIdx

        if (isSig and (len(iGenA) != 2 or len(iGenB) != 4)) or (isTTBar and len(iGenB) < 2):
            print '\n\nEvent %d bizzare error!!! len(iGenA) = %d, len(iGenB) = %d. Quitting!' % (ch.event, len(iGenA), len(iGenB))
            sys.exit()
        if len(iGenA) != len(genAIdx) or len(iGenB) != len(genBIdx):
            print '\n\nnEvent %d bizzare error!!! %d genAIdx, %d iGenA, %d genBIdx, %d iGenB. Quitting!' % (ch.event, len(genAIdx), len(iGenA), len(genBIdx), len(iGenB))
            sys.exit()

        ## Fill TLorentzVectors of b-quarks from 'a' decay
        bVec = []
        if isMC:
            for i in range(len(iGenB)):
                bVec.append(R.TLorentzVector())
                bVec[-1].SetPtEtaPhiM(ch.GenPart_pt[iGenB[i]], ch.GenPart_eta[iGenB[i]], ch.GenPart_phi[iGenB[i]], MASS_B)
            ## Only keep events with all 4 b-quarks passing pT and eta cuts
            if SEL_4B and min( [bVec[0].Pt(), bVec[1].Pt(), bVec[2].Pt(), bVec[3].Pt()] ) < MIN_PT_GEN: continue
            if SEL_4B and max( [abs(bVec[0].Eta()), abs(bVec[1].Eta()), abs(bVec[2].Eta()), abs(bVec[3].Eta())] ) > MAX_ETA_GEN: continue
            nEvtGen4B += 1

        ## Define a list of AK4 jets passing basic quality cuts
        goodJets,candJets = ([],[])
        for iJet in range(nJet):
            if     ch.Jet_pt[iJet]   <= CUTS.MIN_PT_JET:   break
            if abs(ch.Jet_eta[iJet]) >= CUTS.MAX_ETA_JET:  continue
            if ch.Jet_puId[iJet]     <  CUTS.MIN_PUID_JET: continue
            goodJets.append(iJet)
            ## Define a smaller list of potential Higgs candidate AK4 jets, either:
            ##   * candJetsOut: AK4 jets outside the AK8 jet (to be combined with AK8 to form total candidate)
            ##   * candJetsIn:  AK4 jets inside the AK8 jet  (used for di-jet triggering in JetHT dataset)
            if abs(ch.Jet_eta[iJet])      >= 2.4:                          continue 
            if ch.Jet_btagDeepB[iJet]     <= DeepB[CUTS.MIN_TAG_DEEP_JET]: continue
            if ch.Jet_btagDeepFlavB[iJet] <= FlavB[CUTS.MIN_TAG_FLAV_JET]: continue
            candJets.append(iJet)

        # print 'Here3'

        ## Define a list of AK8 jets passing basic quality cuts
        goodFats,candFats = ([],[])
        for iFat in range(nFat):
            if len(goodFats)            == CUTS.MAX_COMB_FAT: break
            if     ch.FatJet_pt[iFat]   <= CUTS.MIN_PT_FAT:   break
            if abs(ch.FatJet_eta[iFat]) >= CUTS.MAX_ETA_FAT: continue
            goodFats.append(iFat)
            ## Test for NaN values
            if math.isnan(ch.FatJet_mass     [iFat]) or \
               math.isnan(ch.FatJet_msoftdrop[iFat]) or \
               math.isnan(ch.FatJet_btagDDBvL[iFat]) or \
               math.isnan(ch.FatJet_btagDeepB[iFat]): continue
            ## Define a smaller list of potential Higgs candidate AK8 jets
            if ch.FatJet_mass[iFat]      <= CUTS.RECO_A_MASS[0]  or \
               ch.FatJet_mass[iFat]      >= CUTS.RECO_H_MASS[1]:          continue
            if ch.FatJet_msoftdrop[iFat] <= CUTS.RECO_A_MSOFT[0] or \
               ch.FatJet_msoftdrop[iFat] >= CUTS.RECO_H_MSOFT[1]:         continue
            if ch.FatJet_msoftdrop[iFat] >= CUTS.RECO_A_MSOFT_VETO[0] and \
               ch.FatJet_msoftdrop[iFat] <= CUTS.RECO_A_MSOFT_VETO[1]:    continue
            if ch.FatJet_btagDDBvL[iFat] <= CUTS.MIN_TAG_DOUB_FAT:        continue
            if ch.FatJet_btagDeepB[iFat] <= DeepB[CUTS.MIN_TAG_DEEP_FAT]: continue
            if ch.FatJet_deepTagMD_H4qvsQCD[iFat] <= CUTS.MIN_TAG_HTO4Q:  continue
            ## Cuts if the fat jet is the full reconstructed Higgs (no AK4 jet)
            if CUTS.NUM_FATS == 1 and CUTS.NUM_JETS_OUT == 0:
                if ch.FatJet_mass[iFat]      <= CUTS.RECO_H_MASS[0] or \
                   ch.FatJet_msoftdrop[iFat] <= CUTS.RECO_H_MSOFT[0]: continue
            candFats.append(iFat)

        # print 'Here4'

        ## Find index of AK8 jet matching 3 or 4 GEN jets
        ## If only 3 match, find index of AK4 jet matching the remaining GEN jet
        fatGenIdx,jetGenIdx = (-99,-99)  ## Index of best-matched AK8 jet [0, nFatJet) and AK4 jet [0 - nJet)
        fatGenVec = R.TLorentzVector()  ## 4-vector of best-matched AK8 jet
        fatGenVecS = R.TLorentzVector()  ## 4-vector of best-matched AK8 jet with SoftDrop mass
        jetGenVec = R.TLorentzVector()  ## 4-vector of best-matched AK4 jet
        fatGenMap = []  ## Map of GEN jets matched to AK8 jet (for signal: b11,b12,b21,b22)
        jetGenMap = []  ## Map of GEN jets matched to AK4 jet (for signal: b11,b12,b21,b22)

        if isMC:
            ## Loop over selected RECO AK8 jets
            for iFat in goodFats:
                iFatVec = R.TLorentzVector()
                iFatVec.SetPtEtaPhiM(ch.FatJet_pt[iFat], ch.FatJet_eta[iFat], ch.FatJet_phi[iFat], ch.FatJet_mass[iFat])
                iFatGenMap = []
                ## GEN b considered "matched" to AK8 jet if dR(jet, b) < 0.8
                for iGen in range(len(bVec)):
                    iFatGenMap.append(     bVec[iGen].DeltaR(iFatVec) < MAX_DR_FAT and \
                                           bVec[iGen].Pt()            > MIN_PT_GEN and \
                                       abs(bVec[iGen].Eta())          < MAX_ETA_GEN )
                ## If AK8 jet matches largest number of b-quarks, considered "the" jet
                if iFatGenMap.count(1) > fatGenMap.count(1):
                    fatGenIdx = iFat
                    fatGenVec = iFatVec
                    fatGenVecS.SetPtEtaPhiM(iFatVec.Pt(), iFatVec.Eta(), iFatVec.Phi(), ch.FatJet_msoftdrop[iFat])
                    fatGenMap = iFatGenMap
                    jetGenIdx = -99
                    jetGenVec = R.TLorentzVector()
                    jetGenMap = []

                    ## Find AK4 jet matched to b-quark (if any) not matched to AK8
                    for iGen in range(len(bVec)):
                        ## Assume this GEN b is not matched to any AK4 jet
                        jetGenMap.append(0)
                        ## Don't consider GEN b quarks already matched to AK8 jet
                        if fatGenMap[iGen] == 1: continue
                        min_dR_jet_b = 99.
                        ## Loop over RECO AK4 jets
                        for iJet in goodJets:
                            iJetVec = R.TLorentzVector()
                            iJetVec.SetPtEtaPhiM(ch.Jet_pt[iJet], ch.Jet_eta[iJet], ch.Jet_phi[iJet], ch.Jet_mass[iJet])
                            ## AK4 jet must be far away from AK8 jet, and closest to b quark
                            if iJetVec.DeltaR(fatGenVec)  < MAX_DR_FAT:   continue
                            if iJetVec.DeltaR(bVec[iGen]) > MAX_DR_JET:   continue
                            if iJetVec.DeltaR(bVec[iGen]) > min_dR_jet_b: continue
                            min_dR_jet_b = iJetVec.DeltaR(bVec[iGen])
                            jetGenIdx = iJet     ## FIXME!!! Doesn't handle signal case with > 1 AK4 jet.  AWB 2020.10.19
                            jetGenVec = iJetVec  ## FIXME!!! Doesn't handle signal case with > 1 AK4 jet.  AWB 2020.10.19
                            jetGenMap[iGen] = 1

                        ## End loop over RECO AK4 jets
                    ## End loop over GEN b-quarks
                ## End conditional for AK8 jet matching the most b-quarks
            ## End loop over AK8 jets
        ## End conditional for signal MC

        ## Only use MC signal events if all 4 b-quarks are matched to 1 AK8 RECO jet, or if
        ## 3 b-quarks are matched to 1 AK8 jet and the other b-quark is matched to an AK4 jet
        if isSig and SEL_4B:
            if fatGenMap.count(1) + jetGenMap.count(1) < 4: continue
            if fatGenMap.count(1) + jetGenMap.count(1) > 4:
                print '\n\nTOTALLY BIZZARE CASE IN EVENT %d!!!  fatGenMap / jetGenMap =' % ch.event
                print fatGenMap
                print jetGenMap
                sys.exit()

        nEvtCatFat1 += 1

        if fatGenMap.count(1) >= 6:
            print '\n\nEvent %d has %d matched b-quarks!!! Here is the list:' % (ch.event, fatGenMap.count(1))
            for iGen in range(len(fatGenMap)):
                if fatGenMap[iGen] != 1: continue
                bIdx = iGenB[iGen]
                genId    = ch.GenPart_pdgId[bIdx]
                bIdxMom  = ch.GenPart_genPartIdxMother[bIdx]
                genMomId = (-9999 if bIdxMom < 0 else ch.GenPart_pdgId[bIdxMom])
                print '*GenPart[%2d] Status = %2d, ID = %2d, Mom[%2d] = %5d, pT = %.3f, eta = %.3f, phi = %.3f' % (bIdx, ch.GenPart_status[bIdx], genId, bIdxMom, genMomId, ch.GenPart_pt[bIdx], ch.GenPart_eta[bIdx], ch.GenPart_phi[bIdx])
            print '\n'
            # sys.exit()




        #####################################################################################
        ###  Check whether different combinations of jets pass the various category cuts  ###
        #####################################################################################

        ## Number of candidate reconstructions
        nCandsSel,nCandsTwoAK4 = (0,0)

        ## Set event weights
        Evt_wgt    = 1.0
        PS_wgt     = ( 1.0 if (PS == 1) else (1.0*ch.GetEntries() / MAX_EVT) )
        ggH_pt_wgt = 1.0
        Xsec_wgt   = 1.0

        ## Apply prescale weight
        Evt_wgt *= PS_wgt

        ## Apply Higgs pT weight for ggH signal
        if isSig and not 'VBF' in SIG_MODEL:
            ## Find "final state" GEN Higgs pT
            for iGen in range(nGen):
                if ch.GenPart_pdgId [iGen] != 25: continue
                if ch.GenPart_status[iGen] != 62: continue
                log2_pt = math.log(ch.GenPart_pt[iGen], 2)
                # log2_pt_bin = h_ggH_pt_wgt.FindBin(log2_pt)
                # log2_pt_bin = max(log2_pt_bin, 1)
                # log2_pt_bin = min(log2_pt_bin, h_ggH_pt_wgt.GetNbinsX())
                # ggH_pt_wgt = h_ggH_pt_wgt.GetBinContent( log2_pt_bin )
                ggH_pt_wgt = max(0.1, 3.9 - 0.4*log2_pt)  ## Approximate functional form for Higgs pT weight
                # if VERBOSE: print '\nFor Higgs pT %.1f (log2 = %.2f, bin %d), weight = %.4f' % (ch.GenPart_pt[iGen], log2_pt, log2_pt_bin, ggH_pt_wgt)
                break
        Evt_wgt *= ggH_pt_wgt


        ## Apply cross-section weight
        if isSig:   ## From theoretial cross section, 100% H --> aa --> 4b branching
            Xsec_wgt = XSEC / N_MC
            if 'HT100' in BASE_STR:  ## Higher-stats HT100 sample
                Xsec_wgt *= 0.10653
        elif isQCD:
            Xsec_wgt = EVT_WGTS.GetWgtLHE(ch.LHE_HT, isQCD_Incl, isQCD_BGen, isQCD_bEnr)
        elif isTTBar:
            Xsec_wgt = 831760.0 / 10244307 ## For TTJets_TuneCP5_13TeV-madgraphMLM-pythia8 (https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO)
        elif isZJets:
            if   ch.LHE_HT >= 400 and ch.LHE_HT < 600: Xsec_wgt = 145400. / 16704355
            elif ch.LHE_HT >= 600 and ch.LHE_HT < 800: Xsec_wgt =  34000. / 14642701
            elif ch.LHE_HT >= 800:                     Xsec_wgt =  18670. / 10561192
        elif isWJets:
            if   ch.LHE_HT >= 400 and ch.LHE_HT < 600: Xsec_wgt = 315600. / 10071273
            elif ch.LHE_HT >= 600 and ch.LHE_HT < 800: Xsec_wgt =  68570. / 15298056
            elif ch.LHE_HT >= 800:                     Xsec_wgt =  34900. / 14627242
        
        Evt_wgt *= Xsec_wgt

        ## Save single per-event value of event weight at this stage
        Evt_wgt_orig = Evt_wgt

        # print 'Here5'

        ## Loop over all Higgs candidate fat jets
        for fatIdx in candFats:
            ## Don't process if we already have enough candidates
            if nCandsSel >= CUTS.MAX_NUM_CANDS: break

            ## Set weights which depend on the selected muon
            WGT      = Evt_wgt_orig
            Lumi_wgt = 1.0
            PU_wgt   = 1.0
            Trg_wgt  = 1.0
            HEM_wgt  = 1.0

            fatVec,fatVecS = (R.TLorentzVector(),R.TLorentzVector())
            fatVec .SetPtEtaPhiM(ch.FatJet_pt[fatIdx], ch.FatJet_eta[fatIdx], ch.FatJet_phi[fatIdx], ch.FatJet_mass[fatIdx])
            fatVecS.SetPtEtaPhiM(ch.FatJet_pt[fatIdx], ch.FatJet_eta[fatIdx], ch.FatJet_phi[fatIdx], ch.FatJet_msoftdrop[fatIdx])

            ## For signal events, AK8 jet must be GEN-matched to at least 3 b-quarks
            if isSig and SEL_4B and fatIdx != fatGenIdx: continue

            ## Indices and 4-vectors of AK4 jets outside or inside the AK8 jet radius
            jetIdxOut,jetIdxInA,jetIdxInB = (-99,-99,-99)
            jetVecOut,jetVecInA,jetVecInB = (R.TLorentzVector(),R.TLorentzVector(),R.TLorentzVector())
            ## 4-vectors of complete Higgs candidate
            higVec,higVecS = (R.TLorentzVector(),R.TLorentzVector())

            ## For now, only a few valid settings
            if CUTS.NUM_JETS_OUT != 0 and CUTS.NUM_JETS_OUT != 1:
                print '\n\n*** CUTS.NUM_JETS_OUT = %d!!! Quitting. ***\n\n' % CUTS.NUM_JETS_OUT
                sys.exit()
            if CUTS.NUM_JETS_IN != 'le1' and CUTS.NUM_JETS_IN != 'ge2' and CUTS.NUM_JETS_IN != 'any':
                print '\n\n*** CUTS.NUM_JETS_IN = %s!!! Quitting. ***\n\n' % CUTS.NUM_JETS_IN
                sys.exit()

            ## Case where no extra AK4 jet is used in Higgs candidate
            if CUTS.NUM_JETS_OUT == 0:
                ## For signal events, AK8 jet must be GEN-matched to all 4 b-quarks
                if isSig and SEL_4B and jetGenIdx >= 0: continue
                ## Higgs candidate is simply AK8 jet
                higVec,higVecS = (fatVec,fatVecS)
            ## End conditional: if CUTS.NUM_JETS_OUT == 0

            # print 'Here6'

            ## Case where AK4 jet is used to account for one b-quark from Higgs
            if CUTS.NUM_JETS_OUT == 1:
                ## AK8 jet must not be too massive
                if fatVec .M() >= CUTS.RECO_A_MASS [1]: continue
                if fatVecS.M() >= CUTS.RECO_A_MSOFT[1]: continue
                ## Find the minimum separation between any qualifying AK4 jet and the candidate AK8 jet
                min_dR_fat_jetOut = CUTS.MAX_H_DR

                ## All possible AK8+AK4 di-jet Higgs candidates
                for jetIdx in candJets:
                    ## For signal events, chosen AK4 jet must match one b-quark from Higgs
                    if isSig and SEL_4B and jetIdx != jetGenIdx: continue
                    ## Define AK4 jet 4-vectors and AK8 + AK4 Higgs candidate 4-vectors
                    iJetVec = R.TLorentzVector()
                    iJetVec.SetPtEtaPhiM(ch.Jet_pt[jetIdx], ch.Jet_eta[jetIdx], ch.Jet_phi[jetIdx], ch.Jet_mass[jetIdx])
                    ## AK4 and AK8 jet must not overlap
                    if iJetVec.DeltaR(fatVec) < MAX_DR_FAT: continue
                    ## Check if this jet has the minimum separation, and update value
                    if iJetVec.DeltaR(fatVec) >= min_dR_fat_jetOut: continue
                    min_dR_fat_jetOut = iJetVec.DeltaR(fatVec)
                    ## Assign this jet as part of the Higgs candidate
                    jetIdxOut = jetIdx
                    jetVecOut = iJetVec
                ## End loop: for jetIdx in candJets

                ## Require a matched AK4 jet
                if jetIdxOut < 0: continue
                ## Construct the AK8+AK4 Higgs candidate
                higVec,higVecS = (fatVec+jetVecOut,fatVecS+jetVecOut)
            ## End conditional: if CUTS.NUM_JETS_OUT == 1

            # print 'Here7'

            ## Higgs candidate must fall in mass range
            if higVec.M()  <= CUTS.RECO_H_MASS[0]        or \
               higVec.M()  >= CUTS.RECO_H_MASS[1]:       continue
            if higVecS.M() <= CUTS.RECO_H_MSOFT[0]       or \
               higVecS.M() >= CUTS.RECO_H_MSOFT[1]:      continue
            if higVecS.M() >= CUTS.RECO_H_MSOFT_VETO[0]  and \
               higVecS.M() <= CUTS.RECO_H_MSOFT_VETO[1]: continue

            # print 'Here8'

            ## Find highst-pT extra jets for ISR and/or VBF system
            isrVec,vbfVec = (R.TLorentzVector(), R.TLorentzVector())
            isrIdx,vbfIdx = (-99,-99)
            for x1 in goodJets:
                if x1 == jetIdxOut: continue
                x1Vec = R.TLorentzVector()
                x1Vec.SetPtEtaPhiM(ch.Jet_pt[x1], ch.Jet_eta[x1], ch.Jet_phi[x1], ch.Jet_mass[x1])
                if x1Vec.DeltaR(fatVec) < MAX_DR_FAT: continue
                ## Highest-pT extra jet assumed to be ISR
                isrIdx = x1
                isrVec = x1Vec
                ## Loop over secondary jets
                for x2 in goodJets:
                    if x2 <= x1 or x2 == jetIdxOut: continue
                    x2Vec = R.TLorentzVector()
                    x2Vec.SetPtEtaPhiM(ch.Jet_pt[x2], ch.Jet_eta[x2], ch.Jet_phi[x2], ch.Jet_mass[x2])
                    if x2Vec.DeltaR(fatVec) < MAX_DR_FAT: continue
                    ## Second-highest-pT extra jet assumed to be part of VBF system
                    vbfIdx = x2
                    vbfVec = x2Vec
                    break
                ## End loop: for x2 in range(x1+1, nJet)
                break
            ## End loop: for x1 in range(nJet)

            # print 'Here9'

            ## Count the number of matched secondary vertices
            nSV_fat,nSV_jet,nSV_ISR = (0,0,0)
            for iSV in range(len(ch.SV_eta)):
                vSV = R.TLorentzVector()
                vSV.SetPtEtaPhiM(ch.SV_pt[iSV], ch.SV_eta[iSV], ch.SV_phi[iSV], ch.SV_mass[iSV])
                if vSV.DeltaR(fatVec)    < 0.8:                    nSV_fat += 1
                if vSV.DeltaR(jetVecOut) < 0.4 and jetIdxOut >= 0: nSV_jet += 1
                if vSV.DeltaR(isrVec)    < 0.4 and isrIdx    >= 0: nSV_ISR += 1
                
            if nSV_fat + nSV_jet < CUTS.MIN_SV_FAT_JET: continue

            # print 'Here10'

            ################################################################################
            ###  At this point, we have *the* leading H --> aa --> 4b candidate AK8 jet  ###
            ################################################################################
            nCandsSel += 1


            ## Loop over all qualifying AK4 jets which overlap the AK8 candidate jet
            for jetIdx in candJets:
                ## Only consider the first two matched AK4 jets
                if jetIdxInA >= 0 and jetIdxInB >= 0: break
                ## Define AK4 jet 4-vectors
                iJetVec = R.TLorentzVector()
                iJetVec.SetPtEtaPhiM(ch.Jet_pt[jetIdx], ch.Jet_eta[jetIdx], ch.Jet_phi[jetIdx], ch.Jet_mass[jetIdx])
                ## Jet must have pT > 140 GeV
                if iJetVec.Pt() <= 140: continue
                ## AK4 and AK8 jet must overlap
                if iJetVec.DeltaR(fatVec) > MAX_DR_FAT: continue
                ## Assign jet as first or second overlapping AK4 jet
                if jetIdxInA < 0:
                    jetIdxInA = jetIdx
                    jetVecInA = iJetVec
                elif jetIdxInB < 0:
                    jetIdxInB = jetIdx
                    jetVecInB = iJetVec
            ## End loop: for jetIdx in candJets


            ##############################
            ## JetHT trigger "region" cuts
            ##############################

            ## Case where >= 2 AK4 jets must overlap the AK8 candidate jet (regions ABC and CX)
            if   CUTS.NUM_JETS_IN == 'ge2':
                if jetIdxInB < 0: continue

            ## Case where <= 1 AK4 jets must overlap the AK8 candidate jet (region AB)
            if CUTS.NUM_JETS_IN == 'le1':
                if jetIdxInB >= 0: break

            ## Case where AK8 jet must not surpass a maximum pT threshold (region CX)
            if fatVec.Pt() > CUTS.MAX_PT_FAT: break


            ###############################
            ## ParkingBPH trigger muon cuts
            ###############################

            ## Track the charges of muons matched to fat jet
            matched_mu_charge = []
            ## Find the minimum dR between signal or ISR jets and triggering muon
            min_dR_mu_Higg,min_dR_mu_ISR = (99.,99.)
            ## If there are multiple muons, pick the one with the maximum effective luminosity
            iMu_sel = -1
            max_eff_lumi,max_PU_wgt        = (-99.,-99.)
            trg_mu_pt,trg_mu_eta,trg_mu_IP = (-99.,-99.,-99.)
            for iMu in range(nMuon):

                muVec = R.TLorentzVector()
                muVec.SetPtEtaPhiM(ch.Muon_pt[iMu], ch.Muon_eta[iMu], ch.Muon_phi[iMu], MASS_MU)
                
                ## Protect against dxyErr = 0
                mu_IP = abs( ch.Muon_dxy[iMu] / max(ch.Muon_dxyErr[iMu], pow(2,-11)) )

                ## Select only "good" (well-reconstructed) muons
                if ch.Muon_softId[iMu]   !=  1:  continue  ## High-efficiecy, low-fake RECO muon ID
                if abs(ch.Muon_ip3d[iMu]) > 0.5: continue  ## Remove highly-displaced muons
                if     muVec.Pt()         < 7.0: continue  ## Minimum trigger pT cut
                if abs(muVec.Eta())       > 2.4: continue  ## Maximum trigger eta cut

                ## ## TODO: study muon-jet dR matching requirement - AWB 2021.03.02
                # ## Require muon to match signal AK8 jet with dR < 0.8, or signal or ISR AK4 jet with dR < 0.4
                # if                muVec.DeltaR(fatVec) > MAX_DR_FAT  and \
                #    (jetIdxOut < 0 or muVec.DeltaR(jetVecOut) > MAX_DR_JET) and \
                #    (isrIdx    < 0 or muVec.DeltaR(isrVec)    > MAX_DR_JET): continue

                ## Track charges of matched muons
                matched_mu_charge.append( ch.Muon_charge[iMu] )

                ## Select only "triggering" muon (based on IP cut)
                if mu_IP < 2.0: continue  ## Minimum trigger IP cut

                ## Access the effective luminosity and pileup weight for this pT/eta/IP bin
                mu_eff_lumi = EVT_WGTS.GetLumi        (muVec.Pt(), mu_IP, muVec.Eta())
                mu_PU_wgt   = EVT_WGTS.GetWgtPileupBPH(muVec.Pt(), mu_IP, muVec.Eta(), ch.PV_npvsGood)
                
                ## If this muon has highest effective luminosity, use it
                if not mu_eff_lumi > max_eff_lumi: continue
                iMu_sel      = iMu
                max_eff_lumi = mu_eff_lumi
                max_PU_wgt   = mu_PU_wgt
                if jetIdxOut < 0: min_dR_mu_Higg =     muVec.DeltaR(fatVec)
                else:             min_dR_mu_Higg = min(muVec.DeltaR(fatVec), muVec.DeltaR(jetVecOut))
                if isrIdx >= 0:   min_dR_mu_ISR  =     muVec.DeltaR(isrVec)
                trg_mu_pt  = muVec.Pt()
                trg_mu_eta = muVec.Eta()
                trg_mu_IP  = mu_IP

            ## End loop over muons

            ## Require a minimum number of muons
            if len(matched_mu_charge) < CUTS.NUM_MUONS: continue
            ## Require triggering muon to come from one of the jets
            if SEL_MU and iMu_sel < 0: continue
            if SEL_MU and max_eff_lumi <= 0.0:
                print '\n\n*** Weird case where Muon[%d] was selected, but max_eff_lumi = %f!!! ***' % (iMu_sel, max_eff_lumi)
                print '    pt = %f, eta = %f, IP = %f\n\n' % (ch.Muon_pt[iMu_sel], ch.Muon_eta[iMu_sel], ch.Muon_dxy[iMu_sel]/ch.Muon_dxyErr[iMu_sel])
                continue


            ## Weight MC events by effective luminosity and pileup
            if isMC:
                Lumi_wgt = 59.83  ## Full 2018 luminosity
                PU_wgt   = 1.0    ## In the absence of any other reweighting scheme

                if SEL_MU:  ## Lumi and PU weights for ParkingBPH
                    Lumi_wgt = max_eff_lumi
                    PU_wgt   = max(max_PU_wgt, 0.001)  ## TODO: check how often the PU weight is tiny - AWB 2021.03.02

                elif CUTS.TRG_JET_HT and not CUTS.TRG_BPH:
                    if 'AB' in CUTS.TRG_JET_HT and ch.HLT_AK8PFJet500:
                        Lumi_wgt = 59.83
                        PU_wgt   = EVT_WGTS.GetWgtPileupJET(ch.Pileup_nTrueInt, True)
                    else:
                        Lumi_wgt = 54.54
                        PU_wgt   = EVT_WGTS.GetWgtPileupJET(ch.Pileup_nTrueInt, False)

                    ## Apply JetHT efficiency scale factors
                    if CUTS.TRG_JET_HT == 'ABC': Trg_wgt = EVT_WGTS.GetTrigSF('ABC', fatVec.Pt())
                    if CUTS.TRG_JET_HT == 'AB':  Trg_wgt = EVT_WGTS.GetTrigSF('AB',  fatVec.Pt())
                    if CUTS.TRG_JET_HT == 'CX':  Trg_wgt = EVT_WGTS.GetTrigSF('CX',  min(ch.Jet_btagDeepB[jetIdxInA], ch.Jet_btagDeepB[jetIdxInB]))

                WGT *= (Lumi_wgt * PU_wgt * Trg_wgt)
                ## Use the first candidate to weight the "whole event" as well
                if nCandsSel == 1:
                    Evt_wgt *= (Lumi_wgt * PU_wgt * Trg_wgt)

            ## End conditional: if isMC


            ## Veto data and weight MC based on HE- veto region
            if fatVec.Eta() < -1.17 and fatVec.Phi() > -1.97 and fatVec.Phi() < -0.47:
                if isData and ch.run > 319077:
                    HEM_wgt = 0.0
                if isMC:
                    if CUTS.TRG_JET_HT != 'CX' and ch.HLT_AK8PFJet500:
                        HEM_wgt = 21.09 / 59.83
                    else:
                        HEM_wgt = 15.80 / 54.54
            ## End conditional: if fatVec.Eta() < -1.17 and fatVec.Phi() > -1.97 and fatVec.Phi() < -0.47
            WGT *= HEM_wgt
            if nCandsSel == 1:
                Evt_wgt *= HEM_wgt


            ## Make sure no events are selected twice
            lheIdx = (0 if isData else math.floor(ch.LHE_HT))
            evtIdx = '%d_%d_%d_%d_%d' % (ch.run, ch.luminosityBlock, ch.event, lheIdx, fatIdx)

            if evtIdx in event_list:
                print '\n\n***** MAJOR ISSUE!!! Event %s appears twice in input file(s)!!!  Skipping. *****\n\n' % evtIdx
                continue
            event_list.append(evtIdx)

            # print 'run = %d, LS = %d, event = %d' % (ch.run, ch.luminosityBlock, ch.event)



            #################################################
            ###  Save NN input variables for Brooks' NNs  ###
            #################################################

            pass_NN = True  ## Track whether event passes NN cuts
            
            ## Store 8 input variable values for this event in proper order

            # NN_input_dict['BMM_2020_06_15'][evtIdx] = [ fatVec.Pt(),
            #                                             abs(fatVec.Eta()),
            #                                             fatVec.Phi(),
            #                                             fatVec.M(),
            #                                             ch.FatJet_btagCSVV2[fatIdx],
            #                                             ch.FatJet_btagDeepB[fatIdx],
            #                                             fatVecS.M(),
            #                                             ch.FatJet_btagDDBvL[fatIdx] ]

            NN_input_dict['BMM_2020_09_29'][evtIdx] = [ fatVec.Pt(),
                                                        abs(fatVec.Eta()),
                                                        fatVec.M(),
                                                        ch.FatJet_btagCSVV2[fatIdx],
                                                        ch.FatJet_btagDeepB[fatIdx],
                                                        fatVecS.M(),
                                                        ch.FatJet_btagDDBvL[fatIdx],
                                                        ch.FatJet_deepTagMD_H4qvsQCD[fatIdx] ]
            
            NN_input_dict['BMM_2020_09_29_btag'][evtIdx] = [ ch.FatJet_btagDeepB[fatIdx],
                                                             ch.FatJet_deepTagMD_H4qvsQCD[fatIdx],
                                                             ch.FatJet_btagDDBvL[fatIdx],
                                                             ch.FatJet_btagCSVV2[fatIdx] ]
            
            NN_input_dict['BMM_2020_09_29_phys'][evtIdx] = [ fatVec.Pt(),
                                                             abs(fatVec.Eta()),
                                                             fatVec.M(),
                                                             fatVecS.M() ]


            ## Plot NN output values
            for NN in NN_list:
                if not DO_NN: continue
                if not evtIdx in NN_output_dict[NN].keys():
                    print 'Event %s not found in %s' % (evtIdx, NN)
                    nLostNN += 1
                    break
                ## Check for mismatches between input and quantities in this event
                mismatch = False
                for iVar in range(len(NN_input_dict[NN][evtIdx])):
                    if NN_input_dict[NN][evtIdx][iVar] != NN_output_dict[NN][evtIdx][iVar]:
                        print '\n\nMismatch between NN input and output!!! Skipping.\n'
                        print NN_input_dict[NN][evtIdx]
                        print NN_output_dict[NN][evtIdx]
                        mismatch = True
                        break
                        # sys.exit()
                if mismatch: continue

                if NN == 'BMM_2020_06_15':  ## Weighted training stored 2nd, after un-weighted training
                    NN_score = NN_output_dict[NN][evtIdx][ len(NN_input_dict[NN][evtIdx]) + 1 ]
                else:
                    NN_score = NN_output_dict[NN][evtIdx][ len(NN_input_dict[NN][evtIdx]) ]

                if len(NN_quant[NN]) != NN_bins[0] - 1:
                    print '\n\n*** Inconsistent binning for %s! %d quantiles, %d histogram bins. Quitting.' % (NN, NN_quant[NN], NN_bins[0])
                    sys.exit()

                ## Only select events with a given minimum NN score
                if NN == 'BMM_2020_09_29' and MIN_NN > 0:
                    if NN_score < NN_quant[NN][int(MIN_NN*10)]:
                        pass_NN = False
                        break

                ## Put NN score into correct quantile
                nBin = NN_bins[0]
                binW = (NN_bins[2] - NN_bins[1]) / nBin
                for iBin in range(nBin-1):
                    if NN_score < NN_quant[NN][iBin]:
                        hst[NN].Fill( NN_bins[1] + binW*(iBin+0.5), WGT )
                        break
                if NN_score >= NN_quant[NN][-1]:
                    hst[NN].Fill( NN_bins[2] - binW*0.5, WGT )

            ## End loop: for NN in NN_list

            ## Only select events with a given minimum NN score
            if not pass_NN: continue


            # ### *********************************************** ###
            # iPrintSig += 1
            # print '\nSignal candidate #%d: run = %d, luminosityBlock = %d, event = %d, PV_npvsGood = %d' % (iPrintSig, ch.event, ch.luminosityBlock, ch.event, ch.PV_npvsGood)
            # print 'FatJet index %d, pt = %.2f, eta = %.2f, mass = %.2f, msoftdrop = %.2f, btagDDBvL = %.3f, btagDeepB = %.3f' % (fatIdx, fatVec.Pt(), fatVec.Eta(), fatVec.M(), fatVecS.M(), ch.FatJet_btagDDBvL[fatIdx], ch.FatJet_btagDeepB[fatIdx])
            # print 'Muon index %d, softId = %d, pt = %.2f, eta = %.2f, dxy/dxyErr = %.2f, ip3d = %.2f' % (iMu_sel, ch.Muon_softId[iMu_sel], trg_mu_pt, trg_mu_eta, trg_mu_IP, ch.Muon_ip3d[iMu_sel])
            # if iPrintSig >= 100: exit()
            # ### *********************************************** ###


            #######################################################
            ###  Fill histograms of Higgs candidate properties  ###
            #######################################################
                
            if isSig and WGT > 25.0:
                print '\n\n\n*** SUPER-BIZZARE EVENT WITH WGT = %.3f (Evt_wgt = %.3f) ***' % (WGT, Evt_wgt)
                print 'PS_wgt = %.3f, Lumi_wgt = %.3f, PU_wgt = %.3f, Xsec_wgt = %.3f, ggH_pt_wgt = %.3f' % (PS_wgt, Lumi_wgt, PU_wgt, Xsec_wgt, ggH_pt_wgt)
                # sys.exit()


            ## Per-candidate weights
            hst['wgt_cand']    .Fill( min( 4.99, max( -4.99, math.log(max(WGT, 0.000001) /PS_wgt, 10) ) ) )
            hst['wgt_cand_wgt'].Fill( min( 4.99, max( -4.99, math.log(max(WGT, 0.000001) /PS_wgt, 10) ) ), WGT )
            
            hst['wgt_lumi']    .Fill( min( 64.99, max( 0.01, Lumi_wgt) ) )
            hst['wgt_lumi_wgt'].Fill( min( 64.99, max( 0.01, Lumi_wgt) ), Lumi_wgt )
            
            hst['wgt_PU']    .Fill( min( 9.99, max( 0.01, PU_wgt) ) )
            hst['wgt_PU_wgt'].Fill( min( 9.99, max( 0.01, PU_wgt) ), PU_wgt )

            hst['wgt_trg']    .Fill( min( 64.99, max( 0.01, Trg_wgt) ) )
            hst['wgt_trg_wgt'].Fill( min( 64.99, max( 0.01, Trg_wgt) ), Trg_wgt )
            
            hst['wgt_HEM']    .Fill( min( 1.049, max( -0.049, HEM_wgt) ) )
            hst['wgt_HEM_wgt'].Fill( min( 1.049, max( -0.049, HEM_wgt) ), HEM_wgt )
            
            ## Track the sign of the Higgs eta value (+/-1)
            etaPM = -1.0 + 2.0*(higVec.Eta() > 0)

            ## Triggering muon matched to candidate or ISR jet
            hst['nMuon'].Fill( len(matched_mu_charge), WGT )
            if max_eff_lumi > 0:
                hst['mu_pt_log']    .Fill( max( log_bins[1] +0.01, min( log_bins[2] -0.01, math.log(trg_mu_pt, 2) ) ), WGT )
                hst['mu_pt']        .Fill( max( ptS_bins[1] +0.01, min( ptS_bins[2] -0.01,        trg_mu_pt ) ), WGT )
                hst['mu_eta']       .Fill( max( eta_bins[1]+0.01, min( eta_bins[2]-0.01,          trg_mu_eta    ) ), WGT )
                hst['mu_IP_log']    .Fill( max( IP_bins[1] +0.01, min( IP_bins[2] -0.01, math.log(trg_mu_IP, 2) ) ), WGT )
                hst['mu_IP']        .Fill( max( 0.01,             min( 19.99,                     trg_mu_IP     ) ), WGT )
                hst['mu_charge']    .Fill( matched_mu_charge[0],   WGT )
                hst['mu_set_charge'].Fill( sum(matched_mu_charge), WGT )

            ## Invariant mass of two highest-pT "extra" AK4 jets in the event
            if vbfIdx >= 0: hst['mass_VBF'].Fill( math.log((isrVec+vbfVec).M(), 2), WGT )

            ## Number of GEN b quarks matched to Higgs candidate
            hst['nGenB_Higg'].Fill( min( 6.0, fatGenMap.count(1) + jetGenMap.count(1) ), WGT )
            hst['nGenB_fatH'].Fill( min( 6.0, fatGenMap.count(1) ), WGT )
            hst['nGenB_jetH'].Fill( min( 6.0, jetGenMap.count(1) ), WGT )
            
            ## Higgs candidate pT, mass, and eta
            hst['pt_log_fatH'].Fill( max(log_bins[1]+0.01, min(log_bins[2]-0.01, math.log(fatVec.Pt(), 2) ) ), WGT )
            hst['pt_log_Higg'].Fill( max(log_bins[1]+0.01, min(log_bins[2]-0.01, math.log(higVec.Pt(), 2) ) ), WGT )
            
            hst['pt_fatH'].Fill( max(ptL_bins[1]+0.01, min(ptL_bins[2]-0.01, fatVec.Pt() ) ), WGT )
            hst['pt_Higg'].Fill( max(ptL_bins[1]+0.01, min(ptL_bins[2]-0.01, higVec.Pt() ) ), WGT )
            
            hst['mass_fatH'] .Fill( max( mass_bins[1]+0.01, min( mass_bins[2]-0.01, fatVec.M()  ) ), WGT )
            hst['msoft_fatH'].Fill( max( mass_bins[1]+0.01, min( mass_bins[2]-0.01, fatVecS.M() ) ), WGT )
            hst['mass_Higg'] .Fill( max( mass_bins[1]+0.01, min( mass_bins[2]-0.01, higVec.M()  ) ), WGT )
            hst['msoft_Higg'].Fill( max( mass_bins[1]+0.01, min( mass_bins[2]-0.01, higVecS.M() ) ), WGT )

            hst['eta_fatH'].Fill( max( eta_bins[1]+0.01, min( eta_bins[2]-0.01, fatVec.Eta()*etaPM ) ), WGT )
            hst['eta_Higg'].Fill( max( eta_bins[1]+0.01, min( eta_bins[2]-0.01, higVec.Eta()*etaPM ) ), WGT )

            if fatVec.Eta() > 0:
                hst['phi_fatH_Pos'].Fill( max( phi_bins[1]+0.01, min( phi_bins[2]-0.01, fatVec.Phi() ) ), WGT )
            else:
                hst['phi_fatH_Neg'].Fill( max( phi_bins[1]+0.01, min( phi_bins[2]-0.01, fatVec.Phi() ) ), WGT )

            
            if jetIdxOut >= 0:
                hst['pt_log_jetH'].Fill( max(  log_bins[1]+0.01, min(  log_bins[2]-0.01, math.log(jetVecOut.Pt(), 2) ) ), WGT )
                hst['pt_jetH']    .Fill( max(  ptL_bins[1]+0.01, min(  ptL_bins[2]-0.01, jetVecOut.Pt() ) ),              WGT )
                hst['mass_jetH']  .Fill( max( mass_bins[1]+0.01, min( mass_bins[2]-0.01, jetVecOut.M() ) ),               WGT )
                hst['eta_jetH']   .Fill( max(  eta_bins[1]+0.01, min(  eta_bins[2]-0.01, jetVecOut.Eta()*etaPM ) ),       WGT )
                
                hst['mass_fatA'] .Fill( max( mass_bins[1]+0.01, min( mass_bins[2]-0.01, fatVec.M()  ) ), WGT )
                hst['msoft_fatA'].Fill( max( mass_bins[1]+0.01, min( mass_bins[2]-0.01, fatVecS.M() ) ), WGT )
                    
                    
                hst['dR_Higg']  .Fill( max(   dR_bins[1]+0.01, min(   dR_bins[2]-0.01,     fatVec.DeltaR  (jetVecOut)  ) ),        WGT )
                hst['dPhi_Higg'].Fill( max( dPhi_bins[1]+0.01, min( dPhi_bins[2]-0.01, abs(fatVec.DeltaPhi(jetVecOut)) ) ),        WGT )
                hst['dEta_Higg'].Fill( max( dEta_bins[1]+0.01, min( dEta_bins[2]-0.01, (fatVec.Eta() - jetVecOut.Eta())*etaPM ) ), WGT )
                    
            if isrIdx >= 0:
                hst['pt_log_ISR'].Fill( max(  log_bins[1]+0.01, min(  log_bins[2]-0.01, math.log(isrVec.Pt(), 2) ) ), WGT )
                hst['pt_ISR']    .Fill( max(  ptL_bins[1]+0.01, min(  ptL_bins[2]-0.01, isrVec.Pt() ) ),              WGT )
                hst['eta_ISR']   .Fill( max(  eta_bins[1]+0.01, min(  eta_bins[2]-0.01, isrVec.Eta()*etaPM ) ),       WGT )

                if isrVec.Eta() > 0:
                    hst['phi_ISR_Pos'].Fill( max( phi_bins[1]+0.01, min( phi_bins[2]-0.01, isrVec.Phi() ) ), WGT )
                else:
                    hst['phi_ISR_Neg'].Fill( max( phi_bins[1]+0.01, min( phi_bins[2]-0.01, isrVec.Phi() ) ), WGT )
                
                hst['dR_ISR']  .Fill( max(   dR_bins[1]+0.01, min(   dR_bins[2]-0.01,     higVec.DeltaR  (isrVec)  ) ),        WGT )
                hst['dPhi_ISR'].Fill( max( dPhi_bins[1]+0.01, min( dPhi_bins[2]-0.01, abs(higVec.DeltaPhi(isrVec)) ) ),        WGT )
                hst['dEta_ISR'].Fill( max( dEta_bins[1]+0.01, min( dEta_bins[2]-0.01, (higVec.Eta() - isrVec.Eta())*etaPM ) ), WGT )

                if min_dR_mu_Higg > 0.8:
                    hst['dR_mu_ISR'] .Fill( max( dRm_bins[1]+0.01, min( dRm_bins[2]-0.01, min_dR_mu_ISR  ) ), WGT )

            if min_dR_mu_ISR > 0.4:
                hst['dR_mu_Higg'].Fill( max( dRm_bins[1]+0.01, min( dRm_bins[2]-0.01, min_dR_mu_Higg ) ), WGT )

            hst['deepB_fatH'] .Fill( max( -0.099, ch.FatJet_btagDeepB[fatIdx] ), WGT )
            hst['doubB_fatH'] .Fill( ch.FatJet_btagDDBvL[fatIdx], WGT )
            hst['bbVsLF_fatH'].Fill( ch.FatJet_deepTagMD_bbvsLight[fatIdx], WGT )
            hst['ZvsQCD_fatH'].Fill( ch.FatJet_deepTagMD_ZbbvsQCD[fatIdx], WGT )
            hst['HvsQCD_fatH'].Fill( ch.FatJet_deepTagMD_HbbvsQCD[fatIdx], WGT )
            hst['XvsQCD_fatH'].Fill( ch.FatJet_deepTagMD_ZHbbvsQCD[fatIdx], WGT )
            hst['fourQ_fatH'] .Fill( ch.FatJet_deepTagMD_H4qvsQCD[fatIdx], WGT )
            hst['nSV_fatH']   .Fill( nSV_fat, WGT )

            if jetIdxOut >= 0:
                hst['QGL_jetH']   .Fill( ch.Jet_qgl[jetIdxOut], WGT )
                hst['jetID_jetH'] .Fill( ch.Jet_jetId[jetIdxOut], WGT )
                hst['jetPU_jetH'] .Fill( ch.Jet_puId[jetIdxOut], WGT )
                hst['deepB_jetH'] .Fill( max( -0.099, ch.Jet_btagDeepB[jetIdxOut] ), WGT )
                hst['flavB_jetH'] .Fill( ch.Jet_btagDeepFlavB[jetIdxOut], WGT )
                hst['nConst_jetH'].Fill( ch.Jet_nConstituents[jetIdxOut], WGT )
                hst['nSV_jetH']   .Fill( nSV_jet, WGT )
                
                hst['deepB_fatA'].Fill( max( -0.099, ch.FatJet_btagDeepB[fatIdx] ), WGT )
                hst['doubB_fatA'].Fill( ch.FatJet_btagDDBvL[fatIdx], WGT )
                hst['fourQ_fatA'].Fill( ch.FatJet_deepTagMD_H4qvsQCD[fatIdx], WGT )
                hst['nSV_fatA']  .Fill( nSV_fat, WGT )
                
            if isrIdx >= 0:
                hst['QGL_ISR']   .Fill( ch.Jet_qgl[isrIdx], WGT )
                hst['jetID_ISR'] .Fill( ch.Jet_jetId[isrIdx], WGT )
                hst['jetPU_ISR'] .Fill( ch.Jet_puId[isrIdx], WGT )
                hst['deepB_ISR'] .Fill( max( -0.099, ch.Jet_btagDeepB[isrIdx] ), WGT )
                hst['flavB_ISR'] .Fill( ch.Jet_btagDeepFlavB[isrIdx], WGT )
                hst['nConst_ISR'].Fill( ch.Jet_nConstituents[isrIdx], WGT )
                hst['nSV_ISR']   .Fill( nSV_ISR, WGT )

            if jetIdxInA >= 0:
                hst['pt_jetInA']    .Fill( max(  ptL_bins[1]+0.01, min(  ptL_bins[2]-0.01, jetVecInA.Pt() ) ), WGT )
                hst['deepB_jetInA'] .Fill( max( -0.099, ch.Jet_btagDeepB[jetIdxInA] ), WGT )
            if jetIdxInB >= 0:
                hst['pt_jetInB']    .Fill( max(  ptL_bins[1]+0.01, min(  ptL_bins[2]-0.01, jetVecInB.Pt() ) ), WGT )
                hst['deepB_jetInB'] .Fill( max( -0.099, ch.Jet_btagDeepB[jetIdxInB] ), WGT )


            hst['nSV_Higg'] .Fill( nSV_fat + nSV_jet, WGT )
            hst['nSV_extra'].Fill( len(ch.SV_eta) - (nSV_fat + nSV_jet + nSV_ISR), WGT )

            ## More properties of AK8 jet and sub-jets
            hst['tau1_fatH'].Fill( max( -0.049, min( 0.649, ch.FatJet_tau1[fatIdx] ) ), WGT )
            hst['tau2_fatH'].Fill( max( -0.049, min( 0.349, ch.FatJet_tau2[fatIdx] ) ), WGT )
            hst['tau3_fatH'].Fill( max( -0.049, min( 0.249, ch.FatJet_tau3[fatIdx] ) ), WGT )
            hst['tau4_fatH'].Fill( max( -0.049, min( 0.149, ch.FatJet_tau4[fatIdx] ) ), WGT )
            
            hst['n2b1_fatH'].Fill( max( -0.099, min( 0.499, ch.FatJet_n2b1[fatIdx] ) ), WGT )
            hst['n3b1_fatH'].Fill( max( -0.099, min(  4.09, ch.FatJet_n3b1[fatIdx] ) ), WGT )

            if ch.FatJet_subJetIdx1[fatIdx] >= 0:
                idx1 = ch.FatJet_subJetIdx1[fatIdx]
                hst['sub1_pt_frac_fatH'].Fill( ch.SubJet_pt[idx1] / ch.FatJet_pt[fatIdx],         WGT )
                hst['sub1_mass_fatH']   .Fill( min( 99.9, ch.SubJet_mass[idx1] ),                 WGT )
                hst['sub1_deepB_fatH']  .Fill( max( -0.09, ch.SubJet_btagDeepB[idx1] ),           WGT )
                hst['sub1_tau1_fatH']   .Fill( max( -0.049, min( 0.649, ch.SubJet_tau1[idx1] ) ), WGT )
                hst['sub1_tau2_fatH']   .Fill( max( -0.049, min( 0.349, ch.SubJet_tau2[idx1] ) ), WGT )
                hst['sub1_tau3_fatH']   .Fill( max( -0.049, min( 0.249, ch.SubJet_tau3[idx1] ) ), WGT )
                hst['sub1_tau4_fatH']   .Fill( max( -0.049, min( 0.149, ch.SubJet_tau4[idx1] ) ), WGT )
                hst['sub1_n2b1_fatH']   .Fill( max(  0.149, min( 0.549, ch.SubJet_n2b1[idx1] ) ), WGT )
                hst['sub1_n3b1_fatH']   .Fill( max( -0.099, min(  4.09, ch.SubJet_n3b1[idx1] ) ), WGT )
            if ch.FatJet_subJetIdx2[fatIdx] >= 0:
                idx2 = ch.FatJet_subJetIdx2[fatIdx]
                hst['sub2_pt_frac_fatH'].Fill( ch.SubJet_pt[idx2] / ch.FatJet_pt[fatIdx],         WGT )
                hst['sub2_mass_fatH']   .Fill( min( 99.9, ch.SubJet_mass[idx2] ),                 WGT )
                hst['sub2_deepB_fatH']  .Fill( max( -0.09, ch.SubJet_btagDeepB[idx2] ),           WGT )
                hst['sub2_tau1_fatH']   .Fill( max( -0.049, min( 0.649, ch.SubJet_tau1[idx2] ) ), WGT )
                hst['sub2_tau2_fatH']   .Fill( max( -0.049, min( 0.349, ch.SubJet_tau2[idx2] ) ), WGT )
                hst['sub2_tau3_fatH']   .Fill( max( -0.049, min( 0.249, ch.SubJet_tau3[idx2] ) ), WGT )
                hst['sub2_tau4_fatH']   .Fill( max( -0.049, min( 0.149, ch.SubJet_tau4[idx2] ) ), WGT )
                hst['sub2_n2b1_fatH']   .Fill( max(  0.149, min( 0.549, ch.SubJet_n2b1[idx2] ) ), WGT )
                hst['sub2_n3b1_fatH']   .Fill( max( -0.099, min(  4.09, ch.SubJet_n3b1[idx2] ) ), WGT )
            if ch.FatJet_subJetIdx1[fatIdx] >= 0 and ch.FatJet_subJetIdx2[fatIdx] >= 0:
                idx1 = ch.FatJet_subJetIdx1[fatIdx]
                idx2 = ch.FatJet_subJetIdx2[fatIdx]
                subVec1,subVec2 = (R.TLorentzVector(),R.TLorentzVector())
                subVec1.SetPtEtaPhiM(ch.SubJet_pt[idx1], ch.SubJet_eta[idx1], ch.SubJet_phi[idx1], ch.SubJet_mass[idx1])
                subVec2.SetPtEtaPhiM(ch.SubJet_pt[idx2], ch.SubJet_eta[idx2], ch.SubJet_phi[idx2], ch.SubJet_mass[idx2])
                hst['sub12_pt_frac_fatH']  .Fill( (subVec1+subVec2).Pt() / ch.FatJet_pt  [fatIdx], WGT )
                hst['sub12_mass_frac_fatH'].Fill( (subVec1+subVec2).M()  / ch.FatJet_mass[fatIdx], WGT )  ## Identical to soft-drop mass
                hst['sub12_dR_fatH']       .Fill( (subVec1+subVec2).DeltaR(fatVec),  WGT )
                hst['sub12_dR']            .Fill(  subVec1         .DeltaR(subVec2), WGT )


        ## End loop: for fatIdx in candFats


        ###################################
        ###  Fill per-event histograms  ###
        ###################################

        ## Fill number of Higgs candidates in each event
        hst['nCandsSel'].Fill(nCandsSel,   Evt_wgt)
        # hst['nCandsSV'] .Fill(nCandsSelSV, Evt_wgt)
        # hst['nCandsHT'] .Fill(nCandsSelHT, Evt_wgt)

        ## For other plots, only consider events with at least one reconstructed Higgs candidate
        if nCandsSel < 1: continue
        nEvtSel += 1

        ## Fill event weight histograms (un-weighted and weighted)
        hst['wgt_evt']    .Fill( min( 4.99, max( -4.99, math.log( max(Evt_wgt, 0.000001) / PS_wgt, 10) ) ) )
        hst['wgt_evt_wgt'].Fill( min( 4.99, max( -4.99, math.log( max(Evt_wgt, 0.000001) / PS_wgt, 10) ) ), Evt_wgt )

        hst['wgt_ggH_pt']    .Fill( ggH_pt_wgt )
        hst['wgt_ggH_pt_wgt'].Fill( ggH_pt_wgt, ggH_pt_wgt )

        hst['wgt_Xsec']    .Fill( min( 1.99, max( -3.99, math.log(Xsec_wgt, 10) ) ) )
        hst['wgt_Xsec_wgt'].Fill( min( 1.99, max( -3.99, math.log(Xsec_wgt, 10) ) ), Xsec_wgt )

        ## Pileup and jets
        hst['nPV_all'] .Fill( min( PV_bins[2]-0.01, max( PV_bins[1]+0.01, ch.PV_npvs     ) ), Evt_wgt )
        hst['nPV_good'].Fill( min( PV_bins[2]-0.01, max( PV_bins[1]+0.01, ch.PV_npvsGood ) ), Evt_wgt )
        hst['nJet']    .Fill(len(goodJets), Evt_wgt)

        ## LHE HT
        if isMC:
            hst['LHE_HT']    .Fill( min( 2499.0, max( 0.001,          ch.LHE_HT     ) ), Evt_wgt )
            hst['LHE_HT_log'].Fill( min( 11.499, max( 3.501, math.log(ch.LHE_HT, 2) ) ), Evt_wgt )

        # ## Fill trigger efficiency histogram
        # if ch.HLT_Mu12_IP6_part0:
        #     hst['trig_eff'].Fill( 1, (Evt_wgt / Lumi_wgt) * 34.79 )
        # elif ch.HLT_Mu9_IP6_part0:
        #     hst['trig_eff'].Fill( 1, (Evt_wgt / Lumi_wgt) * 33.67 )
        # elif ch.HLT_Mu9_IP5_part0:
        #     hst['trig_eff'].Fill( 1, (Evt_wgt / Lumi_wgt) * 20.95 )
        # elif ch.HLT_Mu8_IP5_part0:
        #     hst['trig_eff'].Fill( 1, (Evt_wgt / Lumi_wgt) * 8.26 )
        # elif ch.HLT_Mu7_IP4_part0:
        #     hst['trig_eff'].Fill( 1, (Evt_wgt / Lumi_wgt) * 6.94 )
        # elif ch.HLT_Mu8_IP3_part0:
        #     hst['trig_eff'].Fill( 1, (Evt_wgt / Lumi_wgt) * 1.58 )

        # if (ch.L1_Mu3_Jet120er2p5_dR_Max0p8 and ch.HLT_BTagMu_AK8DiJet170_Mu5_noalgo):
        #     hst['trig_eff'].Fill( 2, (Evt_wgt / Lumi_wgt) * 16.59 )
        # if (ch.L1_Mu3_Jet120er2p5_dR_Max0p4 and ch.HLT_BTagMu_AK4DiJet170_Mu5_noalgo):
        #     hst['trig_eff'].Fill( 3, (Evt_wgt / Lumi_wgt) * 20.23 )
        # if (ch.L1_Mu3_Jet120er2p5_dR_Max0p4 and ch.HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71):
        #     hst['trig_eff'].Fill( 4, (Evt_wgt / Lumi_wgt) * 59.8 )
        # if (ch.L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6 and \
        #     ch.HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71):
        #     hst['trig_eff'].Fill( 5, (Evt_wgt / Lumi_wgt) * 59.8 )
        # if ( (ch.L1_DoubleJet112er2p3_dEta_Max1p6 or ch.L1_DoubleJet150er2p5) and \
        #      ch.HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71 ):
        #     hst['trig_eff'].Fill( 6, (Evt_wgt / Lumi_wgt) * 59.8 )
        # if (ch.HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5):
        #     hst['trig_eff'].Fill( 7, (Evt_wgt / Lumi_wgt) * 59.8 )


    ## End loop over events in chain (iEvt)

    print '\n\n%d events processed, %d with %d fat jets, %d muons, %d trigger, %d GEN 4b, %d cat Fat1, %d selected' % (iEvt, nEvtFatJet, CUTS.NUM_FATS, nEvtMuons, nEvtTrig, nEvtGen4B, nEvtCatFat1, nEvtSel)
    print 'Lost %d events, not found in NN_output_dir' % nLostNN


##########################################
## Save the histograms and NN input values
##########################################

    ## Save NN inputs to .pkl file, only if we're processing all events
    if MAX_EVT < 0:
        for NN in NN_list:
            with open('data/NN_inputs/%s_%s.pkl' % (BASE_STR, NN), 'wb') as NN_file:
                pickle.dump(NN_input_dict[NN], NN_file, protocol=2)
        print 'Saved NN inputs to .pkl files in data/NN_inputs/'
    else:
        print '\n*** DID NOT SAVE NN INPUTS!!!  PROCESSED ONLY %d EVENTS. ***' % MAX_EVT


    ## Navigate to output ROOT file to save histograms
    out_file.cd()

    c0 = R.TCanvas('c0')
    c0.cd()

    ## Figure out if PU re-weighting normalizes to 1.0, as expected
    PU_wgt_corr = 1.0
    if hst['wgt_PU_wgt'].Integral() != hst['wgt_PU'].Integral():
        print '\n*** PU re-weighting not quite normalized!!!'
        PU_wgt_corr = hst['wgt_PU'].Integral() / hst['wgt_PU_wgt'].Integral()
        # print 'Ratio = %.4f, will correct all plots.' % PU_wgt_corr
        print '\n\n*** Ratio = %.4f, in future should correct all plots: TODO! ***\n\n' % PU_wgt_corr

    ## Common formatting for histograms
    for key in hst.keys():
        hst[key].SetLineWidth(2)
        ## Correct PU re-weighting normalization
        # if not key.startswith('wgt_'):
        #     hst[key].Scale(PU_wgt_corr)
    
    hst['trig_eff'].SetLineColor(R.kBlack)
    hst['trig_eff'].Draw('histtext')
    hst['trig_eff'].GetXaxis().SetBinLabel(2, 'B Parking')
    hst['trig_eff'].GetXaxis().SetBinLabel(3, 'AK8 Mu DiJet')
    hst['trig_eff'].GetXaxis().SetBinLabel(4, 'AK4 Mu DiJet')
    hst['trig_eff'].GetXaxis().SetBinLabel(5, 'AK4 Mu DiPFJet')
    hst['trig_eff'].GetXaxis().SetBinLabel(6, 'Mu12 DiBtag')
    hst['trig_eff'].GetXaxis().SetBinLabel(7, 'DiBtag dEta')
    hst['trig_eff'].GetXaxis().SetBinLabel(8, 'HT QuadJet')
    c0.SaveAs(png_dir+'h_trig_eff.png')

    hst['nCandsSel'].SetLineColor(R.kBlack)
    hst['nCandsSel'].Draw('hist')
    # hst['nCandsSV'].SetLineColor(R.kBlue)
    # hst['nCandsSV'].Draw('histsame')
    # hst['nCandsHT'].SetLineColor(R.kRed)
    # hst['nCandsHT'].Draw('histsame')
    # hst['nCandsSel'].GetYaxis().SetRangeUser(0.01, 200)
    # hst['nCandsSel'].GetYaxis().SetRangeUser(0, hst['nCandsSel'].GetBinContent(2)*1.2)
    # c0.SetLogy()
    c0.SaveAs(png_dir+'h_nCandsSel.png')
    c0.SetLogy(0)

    # ## Full event variables
    # hst['evt_HT'].SetLineColor(R.kRed)
    # hst['evt_HT'].Draw('hist')
    # hst['evt_MT'].SetLineColor(R.kBlue)
    # hst['evt_MT'].Draw('histsame')
    # hst['evt_mass'].SetLineColor(R.kBlack)
    # hst['evt_mass'].Draw('histsame')
    # c0.SaveAs(png_dir+'h_evt_HT_MT_mass.png')

    ## LHE HT for MC
    if isMC:
        hst['LHE_HT'].SetLineColor(R.kBlack)
        c0.SetLogy()
        hst['LHE_HT'].Draw('hist')
        c0.SaveAs(png_dir+'h_LHE_HT.png')
        c0.SetLogy(0)

        hst['LHE_HT_log'].SetLineColor(R.kBlack)
        hst['LHE_HT_log'].Draw('hist')
        c0.SaveAs(png_dir+'h_LHE_HT_log.png')

    ## Pileup distributions, overlaid
    hst['nPV_good'].SetLineColor(R.kBlue)
    hst['nPV_good'].Draw('hist')
    hst['nPV_all'].SetLineColor(R.kBlack)
    hst['nPV_all'].Draw('histsame')
    c0.SaveAs(png_dir+'h_nPV_all_good.png')

    ## Jet and muon multiplicities
    for var in ['nJet', 'nMuon']:
        hst[var].SetLineColor(R.kBlack)
        hst[var].Draw('hist')
        c0.SaveAs(png_dir+'h_%s.png' % var)

    ## Common kinematic variables for muons
    for var in ['pt', 'eta', 'IP', 'charge', 'set_charge']:
        hst['mu_%s' % var].SetLineColor(R.kBlack)
        hst['mu_%s' % var].Draw('hist')
        c0.SaveAs(png_dir+'h_mu_%s.png' % var)

    ## GEN matched b quarks for jets, overlaid
    hst['nGenB_Higg'].SetLineColor(R.kBlack)
    hst['nGenB_fatH'].SetLineColor(R.kBlue)
    hst['nGenB_jetH'].SetLineColor(R.kRed)
    hst['nGenB_Higg'].Draw('hist')
    hst['nGenB_fatH'].Draw('histsame')
    hst['nGenB_jetH'].Draw('histsame')
    c0.SaveAs(png_dir+'h_nGenB_Higg_fatH_jetH.png')

    ## Common kinematic variables for jets, overlaid
    for var in ['pt', 'eta', 'mass', 'msoft']:
        hst['%s_fatH' % var].SetLineColor(R.kBlue)
        hst['%s_Higg' % var].SetLineColor(R.kBlack)
        hst['%s_Higg' % var].Draw('hist')
        if CUTS.NUM_JETS_OUT > 0:
            hst['%s_fatH' % var].Draw('histsame')
            if var == 'pt' or var == 'eta':
                hst['%s_jetH' % var].Draw('histsame')
                hst['%s_ISR'  % var].SetLineColor(R.kMagenta)
                # hst['%s_ISR' % var].Draw('histsame')
            if var == 'mass' or var == 'msoft':
                hst['%s_jetH' % var].SetLineColor(R.kRed)
                hst['%s_fatA' % var].SetLineColor(R.kTeal)
                hst['%s_fatA' % var].Draw('histsame')
        c0.SaveAs(png_dir+'h_%s_fatH_jetH_Higg.png' % var)

    ## VBF system mass
    hst['mass_VBF'].SetLineColor(R.kBlack)
    hst['mass_VBF'].Draw('hist')
    c0.SaveAs(png_dir+'h_mass_VBF.png')

    ## Common separation variables for jets, overlaid
    for var in ['dR', 'dEta', 'dPhi', 'dR_mu']:
        hst['%s_ISR'  % var].SetLineColor(R.kMagenta)
        hst['%s_Higg' % var].SetLineColor(R.kRed)

        if var == 'dR_mu':
            hst['%s_Higg' % var].Draw('hist')
            hst['%s_ISR' % var].Draw('histsame')
        else:
            hst['%s_ISR' % var].Draw('hist')
            if CUTS.NUM_JETS_OUT > 0:
                hst['%s_Higg' % var].Draw('histsame')
        c0.SaveAs(png_dir+'h_%s_Higg_ISR.png' % var)

    ## ID variables for jets, overlaid
    for var in ['QGL', 'jetID', 'jetPU', 'nConst']:
        hst['%s_ISR' % var].SetLineColor(R.kMagenta)
        # hst['%s_ISR' % var].Draw('hist')
        if CUTS.NUM_JETS_OUT > 0:
            hst['%s_jetH' % var].SetLineColor(R.kRed)
            hst['%s_jetH' % var].Draw('hist')
            c0.SaveAs(png_dir+'h_%s_jetH_ISR.png' % var)

    ## b-tagging variables for jets, overlaid
    for var in ['deepB', 'flavB']:
        hst['%s_ISR'  % var].SetLineColor(R.kMagenta)
        hst['%s_jetH' % var].SetLineColor(R.kRed)
        if var == 'deepB':
            hst['%s_fatH' % var].SetLineColor(R.kBlue)
            hst['%s_fatA' % var].SetLineColor(R.kTeal)
        if var == 'deepB':
            hst['%s_fatH' % var].Draw('hist')
            # hst['%s_ISR'  % var].Draw('hist')
            if CUTS.NUM_JETS_OUT > 0:
                hst['%s_jetH' % var].Draw('histsame')
                hst['%s_fatA' % var].Draw('histsame')
            c0.SaveAs(png_dir+'h_%s_jetH_fatH_ISR.png' % var)
        if var == 'flavB':
            # hst['%s_ISR'  % var].Draw('hist')
            if CUTS.NUM_JETS_OUT > 0:
                hst['%s_jetH' % var].Draw('hist')
                c0.SaveAs(png_dir+'h_%s_jetH_ISR.png' % var)

    ## Tagging variables for fat jets, overlaid
    for var in ['doubB', 'fourQ']:
        hst['%s_fatH' % var].SetLineColor(R.kBlue)
        hst['%s_fatH' % var].Draw('hist')
        if CUTS.NUM_JETS_OUT > 0:
            hst['%s_fatA' % var].SetLineColor(R.kTeal)
            hst['%s_fatA' % var].Draw('histsame')
        c0.SaveAs(png_dir+'h_%s_fatH.png' % var)

    ## Other fat jet quantities
    for var in ['tau1', 'tau2', 'tau3', 'tau4', 'n2b1', 'n3b1']:
        hst['%s_fatH' % var].SetLineColor(R.kBlue)
        hst['%s_fatH' % var].Draw('hist')
        c0.SaveAs(png_dir+'h_%s_fatH.png' % var)

    ## Fat jet sub-jet quantities
    for sub in ['sub1', 'sub2']:
        for var in ['pt_frac', 'mass', 'deepB', 'n2b1', 'n3b1',
                    'tau1', 'tau2', 'tau3', 'tau4']:
            hst['%s_%s_fatH' % (sub, var)].SetLineColor(R.kBlue)
            hst['%s_%s_fatH' % (sub, var)].Draw('hist')
            c0.SaveAs(png_dir+'h_%s_%s_fatH.png' % (sub, var))


    ## Secondary vertices associated with jets
    hst['nSV_ISR'].SetLineColor(R.kMagenta)
    # hst['nSV_ISR'].Draw('hist')
    hst['nSV_extra'].SetLineColor(R.kGreen)
    hst['nSV_extra'].Draw('hist')
    if CUTS.NUM_JETS_OUT > 0:
        hst['nSV_jetH'].SetLineColor(R.kRed)
        hst['nSV_jetH'].Draw('histsame')
    hst['nSV_fatH'].SetLineColor(R.kBlue)
    hst['nSV_fatH'].Draw('histsame')
    hst['nSV_Higg'].SetLineColor(R.kBlack)
    hst['nSV_Higg'].Draw('histsame')
    c0.SaveAs(png_dir+'h_nSV_ISR_jetH_Higg_extra.png')

    for NN in NN_list:
        hst[NN].SetLineColor(R.kBlue)
        hst[NN].Draw('hist')
        c0.SaveAs(png_dir+'h_NN_%s.png' % NN)


    ## Write out all histograms, after format changes
    for key in hst.keys():
        hst[key].Write()

    ## Delete output file from local memory
    del out_file


## Define 'main' function as primary executable
if __name__ == '__main__':
    main()
