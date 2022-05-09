#define read_d4pd_cxx
//https://svnweb.cern.ch/trac/atlasinst/browser/Institutes/Pavia/twolep/trunk/charginob/fillObjs_d3pd.C
// The class definition in read_d4pd.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.
// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("read_d4pd.C")
// Root > T->Process("read_d4pd.C","some options")
// Root > T->Process("read_d4pd.C+")
 
#include "read_d4pd.h"
#include <TH2.h>
#include <TStyle.h> 
#include <iostream>
#include <fstream>
#include <TMath.h>
#include "myown_util/mt2_bisect.h"
#include "myown_util/mctlib.h"


void read_d4pd::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
} 

void read_d4pd::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  std::cout<<"Slave begin"<<std::endl;
  TString option = GetOption();
  
  ReadJobConfiguration("/Users/sdarmora/ST_3_20/run/config_info.dat"); 
  // >> dbg >> isdata >> ismc >> isfullsim >> which_syst;

  
  //  debug=false; //true;
  //m_isDATA=false;; //true;
  //m_isMC=true;
  //m_isFS=true;
  whichsyste = SystErr::NONE;
  whichmet=SUSYMet::Default; 
  jetTagger = SUSYBTagger::MV1;
  m_ismc12b=true;//change 
  m_istriglep=true;
  
  
  std::cout<< "configuration flags: "<< "debug" << '\t' << "m_isDATA" <<'\t'  << "m_isMC" <<'\t' << "m_isFASTSim"<< '\t' << "m_stream" << '\t' << "Ptcut" << '\t' << "CreateConfigFile pileup" <<endl;
  std::cout<< "configuration flags: "<< debug <<"               "<< m_isDATA <<"            "<< m_isMC <<"         "<< m_isFS<<"              "<< m_stream <<"             "<< m_isPtcut  <<"               " << CreateConfigFile <<'\n';
  
  
  
  
  BookHistos();
  susyobj.initialize(m_isDATA,m_isFS, m_ismc12b, m_istriglep); 
  susyobj.SetJetCalib(true);
  
  m_fakemet.initialize("fest_periodF_v1.root");
  // GRL initialization 
  m_eventInfo.SetGlobalEventInfo("/Users/sdarmora/ST_3_20/run/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.dat");
  
  
  // MC info initialization 
  // weight from cross section eff and number of events. // to be checked. 
  
    ReadWeights("/Users/sdarmora/ST_3_20/run/mc_info_nov_2014.dat"); 
 
  
    bool isJVF = true;
    BTagComputation = new BTagCalib("MV1","/Users/sdarmora/ST_3_20/run/files_for_run/BTagCalibration.env","/Users/sdarmora/ST_3_20/run/files_for_run/","0_7892",isJVF,0.7892);      
    
    
    
    doPUreweight=true;
    if (doPUreweight) {  
      //CreateConfigFile = false; //true; //false;
      //pileup   = new Root::TPileupReweighting("mytool");
      //pileup.Initialize();
      if(CreateConfigFile) pileup.UsePeriodConfig("MC12a"); //test per generare ConfigFile
      if(!CreateConfigFile)   pileup.AddConfigFile("/Users/sdarmora/ST_3_20/run/PileUpFiles/ConfFile_xxxyyy_p1328.root");
      pileup.SetDataScaleFactors(0.9);
      //pileup.SetDataScaleFactors(1);
      if(!CreateConfigFile)  pileup.AddLumiCalcFile("/Users/sdarmora/ST_3_20/run/PileUpFiles/ilumicalc_histograms_None_200842-215643.root");
      if(!CreateConfigFile)   pileup.SetUnrepresentedDataAction(2);
      pileup.initialize();
    }
    
    
    // ntuple for MVA
    
    theNtuple = new TNtuple("susy","susy","RunNumber:EventNumber:weight:met:ptJet1:ptJet2:ptEle1:ptEle2:ptMuo1:ptMuo2:InvM_ll:meff:met_over_meff:mT2:mctcorrll:bJet1:bJet2:DPhiJets:DThetaJets:DPhiJet1Met:DPhiJet2Met:DPhiLepton:DThetaLepton:DPhiLep1Met:DPhiLep2Met:DPhiL1J1:DPhiL1J2:DPhiL2J1:DPhiL2J2:DThetaL1J1:DThetaL1J2:DThetaL2J1:DThetaL2J2:Pbll:tightEle1:tightEle2:tightMuo1:tightMuo2:BWeight:DEtaLepton:had_mT2:DPhiMetPbll:DRll:btag1:btag2:M1_looseIso:M2_looseIso:E1_looseIso:E2_looseIso:M1_mediumIso:M2_mediumIso:E1_mediumIso:E2_mediumIso:M1_tightIso:M2_tightIso:E1_tightIso:E2_tightIso"); // for TMVA purposes       
    fOutput->Add(theNtuple);
    
    
    
    
    std::cout<<"Slave begin end" << std::endl;
}

Bool_t read_d4pd::Process(Long64_t entry)
{
  //std::cout<<"Process entry  "<<std::endl;
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either read_d4pd::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  
  
  float events;
 
  m_passGRL=false; 
  m_passTrigger_el= false;
  m_passTrigger_mu= false;
  m_passTrigger_met= false;
  m_GoodVertex= false; 
  m_GoodElectron= false; 
  m_GoodMuon= false; 
  m_GoodJet=false;
  m_BadJetEvent=false;
  m_BadJetEvent1=false; 
  //m_isMC=false;
  
  Jets.clear();
  Muons.clear(); 
  Electrons.clear(); 
  
  //Reset intern objects
  susyobj.Reset();
  
  float thr_el_pt=25000;  // MeV 
  float thr_el_pt_0=20000;  // MeV 
  float thr_mu_pt=18000;  // MeV 
  float thr_mu_pt_0=10000;  // MeV
  float thr_jet_pt=25000;// MeV 
  float thr_etmiss=40; // GeV 
  float eta_jet=2.5; 
  //std::cout<<m_isDATA<<"  "<<m_isMC<<endl;  
  
  w=1;
  double mcevtw=1;
  
  b_RunNumber->GetEntry(entry);  
  b_averageIntPerXing->GetEntry(entry);
  b_EventNumber->GetEntry(entry);

  Int_t runnumber=0; 
  if (m_isDATA ) {
    runnumber=RunNumber;
  } 
  else {
    b_mc_channel_number->GetEntry(entry);
    runnumber=mc_channel_number;
    b_mcevt_weight->GetEntry(entry);
    mcevtw=(mcevt_weight->at(0)).at(0);
   
    int size_weight=0;
    if(mcevt_weight){
      size_weight = mcevt_weight->at(0).size();
    }    
 
    if (doPUreweight) {
      if(CreateConfigFile && size_weight!=0) pileup.Fill(RunNumber,mc_channel_number,(mcevt_weight->at(0)).at(0),averageIntPerXing);
      if(CreateConfigFile) return kTRUE;
     
    }
  }


  b_isSimulation->GetEntry(entry);  

  //  std::cout <<"isSimulation: " << isSimulation <<'\n';
  //std::cout <<"RunNumber: " << RunNumber <<'\n';
  //std::cout <<"mc_cha_number: " <<   mc_channel_number <<'\n';
  //std::cout <<"isaMCRun(): " << isaMCRun(RunNumber) <<'\n';
  //std::cout<<"is_DATA"<<"  "<<m_isDATA<<endl;
  
  // isMC check or GRL flag  
  if (isaMCRun(runnumber) || isSimulation) { 
    m_isMC=true;
    if (isaMCRun(runnumber) != isSimulation)  std::cout <<"(isaMCRun(runnumber) != isSimulation) is list updated?  " <<endl; 
    
    double myPileupWeight=1.;
    
    if (doPUreweight) myPileupWeight = pileup.GetCombinedWeight(RunNumber,mc_channel_number,averageIntPerXing);
    if(mcevtw!=0 ){
      
      w=Weight(runnumber)*mcevtw*myPileupWeight;
     
      
    }
    
    
    if(debug) {
      if (w==0 || entry%10000==0)  std::cout << " RunNumber,mc_channel_number,averageIntPerXing "<< RunNumber <<" " <<mc_channel_number<<" "<<averageIntPerXing << " myPileupWeight " << myPileupWeight<<'\n';
    }
  }
   else{ 

     if (isaGoodLumiblock(entry)) m_passGRL=true;
       
    if (debug) std::cout << " Data.   GRL pass: " << m_passGRL <<endl;
  }
  
  
  
  //enum cutflow {all=0,grl,tri_el,tri_mu,gvtx}; 
  // ***all events 
  
  fillcutFlow(0.0); 
  h_totevent->Fill(runnumber); 
  
  
  //    if (!m_isMC && !m_passGRL) return  kTRUE; 
  // *** GRL 
  fillcutFlow(1.0); 

  //std::cout << "RunNumber: " << RunNumber << "Weight: " << Weight(RunNumber) << "entry: " <<entry <<'\n';
   
  
  // larError  or core flag 
  b_larError->GetEntry(entry); 
  b_coreFlags->GetEntry(entry);
  if(!m_isMC && (coreFlags&0x40000) != 0) return kTRUE; // This is an incomplete event remove from analysis }
 
  fillcutFlow(2.0);
  if (!m_isMC && larError!=0) return kTRUE; 
  fillcutFlow(3.0); 

  
    
  //  
  // -------------------Trigger 
  // 
  
  
     m_stream=4;

  bool trigger_pass=trigger(entry,m_stream);
  
  
  if (m_passTrigger_mu) {
    
    fillcutFlow(4.0); 
  }
  
  if (m_passTrigger_el) { 
    
    fillcutFlow(5.0); 
  } 
  

  if ( !trigger_pass ) return kTRUE; 
    fillcutFlow(6.0); 
    h_event_yield->Fill(runnumber); 
   
  // ------------------leptons selection   
 
  // -------------------muons --------------------------------------------------------
  b_mu_staco_pt->GetEntry(entry);
  b_mu_staco_eta->GetEntry(entry);
  b_mu_staco_phi->GetEntry(entry);
  b_mu_staco_E->GetEntry(entry); 
  b_mu_staco_charge->GetEntry(entry); 
  b_mu_staco_isCombinedMuon->GetEntry(entry); 
  b_mu_staco_z0_exPV->GetEntry(entry);
  b_mu_staco_d0_exPV->GetEntry(entry);
  b_mu_staco_ptcone20->GetEntry(entry); 
  //b_trig_L1_mu_eta->GetEntry(entry); 
  //b_trig_L1_mu_phi->GetEntry(entry); 
  b_mu_staco_id_theta_exPV->GetEntry(entry);
  b_mu_staco_me_theta_exPV->GetEntry(entry);
  b_mu_staco_id_theta->GetEntry(entry);
  b_mu_staco_isSegmentTaggedMuon->GetEntry(entry);
  b_mu_staco_loose->GetEntry(entry);
  b_mu_staco_expectBLayerHit->GetEntry(entry);
  b_mu_staco_nBLHits->GetEntry(entry);
  b_mu_staco_nPixHits->GetEntry(entry);
  b_mu_staco_nPixelDeadSensors->GetEntry(entry);
  b_mu_staco_nPixHoles->GetEntry(entry);
  b_mu_staco_nSCTHits->GetEntry(entry);
  b_mu_staco_nSCTDeadSensors->GetEntry(entry);
  b_mu_staco_nSCTHoles->GetEntry(entry);
  b_mu_staco_nTRTHits->GetEntry(entry);
  b_mu_staco_nTRTOutliers->GetEntry(entry);
  // b_mu_staco_cov_qoverp_exPV->GetEntry(entry);
  // b_mu_staco_qoverp_exPV->GetEntry(entry);
  b_mu_staco_me_qoverp_exPV->GetEntry(entry); 
  b_mu_staco_id_qoverp_exPV->GetEntry(entry); 
  b_trig_EF_trigmuonef_track_n->GetEntry(entry);
  
  //b_trig_EF_trigmuonef_track_CB_pt->GetEntry(entry);   
  b_trig_EF_trigmuonef_track_CB_eta->GetEntry(entry);   
  b_trig_EF_trigmuonef_track_CB_phi->GetEntry(entry);   
  b_trig_EF_trigmuonef_track_CB_hasCB->GetEntry(entry);   
  //b_trig_EF_trigmuonef_EF_mu24->GetEntry(entry);   
  //b_mu_staco_MET_Egamma10NoTau_STVF_wet->GetEntry(entry);
  b_mu_staco_trackIPEstimate_d0_unbiasedpvunbiased->GetEntry(entry);
  b_mu_staco_trackIPEstimate_z0_unbiasedpvunbiased->GetEntry(entry);
  b_mu_staco_trackIPEstimate_sigd0_unbiasedpvunbiased->GetEntry(entry);
  b_mu_staco_ptcone30_trkelstyle->GetEntry(entry);
  b_vx_nTracks->GetEntry(entry);
  b_mu_staco_etcone30->GetEntry(entry);


  bool muon_trigger=false;  
  vector<int> v_mu_idx;  //vectors of muon indexes  used by MET re-computation macro
  vector<int> mu_id;
  
  for (unsigned int iMu=0; iMu<mu_staco_pt->size(); iMu++) {
    bool is_cosmic=false;
    bool muon_signal=false;
    bool muon_mediumIso=false;
    bool muon_looseIso=false;
    bool muon_tightIso=false;
    bool muon_baseline=false;
 
    // all  quality definition  deemed to SUSYTOOLS  
    bool muon_id = susyobj.FillMuon(iMu, mu_staco_pt->at(iMu), mu_staco_eta->at(iMu), mu_staco_phi->at(iMu),
				    mu_staco_me_qoverp_exPV->at(iMu), mu_staco_id_qoverp_exPV->at(iMu),
				    mu_staco_me_theta_exPV->at(iMu), mu_staco_id_theta_exPV->at(iMu), 
				    mu_staco_id_theta->at(iMu), mu_staco_charge->at(iMu),
				    mu_staco_isCombinedMuon->at(iMu), mu_staco_isSegmentTaggedMuon->at(iMu),
				    mu_staco_loose->at(iMu), 
				    //mu_staco_expectBLayerHit->at(iMu), mu_staco_nBLHits->at(iMu),
				    mu_staco_nPixHits->at(iMu), mu_staco_nPixelDeadSensors->at(iMu),
				    mu_staco_nPixHoles->at(iMu), mu_staco_nSCTHits->at(iMu),
				    mu_staco_nSCTDeadSensors->at(iMu),mu_staco_nSCTHoles->at(iMu),
				    mu_staco_nTRTHits->at(iMu), mu_staco_nTRTOutliers->at(iMu),
				    //				    m_isDATA,
				    6000.,2.5);
    
    Int_t classific=1000; /// bit mask:  SIGNALMUON  baseline BADMUON COSMICMUON TRIGGERMUON
    if (muon_id )  { // baseline  muon definition 
      
      v_mu_idx.push_back(iMu); //? or *mu_staco_MET_Egamma10NoTau_STVF_wet;  //  to be used for METUtils 
      //      Int_t classific=1000; /// bit mask:  SIGNALMUON  baseline BADMUON COSMICMUON TRIGGERMUON  
      if (susyobj.IsSignalMuon(iMu, 
			       mu_staco_ptcone20->at(iMu), 
			       10000.,
			       1800.) ) { classific+=10000;muon_signal=true;}
      
      if(susyobj.IsSignalMuonExp(iMu,
                                 vx_nTracks,
                                 mu_staco_ptcone30_trkelstyle->at(iMu),
                                 mu_staco_etcone30->at(iMu),
                                 mu_staco_trackIPEstimate_d0_unbiasedpvunbiased->at(iMu),
                                 mu_staco_trackIPEstimate_z0_unbiasedpvunbiased->at(iMu),
                                 mu_staco_trackIPEstimate_sigd0_unbiasedpvunbiased->at(iMu),
                                 SignalIsoExp::LooseIso )) muon_looseIso=true;
      
      
      if ( susyobj.IsSignalMuonExp(iMu,
				   vx_nTracks,
				   mu_staco_ptcone30_trkelstyle->at(iMu),
				   mu_staco_etcone30->at(iMu),
				   mu_staco_trackIPEstimate_d0_unbiasedpvunbiased->at(iMu),
				   mu_staco_trackIPEstimate_z0_unbiasedpvunbiased->at(iMu),
				   mu_staco_trackIPEstimate_sigd0_unbiasedpvunbiased->at(iMu),
				   SignalIsoExp::MediumIso ) ) muon_mediumIso=true;
      
      if ( susyobj.IsSignalMuonExp(iMu,
                                   vx_nTracks,
                                   mu_staco_ptcone30_trkelstyle->at(iMu),
                                   mu_staco_etcone30->at(iMu),
                                   mu_staco_trackIPEstimate_d0_unbiasedpvunbiased->at(iMu),
                                   mu_staco_trackIPEstimate_z0_unbiasedpvunbiased->at(iMu),
                                   mu_staco_trackIPEstimate_sigd0_unbiasedpvunbiased->at(iMu),
                                   SignalIsoExp::TightIso ) ) muon_tightIso=true;
      
      
      
      if (susyobj.IsCosmicMuon(mu_staco_z0_exPV->at(iMu), mu_staco_d0_exPV->at(iMu),1.0,0.2) ) { is_cosmic=true; classific+=10;}
      
      //const vector<int> * trig_EF_trigmuonef_signature; 
      int EFindex(-1); 
      int EFtrackindex(-1);
          
      Muon mu;
	mu.SetPtEtaPhiE(susyobj.GetMuonTLV(iMu).Pt(),  susyobj.GetMuonTLV(iMu).Eta(), susyobj.GetMuonTLV(iMu).Phi(), susyobj.GetMuonTLV(iMu).E());//Pt/E/1000                                                                      

      mu.isSignal= muon_signal;
      mu.isLooseIso=muon_looseIso;
      mu.isMediumIso=muon_mediumIso;
      mu.isTightIso=muon_tightIso;
      mu.isCosmic=is_cosmic;  
      mu.isBaseline=muon_id;
      mu.charge=mu_staco_charge->at(iMu);
      mu.truthIndex=iMu;
      mu.classification=classific;    
      Muons.push_back(mu);  
      m_GoodMuon=true; // one good baseline muon 
      
      
    }//closes for muon_id
  }
  
  
  if (Muons.size()>0 )  std::sort(Muons.begin(), Muons.end());

 
    
  if (m_GoodMuon) {
    fillcutFlow(9.0); 
    // h_mu_n->Fill(Muons.size(),w); 
  }
  
 
  //--------------------------------------- electrons 
  b_el_author->GetEntry(entry);
  b_el_ptcone20->GetEntry(entry);
  b_el_ptcone30->GetEntry(entry);
  b_el_cl_eta->GetEntry(entry);
  b_el_eta->GetEntry(entry);
  b_el_phi->GetEntry(entry);
  b_el_charge->GetEntry(entry);
  b_trig_EF_el_eta->GetEntry(entry); 
  b_trig_EF_el_phi->GetEntry(entry); 
  // new for SO 
  b_el_mediumPP->GetEntry(entry); 
  b_el_tightPP->GetEntry(entry); 
  b_el_cl_phi->GetEntry(entry);
  b_el_cl_E->GetEntry(entry);
  b_el_tracketa->GetEntry(entry);
  b_el_trackphi->GetEntry(entry);
  b_el_OQ->GetEntry(entry);
  b_el_nPixHits->GetEntry(entry); 
  b_el_nSCTHits->GetEntry(entry); 
  b_el_MET_Egamma10NoTau_wet->GetEntry(entry);  
  b_el_trackd0pv->GetEntry(entry);
  b_el_trackz0pv->GetEntry(entry);
  
  // isMediumPP and isTightPP
  // b_el_etas2->GetEntry(entry);
  // b_el_f3->GetEntry(entry);
  //b_el_Ethad->GetEntry(entry);
  //b_el_Ethad1->GetEntry(entry);
  //b_el_reta->GetEntry(entry);
  //b_el_weta2->GetEntry(entry);
  //b_el_f1->GetEntry(entry);
  //b_el_wstot->GetEntry(entry);
  //b_el_emaxs1->GetEntry(entry);
  //b_el_Emax2->GetEntry(entry);
  //b_el_deltaeta1->GetEntry(entry);
  //b_el_TRTHighTOutliersRatio->GetEntry(entry);
  //b_el_nTRTOutliers->GetEntry(entry);
  //b_el_nTRTHits->GetEntry(entry);
  //b_el_nSiHits->GetEntry(entry);
  //b_el_nSCTOutliers->GetEntry(entry);
  //b_el_nPixelOutliers->GetEntry(entry);
  //b_el_nPixHits->GetEntry(entry);
  //b_el_nBLHits->GetEntry(entry);
  //b_el_nBLayerOutliers->GetEntry(entry);
  b_el_expectHitInBLayer->GetEntry(entry);
  b_el_trackd0_physics->GetEntry(entry);
  b_el_deltaphi2->GetEntry(entry);
  b_el_trackqoverp->GetEntry(entry);
  //b_el_isEM->GetEntry(entry);
  b_el_topoEtcone30_corrected->GetEntry(entry);
  b_el_trackIPEstimate_d0_unbiasedpvunbiased->GetEntry(entry);
  b_el_trackIPEstimate_z0_unbiasedpvunbiased->GetEntry(entry);
  b_el_trackIPEstimate_sigd0_unbiasedpvunbiased->GetEntry(entry);


  double etcut =7000.;
  // //Fill the electrons 
  bool electron_trigger=false;
  vector<int> v_el_idx;  //vectors of electrons index by MET re-computation macro
  
  for (unsigned int iEl=0; iEl<el_author->size(); iEl++) {
    bool electron_signal=false;
    
    bool electron_mediumIso=false;
    bool electron_looseIso=false;
    bool electron_tightIso=false;
    
    bool   ismediumPP= el_mediumPP->at(iEl); 
    bool  istightPP = el_tightPP->at(iEl);
    
    // all definitions  are in  SUSYTools
    bool electron_id = susyobj.FillElectron(iEl,
					    el_eta->at(iEl),
					    el_phi->at(iEl),					    
					    el_cl_eta->at(iEl),
					    el_cl_phi->at(iEl),
					    el_cl_E->at(iEl),
					    el_tracketa->at(iEl),
					    el_trackphi->at(iEl),
					    el_author->at(iEl),
					    ismediumPP,
					    //el_mediumPP->at(iEl),
					    el_OQ->at(iEl),
					    el_nPixHits->at(iEl),
					    el_nSCTHits->at(iEl), 
					    el_MET_Egamma10NoTau_wet->at(iEl).at(0), 
					    // m_isDATA,
					    etcut,
					    2.47,
					    whichsyste);
    
    //if ( ismediumPP!= el_mediumPP->at(iEl) )  std::cout <<" medium/tight " << ismediumPP <<" : " <<  el_mediumPP->at(iEl)<<" : "  << istightPP <<" : " <<  el_tightPP->at(iEl) <<'\n';
    if (electron_id) {
      if ( (el_MET_Egamma10NoTau_wet->at(iEl)).at(0) != 0 ) v_el_idx.push_back(iEl); // used  MET re-computation macro
      Int_t classific=1000; /// bit mask: BBBBB  SIGNAL  BASELINE BAD(CRACK) COSMIC(EMPTY) TRIGGER  
      if (susyobj.IsSignalElectron(iEl,
				   istightPP,
				   //el_tightPP->at(iEl),
				   el_ptcone20->at(iEl),
				   el_trackd0pv->at(iEl),
				   el_trackz0pv->at(iEl),10000,.1,1,2) && !susyobj.IsInCrack(el_cl_eta->at(iEl))) { classific+=10000; electron_signal=true;}
      // el_trackz0pv->at(iEl),10000,.1,1,2)) { classific+=10000; electron_signal=true;}
      if (susyobj.IsInCrack(el_cl_eta->at(iEl)) ) {  classific+=100;}
      
            
      if(susyobj.IsSignalElectronExp(iEl,                                                                     
				     //el_tightPP->at(iEl),                                                    
				     istightPP,
				     vx_nTracks,                                                              
				     el_ptcone30->at(iEl),                                                    
				     el_topoEtcone30_corrected->at(iEl),                                      
				     el_trackIPEstimate_d0_unbiasedpvunbiased->at(iEl),                       
				     el_trackIPEstimate_z0_unbiasedpvunbiased->at(iEl),                       
				     el_trackIPEstimate_sigd0_unbiasedpvunbiased->at(iEl),                    
				     SignalIsoExp::LooseIso)) {electron_looseIso=true;}                     
      
	
      if(susyobj.IsSignalElectronExp(iEl,                                                                     
                                     //el_tightPP->at(iEl),                                                    
				     istightPP,
                                     vx_nTracks,                                                              
                                     el_ptcone30->at(iEl),                                                    
                                     el_topoEtcone30_corrected->at(iEl),                                      
                                     el_trackIPEstimate_d0_unbiasedpvunbiased->at(iEl),   
                                     el_trackIPEstimate_z0_unbiasedpvunbiased->at(iEl),                       
                                     el_trackIPEstimate_sigd0_unbiasedpvunbiased->at(iEl),                    
                                     SignalIsoExp::MediumIso)) {electron_mediumIso=true;}                       
      
      if(susyobj.IsSignalElectronExp(iEl,
                                     //el_tightPP->at(iEl),
				     istightPP,
                                     vx_nTracks,
                                     el_ptcone30->at(iEl),
                                     el_topoEtcone30_corrected->at(iEl),
                                     el_trackIPEstimate_d0_unbiasedpvunbiased->at(iEl),
                                     el_trackIPEstimate_z0_unbiasedpvunbiased->at(iEl),
                                     el_trackIPEstimate_sigd0_unbiasedpvunbiased->at(iEl),
                                     SignalIsoExp::TightIso)) {electron_tightIso=true;}
      
      
      // trigger matching 
      for (unsigned int j=0; j<trig_EF_el_eta->size(); j++) {
	TLorentzVector trig_el; 
	TLorentzVector el=susyobj.GetElecTLV(iEl);
	trig_el.SetPtEtaPhiE(1.0,trig_EF_el_eta->at(j),trig_EF_el_phi->at(j),1.0); 
	float dr=trig_el.DeltaR(el); 
	if (iEl<3) h_el_trig_dr[iEl]->Fill(dr,w); 
	if (dr<0.1)  { electron_trigger=true; classific+=1;}
      } 
      
      
      Electron ele; 
        ele.SetPtEtaPhiE(susyobj.GetElecTLV(iEl).Pt(), susyobj.GetElecTLV(iEl).Eta(), susyobj.GetElecTLV(iEl).Phi(), susyobj.GetElecTLV(iEl).E() );//E,Pt/1000                                                                     
  
      ele.isLooseIso=electron_looseIso;     
      ele.isMediumIso=electron_mediumIso;
      ele.isTightIso=electron_tightIso;
      ele.isSignal=electron_signal;
      ele.isBaseline=electron_id;
      ele.isTight=istightPP;
      ele.charge=el_charge->at(iEl);
      ele.truthIndex=iEl; 
      ele.classification=classific;
      Electrons.push_back(ele);
      m_GoodElectron=true;
    }
    
  }
  
  if ( Electrons.size()>0) std::sort(Electrons.begin(), Electrons.end()); 
  
 
  
  if (m_GoodElectron) { 
   fillcutFlow(10.0); 
   //h_el_n->Fill(Electrons.size(),w); 
  }



  
  // // *** one ele or one muon matched with trigger 
  // if (!(m_GoodMuon || m_GoodElectron)) return kTRUE;  
  // fillcutFlow(11.0); 
  if (!(m_GoodMuon || m_GoodElectron)){
    fillcutFlow(11.0);
      }
  //std::cout << "runnumber: " << RunNumber << "Weight: " << Weight(RunNumber) <<'\n';
  //return kTRUE; 
  
  //
  // Jets 
  //
  b_jet_AntiKt4LCTopo_E->GetEntry(entry);
  b_jet_AntiKt4LCTopo_pt->GetEntry(entry);
  b_jet_AntiKt4LCTopo_eta->GetEntry(entry);
  b_jet_AntiKt4LCTopo_phi->GetEntry(entry);
  //b_jet_AntiKt4LCTopo_n90->GetEntry(entry);
  b_jet_AntiKt4LCTopo_Timing->GetEntry(entry);
  b_jet_AntiKt4LCTopo_LArQuality->GetEntry(entry);
  b_jet_AntiKt4LCTopo_sumPtTrk->GetEntry(entry);
  b_jet_AntiKt4LCTopo_HECQuality->GetEntry(entry);
  b_jet_AntiKt4LCTopo_NegativeE->GetEntry(entry);
  b_jet_AntiKt4LCTopo_fracSamplingMax->GetEntry(entry);
  b_jet_AntiKt4LCTopo_hecf->GetEntry(entry);
  b_jet_AntiKt4LCTopo_emfrac->GetEntry(entry);
  b_jet_AntiKt4LCTopo_AverageLArQF->GetEntry(entry);
  b_jet_AntiKt4LCTopo_SamplingMax->GetEntry(entry);
  b_jet_AntiKt4LCTopo_NegativeE->GetEntry(entry);
  //b_jet_AntiKt4LCTopo_emscale_E->GetEntry(entry);
  //b_jet_AntiKt4LCTopo_emscale_pt->GetEntry(entry);
  //b_jet_AntiKt4LCTopo_emscale_eta->GetEntry(entry);
  //b_jet_AntiKt4LCTopo_emscale_phi->GetEntry(entry);
  b_jet_AntiKt4LCTopo_jvtxf->GetEntry(entry);
  // b_jet_AntiKt4LCTopo_flavor_weight_TrackCounting2D->GetEntry(entry);  
  //b_jet_AntiKt4LCTopo_flavor_weight_JetProb->GetEntry(entry);
  b_jet_AntiKt4LCTopo_flavor_weight_SV0->GetEntry(entry);
  b_jet_AntiKt4LCTopo_flavor_weight_MV1->GetEntry(entry);
  // b_averageIntPerXing->GetEntry(entry);
  // new
  b_jet_AntiKt4LCTopo_BCH_CORR_DOTX->GetEntry(entry); 
  b_jet_AntiKt4LCTopo_BCH_CORR_CELL->GetEntry(entry);
  b_jet_AntiKt4LCTopo_BCH_CORR_JET->GetEntry(entry);
  if (!m_isDATA)  b_jet_AntiKt4LCTopo_flavor_truth_label->GetEntry(entry);
  // b_jet_AntiKt4LCTopo_EtaOrigin->GetEntry(entry);
  //b_jet_AntiKt4LCTopo_PhiOrigin->GetEntry(entry);
  //b_jet_AntiKt4LCTopo_MOrigin->GetEntry(entry);
  b_jet_AntiKt4LCTopo_constscale_eta->GetEntry(entry);
  b_jet_AntiKt4LCTopo_constscale_phi->GetEntry(entry);
  b_jet_AntiKt4LCTopo_constscale_E->GetEntry(entry);
  b_jet_AntiKt4LCTopo_constscale_m->GetEntry(entry);
  b_jet_AntiKt4LCTopo_ActiveAreaPx->GetEntry(entry);
  b_jet_AntiKt4LCTopo_ActiveAreaPy->GetEntry(entry);
  b_jet_AntiKt4LCTopo_ActiveAreaPz->GetEntry(entry);
  b_jet_AntiKt4LCTopo_ActiveAreaE->GetEntry(entry);
  b_Eventshape_rhoKt4LC->GetEntry(entry);
  //  b_vx_nTracks->GetEntry(entry);
  b_jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN->GetEntry(entry);
    
  /*
  for(int iJet=0; iJet<(int)jet_AntiKt4LCTopo_n; iJet++){
    susyobj.SetJetTLV(iJet,jet_AntiKt4LCTopo_pt->at(iJet),jet_AntiKt4LCTopo_eta->at(iJet), jet_AntiKt4LCTopo_phi->at(iJet), jet_AntiKt4LCTopo_E->at(iJet));
  }
  */

 

  for (unsigned int iJet=0; iJet<jet_AntiKt4LCTopo_pt->size(); iJet++) {
    
    //--- flag as good the jets that satisfy :
    // - Pt>20GeV
    // - |eta|<10
    // - passes JetID::VeryLooseBad      

    int local_truth_label=0; 
    if (!m_isDATA)  local_truth_label= jet_AntiKt4LCTopo_flavor_truth_label->at(iJet);
    bool is_jet = false;
    if ( is_jet=susyobj.FillJet(iJet,jet_AntiKt4LCTopo_pt->at(iJet),jet_AntiKt4LCTopo_eta->at(iJet),
			 //    bool  jet_all = susyobj.FillJet(iJet,jet_AntiKt4LC1Topo_pt->at(iJet),jet_AntiKt4LCTopo_eta->at(iJet), 
			 jet_AntiKt4LCTopo_phi->at(iJet), jet_AntiKt4LCTopo_E->at(iJet),
			 jet_AntiKt4LCTopo_constscale_eta->at(iJet), jet_AntiKt4LCTopo_constscale_phi->at(iJet),
			 jet_AntiKt4LCTopo_constscale_E->at(iJet), jet_AntiKt4LCTopo_constscale_m->at(iJet),
			 jet_AntiKt4LCTopo_ActiveAreaPx->at(iJet), jet_AntiKt4LCTopo_ActiveAreaPy->at(iJet),
			 jet_AntiKt4LCTopo_ActiveAreaPz->at(iJet), jet_AntiKt4LCTopo_ActiveAreaE->at(iJet),
			 Eventshape_rhoKt4LC,
			 averageIntPerXing,
				vx_nTracks));
  
    {
      bool Jet_cal= susyobj.ApplyJetSystematics(iJet, 

					       jet_AntiKt4LCTopo_constscale_eta->at(iJet),
					       local_truth_label, 
					       averageIntPerXing, 
					       vx_nTracks, 
					      	whichsyste);
      
    
      bool jet_id=susyobj.IsGoodJet(iJet, 
				    jet_AntiKt4LCTopo_constscale_eta->at(iJet),
				    jet_AntiKt4LCTopo_emfrac->at(iJet), 
				    jet_AntiKt4LCTopo_hecf->at(iJet), 
				    jet_AntiKt4LCTopo_LArQuality->at(iJet),
				    jet_AntiKt4LCTopo_HECQuality->at(iJet), 
				    jet_AntiKt4LCTopo_AverageLArQF->at(iJet), 
				    jet_AntiKt4LCTopo_Timing->at(iJet), 
				    jet_AntiKt4LCTopo_sumPtTrk->at(iJet), 
				    jet_AntiKt4LCTopo_fracSamplingMax->at(iJet),
				    jet_AntiKt4LCTopo_SamplingMax->at(iJet), 
				    jet_AntiKt4LCTopo_NegativeE->at(iJet), 
				    RunNumber,        
				    20000.,
				    10,
				    JetID::VeryLooseBad);
      
      
      Int_t classific=1000;  // SIGNAL BASELINE BAD (LARHOLE) 
      
            bool jet_signal=false;
	    bool jet_bad=false;
	    if(!jet_id&&susyobj.GetJetTLV(iJet).Pt() >20000)jet_bad=true;
	    //if (TMath::Abs(susyobj.GetJetTLV(iJet).Eta()) < 2.5) jet_signal=true;     
	     //        if (susyobj.GetJetTLV(iJet).Pt()<50000. &&  TMath::Abs(susyobj.GetJetTLV(iJet).Eta()) <2.4  &&  jet_AntiKt4LCTopo_jvtxf->at(iJet)>0.5) jet_signal=true; 
	
	    if (TMath::Abs(susyobj.GetJetTLV(iJet).Eta()) < 2.5){
	    if(susyobj.GetJetTLV(iJet).Pt()>50000. || (TMath::Abs(susyobj.GetJetTLV(iJet).Eta()) > 2.4 && TMath::Abs(susyobj.GetJetTLV(iJet).Eta()) < 2.5 )|| TMath::Abs(susyobj.GetJetTLV(iJet).Eta()) <2.4  &&  jet_AntiKt4LCTopo_jvtxf->at(iJet)>0.5) jet_signal=true;
	    }
	    
	  
	      Jet je; 
	      je.SetPtEtaPhiE(susyobj.GetJetTLV(iJet).Pt(),susyobj.GetJetTLV(iJet).Eta(),susyobj.GetJetTLV(iJet).Phi(),susyobj.GetJetTLV(iJet).E());//Pt/E/1000
	      je.isBaseline=  Jet_cal &&  susyobj.GetJetTLV(iJet).Pt() >20000 && TMath::Abs(susyobj.GetJetTLV(iJet).Eta()) < 2.8;
	      je.isSignal=jet_id &&Jet_cal && jet_signal; 
	      // je.isBad=!jet_id;
	      je.isBad=!jet_id&&susyobj.GetJetTLV(iJet).Pt() >20000;
	      // je.isBad=jet_bad;
	      je.truthIndex=iJet; // index of object in the original tree
	      je.isBJet= susyobj.IsBJet(jet_AntiKt4LCTopo_flavor_weight_MV1->at(iJet),0.7892) && susyobj.GetJetTLV(iJet).Pt() >20000 && TMath::Abs(susyobj.GetJetTLV(iJet).Eta()) < 2.5; 
	      je.classification=classific;
	      if(susyobj.GetJetTLV(iJet).Pt()>20000.)Jets.push_back(je);   
	      //  Jets.push_back(je);
	      m_GoodJet=true;  
    }
  } // loop 
 
 
  if (Jets.size()>0) std::sort(Jets.begin(), Jets.end());

  
  
  //std::cout << "runnumber: " << RunNumber << "Weight: " << Weight(RunNumber) <<'\n';

 

 


  //------------- resolve overlaps  between dif. particles  ------------------------- 
  // overlap flags: 
  // -1 no overlap
  // 0  it overlap with one other objet  but it is keep as the good.  
  // 1 it will  be not used   
  // "ele"  from a muon  
  // 2  "ele"  from an ele
  // 3 "jet" from an ele. 
  //    "jet from a muon ?? 
  //---------------------------------------------------------
      
  
 
 
  int ms=Muons.size();
  int es=Electrons.size();
  int js=Jets.size();
  
  if (debug) std::cout << " * resolve Overlap  : muon,electron,jet sizes: " <<ms <<"; "<<es<< "; " <<js<<'\n';  

    
  
   for (unsigned int i=0;i<(int)Electrons.size(); i++) {
    if (!(Electrons.at(i)).isBaseline) continue;
    for (unsigned int j=0; j<(int)Jets.size() ; ) {
      // if (!(Jets.at(j)).isBaseline) continue;
      float dr= (Jets.at(j)).DeltaR(Electrons.at(i));
      if (i<3) h_el_jet_dr[i]->Fill(dr,w);
      if ( dr< 0.2 ) {
	Jets.erase(Jets.begin()+j);
	//(Jets.at(j)).overlap=1;
	(Electrons.at(i)).overlap=0;
      } else {
	j++;
      }
    }   
   }
  


  
   
   for (unsigned int j=0;j<Jets.size(); j++) {
    if (!(Jets.at(j)).isBaseline) continue;
    for (unsigned int i=0;i<Electrons.size(); ) {
      if (!(Electrons.at(i)).isBaseline) continue;
      float dr= (Electrons.at(i)).DeltaR(Jets.at(j));
      if (i<3) h_el_el_dr[i]->Fill(dr,w);
      if ( dr< 0.4 ) {
	Electrons.erase(Electrons.begin()+i);
	(Jets.at(j)).overlap=0;
      } else {
	i++;
      }
    }
  }
  
     
  
  for (unsigned int j=0;j<Jets.size(); j++) {
    if (!(Jets.at(j)).isBaseline) continue;
    for (unsigned int i=0;i<Muons.size(); ) {
      if (!(Muons.at(i)).isBaseline) continue;
      float dr= (Muons.at(i)).DeltaR(Jets.at(j));
      if (i<3) h_el_mu_dr[i]->Fill(dr,w);
      if ( dr< 0.4  ) {
	Muons.erase(Muons.begin()+i);
	(Jets.at(j)).overlap=0;
      } else {
	i++;
      }
    }  
  }
  




  if (debug && (Muons.size()!=ms || Electrons.size()!=es || Jets.size()!=js))  {
    std::cout << " Overlap  removed   : muon,electron,jet sizes: "<< ms <<":"<<es<<":"<<js  <<" | " 
	      << Muons.size()  <<"; "<<Electrons.size()<< "; " <<Jets.size()<<'\n';   
  }
  

  
  //   std::cout << "RunNumber: " << RunNumber << "Weight: " << Weight(RunNumber) << "entry: " <<entry <<'\n';
  //  return kTRUE;
  // 
  //----------- recompute MET -------------------------------------------- 
  //
  b_el_MET_Egamma10NoTau_wpx->GetEntry(entry);
  b_el_MET_Egamma10NoTau_wpy->GetEntry(entry);
  //b_el_MET_Egamma10NoTau_wet->GetEntry(entry);
  b_el_MET_Egamma10NoTau_statusWord->GetEntry(entry);
  b_jet_AntiKt4LCTopo_MET_Egamma10NoTau_wet->GetEntry(entry);  
  b_jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpx->GetEntry(entry);  
  b_jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpy->GetEntry(entry);  
  b_jet_AntiKt4LCTopo_MET_Egamma10NoTau_statusWord->GetEntry(entry);  
  b_MET_Egamma10NoTau_CellOut_etx->GetEntry(entry); 
  b_MET_Egamma10NoTau_CellOut_ety->GetEntry(entry);
  b_MET_Egamma10NoTau_CellOut_sumet->GetEntry(entry); 
  b_MET_Egamma10NoTau_CellOut_Eflow_STVF_etx->GetEntry(entry); 
  b_MET_Egamma10NoTau_CellOut_Eflow_STVF_ety->GetEntry(entry); 
  b_MET_Egamma10NoTau_CellOut_Eflow_STVF_sumet->GetEntry(entry); 
  b_MET_Egamma10NoTau_RefGamma_etx->GetEntry(entry); 
  b_MET_Egamma10NoTau_RefGamma_ety->GetEntry(entry); 
  b_MET_Egamma10NoTau_RefGamma_sumet->GetEntry(entry); 
  b_mu_staco_ms_qoverp->GetEntry(entry); 
  b_mu_staco_ms_theta->GetEntry(entry); 
  b_mu_staco_ms_phi->GetEntry(entry); 
  //  b_mu_staco_charge->GetEntry(entry);
 
  
  vector<int> muon_index_noiso;
  vector<int> elec_index_noiso;
 
  //for(int i=0; i<(int)el_n; i++){
  // if((el_MET_Egamma10NoTau_wet->at(i)).at(0)!=0)  elec_index_noiso.push_back(i);
  //}
  //  for(int i=0; i<(int)mu_staco_n; i++){
  //if((Muons.at(i)).isBaseline)  muon_index_noiso.push_back(i);
  //}

  bool MetFix = false;
  bool doMuonElossCorrection = false;
  
  std::vector<float>* mu_staco_energyLossPar=0;
    
  TVector2 met_vec; 
  met_vec = susyobj.GetMET(
			   //jet_AntiKt4LCTopo_pt,
			   jet_AntiKt4LCTopo_MET_Egamma10NoTau_wet,
			   jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpx,
			   jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpy,
			   jet_AntiKt4LCTopo_MET_Egamma10NoTau_statusWord,
			   v_el_idx,
			   //elec_index_noiso,
			   el_MET_Egamma10NoTau_wet,
			   el_MET_Egamma10NoTau_wpx,
			   el_MET_Egamma10NoTau_wpy,
			   el_MET_Egamma10NoTau_statusWord,
			   MET_Egamma10NoTau_CellOut_etx, 
			   MET_Egamma10NoTau_CellOut_ety, 
			   MET_Egamma10NoTau_CellOut_sumet,
			   MET_Egamma10NoTau_CellOut_Eflow_STVF_etx, 
			   MET_Egamma10NoTau_CellOut_Eflow_STVF_ety, 
			   MET_Egamma10NoTau_CellOut_Eflow_STVF_sumet,
			   MET_Egamma10NoTau_RefGamma_etx,
			   MET_Egamma10NoTau_RefGamma_ety,
			   MET_Egamma10NoTau_RefGamma_sumet,
			   v_mu_idx,
			   //muon_index_noiso,
			   mu_staco_ms_qoverp, 
			   mu_staco_ms_theta, 
			   mu_staco_ms_phi, 
			   mu_staco_charge,
			   mu_staco_energyLossPar,
			   averageIntPerXing,
			   whichmet, 
			   whichsyste, 
			   doMuonElossCorrection,
			   MetFix);
  // false,
  //		   false);
  
  
  
  h_etmiss_mod->Fill( met_vec.Mod()/1000,w);
  h_etmiss_phi->Fill(met_vec.Phi_mpi_pi(met_vec.Phi()),w);

  
  double NewMET_x;
  double NewMET_y;
  double NewMET_phi;

  NewMET_x = met_vec.X();
  NewMET_x = met_vec.Y();
  NewMET_phi = met_vec.Phi();
 
  float NewMET = sqrt(pow(NewMET_x,2)+pow(NewMET_y,2));
  
  TLorentzVector etm;

  //    etm.SetPtEtaPhiM(NewMET/1000,0.,NewMET_phi,0.);
  //  etm.SetPtEtaPhiE( met_vec.Mod()/1000,0.0, met_vec.Phi_mpi_pi(met_vec.Phi()), met_vec.Mod()/1000);
  etm.SetPtEtaPhiE( met_vec.Mod(),0.0, met_vec.Phi(), met_vec.Mod());
 
 
 
  //
  // SMART veto 
  //
          
  bool inlarhole = false;
  bool inlarhole1 = false;
  


  for(unsigned int iJet=0; iJet<jet_AntiKt4LCTopo_pt->size(); iJet++){
    if(susyobj.GetJetTLV(iJet).Pt()>20000.&& susyobj.GetJetTLV(iJet).E()>0.){
    if(m_fakemet.isBad(jet_AntiKt4LCTopo_pt->at(iJet),
		       jet_AntiKt4LCTopo_BCH_CORR_JET->at(iJet),
		       jet_AntiKt4LCTopo_BCH_CORR_CELL->at(iJet),
		       jet_AntiKt4LCTopo_BCH_CORR_DOTX->at(iJet),
		       jet_AntiKt4LCTopo_phi->at(iJet), met_vec.Px(),met_vec.Py(),10000.,10.,-1.,-1.)) inlarhole = true;
    inlarhole1 = susyobj.GetJetTLV(iJet).Pt()>40000.&& jet_AntiKt4LCTopo_BCH_CORR_JET->at(iJet)>0.05 && TMath::Abs(susyobj.GetJetTLV(iJet).DeltaPhi(etm))<0.3;
    if (inlarhole1) return kTRUE;
    }
  }


  
  
  fillcutFlow(14.0);
  


  //
  // Jet Cleaning 
  //
  
    
  for (unsigned int iJet=0;iJet<Jets.size(); iJet++) {
    if ((Jets.at(iJet)).isBad) m_BadJetEvent=true;
    if (m_BadJetEvent) return kTRUE;
    
  }
  	  
 
  // *** EVENT with BAD Jets 
  fillcutFlow(12.0); 


    
  float  BTag_SF =1;
  
  vector<float>  pt_last_jets;
  vector<float>  eta_last_jets;
  vector<float>  val_btag;
  vector<int>    truth_label_last_jets;
  if(!m_isDATA){
  for(int i=0; i<Jets.size(); i++){
    if((Jets.at(i)).Pt()>20000. && TMath::Abs((Jets.at(i)).Eta())<2.5){
    pt_last_jets.push_back((Jets.at(i)).Pt());
    //      eta_last_jets.push_back(jet_AntiKt4LCTopo_eta->at((Jets.at(i)).truthIndex));                                                                                                                                                   
    eta_last_jets.push_back((Jets.at(i)).Eta());
    val_btag.push_back(jet_AntiKt4LCTopo_flavor_weight_MV1->at((Jets.at(i)).truthIndex));
    truth_label_last_jets.push_back(jet_AntiKt4LCTopo_flavor_truth_label->at((Jets.at(i)).truthIndex));
  }
  }
  std::pair<vector<float>,vector<float> > BTag_wei;
  BTag_wei =  BTagComputation->BTagCalibrationFunction(pt_last_jets,eta_last_jets,val_btag,truth_label_last_jets);
  BTag_SF = BTag_wei.first.at(0);
  }
  
  
  
  //
  // Primary vertex
  //
 
  // count vertexes 
  int ngoodvx=0;
  for(Int_t i=0; i<vx_nTracks->size(); i++){  
    if(vx_nTracks->at(i) >4 ){ 
      ngoodvx++; 
    }
  }
  h_vertex_n->Fill(ngoodvx,w);
  m_GoodVertex=0;
  m_GoodVertex=susyobj.IsGoodVertex(vx_nTracks); // good primary vertex
  if (!m_GoodVertex) return kTRUE; 
  fillcutFlow(7.0); 
   


  //if dr<0.01 veto the event                                                                                                                                                                                          
  for (unsigned int j=0;j<Electrons.size(); j++) {                                                                                                                                                            
    if (!(Electrons.at(j)).isBaseline) continue;                                                                                                                                                                                          
    for (unsigned int i=0;i<Muons.size(); i++) {                                                                                                                                                                                          
      if (!(Muons.at(i)).isBaseline) continue;                                                                                                                                                                                           
      float dr= (Muons.at(i)).DeltaR(Electrons.at(j));                                                                                                                                                                                  
      if ( dr< 0.01 )return kTRUE;                                                                                                                                                                                                     
    }                                                                                                                                                                                                                                     
  }     

  
  //if dr<0.05; remove the electron with the lowest pt                                                                                                                                                                                  
  for (unsigned int j=0;j<Electrons.size();j++) {                                                                                                                                                                                        
    if (!(Electrons.at(j)).isBaseline) continue;                                                                                                                                                                                           
    for (unsigned int i=0;i<Electrons.size();) {                                                                                                                                                                                       
      if (!(Electrons.at(i)).isBaseline) continue;  
      float dr= (Electrons.at(j)).DeltaR(Electrons.at(i));   
      if ( dr< 0.05 && i>j) { 
	Electrons.erase(Electrons.begin()+i); 
	(Electrons.at(j)).overlap=0;
	
      }else {
	i++;
      }
    }  
  }                                                                                                                                                                                                                                        
  
  
  //                                                                                                                                                              
  // cosmic veto                                                                                                                                                                                                                          
  //                                                                                                                                                                                                                                       
 
 
  
  bool cosmic_veto=false;
  for (unsigned int i=0;i<Muons.size(); i++) {
    if((Muons.at(i)).isCosmic && (Muons.at(i)).isBaseline) cosmic_veto=true;
    if (cosmic_veto) return kTRUE;
  }
  
  fillcutFlow(8.0);
  
     
  
 
  // Energy averaged time
  
 
  double sume=0; double sumt=0;double ratio=0;
  
  
  for(int iJet=0; iJet<Jets.size(); iJet++){
    if((Jets.at(iJet)).Pt()>20000){//change 20
 
	sumt +=    jet_AntiKt4LCTopo_Timing->at((Jets.at(iJet)).truthIndex)*(Jets.at(iJet)).E();
	sume += (Jets.at(iJet)).E();
      }      
    }
  
  if(Jets.size()!=0) {
    if(sumt/sume>5) return kTRUE;
  }
  
  

 
  fillcutFlow(38.0);
  
  
  
  // Bad Ledaing Jet cut
  
  int Jets_size;
  if (Jets.size()<2)Jets_size = Jets.size();
  if (Jets.size()>=2)Jets_size = 2;
  int BadLeadingJet = 0;
  
 

  for(int iJet=0; iJet<Jets_size; iJet++){

    if((Jets.at(iJet)).Pt() >100000. && TMath::Abs((Jets.at(iJet)).Eta()) < 2 ) {//change 100  
if( (jet_AntiKt4LCTopo_sumPtTrk->at((Jets.at(iJet)).truthIndex)/(1.0*(Jets.at(iJet)).Pt()))<0.02 || (jet_AntiKt4LCTopo_sumPtTrk->at((Jets.at(iJet)).truthIndex)/(1.0*(Jets.at(iJet)).Pt()))<0.05 && jet_AntiKt4LCTopo_emfrac->at((Jets.at(iJet)).truthIndex)>0.9)BadLeadingJet++;
    }
  }



    if(BadLeadingJet>0)return kTRUE;
    
    fillcutFlow(39.0);
   
   
    // if (Jets.size() <2) return kTRUE;
    // fillcutFlow(15.0); 
    
    //
    // B tagg
    //
  //  if ( ! ((Jets.at(0)).isBJet ||  (Jets.at(1)).isBJet )  ) return kTRUE; 
  //fillcutFlow(16.0); 
    
    
  // // 
  // // --- multiplicities  of leptons and jets  
  // //  
  /*
    for (unsigned int iMu=0; iMu<Muons.size(); iMu++) {
    if(m_GoodMuon) mu_id.push_back(iMu);
    
    }
  */
   
  bool recoSF=true;
  bool idSF=true;
  bool triggerSF=false;

  //  double MuonWeight=1;
  //double ElectronWeight=1;

  bool El1_tight = false;
  bool El2_tight = false;
  bool Mu1_tight = false;
  bool Mu2_tight = false;

  bool Mu2_tightIso=false;
  bool Mu1_tightIso=false;
  bool El1_tightIso=false;
  bool El2_tightIso=false;
  bool El1_mediumIso=false;
  bool El2_mediumIso=false;
  bool Mu1_mediumIso=false;
  bool Mu2_mediumIso=false;
  bool El1_looseIso=false;
  bool El2_looseIso=false;
  bool Mu1_looseIso=false;
  bool Mu2_looseIso=false;



  unsigned int n_bmuon=0;
  unsigned int n_smuon=0;
  unsigned int n_tightIsoMuon=0;
  unsigned int n_mediumIsoMuon=0;
  unsigned int n_looseIsoMuon=0;
  unsigned int n_tightIsoEle=0;
  unsigned int n_mediumIsoEle=0;
  unsigned int n_looseIsoEle=0;
  SystErr::Syste whichsyste_muSF = SystErr::NONE;

  for (unsigned int i=0;i<Muons.size(); i++) {
    
    if ((Muons.at(i)).isBaseline){
      
      n_bmuon++; // baseline 
      if ((Muons.at(i)).isSignal) {
	
	if(i==0)Mu1_tight =true;
	if(i==1)Mu2_tight =true;
	float dr= (Muons.at(n_smuon)).DeltaR(etm);
	if (i<n_smuon) {  
	  h_mu_pt[n_smuon]->Fill((Muons.at(n_smuon)).Pt(),w);
	  h_mu_eta[n_smuon]->Fill((Muons.at(n_smuon)).Eta(),w);
	  h_mu_phi[n_smuon]->Fill((Muons.at(n_smuon)).Phi(),w);
	  unsigned int ind=(Muons.at(n_smuon)).truthIndex;
	  h_mu_ptcone20[n_smuon]->Fill(mu_staco_ptcone20->at(ind),w);	      
	  h_mu_d0_exPV[n_smuon]->Fill(mu_staco_d0_exPV->at(ind),w); 
	  h_mu_z0_exPV[n_smuon]->Fill(mu_staco_z0_exPV->at(ind),w); 	 
	}
	n_smuon++;              // signal
      }//closes for eignal mu

      if((Muons.at(i)).isTightIso) {
        n_tightIsoMuon++;
	if(i==0)Mu1_tightIso=true;
        if(i==1)Mu2_tightIso=true;
	
      }//closes for tight      



      if ((Muons.at(i)).isMediumIso) {
	if(i==0)Mu1_mediumIso=true;
	if(i==1)Mu2_mediumIso=true;
	n_mediumIsoMuon++;
      }//closes for medium    
      
         
      if ((Muons.at(i)).isLooseIso) {
	n_looseIsoMuon++;
	if(i==0)Mu1_looseIso=true;
	if(i==1)Mu2_looseIso=true;
      }//closes for loose  
      
      
      
    }//closes for baseline mu
  }
  
  h_mu_or_n->Fill(n_bmuon,w);  
 
  
  SystErr::Syste whichsyste_eleSF = SystErr::NONE;
  unsigned int n_belec=0;
  unsigned int n_selec=0;
  for (unsigned int i=0;i<Electrons.size(); i++) {
    
    if ((Electrons.at(i)).isBaseline){
      n_belec++; // no. baseline electrons   
      
      if ((Electrons.at(i)).isSignal) { 

	    if(i==0)El1_tight =true;
	    if(i==1)El2_tight =true;
	    
	    float dr= (Electrons.at(n_selec)).DeltaR(etm);
	if (n_selec<3) {
	  h_el_etm_dr[n_selec]->Fill(dr,w);
	  h_el_pt[n_selec]->Fill((Electrons.at(n_selec)).Pt(),w); 
	  h_el_eta[n_selec]->Fill((Electrons.at(n_selec)).Eta(),w); 
	  h_el_phi[n_selec]->Fill((Electrons.at(n_selec)).Phi(),w); 
	  unsigned int ind=(Electrons.at(n_selec)).truthIndex;
	  h_el_trackd0[n_selec]->Fill( el_trackd0pv->at(ind),w);	
	}
	n_selec++;   //  signal electrons
           }//closes for signal el

   
      if((Electrons.at(i)).isTightIso){
        if(i==0)El1_tightIso=true;
        if(i==1)El2_tightIso=true;
        n_tightIsoEle++;
      }//closes for tight   
      
      if((Electrons.at(i)).isMediumIso){
        if(i==0)El1_mediumIso=true;
        if(i==1)El2_mediumIso=true;
        n_mediumIsoEle++;
      }//closes for medium      
      
      if ((Electrons.at(i)).isLooseIso){
        if(i==0)El1_looseIso=true;
        if(i==1)El2_looseIso=true;
        n_looseIsoEle++;
      }//closes for loose
      
    }//closes for baseline el 
  }
  
  
  h_el_or_n->Fill(n_belec,w);

    bool btag1=false;
   bool btag2=false; 


  int  n_sjet1 = -1;
  int  n_sjet2 = -1;
  int  n_bjet1 = -1;
  int  n_bjet2 = -1;

 
  unsigned int n_bjet=0;  
  unsigned int n_jet=0;  
  for (unsigned int i=0; i<Jets.size();i++) { 
    if ((Jets.at(i)).isBJet){
      if(i==0)btag1=true;
      if(i==1)btag2=true;  
      if(n_bjet==0) n_bjet1=i;
      if(n_bjet==1) n_bjet2=i;
      
      n_bjet++;
    }
    
    if ((Jets.at(i)).isSignal )
      {
	if(n_jet==0) n_sjet1=i;
	if(n_jet==1) n_sjet2=i;
	n_jet++;	
      }
    
  }
  

  h_jet_b_or_n->Fill(n_bjet,w);
  h_jet_or_n->Fill(n_jet,w);
  


 
  
  //
  // select events with at least two baseline leptons 
  //
  
  if ( !(n_bmuon >= 2  || n_belec >=2  || n_bmuon+n_belec >= 2 )   ) return kTRUE; 

 
  fillcutFlow(17.0);
  
  // 
  //  exactly two lepton selection... 
  // 

   bool one_base_mu = (n_bmuon==1) ? true:false;
  bool one_base_el = (n_belec==1)? true:false;
  bool two_base_leptons =  (n_belec==1 &&n_bmuon==1) ?  true:false; 
  bool two_base_el = (n_belec==2) ?  true:false;
  bool two_base_mu = (n_bmuon==2 ) ?  true:false;
   bool one_tight_mu =(n_smuon==1)? true:false;
  bool one_tight_el =(n_selec==1) ? true:false;
  bool two_tight_leptons= (n_selec==1&& n_smuon==1) ? true:false;  
  bool two_tight_el = (n_selec==2) ?  true:false;
  bool two_tight_mu = (n_smuon==2 ) ?  true:false;  
  bool sameflavour = false; 
  bool diffflavour = false; 
  bool samecharge = false;
  bool oppositecharge = false; 
  
  float InvM_ll=0; 
  
  if((n_bmuon + n_belec)!=2)return kTRUE;  
  fillcutFlow(63.0);

 
  if ( two_base_mu)  {  // 2 muons  
    
    double MuonWeight=1;
    double ElectronWeight=1;
    for (unsigned int i=0;i<Muons.size(); i++) {
	double MuonWeight_1 =1;
	if (!m_isDATA  && m_GoodMuon){
	  MuonWeight*=susyobj.GetSignalMuonSF(i ,whichsyste_muSF);
      }
    }

    w=w*MuonWeight*ElectronWeight;   
    h_btag_weight->Fill(BTag_SF,w);
    fillcutFlow(18.0); 
    sameflavour=true;
    if ( (Muons.at(0)).charge   *  (Muons.at(1)).charge  < 0  ) {
      oppositecharge=true;
      fillcutFlow(53.0);
      
    }//closes for opp charge for basqe
    else {
      
      samecharge=true;
      fillcutFlow(54.0);
    }//closes for same charge for base
  }//closes for two_base_mu
  
    
    
    if ( two_tight_mu )  {  // 2 muons  
      sameflavour=true;
      fillcutFlow(19.0);      
      if ( (Muons.at(0)).charge   *  (Muons.at(1)).charge  < 0  ) {
	oppositecharge=true; 
	fillcutFlow(21.0); 
      }//closes for same charge for tight
      else {
	
	samecharge=true;
	fillcutFlow(55.0);
      }//closes for same charge for tight
      
    }//for two tight
  
      
  
      if (two_base_el) { // two baseline electrons  
      
      double MuonWeight=1;
      double ElectronWeight=1;
      for (unsigned int i=0;i<Electrons.size(); i++) { 
	if(!m_isDATA && m_GoodElectron){
	  double ElectronWeight_1 =1;
	  ElectronWeight_1 =susyobj.GetSignalElecSF( el_cl_eta->at((Electrons.at(i)).truthIndex), susyobj.GetElecTLV((Electrons.at(i)).truthIndex).Pt(),recoSF ,idSF, triggerSF,200841,whichsyste_eleSF);
	  ElectronWeight*=ElectronWeight_1;
	}
      }

     
      w=w*ElectronWeight*MuonWeight;

	fillcutFlow(22.0); 
	sameflavour=true;
	if ( (Electrons.at(0)).charge * (Electrons.at(1)).charge < 0  ) {
	  oppositecharge=true;
	  fillcutFlow(56.0);
	}//closes for opp charge for base
	else {
	  samecharge=true;
	  fillcutFlow(57.0);
	}//closes for same charge for base el  
    }//closes for two_base_el
    
    
    
    if (two_tight_el) { // two tight  electrons 
     
      fillcutFlow(23.0); 
      sameflavour=true;
      if ( (Electrons.at(0)).charge * (Electrons.at(1)).charge < 0  ) {
	oppositecharge=true;  
	fillcutFlow(25.0); 
      }//closes for opp charge for tight
      
      else {
        samecharge=true;
        fillcutFlow(58.0);
      }//closes for same charge for tight
      
    }//two_tight_el

  
  
 
    
    if (one_base_el && one_base_mu) { 

     
      double MuonWeight=1;
      double ElectronWeight=1;
      for (unsigned int i=0;i<Muons.size(); i++) {
        double MuonWeight_1 =1;
        if (!m_isDATA){
          MuonWeight*=susyobj.GetSignalMuonSF(i);
	}
      }
      for (unsigned int i=0;i<Electrons.size(); i++) {
	  if(!m_isDATA && m_GoodElectron){
	    double ElectronWeight_1 =1;
	    ElectronWeight_1 =susyobj.GetSignalElecSF( el_cl_eta->at((Electrons.at(i)).truthIndex ), susyobj.GetElecTLV((Electrons.at(i)).truthIndex ).Pt(),recoSF ,idSF, triggerSF,200841,whichsyste_eleSF);
	    ElectronWeight*=ElectronWeight_1;
	  }
	}
      
      w=w*MuonWeight*ElectronWeight;
	fillcutFlow(26.0); 
	diffflavour=true;
	if ( (Electrons.at(0 )).charge * (Muons.at(0)).charge <0 ) {
	  oppositecharge=true;
	  fillcutFlow(60.0);
	}//closes for opp charge for base
	
	else {
	  samecharge=true;
	}//closes for same charge for base
    }//closes for base lepton
    
    
    
    if (one_tight_mu && one_tight_el){

      fillcutFlow(27.0); 
      diffflavour=true;
      if ( (Electrons.at(0 )).charge * (Muons.at(0)).charge <0 ) {
	oppositecharge=true;
	fillcutFlow(28.0); 
      }//closes for opp charge
      
      else {
        samecharge=true;
        fillcutFlow(61.0);
      }//closes for same charge for tight
    }//closes for two_tight_leptons
    


    if(m_isPtcut){    
      
      if(two_base_mu){
	if( (Muons.at(0) + Muons.at(1)).Pt()  >69999)
	  return kTRUE;
      }
      if(two_base_el){
	if( (Electrons.at(0) + Electrons.at(1)).Pt()  >69999 )      
	  return kTRUE;
      }
      if(two_base_leptons){
	if(  (Electrons.at(0)  + Muons.at(0)).Pt() >69999)
	  return kTRUE;
      }
      
    }
    
//////////////////////////////////////////////////////////////////////////////////////////////////////     
//ntuple for TMVA 
   
  bool fillNtuple=true;
  if(fillNtuple){
    //if   (n_jet >=1){
    float had_mT2=0,mt2=0,minv=0,etmiss=0;
      float MVA_weight=-9999,ptEle1=0 , ptEle2=0, ptMuo1=0, ptMuo2=0,ptJet1=0,ptJet2=0,DThetaLepton=-9999,DPhiLepton=-9999,DEtaLepton=-9999,DPhiLep1Met=-9999,DPhiLep2Met=-9999, DPhiL1J1=-9999, DPhiL2J1=-9999, DPhiL1J2=-9999,DPhiL2J2=-9999, DRll=-9999;
      float DPhiJet1Met=-9999,DPhiJet2Met=-9999,DPhiJets=-9999,DThetaJets=-9999,DThetaL1J1=-9999,DThetaL2J1=-9999,DThetaL1J2=-9999,DThetaL2J2=-9999,mctcorrll=-9999,Pbll=-9999,DPhiJetMet=-9999,DPhiMetPbll=-9999;
      float asym;
float etm_dphi1=-9999,etm_dphi2=-9999;
// bool M1_mediumIso,M2_mediumIso;
      float Filler[57]={57*0};
          double myFinalWeight=1.;
 
      double bJet1=1;
      double bJet2=2;
  
      
      if   (n_jet >=1)
	{
	  ptJet1=(Jets.at(0)).Pt();
	  btag1= (Jets.at(0)).isBJet;
	  //	  bJet1=jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN->at(0);
	  bJet1=jet_AntiKt4LCTopo_flavor_weight_MV1->at(0);
	  DPhiJet1Met=TMath::Abs((Jets.at(0)).DeltaPhi(etm));
	}
      
    
      if   (n_jet >=2)
        {
	  ptJet2=(Jets.at(1)).Pt();
	  btag2= (Jets.at(1)).isBJet;
	  //	   bJet2=jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN->at(1);
	   bJet2=jet_AntiKt4LCTopo_flavor_weight_MV1->at(1);
	   DPhiJet2Met=TMath::Abs((Jets.at(1)).DeltaPhi(etm));
	   DPhiJets=TMath::Abs((Jets.at(0)).DeltaPhi(Jets.at(1)));
	   DThetaJets=TMath::Abs((Jets.at(0)).Theta()-(Jets.at(1)).Theta());
	}
    

      //if ( two_tight_el ) {
	     if(two_base_el){
	 DRll= (Electrons.at(0)).DeltaR(Electrons.at(1));
	asym=  ASYM(Electrons.at(0),Electrons.at(1),etm);
	etmiss = met_vec.Mod();
	minv=(Electrons.at(0)+Electrons.at(1)).M();
	DThetaLepton= TMath::Abs((Electrons.at(0)).Theta()-(Electrons.at(1)).Theta());
	DPhiLepton = TMath::Abs((Electrons.at(0)).DeltaPhi(Electrons.at(1)));
	DPhiLep1Met=TMath::Abs((Electrons.at(0)).DeltaPhi(etm));
	DPhiLep2Met=TMath::Abs((Electrons.at(1)).DeltaPhi(etm));
	DEtaLepton = TMath::Abs((Electrons.at(0)).Eta()-(Electrons.at(1)).Eta());
	mt2=MT2(Electrons.at(0),Electrons.at(1),etm);
	ptEle1=Electrons.at(0).Pt() * (Electrons.at(0)).charge;
	ptEle2=Electrons.at(1).Pt() * (Electrons.at(1)).charge;   
       
	TLorentzVector vds;
	 if (n_jet >=2) vds= Jets.at(0) + Jets.at(1);
	double v1t[4] = {(Electrons.at(0)).E(),(Electrons.at(0)).Px(),(Electrons.at(0)).Py(),(Electrons.at(0)).Pz()};
	double v2t[4] = {(Electrons.at(1)).E(),(Electrons.at(1)).Px(),(Electrons.at(1)).Py(),(Electrons.at(1)).Pz()};
	double vdst[4] = {0,0,0,0};//{vds.E(),vds.Px(),vds.Py(),vds.Pz()};                                                                                                                                                                       
	double ptmt[2] = {etm.Px(),etm.Py()};
	double ecm = 8000.0;
	double mxlo = 0.0;
	mctlib::mctlib mc;
	mctcorrll = mc.mt2(v1t,v2t,vdst,ptmt,ecm,mxlo);


	TLorentzVector Pbll_vec = etm + (Electrons.at(0)) + (Electrons.at(1));
       	TLorentzVector Pbll_vec1=etm + (Electrons.at(0)).Pt() + (Electrons.at(1)).Pt();

	 Pbll= TMath::Sqrt(std::pow((Electrons.at(0)).Py()+ (Electrons.at(1)).Py()+met_vec.Py(),2)+std::pow((Electrons.at(0)).Px()+ (Electrons.at(1)).Px()+met_vec.Px(),2)); 
	double Pbll_X=Pbll_vec.Px();
	double Pbll_Y=Pbll_vec.Py();
	    
	TLorentzVector Pbll_XY;
	Pbll_XY.SetPxPyPzE(Pbll_X,Pbll_Y,0.,0.);
	DPhiMetPbll=TMath::Abs(etm.DeltaPhi(Pbll_XY));
     

	if   (n_jet >=1)
	  {
	    DPhiL1J1=TMath::Abs((Electrons.at(0)).DeltaPhi(Jets.at(0)));
            DPhiL2J1=TMath::Abs((Electrons.at(1)).DeltaPhi(Jets.at(0)));
            DThetaL1J1=TMath::Abs((Electrons.at(0)).Theta()-(Jets.at(0)).Theta());
            DThetaL2J1=TMath::Abs((Electrons.at(1)).Theta()-(Jets.at(0)).Theta());
	  }


	if(n_jet>1){
	  had_mT2=had_MT2(Jets.at(n_sjet1),Jets.at(n_sjet2),Pbll_vec);
	}

	    if   (n_jet >=2){


	      DPhiL1J2=TMath::Abs((Electrons.at(0)).DeltaPhi(Jets.at(1)));
	      DPhiL2J2=TMath::Abs((Electrons.at(1)).DeltaPhi(Jets.at(1)));
	      DThetaL1J2=TMath::Abs((Electrons.at(0)).Theta()-(Jets.at(1)).Theta());
	      DThetaL2J2=TMath::Abs((Electrons.at(1)).Theta()-(Jets.at(1)).Theta());
	    
	  }//ends for n_jet >=2                                                                                                                                                                                                            
      
      }//ends for baseline/tight el                             
	  
	  
	     //   if (two_tight_mu) {
	     if(two_base_mu){
	       DRll= (Muons.at(0)).DeltaR(Muons.at(1));
	       etmiss = met_vec.Mod();
	       minv=(Muons.at(0)+Muons.at(1)).M();
	       DThetaLepton = TMath::Abs((Muons.at(0)).Theta()-(Muons.at(1)).Theta());
	       DPhiLepton= TMath::Abs((Muons.at(0)).DeltaPhi(Muons.at(1)));
	       DPhiLep1Met=TMath::Abs((Muons.at(0)).DeltaPhi(etm));
	       DPhiLep2Met=TMath::Abs((Muons.at(1)).DeltaPhi(etm));
	       DEtaLepton = TMath::Abs((Muons.at(0)).Eta()-(Muons.at(1)).Eta());
	       mt2=MT2(Muons.at(0),Muons.at(1),etm);
	       ptMuo1=Muons.at(0).Pt() * Muons.at(0).charge;
	       ptMuo2=Muons.at(1).Pt() * Muons.at(1).charge;

	       TLorentzVector vds;
	       if (n_jet >=2) vds= Jets.at(0) + Jets.at(1);
	       double v1t[4] = {(Muons.at(0)).E(),(Muons.at(0)).Px(),(Muons.at(0)).Py(),(Muons.at(0)).Pz()};
	       double v2t[4] = {(Muons.at(1)).E(),(Muons.at(1)).Px(),(Muons.at(1)).Py(),(Muons.at(1)).Pz()};
	       double vdst[4] = {0,0,0,0};//{vds.E(),vds.Px(),vds.Py(),vds.Pz()};                                                                                                                                                                       
	       double ptmt[2] = {etm.Px(),etm.Py()};
	       double ecm = 8000.0;
	       double mxlo = 0.0;
	       mctlib::mctlib mc;
	       mctcorrll = mc.mt2(v1t,v2t,vdst,ptmt,ecm,mxlo);

	   
	       TLorentzVector Pbll_vec = etm + (Muons.at(0)) + (Muons.at(1));
	       TLorentzVector Pbll_vec1=etm + (Muons.at(0)).Pt() + (Muons.at(1)).Pt();
	       double Pbll_X=Pbll_vec.Px();
	       double Pbll_Y=Pbll_vec.Py();
	        Pbll= TMath::Sqrt(std::pow((Muons.at(0)).Py()+ (Muons.at(1)).Py()+met_vec.Py(),2)+std::pow((Muons.at(0)).Px()+ (Muons.at(1)).Px()+met_vec.Px(),2)); 
	       //	              Pbll=  TMath::Sqrt(std::pow(Pbll_vec.Y(),2)+std::pow(Pbll_vec.X(),2));   
	       TLorentzVector Pbll_XY;
	       Pbll_XY.SetPxPyPzE(Pbll_X,Pbll_Y,0.,0.);
	       DPhiMetPbll=TMath::Abs(etm.DeltaPhi(Pbll_XY));

	    
	      
	       if   (n_jet >=1)
		 {
		   
		   DPhiL1J1=TMath::Abs((Muons.at(0)).DeltaPhi(Jets.at(0)));
		   DPhiL2J1=TMath::Abs((Muons.at(1)).DeltaPhi(Jets.at(0)));
		   DThetaL1J1=TMath::Abs((Muons.at(0)).Theta()-(Jets.at(0)).Theta());
		   DThetaL2J1=TMath::Abs((Muons.at(1)).Theta()-(Jets.at(0)).Theta());
		 }
	       
	       if(n_jet>1){
		 had_mT2=had_MT2(Jets.at(n_sjet1),Jets.at(n_sjet2),Pbll_vec);
		   }




	       if   (n_jet >=2)
                 {
		   
		   DPhiL1J2=TMath::Abs((Muons.at(0)).DeltaPhi(Jets.at(1)));
		   DPhiL2J2=TMath::Abs((Muons.at(1)).DeltaPhi(Jets.at(1)));
		   DThetaL1J2=TMath::Abs((Muons.at(0)).Theta()-(Jets.at(1)).Theta());
		   DThetaL2J2=TMath::Abs((Muons.at(1)).Theta()-(Jets.at(1)).Theta());
		 }//ends for n_jet >=2 
	     }//ends for two_tight_mu                                                                                                                                                                                                             


	     
	      //if (two_tight_leptons && diffflavour) {
	     if(two_base_leptons && diffflavour) {
	       DRll= (Electrons.at(0)).DeltaR(Muons.at(0));	       
	       etmiss = met_vec.Mod();
	       minv=(Electrons.at(0)+Muons.at(0)).M();
	       DPhiLepton= TMath::Abs((Muons.at(0)).DeltaPhi(Electrons.at(0)));
	       
	       if  ( (Electrons.at(0)).Pt() > (Muons.at(0)).Pt() )
		 {
		   
		   DPhiLep1Met=TMath::Abs((Electrons.at(0)).DeltaPhi(etm));
		   DPhiLep2Met=TMath::Abs((Muons.at(0)).DeltaPhi(etm));	    
		 }
	       
	       else {
		 
		 DPhiLep1Met=TMath::Abs((Muons.at(0)).DeltaPhi(etm));
		 DPhiLep2Met=TMath::Abs((Electrons.at(0)).DeltaPhi(etm));
	       }
	       
	       
	       DThetaLepton = TMath::Abs((Electrons.at(0)).Theta()-(Muons.at(0)).Theta());
	       DEtaLepton = TMath::Abs((Electrons.at(0)).Eta()-(Muons.at(0)).Eta());
	       mt2=MT2(Electrons.at(0),Muons.at(0),etm);
	       ptEle1=Electrons.at(0).Pt()* Electrons.at(0).charge;
	       ptMuo1=Muons.at(0).Pt()* Muons.at(0).charge;
	    
	      
	       TLorentzVector vds;
	       if (n_jet >=2) vds= Jets.at(0) + Jets.at(1);
	       double v1t[4] = {(Electrons.at(0)).E(),(Electrons.at(0)).Px(),(Electrons.at(0)).Py(),(Electrons.at(0)).Pz()};
	       double v2t[4] = {(Muons.at(0)).E(),(Muons.at(0)).Px(),(Muons.at(0)).Py(),(Muons.at(0)).Pz()};
	       double vdst[4] = {0,0,0,0};//{vds.E(),vds.Px(),vds.Py(),vds.Pz()};                                                                                                                                                                       
	       double ptmt[2] = {etm.Px(),etm.Py()};
	       double ecm = 8000.0;
	       double mxlo = 0.0;
	       mctlib::mctlib mc;
	       mctcorrll = mc.mt2(v1t,v2t,vdst,ptmt,ecm,mxlo);

	     
	       TLorentzVector Pbll_vec = etm + (Electrons.at(0)) + (Muons.at(0));
	       TLorentzVector Pbll_vec1=etm + (Electrons.at(0)).Pt() + (Muons.at(0)).Pt();
	       double Pbll_X=Pbll_vec.Px();
	       double Pbll_Y=Pbll_vec.Py();

	        Pbll= TMath::Sqrt(std::pow((Electrons.at(0)).Py()+ (Muons.at(0)).Py()+met_vec.Py(),2)+std::pow((Electrons.at(0)).Px()+ (Muons.at(0)).Px()+met_vec.Px(),2)); 
	      
	       TLorentzVector Pbll_XY;
	       Pbll_XY.SetPxPyPzE(Pbll_X,Pbll_Y,0.,0.);
	       DPhiMetPbll=TMath::Abs(etm.DeltaPhi(Pbll_XY));


	     
               if   (n_jet >=1)
		 {
		   if  ( (Electrons.at(0)).Pt() > (Muons.at(0)).Pt() )
                     {

                       DPhiL1J1=TMath::Abs((Electrons.at(0)).DeltaPhi(Jets.at(0)));
                       DPhiL2J1=TMath::Abs((Muons.at(0)).DeltaPhi(Jets.at(0)));
                       DThetaL1J1=TMath::Abs((Electrons.at(0)).Theta()-(Jets.at(0)).Theta());
                       DThetaL2J1=TMath::Abs((Muons.at(0)).Theta()-(Jets.at(0)).Theta());
		     }
		   else {
		  
		     DPhiL1J1=TMath::Abs((Muons.at(0)).DeltaPhi(Jets.at(0)));
                     DPhiL2J1=TMath::Abs((Electrons.at(0)).DeltaPhi(Jets.at(0)));
                     DThetaL1J1=TMath::Abs((Muons.at(0)).Theta()-(Jets.at(0)).Theta());
                     DThetaL2J1=TMath::Abs((Electrons.at(0)).Theta()-(Jets.at(0)).Theta());

		   }

		   }//closes for (n_jet >=1)
	 

	       if(n_jet>1){
		 had_mT2=had_MT2(Jets.at(n_sjet1),Jets.at(n_sjet2),Pbll_vec);
		   }



	       if   (n_jet >=2)
		 {
		   
		   if  ( (Electrons.at(0)).Pt() > (Muons.at(0)).Pt() )
		     {
		       DPhiL1J2=TMath::Abs((Electrons.at(0)).DeltaPhi(Jets.at(1)));
		       DPhiL2J2=TMath::Abs((Muons.at(0)).DeltaPhi(Jets.at(1)));
		       DThetaL1J2=TMath::Abs((Electrons.at(0)).Theta()-(Jets.at(1)).Theta());
		       DThetaL2J2=TMath::Abs((Muons.at(0)).Theta()-(Jets.at(1)).Theta());
		     }
		   else {
		     
		     DPhiL1J2=TMath::Abs((Muons.at(0)).DeltaPhi(Jets.at(1)));
		     DPhiL2J2=TMath::Abs((Electrons.at(0)).DeltaPhi(Jets.at(1)));
		     DThetaL1J2=TMath::Abs((Muons.at(0)).Theta()-(Jets.at(1)).Theta());
		     DThetaL2J2=TMath::Abs((Electrons.at(0)).Theta()-(Jets.at(1)).Theta());
		   }
		   
		 }//ends for n_jet >=2
	   
	     }//ends for (two baseline/tight_leptons && diffflavour)       
	     
	     myFinalWeight=w;       
	   
	    
	     float meff=0.;
	     meff = etmiss+ptJet1+ptJet2+fabs(ptEle1)+fabs(ptEle2)+fabs(ptMuo1)+fabs(ptMuo2);
	     float met_over_meff=0;
	     if(meff>0) met_over_meff=etmiss/meff;
	     
	   
	     Filler[0]=RunNumber;
	     Filler[1]=EventNumber;
	     Filler[2]=myFinalWeight;
	     Filler[3]=etmiss;
	     Filler[4]=ptJet1;
	     Filler[5]=ptJet2;
	     Filler[6]=ptEle1;
	     Filler[7]=ptEle2;
	     Filler[8]=ptMuo1;
	     Filler[9]=ptMuo2;
	     Filler[10]=minv;
	     Filler[11]=meff;
	     Filler[12]=met_over_meff;
	     Filler[13]=mt2;
	     Filler[14]=mctcorrll;
	     Filler[15]=bJet1;
	     Filler[16]=bJet2;
	     Filler[17]=DPhiJets;
	     Filler[18]=DThetaJets;
	     Filler[19]=DPhiJet1Met;
	     Filler[20]=DPhiJet2Met;
	     Filler[21]=DPhiLepton;
	     Filler[22]=DThetaLepton;
	     Filler[23]=DPhiLep1Met;
	     Filler[24]=DPhiLep2Met;
	     Filler[25]=DPhiL1J1;
	     Filler[26]=DPhiL1J2;
	     Filler[27]=DPhiL2J1;
	     Filler[28]=DPhiL2J2;
	     Filler[29]=DThetaL1J1;
	     Filler[30]=DThetaL1J2;
	     Filler[31]=DThetaL2J1;
	     Filler[32]=DThetaL2J2;
	     Filler[33]=Pbll; 
	     Filler[34]=El1_tight;
	     Filler[35]=El2_tight;
	     Filler[36]=Mu1_tight;
	     Filler[37]=Mu2_tight;
	     Filler[38]=BTag_SF;  
	     Filler[39]=DEtaLepton;
	     Filler[40]=had_mT2;      
             Filler[41]=DPhiMetPbll;	     
	     Filler[42]=DRll;
	     Filler[43]=btag1;
	     Filler[44]=btag2;	   
	     Filler[45]=Mu1_looseIso;
             Filler[46]=Mu2_looseIso;
             Filler[47]=El1_looseIso;
             Filler[48]=El2_looseIso;
	     Filler[49]=Mu1_mediumIso;                                                                
             Filler[50]=Mu2_mediumIso;   
             Filler[51]=El1_mediumIso;                                                                      
             Filler[52]=El2_mediumIso;                                                                      
	     Filler[53]=Mu1_tightIso;	     
	     Filler[54]=Mu2_tightIso;
	     Filler[55]=El1_tightIso;
	     Filler[56]=El2_tightIso;
	     if (oppositecharge){
	       
	       theNtuple->Fill(Filler);
	     }
	     //    }  
  }//ends for if(fillNtuple)    
  
  

  //////////////////////////////////////////////////////////////////////////////////////////////////////

 
  //if(n_smuon + n_selec==2)fillcutFlow(64.0);
  
  float had_mT2=0,DPhiMetPbll=-9999,mctcorrll=-9999, Pbll=-9999,DRll=-9999;  
  TLorentzVector vds;
  if(n_jet >=2) vds= Jets.at(0) + Jets.at(1);
  double vdst[4] = {0,0,0,0};//{vds.E(),vds.Px(),vds.Py(),vds.Pz()};
  double ptmt[2] = {etm.Px(),etm.Py()};
  double ecm = 8000.0;
  double mxlo = 0.0;
  mctlib::mctlib mc;
  TLorentzVector Pbll_XY;
  
    if (two_base_mu  && sameflavour) {  
      //       if (two_tight_mu  && sameflavour) {
      fillcutFlow(65.0);      
      float asym=  ASYM(Muons.at(0),Muons.at(1),etm);
      double minv=(Muons.at(0)+Muons.at(1)).M();
      double mt2=MT2(Muons.at(0),Muons.at(1),etm);
      double v1t[4] = {(Muons.at(0)).E(),(Muons.at(0)).Px(),(Muons.at(0)).Py(),(Muons.at(0)).Pz()};
      double v2t[4] = {(Muons.at(1)).E(),(Muons.at(1)).Px(),(Muons.at(1)).Py(),(Muons.at(1)).Pz()};
      mctcorrll = mc.mt2(v1t,v2t,vdst,ptmt,ecm,mxlo);

      TLorentzVector Pbll_vec = etm + (Muons.at(0)) + (Muons.at(1));
      TLorentzVector Pbll_vec1=etm + (Muons.at(0)).Pt() + (Muons.at(1)).Pt();
      double Pbll_X=Pbll_vec.Px();
      double Pbll_Y=Pbll_vec.Py();
      Pbll= TMath::Sqrt(std::pow((Muons.at(0)).Py()+ (Muons.at(1)).Py()+met_vec.Py(),2)+std::pow((Muons.at(0)).Px()+ (Muons.at(1)).Px()+met_vec.Px(),2));
      TLorentzVector Pbll_XY;
      Pbll_XY.SetPxPyPzE(Pbll_X,Pbll_Y,0.,0.);
      DPhiMetPbll=TMath::Abs(etm.DeltaPhi(Pbll_XY));
  
    if (oppositecharge){
      fillcutFlow(66.0);
      // if ( minv >20 ) {
	
        fillcutFlow(20.0);
        
	//if ((Muons.at(0)).Pt() > 25){                                         
	  
	  fillcutFlow(31.0);

	    if   (n_jet > 1){
	      had_mT2=had_MT2(Jets.at(0),Jets.at(1),Pbll_vec);
	      DRll= (Muons.at(0)).DeltaR(Muons.at(1));    
	      fillcutFlow(32.0);
	           if (met_vec.Mod()/1000 > 120){

		double a=(Muons.at(0)).Pt()+(Muons.at(1)).Pt();
                double b= ((Muons.at(0)).Pt()+(Muons.at(1)).Pt())/((Jets.at(0)).Pt()+(Jets.at(1)).Pt());
              h_mu_os_minv->Fill(minv,w);
              h_mu_os_dphi->Fill(TMath::Abs((Muons.at(0)).DeltaPhi(Muons.at(1))),w);
              h_mu_os_dtheta->Fill(TMath::Abs((Muons.at(0)).Theta()-(Muons.at(1)).Theta()),w);
              h_mu_os_deta->Fill(TMath::Abs((Muons.at(0)).Eta()-(Muons.at(1)).Eta()),w);
              h_mu_etm_dphi[0]->Fill(TMath::Abs((Muons.at(0)).DeltaPhi(etm)),w);
              h_mu_etm_dphi[1]->Fill(TMath::Abs((Muons.at(1)).DeltaPhi(etm)),w);
              h_mu_jet_dphi[0]->Fill(TMath::Abs((Muons.at(0)).DeltaPhi(Jets.at(0))),w);
              h_mu_jet_dphi[1]->Fill(TMath::Abs((Muons.at(1)).DeltaPhi(Jets.at(0))),w);
              h_mu_os_asym->Fill(asym,w);
              h_mu_mt2->Fill(mt2,w);
              h_etmiss_mu->Fill( met_vec.Mod()/1000,w);
              h_mu_os_pt[0]->Fill((Muons.at(0).Pt()),w);
              h_mu_os_pt[1]->Fill((Muons.at(1).Pt()),w);
              h_mu_os_Pbll->Fill(Pbll,w);
              h_mu_os_had_mt2->Fill(had_mT2,w);
	      h_mu_os_DRll->Fill(DRll,w);
              h_mu_os_PtLep->Fill(a,w);
              h_mu_os_PtRatio->Fill(b,w);
	      h_mu_os_DPhiMetPbll->Fill(DPhiMetPbll,w);

	     } 
	    	    }//closes for (n_jet > 1)
	    //}//closes for ((Muons.at(0)).Pt() > 25)
	    // }//closes for minv
    } else { // same charge                                                                               
      
      h_mu_ss_minv->Fill(minv,w);
      h_mu_ss_asym->Fill(asym,w);
    }//closes for same charge
    
	 }//closes for (two_tight_mu )
	 
 
  
   
  if ( two_base_el && sameflavour ){
    //if ( two_tight_el && sameflavour ){
    fillcutFlow(67.0);
    
    float asym=  ASYM(Electrons.at(0),Electrons.at(1),etm);
    double minv=(Electrons.at(0)+Electrons.at(1)).M();
    double mt2=MT2(Electrons.at(0),Electrons.at(1),etm);

    double v1t[4] = {(Electrons.at(0)).E(),(Electrons.at(0)).Px(),(Electrons.at(0)).Py(),(Electrons.at(0)).Pz()};
    double v2t[4] = {(Electrons.at(1)).E(),(Electrons.at(1)).Px(),(Electrons.at(1)).Py(),(Electrons.at(1)).Pz()};
     mctcorrll = mc.mt2(v1t,v2t,vdst,ptmt,ecm,mxlo);
 
    TLorentzVector Pbll_vec = etm + (Electrons.at(0)) + (Electrons.at(1));
    TLorentzVector Pbll_vec1=etm + (Electrons.at(0)).Pt() + (Electrons.at(1)).Pt();

    Pbll= TMath::Sqrt(std::pow((Electrons.at(0)).Py()+ (Electrons.at(1)).Py()+met_vec.Py(),2)+std::pow((Electrons.at(0)).Px()+ (Electrons.at(1)).Px()+met_vec.Px(),2));
    double Pbll_X=Pbll_vec.Px();
    double Pbll_Y=Pbll_vec.Py();

    TLorentzVector Pbll_XY;
    Pbll_XY.SetPxPyPzE(Pbll_X,Pbll_Y,0.,0.);
    DPhiMetPbll=TMath::Abs(etm.DeltaPhi(Pbll_XY));
    
    if (oppositecharge){
      fillcutFlow(68.0);

      //      if ( minv >20 ) {
	
	fillcutFlow(33.0);
	//if ((Electrons.at(0)).Pt() > 25 ){ 
	  fillcutFlow(34.0);                                                 
	  if   (n_jet >1){                                                                                
	    had_mT2=had_MT2(Jets.at(0),Jets.at(1),Pbll_vec);	
	    DRll= (Electrons.at(0)).DeltaR(Electrons.at(1));  
	      fillcutFlow(24.0);

	           if (met_vec.Mod()/1000 > 120){
		double a=(Electrons.at(0)).Pt()+(Electrons.at(1)).Pt();
		double b= ((Electrons.at(0)).Pt()+(Electrons.at(1)).Pt())/((Jets.at(0)).Pt()+(Jets.at(1)).Pt());
              h_el_os_minv->Fill(minv,w);
              h_el_os_dphi->Fill(TMath::Abs((Electrons.at(0)).DeltaPhi(Electrons.at(1))),w);
              h_el_os_deta->Fill(TMath::Abs((Electrons.at(0)).Eta()-(Electrons.at(1)).Eta()),w);
              h_el_os_dtheta->Fill(TMath::Abs((Electrons.at(0)).Theta()-(Electrons.at(1)).Theta()),w);
              h_el_etm_dphi[0]->Fill(TMath::Abs((Electrons.at(0)).DeltaPhi(etm)),w);
              h_el_etm_dphi[1]->Fill(TMath::Abs((Electrons.at(1)).DeltaPhi(etm)),w);
              h_el_jet_dphi[0]->Fill(TMath::Abs((Electrons.at(0)).DeltaPhi(Jets.at(0))),w);
              h_el_jet_dphi[1]->Fill(TMath::Abs((Electrons.at(1)).DeltaPhi(Jets.at(0))),w);
              h_el_os_asym->Fill(asym,w);
              h_el_mt2->Fill(mt2,w);
              h_etmiss_e->Fill( met_vec.Mod()/1000,w);
              h_el_os_pt[0]->Fill((Electrons.at(0).Pt()),w);
              h_el_os_pt[1]->Fill((Electrons.at(1).Pt()),w);
              h_el_os_Pbll->Fill(Pbll,w);
              h_el_os_had_mt2->Fill(had_mT2,w);
	      h_el_os_DRll->Fill(DRll,w);
	      h_el_os_PtLep->Fill(a,w);
	      h_el_os_PtRatio->Fill(b,w);	   
	      h_el_os_DPhiMetPbll->Fill(DPhiMetPbll,w);	      
	       }	    
	    //	    }//closes for El.Pt()>25
	}//closes for (n_jet >=2)
	  // }//closes for minv
      
    }else {    // same sign                                                                             
      samecharge=true;
      h_el_ss_minv->Fill(minv,w);
      h_el_ss_asym->Fill(asym,w);
      }
  }//for tight el
  
  
  
  
  
  if (one_base_mu && one_base_el){ 
    //         if (one_tight_mu && one_tight_el){
    fillcutFlow(69.0);
    
    float asym= ASYM(Electrons.at(0),Muons.at(0),etm);
    double minv=(Electrons.at(0)+Muons.at(0)).M();
    double mt2=MT2(Electrons.at(0),Muons.at(0),etm);
    
    double v1t[4] = {(Electrons.at(0)).E(),(Electrons.at(0)).Px(),(Electrons.at(0)).Py(),(Electrons.at(0)).Pz()};
    double v2t[4] = {(Muons.at(0)).E(),(Muons.at(0)).Px(),(Muons.at(0)).Py(),(Muons.at(0)).Pz()};
    mctcorrll = mc.mt2(v1t,v2t,vdst,ptmt,ecm,mxlo);
    double mt2_c = mc.mt2(v1t,v2t,vdst,ptmt,ecm,mxlo);
    
    
    TLorentzVector Pbll_vec = etm + (Electrons.at(0)) + (Muons.at(0));
    TLorentzVector Pbll_vec1=etm + (Electrons.at(0)).Pt() + (Muons.at(0)).Pt();
    double Pbll_X=Pbll_vec.Px();
    double Pbll_Y=Pbll_vec.Py();
    
    Pbll= TMath::Sqrt(std::pow((Electrons.at(0)).Py()+ (Muons.at(0)).Py()+met_vec.Py(),2)+std::pow((Electrons.at(0)).Px()+ (Muons.at(0)).Px()+met_vec.Px(),2));
    
    TLorentzVector Pbll_XY;
    Pbll_XY.SetPxPyPzE(Pbll_X,Pbll_Y,0.,0.);
    DPhiMetPbll=TMath::Abs(etm.DeltaPhi(Pbll_XY));
    
    
    if (oppositecharge ) {     // opposite sign                                                                 
      fillcutFlow(70.0);
      //      if (minv  >  20){                                                                                           
      
      fillcutFlow(35.0);
      //	if (((Electrons.at(0)).Pt() > 25) || ((Muons.at(0)).Pt() > 25)){                                          
      fillcutFlow(36.0);
      if   (n_jet >1){                                                                                         
	had_mT2=had_MT2(Jets.at(0),Jets.at(1),Pbll_vec);	
	DRll= (Muons.at(0)).DeltaR(Electrons.at(0));
	fillcutFlow(37.0);
	
	if (met_vec.Mod()/1000 > 120){
	double a=(Electrons.at(0)).Pt()+(Muons.at(0)).Pt();
	double b= ((Electrons.at(0)).Pt()+(Muons.at(0)).Pt())/((Jets.at(0)).Pt()+(Jets.at(1)).Pt());
	
	h_etmiss_emu->Fill( met_vec.Mod()/1000,w);
	
	h_elmu_os_minv->Fill(minv,w);
	h_elmu_os_dphi->Fill(TMath::Abs((Electrons.at(0)).DeltaPhi(Muons.at(0))),w);
	h_elmu_os_dtheta->Fill(TMath::Abs((Electrons.at(0)).Theta()-(Muons.at(0)).Theta()),w);
	h_elmu_os_deta->Fill(TMath::Abs((Electrons.at(0)).Eta()-(Muons.at(0)).Eta()),w);
	h_elmu_os_Pbll->Fill(Pbll,w);
	h_elmu_os_had_mt2->Fill(had_mT2,w);
	if  ( (Electrons.at(0)).Pt() > (Muons.at(0)).Pt() )
	  {
	    h_elmu_etm_dphi[0]->Fill(TMath::Abs((Electrons.at(0)).DeltaPhi(etm)),w);
	    h_elmu_etm_dphi[1]->Fill(TMath::Abs((Muons.at(0)).DeltaPhi(etm)),w);
	    h_elmu_os_pt[0]->Fill((Electrons.at(0).Pt()),w);
	    h_elmu_os_pt[1]->Fill((Muons.at(0).Pt()),w);
	    h_elmu_jet_dphi[0]->Fill(TMath::Abs((Electrons.at(0)).DeltaPhi(Jets.at(0))),w);
	    h_elmu_jet_dphi[1]->Fill(TMath::Abs((Muons.at(0)).DeltaPhi(Jets.at(0))),w);
	  }
	else {
	  h_elmu_etm_dphi[0]->Fill(TMath::Abs((Muons.at(0)).DeltaPhi(etm)),w);
	  h_elmu_etm_dphi[1]->Fill(TMath::Abs((Electrons.at(0)).DeltaPhi(etm)),w);
	  h_elmu_jet_dphi[0]->Fill(TMath::Abs((Muons.at(0)).DeltaPhi(Jets.at(0))),w);
	  h_elmu_jet_dphi[1]->Fill(TMath::Abs((Electrons.at(0)).DeltaPhi(Jets.at(0))),w);
	  h_elmu_os_pt[0]->Fill((Muons.at(0).Pt()),w);
	  h_elmu_os_pt[1]->Fill((Electrons.at(0).Pt()),w);
	}
	
	
	h_elmu_mt2->Fill(mt2,w);
	h_elmu_mt2_cor->Fill(mt2_c,w);
	h_elmu_os_asym->Fill(asym,w);
	h_elmu_os_DRll->Fill(DRll,w);
	h_elmu_os_PtLep->Fill(a,w);
	h_elmu_os_PtRatio->Fill(b,w);
	h_elmu_os_DPhiMetPbll->Fill(DPhiMetPbll,w);
	
	 }//closes for MET>120
      }//closes for (n_jet >1)
      //}//closes for Pt()>25
      // }//closes for minv
      
    }else {
      
      // same sign                                                                                              
      samecharge=true;
      // h_cutFlow->Fill(14.0);                                                                                 
      //h_cutFlow_w->Fill(14.0,w);                                                                              
      //h_elmu_ss_minv->Fill(minv,w);                                                                           
      h_elmu_ss_asym->Fill(asym,w);
      
    }
  }
  
 
   

  //return kTRUE;      
  

  if ( two_base_leptons && oppositecharge) {
    fillcutFlow(29.0);
    h_etmiss->Fill( met_vec.Mod()/1000,w); 
    
    if (Jets.size() <2)  return kTRUE;
    fillcutFlow(15.0); 
    
    // jet molteplicity  
    // >2 
    h_jet_n->Fill(Jets.size(),w);
    float dphi= TMath::Abs((Jets.at(0)).DeltaPhi(Jets.at(1)));
    h_jet_jet_dphi->Fill(dphi,w);
    
    float ht=0;
    for (unsigned int i=0;i<Jets.size(); i++) { 
      unsigned int ind=(Jets.at(i)).truthIndex;
      ht+=jet_AntiKt4LCTopo_pt->at(ind) /1000; 
      if (i<3 ) {
  	float dr= (Jets.at(i)).DeltaR(etm);
  	h_jet_etm_dr[i]->Fill(dr,w);
  	//h_jet_pt[i]->Fill(jet_AntiKt4LCTopo_pt->at(ind) /1000.,w); 
  	//h_jet_eta[i]->Fill(jet_AntiKt4LCTopo_emscale_eta->at(ind),w ); 
	//	h_jet_phi[i]->Fill(jet_AntiKt4LCTopo_emscale_phi->at(ind),w ); 
	//	h_jet_m[i]->Fill((Jets.at(i)).M(),w);
  	//h_jet_ftag_sv0[i]->Fill(jet_AntiKt4LCTopo_flavor_weight_SV0->at(ind),w );
  	//h_jet_ftag_JProb[i]->Fill(jet_AntiKt4LCTopo_flavor_weight_MV1->at(ind),w );
  	//	h_jet_ftag_ntrack[i]->Fill(jet_AntiKt4LCTopo_flavor_weight_TrackCounting2D->at(ind),w );	
      } 
    }
    

   
    h_jet_ht->Fill(ht,w); 
    
    // compute mt2... 
    //pa, pb = {mass, px, py}
    //pmiss  = {NULL, pxmiss, pymiss}
    //mn     = invisible particle mass
    double pa[3] = { 0, 0, 0 };
    double pb[3] = { 0, 0, 0 };
    double pmiss[3] = { 0, 0,0 };
    double mn    = 1.;
    pa[0]=(Jets.at(0)).M();
    pa[1]=(Jets.at(0)).X();
    pa[2]=(Jets.at(0)).Y();
    pb[0]=(Jets.at(1)).M();
    pb[1]=(Jets.at(1)).X();
    pb[2]=(Jets.at(1)).Y();
    pmiss[1]=etm.X();
    pmiss[2]=etm.Y();
    
    mt2_bisect::mt2 mt2_jet;
    mt2_jet.set_momenta(pa,pb,pmiss);
    mt2_jet.set_mn(mn);    
    
    h_jet_mt2->Fill( mt2_jet.get_mt2(),w); 

    
    }
    

  ////
  
  
  



  //  Fillcutflow(30.0);
  
  





  
  
    //Std::cout<<"Process exit  "<<std::endl;
    return kTRUE;
}

void read_d4pd::SlaveTerminate()
{

  std::cout<<"slave terminate   "<<std::endl;
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  if (doPUreweight) {
    if(CreateConfigFile) pileup.WriteToFile("/Users/sdarmora/ST_3_20/run/PileUpFiles/ConfFile_xxxyyy_p1328.root");
  }
  
  // save the  Event yield 
  // m_eventInfo.printTotEvt(); 

  // m_totevt= m_eventInfo.getTotEvt();


  // std::map< int,int>::const_iterator itr; 
  // unsigned int in=0;
  // for (itr=m_totevt.begin() ;itr !=m_totevt.end(); itr++) {
  //   in++;
  //   g_event_yield->SetPoint(in, (*itr).first, (*itr).second);  
  //   //std::cout << "Terminate: Run and tot events in GRL " << (*itr).first <<'\t' << (*itr).second <<'\n';
  // } 


}


void read_d4pd::Terminate()
{
  
  std::cout<<"terminate   "<<std::endl;
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.


  std::cout<<"hello"<<endl;

  // save the  Event yield 
  // m_eventInfo.printTotEvt(); 
 
  //std::map <int,int> totevt= m_eventInfo.getTotEvt();
  //TGraph *tg =new TGraph();
 
  //std::map< int,int>::const_iterator itr; 
  //unsigned int in=0;
  //for (itr=totevt.begin() ;itr !=totevt.end(); itr++) {
  // in++;
  // tg->SetPoint(in, (*itr).first, (*itr).second);  
  // //std::cout << "Terminate: Run and tot events in GRL " << (*itr).first <<'\t' << (*itr).second <<'\n';
  //} 
  //tg->SetName("h_evt_yield_run"); 

  //SaveHistos(); 
 
  TFile *fileOut = new TFile( "histograms.root","RECREATE");
  fileOut->cd();
  
  theNtuple = dynamic_cast <TNtuple *>(fOutput->FindObject("susy")); 
  theNtuple->Write();


  h_btag_weight = dynamic_cast <TH1F *>(fOutput->FindObject(Form("h_btag_weight")));  
  h_btag_weight->Write();


  h_el_trigg_freq = dynamic_cast <TH1D *>(fOutput->FindObject(Form("h_el_trigg_freq")));
  h_el_trigg_freq->Write();

  h_mu_trigg_freq = dynamic_cast <TH1D *>(fOutput->FindObject(Form("h_mu_trigg_freq")));
  h_mu_trigg_freq->Write();

  
  h_event_yield = dynamic_cast <TH1F *>(fOutput->FindObject(Form("h_event_yield"))); 
  h_event_yield->Write(); 
  
  h_totevent = dynamic_cast <TH1D *>(fOutput->FindObject(Form("h_totevent"))); 
  h_totevent->Write(); 

  h_cutFlow = dynamic_cast<TH1D *>(fOutput->FindObject(Form("h_cutFlow")));
  h_cutFlow->Write();
  
  h_cutFlow_w = dynamic_cast<TH1D *>(fOutput->FindObject(Form("h_cutFlow_w")));
  h_cutFlow_w->Write();

  h_el_n = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_n")));
  h_el_n->Write();
 
 h_mu_n = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_n")));
 h_mu_n->Write();


 h_elmu_n = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_n")));
 h_elmu_n->Write();
  
 h_jet_n = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_n")));
 h_jet_n->Write();
 
 h_el_or_n = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_or_n")));
 h_el_or_n->Write();
 
 h_mu_or_n = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_or_n")));
 h_mu_or_n->Write();
  
 h_jet_or_n = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_or_n")));
 h_jet_or_n->Write();

 h_jet_b_or_n = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_b_or_n")));
 h_jet_b_or_n->Write(); 

 h_jet_jet_dphi = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_jet_dphi")));
 h_jet_jet_dphi->Write();

 h_jet_ht = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_ht")));
 h_jet_ht->Write();

 h_vertex_n = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_vertex_n")));
 h_vertex_n->Write();
 
 h_etmiss_mod = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_mod")));
 h_etmiss_mod->Write();

 h_etmiss = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss")));
 h_etmiss->Write();


 h_etmiss_phi = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_phi")));
 h_etmiss_phi->Write();

 h_etmiss_e = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_e")));
 h_etmiss_e->Write();

 h_etmiss_emu = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_emu")));
 h_etmiss_emu->Write();

 h_etmiss_mu = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_mu")));
 h_etmiss_mu->Write();

 


 //lepton dphi 
 h_mu_os_dphi = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_dphi")));
 h_mu_os_dphi->Write();
 h_el_os_dphi = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_dphi")));
 h_el_os_dphi->Write();
 h_elmu_os_dphi = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_dphi")));
 h_elmu_os_dphi->Write();

 //lepton dtheta

 h_mu_os_dtheta = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_dtheta")));
 h_mu_os_dtheta->Write();
 h_el_os_dtheta = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_dtheta")));
 h_el_os_dtheta->Write();
 h_elmu_os_dtheta = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_dtheta")));
 h_elmu_os_dtheta->Write();
 

 //lepton deta

 h_mu_os_deta = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_deta")));
 h_mu_os_deta->Write();
 h_el_os_deta = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_deta")));
 h_el_os_deta->Write();
 h_elmu_os_deta = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_deta")));
 h_elmu_os_deta->Write();

 //Pblll only

 h_mu_os_Pbll = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_Pbll")));
 h_mu_os_Pbll->Write();
 h_el_os_Pbll = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_Pbll")));
 h_el_os_Pbll->Write();
 h_elmu_os_Pbll = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_Pbll")));
 h_elmu_os_Pbll->Write();


 h_el_os_DPhiMetPbll = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_DPhiMetPbll")));
 h_el_os_DPhiMetPbll->Write();
 h_mu_os_DPhiMetPbll = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_DPhiMetPbll")));
 h_mu_os_DPhiMetPbll->Write();
 h_elmu_os_DPhiMetPbll = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_DPhiMetPbll")));
 h_elmu_os_DPhiMetPbll->Write();




 //hadronic mt2

 h_mu_os_had_mt2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_had_mt2")));
 h_mu_os_had_mt2->Write();
 h_el_os_had_mt2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_had_mt2")));
 h_el_os_had_mt2->Write();
 h_elmu_os_had_mt2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_had_mt2")));
 h_elmu_os_had_mt2->Write();




 //delta R betwenn leptons
 h_el_os_DRll = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_DRll")));
 h_el_os_DRll->Write();
 h_mu_os_DRll = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_DRll")));
 h_mu_os_DRll->Write();
 h_elmu_os_DRll = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_DRll")));
 h_elmu_os_DRll->Write();

 h_el_os_PtLep = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_PtLep")));
 h_el_os_PtLep->Write();
 h_mu_os_PtLep = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_PtLep")));
 h_mu_os_PtLep->Write();
 h_elmu_os_PtLep = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_PtLep")));
 h_elmu_os_PtLep->Write();

 h_el_os_PtRatio = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_PtRatio")));
 h_el_os_PtRatio->Write();
 h_mu_os_PtRatio = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_PtRatio")));
 h_mu_os_PtRatio->Write();
 h_elmu_os_PtRatio = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_PtRatio")));
 h_elmu_os_PtRatio->Write();

 // mt2 
 h_mu_mt2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_mt2")));
 h_mu_mt2->Write();
 
 h_el_mt2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_mt2")));
 h_el_mt2->Write();

  h_elmu_mt2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_mt2")));
 h_elmu_mt2->Write();

 h_elmu_mt2_cor = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_mt2_cor")));
 h_elmu_mt2_cor->Write();

 h_jet_mt2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_mt2")));
 h_jet_mt2->Write();

 h_jet_elmu_mt2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_elmu_mt2")));
h_jet_elmu_mt2->Write();

 h_el_mt2_bl = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_mt2_bl")));
 h_el_mt2_bl->Write();

 h_elmu_mt2_bl = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_mt2_bl")));
 h_elmu_mt2_bl->Write();

 h_mu_mt2_bl = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_mt2_bl")));
 h_mu_mt2_bl->Write();


 // minv 
 h_mu_ss_minv = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_ss_minv")));
 h_mu_ss_minv->Write();
 h_mu_os_minv = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_minv")));
 h_mu_os_minv->Write();
 
 h_el_ss_minv = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_ss_minv")));
 h_el_ss_minv->Write();
 h_el_os_minv = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_minv")));
 h_el_os_minv->Write();
 
 h_elmu1_jet_minv = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu1_jet_minv")));
 h_elmu1_jet_minv->Write();
 h_elmu2_jet_minv = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu2_jet_minv")));
 h_elmu2_jet_minv->Write();
 
 h_elmu_os_minv = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_minv")));
 h_elmu_os_minv->Write(); 
 // asym 
 h_mu_ss_asym = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_ss_asym")));
 h_mu_ss_asym->Write();
 h_mu_os_asym = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_asym")));
 h_mu_os_asym->Write();
 
 h_el_ss_asym = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_ss_asym")));
 h_el_ss_asym->Write();
 h_el_os_asym = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_asym")));
 h_el_os_asym->Write();
 
 h_elmu_ss_asym = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_ss_asym")));
 h_elmu_ss_asym->Write();
 h_elmu_os_asym = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_asym")));
 h_elmu_os_asym->Write(); 


 //histo for DF and SF (pt,eta)
 

  char hname[30]; 
 for (int i=0; i<3; i++) {
   
   // electrons 
   sprintf(hname,"h_el_pt_%i",i);
   h_el_pt[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_pt[i]->Write();
   
   sprintf(hname,"h_el_eta_%i",i);      
   h_el_eta[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_eta[i]->Write();
   
   sprintf(hname,"h_el_phi_%i",i);      
   h_el_phi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_phi[i]->Write();  
   
  
   if (i<2) {
     sprintf(hname,"h_el_ss_pt_%i",i);
     h_el_ss_pt[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_ss_pt[i]->Write();
     
     sprintf(hname,"h_el_ss_eta_%i",i); 
     h_el_ss_eta[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_ss_eta[i]->Write();
     
     sprintf(hname,"h_el_ss_phi_%i",i); 
     h_el_ss_phi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_ss_phi[i]->Write();  

     sprintf(hname,"h_el_os_pt_%i",i);
     h_el_os_pt[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_os_pt[i]->Write();
     
     sprintf(hname,"h_el_os_eta_%i",i); 
     h_el_os_eta[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_os_eta[i]->Write();
     
     sprintf(hname,"h_el_os_phi_%i",i); 
     h_el_os_phi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_os_phi[i]->Write();  
     //  muons 
     sprintf(hname,"h_mu_ss_pt_%i",i);
     h_mu_ss_pt[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_ss_pt[i]->Write();
     
     sprintf(hname,"h_mu_ss_eta_%i",i); 
     h_mu_ss_eta[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_ss_eta[i]->Write();
     
     sprintf(hname,"h_mu_ss_phi_%i",i); 
     h_mu_ss_phi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_ss_phi[i]->Write();  

     sprintf(hname,"h_mu_os_pt_%i",i);
     h_mu_os_pt[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_os_pt[i]->Write();
     
     sprintf(hname,"h_mu_os_eta_%i",i); 
     h_mu_os_eta[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_os_eta[i]->Write();
     
     sprintf(hname,"h_mu_os_phi_%i",i); 
     h_mu_os_phi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_os_phi[i]->Write();  
     
     // lepton jet dphi 
     sprintf(hname,"h_mu_jet_dphi_%i",i);
     h_mu_jet_dphi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_jet_dphi[i]->Write();
     
     sprintf(hname,"h_el_jet_dphi_%i",i);
     h_el_jet_dphi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_jet_dphi[i]->Write();

     sprintf(hname,"h_elmu_jet_dphi_%i",i);
     h_elmu_jet_dphi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_elmu_jet_dphi[i]->Write();    

     // lepton met dphi 
     sprintf(hname,"h_mu_etm_dphi_%i",i);
     h_mu_etm_dphi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_etm_dphi[i]->Write();
     
     sprintf(hname,"h_el_etm_dphi_%i",i);
     h_el_etm_dphi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_etm_dphi[i]->Write();

     sprintf(hname,"h_elmu_etm_dphi_%i",i);
     h_elmu_etm_dphi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_elmu_etm_dphi[i]->Write();     


     sprintf(hname,"h_elmu_os_pt_%i",i);
     h_elmu_os_pt[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_elmu_os_pt[i]->Write();

     
   }


   sprintf(hname,"h_el_trackd0_%i",i);  
   h_el_trackd0[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_trackd0[i]->Write();
   
   sprintf(hname,"h_el_etcone20_%i",i); 
   h_el_etcone20[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_etcone20[i]->Write();
   
   sprintf(hname,"h_el_etcone30_%i",i); 
   h_el_etcone30[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_etcone30[i]->Write();
   
 



   // muons 
   sprintf(hname,"h_mu_pt_%i",i); 
   h_mu_pt[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_pt[i]->Write();
  
   sprintf(hname,"h_mu_eta_%i",i); 
   h_mu_eta[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_eta[i]->Write();
  
   sprintf(hname,"h_mu_phi_%i",i); 
   h_mu_phi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_phi[i]->Write();

   sprintf(hname,"h_mu_etcone20_%i",i); 
   h_mu_etcone20[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_etcone20[i]->Write();

   sprintf(hname,"h_mu_etcone30_%i",i); 
   h_mu_etcone30[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_etcone30[i]->Write();

   sprintf(hname,"h_mu_ptcone20_%i",i); 
   h_mu_ptcone20[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_ptcone20[i]->Write();

   sprintf(hname,"h_mu_ptcone30_%i",i); 
   h_mu_ptcone30[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_ptcone30[i]->Write();

   sprintf(hname,"h_mu_d0_exPV_%i",i);  
   h_mu_d0_exPV[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_d0_exPV[i]->Write();

   sprintf(hname,"h_mu_z0_exPV_%i",i);  
   h_mu_z0_exPV[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_z0_exPV[i]->Write();
   
  
   
   // jets 
   sprintf(hname,"h_jet_pt_%i",i);
   h_jet_pt[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_jet_pt[i]->Write();
  
   sprintf(hname,"h_jet_eta_%i",i);     
   h_jet_eta[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_jet_eta[i]->Write();

   sprintf(hname,"h_jet_phi_%i",i);     
   h_jet_phi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_jet_phi[i]->Write();  

   sprintf(hname,"h_jet_m_%i",i);       
   h_jet_m[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_jet_m[i]->Write();

   sprintf(hname,"h_jet_ftag_sv0_%i",i);        
   h_jet_ftag_sv0[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_jet_ftag_sv0[i]->Write();
   
   sprintf(hname,"h_jet_ftag_JProb_%i",i);      
   h_jet_ftag_JProb[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_jet_ftag_JProb[i]->Write();
   
   sprintf(hname,"h_jet_ftag_ntrack_%i",i);     
   h_jet_ftag_ntrack[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_jet_ftag_ntrack[i]->Write();
   

   // dr  

   sprintf(hname,"h_el_mu_dr_%i",i);
   h_el_mu_dr[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_mu_dr[i]->Write();

   sprintf(hname,"h_el_jet_dr_%i",i);
   h_el_jet_dr[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_jet_dr[i]->Write();

   sprintf(hname,"h_el_el_dr_%i",i);
   h_el_el_dr[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_el_dr[i]->Write();

   sprintf(hname,"h_el_trig_dr_%i",i);
   h_el_trig_dr[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_trig_dr[i]->Write();

   sprintf(hname,"h_mu_trig_dr_%i",i);
   h_mu_trig_dr[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_trig_dr[i]->Write();

   sprintf(hname,"h_mu_etm_dr_%i",i);
   h_mu_etm_dr[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_etm_dr[i]->Write();

   sprintf(hname,"h_el_etm_dr_%i",i);
   h_el_etm_dr[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_etm_dr[i]->Write();

   sprintf(hname,"h_jet_etm_dr_%i",i);
   h_jet_etm_dr[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_jet_etm_dr[i]->Write();
 }


 //*************************************************************************************************START C1
 //************************************************************************************************




 h_etmiss_e_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_e_C1")));
 h_etmiss_e_C1->Write();

 h_etmiss_emu_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_emu_C1")));
 h_etmiss_emu_C1->Write();

 h_etmiss_mu_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_mu_C1")));
 h_etmiss_mu_C1->Write();

 h_jet_mt2_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_mt2_C1")));
 h_jet_mt2_C1->Write();



 //lepton dphi 
 h_mu_os_dphi_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_dphi_C1")));
 h_mu_os_dphi_C1->Write();

 h_el_os_dphi_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_dphi_C1")));
 h_el_os_dphi_C1->Write();

 h_elmu_os_dphi_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_dphi_C1")));
 h_elmu_os_dphi_C1->Write();
 //lepton dtheta


 h_mu_os_dtheta_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_dtheta_C1")));
 h_mu_os_dtheta_C1->Write();
 h_el_os_dtheta_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_dtheta_C1")));
 h_el_os_dtheta_C1->Write();
 h_elmu_os_dtheta_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_dtheta_C1")));
 h_elmu_os_dtheta_C1->Write();


 h_mu_os_deta_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_deta_C1")));
 h_mu_os_deta_C1->Write();
 h_el_os_deta_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_deta_C1")));
 h_el_os_deta_C1->Write();
 h_elmu_os_deta_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_deta_C1")));
 h_elmu_os_deta_C1->Write();


 // mt2 
 h_mu_mt2_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_mt2_C1")));
 h_mu_mt2_C1->Write();
 
 h_el_mt2_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_mt2_C1")));
 h_el_mt2_C1->Write();

  h_elmu_mt2_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_mt2_C1")));
 h_elmu_mt2_C1->Write();


 h_jet_mt2_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_mt2_C1")));
 h_jet_mt2_C1->Write();



 // minv 
 h_mu_os_minv_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_minv_C1")));
 h_mu_os_minv_C1->Write();
 
 h_el_os_minv_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_minv_C1")));
 h_el_os_minv_C1->Write();
  
 h_elmu_os_minv_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_minv_C1")));
 h_elmu_os_minv_C1->Write(); 
 // asym 
 h_mu_os_asym_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_asym_C1")));
 h_mu_os_asym_C1->Write();
 

 h_el_os_asym_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_asym_C1")));
 h_el_os_asym_C1->Write();
 
 h_elmu_os_asym_C1 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_asym_C1")));
 h_elmu_os_asym_C1->Write(); 

 
 //histo for DF and SF (pt,eta)
 //char hname[30];
 for (int i=0; i<2; i++) {
                                                                                                                 
   sprintf(hname,"h_el_os_pt_C1_%i",i);
   h_el_os_pt_C1[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_os_pt_C1[i]->Write();

   sprintf(hname,"h_el_os_eta_C1_%i",i);
   h_el_os_eta_C1[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_os_eta_C1[i]->Write();

   sprintf(hname,"h_mu_os_pt_C1_%i",i);
   h_mu_os_pt_C1[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_os_pt_C1[i]->Write();

   sprintf(hname,"h_mu_os_eta_C1_%i",i);
   h_mu_os_eta_C1[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_os_eta_C1[i]->Write();


   sprintf(hname,"h_elmu_os_pt_C1_%i",i);
   h_elmu_os_pt_C1[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_elmu_os_pt_C1[i]->Write();

   sprintf(hname,"h_elmu_os_eta_C1_%i",i);
   h_elmu_os_eta_C1[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_elmu_os_eta_C1[i]->Write();
   
 }
 
 // char hname[30]; 
 for (int i=0; i<2; i++) {
     
    
 
     
     // lepton jet dphi 
     sprintf(hname,"h_mu_jet_dphi_C1_%i",i);
     h_mu_jet_dphi_C1[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_jet_dphi_C1[i]->Write();
     
     sprintf(hname,"h_el_jet_dphi_C1_%i",i);
     h_el_jet_dphi_C1[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_jet_dphi_C1[i]->Write();

     sprintf(hname,"h_elmu_jet_dphi_C1_%i",i); //CHANGE NAME
     h_elmu_jet_dphi_C1[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_elmu_jet_dphi_C1[i]->Write();    

     // lepton met dphi 
     sprintf(hname,"h_mu_etm_dphi_C1_%i",i); //CHANGE NAME
     h_mu_etm_dphi_C1[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_etm_dphi_C1[i]->Write();
     
     sprintf(hname,"h_el_etm_dphi_C1_%i",i);
     h_el_etm_dphi_C1[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_etm_dphi_C1[i]->Write();

     sprintf(hname,"h_elmu_etm_dphi_C1_%i",i);
     h_elmu_etm_dphi_C1[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_elmu_etm_dphi_C1[i]->Write();     
     
 

 }


 //*******************************************************************************end C1

 //************************************************************************************************* START C2
 //************************************************************************************************

h_etmiss_e_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_e_C2")));
 h_etmiss_e_C2->Write();

 h_etmiss_emu_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_emu_C2")));
 h_etmiss_emu_C2->Write();

 h_etmiss_mu_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_mu_C2")));
 h_etmiss_mu_C2->Write();

 h_jet_mt2_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_mt2_C2")));
 h_jet_mt2_C2->Write();



 //lepton dphi 
 h_mu_os_dphi_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_dphi_C2")));
 h_mu_os_dphi_C2->Write();
 h_el_os_dphi_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_dphi_C2")));
 h_el_os_dphi_C2->Write();
 h_elmu_os_dphi_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_dphi_C2")));
 h_elmu_os_dphi_C2->Write();

 //lepton dtheta

 h_mu_os_dtheta_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_dtheta_C2")));
 h_mu_os_dtheta_C2->Write();
 h_el_os_dtheta_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_dtheta_C2")));
 h_el_os_dtheta_C2->Write();
 h_elmu_os_dtheta_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_dtheta_C2")));
 h_elmu_os_dtheta_C2->Write();
 
//lepton deta

 h_mu_os_deta_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_deta_C2")));
 h_mu_os_deta_C2->Write();
 h_el_os_deta_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_deta_C2")));
 h_el_os_deta_C2->Write();
 h_elmu_os_deta_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_deta_C2")));
 h_elmu_os_deta_C2->Write();




 // mt2 
 h_mu_mt2_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_mt2_C2")));
 h_mu_mt2_C2->Write();
 
 h_el_mt2_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_mt2_C2")));
 h_el_mt2_C2->Write();
  h_elmu_mt2_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_mt2_C2")));
 h_elmu_mt2_C2->Write();




 // minv 
 h_mu_os_minv_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_minv_C2")));
 h_mu_os_minv_C2->Write();
 
 h_el_os_minv_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_minv_C2")));
 h_el_os_minv_C2->Write();
 

 
 h_elmu_os_minv_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_minv_C2")));
 h_elmu_os_minv_C2->Write(); 
 
// asym 
 h_mu_os_asym_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_asym_C2")));
 h_mu_os_asym_C2->Write();
 

 h_el_os_asym_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_asym_C2")));
 h_el_os_asym_C2->Write();
 
 h_elmu_os_asym_C2 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_asym_C2")));
 h_elmu_os_asym_C2->Write(); 

 
 //histo for DF and SF (pt,eta)
 //char hname[30];
 for (int i=0; i<2; i++) {
                                                                                                                 
   sprintf(hname,"h_el_os_pt_C2_%i",i);
   h_el_os_pt_C2[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_os_pt_C2[i]->Write();

   sprintf(hname,"h_el_os_eta_C2_%i",i);
   h_el_os_eta_C2[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_os_eta_C2[i]->Write();

   sprintf(hname,"h_mu_os_pt_C2_%i",i);
   h_mu_os_pt_C2[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_os_pt_C2[i]->Write();

   sprintf(hname,"h_mu_os_eta_C2_%i",i);
   h_mu_os_eta_C2[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_os_eta_C2[i]->Write();


   sprintf(hname,"h_elmu_os_pt_C2_%i",i);
   h_elmu_os_pt_C2[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_elmu_os_pt_C2[i]->Write();

   sprintf(hname,"h_elmu_os_eta_C2_%i",i);
   h_elmu_os_eta_C2[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_elmu_os_eta_C2[i]->Write();
   
 }

 // char hname[30]; 
 for (int i=0; i<2; i++) {
     
    

     
     // lepton jet dphi 
     sprintf(hname,"h_mu_jet_dphi_C2_%i",i);
     h_mu_jet_dphi_C2[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_jet_dphi_C2[i]->Write();
     
     sprintf(hname,"h_el_jet_dphi_C2_%i",i);
     h_el_jet_dphi_C2[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_jet_dphi_C2[i]->Write();

     sprintf(hname,"h_elmu_jet_dphi_C2_%i",i);
     h_elmu_jet_dphi_C2[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_elmu_jet_dphi_C2[i]->Write();    

     // lepton met dphi 
     sprintf(hname,"h_mu_etm_dphi_C2_%i",i);
     h_mu_etm_dphi_C2[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_etm_dphi_C2[i]->Write();
     
     sprintf(hname,"h_el_etm_dphi_C2_%i",i);
     h_el_etm_dphi_C2[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_etm_dphi_C2[i]->Write();

     sprintf(hname,"h_elmu_etm_dphi_C2_%i",i);
     h_elmu_etm_dphi_C2[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_elmu_etm_dphi_C2[i]->Write();     
     
  


 } 


 //***********************************************************************************************END C2
 //***************************************************************************************************

 //************************************************************************************************* START C3
 //************************************************************************************************

h_etmiss_e_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_e_C3")));
 h_etmiss_e_C3->Write();

 h_etmiss_emu_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_emu_C3")));
 h_etmiss_emu_C3->Write();

 h_etmiss_mu_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_mu_C3")));
 h_etmiss_mu_C3->Write();


 h_jet_mt2_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_mt2_C3")));
 h_jet_mt2_C3->Write();


 //lepton dphi 
 h_mu_os_dphi_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_dphi_C3")));
 h_mu_os_dphi_C3->Write();
 h_el_os_dphi_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_dphi_C3")));
 h_el_os_dphi_C3->Write();
 h_elmu_os_dphi_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_dphi_C3")));
 h_elmu_os_dphi_C3->Write();

 //lepton dtheta

 h_mu_os_dtheta_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_dtheta_C3")));
 h_mu_os_dtheta_C3->Write();
 h_el_os_dtheta_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_dtheta_C3")));
 h_el_os_dtheta_C3->Write();
 h_elmu_os_dtheta_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_dtheta_C3")));
 h_elmu_os_dtheta_C3->Write();

 //lepton deta

 h_mu_os_deta_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_deta_C3")));
 h_mu_os_deta_C3->Write();
 h_el_os_deta_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_deta_C3")));
 h_el_os_deta_C3->Write();
 h_elmu_os_deta_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_deta_C3")));
 h_elmu_os_deta_C3->Write();



 // mt2 
 h_mu_mt2_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_mt2_C3")));
 h_mu_mt2_C3->Write();
 
 h_el_mt2_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_mt2_C3")));
 h_el_mt2_C3->Write();
  h_elmu_mt2_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_mt2_C3")));
 h_elmu_mt2_C3->Write();




 // minv 
 h_mu_os_minv_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_minv_C3")));
 h_mu_os_minv_C3->Write();
 
 h_el_os_minv_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_minv_C3")));
 h_el_os_minv_C3->Write();
 

 
 h_elmu_os_minv_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_minv_C3")));
 h_elmu_os_minv_C3->Write(); 
 
// asym 
 h_mu_os_asym_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_asym_C3")));
 h_mu_os_asym_C3->Write();
 

 h_el_os_asym_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_asym_C3")));
 h_el_os_asym_C3->Write();
 
 h_elmu_os_asym_C3 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_asym_C3")));
 h_elmu_os_asym_C3->Write(); 

 
 //histo for DF and SF (pt,eta)
 //char hname[30];
 for (int i=0; i<2; i++) {
                                                                                                                 
   sprintf(hname,"h_el_os_pt_C3_%i",i);
   h_el_os_pt_C3[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_os_pt_C3[i]->Write();

   sprintf(hname,"h_el_os_eta_C3_%i",i);
   h_el_os_eta_C3[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_os_eta_C3[i]->Write();

   sprintf(hname,"h_mu_os_pt_C3_%i",i);
   h_mu_os_pt_C3[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_os_pt_C3[i]->Write();

   sprintf(hname,"h_mu_os_eta_C3_%i",i);
   h_mu_os_eta_C3[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_os_eta_C3[i]->Write();


   sprintf(hname,"h_elmu_os_pt_C3_%i",i);
   h_elmu_os_pt_C3[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_elmu_os_pt_C3[i]->Write();

   sprintf(hname,"h_elmu_os_eta_C3_%i",i);
   h_elmu_os_eta_C3[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_elmu_os_eta_C3[i]->Write();
   
 }

 // char hname[30]; 
 for (int i=0; i<2; i++) {
     
    
 
     
     // lepton jet dphi 
     sprintf(hname,"h_mu_jet_dphi_C3_%i",i);
     h_mu_jet_dphi_C3[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_jet_dphi_C3[i]->Write();
     
     sprintf(hname,"h_el_jet_dphi_C3_%i",i);
     h_el_jet_dphi_C3[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_jet_dphi_C3[i]->Write();

     sprintf(hname,"h_elmu_jet_dphi_C3_%i",i);
     h_elmu_jet_dphi_C3[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_elmu_jet_dphi_C3[i]->Write();    

     // lepton met dphi 
     sprintf(hname,"h_mu_etm_dphi_C3_%i",i);
     h_mu_etm_dphi_C3[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_etm_dphi_C3[i]->Write();
     
     sprintf(hname,"h_el_etm_dphi_C3_%i",i);
     h_el_etm_dphi_C3[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_etm_dphi_C3[i]->Write();

     sprintf(hname,"h_elmu_etm_dphi_C3_%i",i);
     h_elmu_etm_dphi_C3[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_elmu_etm_dphi_C3[i]->Write();     
     
  


 } 


 //***********************************************************************************************END C3
 //***************************************************************************************************

 //************************************************************************************************* START C4
 //************************************************************************************************

h_etmiss_e_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_e_C4")));
 h_etmiss_e_C4->Write();

 h_etmiss_emu_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_emu_C4")));
 h_etmiss_emu_C4->Write();

 h_etmiss_mu_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_mu_C4")));
 h_etmiss_mu_C4->Write();

 h_jet_mt2_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_mt2_C4")));
 h_jet_mt2_C4->Write();



 //lepton dphi 
 h_mu_os_dphi_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_dphi_C4")));
 h_mu_os_dphi_C4->Write();
 h_el_os_dphi_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_dphi_C4")));
 h_el_os_dphi_C4->Write();
 h_elmu_os_dphi_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_dphi_C4")));
 h_elmu_os_dphi_C4->Write();

 //lepton dtheta

 h_mu_os_dtheta_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_dtheta_C4")));
 h_mu_os_dtheta_C4->Write();
 h_el_os_dtheta_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_dtheta_C4")));
 h_el_os_dtheta_C4->Write();
 h_elmu_os_dtheta_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_dtheta_C4")));
 h_elmu_os_dtheta_C4->Write();

 //lepton deta

 h_mu_os_deta_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_deta_C4")));
 h_mu_os_deta_C4->Write();
 h_el_os_deta_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_deta_C4")));
 h_el_os_deta_C4->Write();
 h_elmu_os_deta_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_deta_C4")));
 h_elmu_os_deta_C4->Write();





 // mt2 
 h_mu_mt2_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_mt2_C4")));
 h_mu_mt2_C4->Write();
 
 h_el_mt2_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_mt2_C4")));
 h_el_mt2_C4->Write();
  h_elmu_mt2_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_mt2_C4")));
 h_elmu_mt2_C4->Write();




 // minv 
 h_mu_os_minv_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_minv_C4")));
 h_mu_os_minv_C4->Write();
 
 h_el_os_minv_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_minv_C4")));
 h_el_os_minv_C4->Write();
 

 
 h_elmu_os_minv_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_minv_C4")));
 h_elmu_os_minv_C4->Write(); 
 
// asym 
 h_mu_os_asym_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_mu_os_asym_C4")));
 h_mu_os_asym_C4->Write();
 

 h_el_os_asym_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_el_os_asym_C4")));
 h_el_os_asym_C4->Write();
 
 h_elmu_os_asym_C4 = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_elmu_os_asym_C4")));
 h_elmu_os_asym_C4->Write(); 

 
 //histo for DF and SF (pt,eta)
 //char hname[30];
 for (int i=0; i<2; i++) {
                                                                                                                 
   sprintf(hname,"h_el_os_pt_C4_%i",i);
   h_el_os_pt_C4[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_os_pt_C4[i]->Write();

   sprintf(hname,"h_el_os_eta_C4_%i",i);
   h_el_os_eta_C4[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_el_os_eta_C4[i]->Write();

   sprintf(hname,"h_mu_os_pt_C4_%i",i);
   h_mu_os_pt_C4[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_os_pt_C4[i]->Write();

   sprintf(hname,"h_mu_os_eta_C4_%i",i);
   h_mu_os_eta_C4[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_mu_os_eta_C4[i]->Write();


   sprintf(hname,"h_elmu_os_pt_C4_%i",i);
   h_elmu_os_pt_C4[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_elmu_os_pt_C4[i]->Write();

   sprintf(hname,"h_elmu_os_eta_C4_%i",i);
   h_elmu_os_eta_C4[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
   h_elmu_os_eta_C4[i]->Write();
   
 }

 // char hname[30]; 
 for (int i=0; i<2; i++) {
     
      
     // lepton jet dphi 
     sprintf(hname,"h_mu_jet_dphi_C4_%i",i);
     h_mu_jet_dphi_C4[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_jet_dphi_C4[i]->Write();
     
     sprintf(hname,"h_el_jet_dphi_C4_%i",i);
     h_el_jet_dphi_C4[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_jet_dphi_C4[i]->Write();

     sprintf(hname,"h_elmu_jet_dphi_C4_%i",i);
     h_elmu_jet_dphi_C4[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_elmu_jet_dphi_C4[i]->Write();    

     // lepton met dphi 
     sprintf(hname,"h_mu_etm_dphi_C4_%i",i);
     h_mu_etm_dphi_C4[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_mu_etm_dphi_C4[i]->Write();
     
     sprintf(hname,"h_el_etm_dphi_C4_%i",i);
     h_el_etm_dphi_C4[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_el_etm_dphi_C4[i]->Write();

     sprintf(hname,"h_elmu_etm_dphi_C4_%i",i);
     h_elmu_etm_dphi_C4[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
     h_elmu_etm_dphi_C4[i]->Write();     
     
  


 } 



 //***********************************************************************************************END C4
 //***************************************************************************************************


 
 
 fileOut->Close();



}




void read_d4pd::BookHistos()
{


  h_event_yield = new TH1F("h_event_yield", "h_event_yield",90000,100000,190000); 
  fOutput->Add(h_event_yield); 
  
  h_totevent = new TH1D("h_totevent", "h_totevent",90000,100000,190000); 
  fOutput->Add(h_totevent); 

  h_btag_weight=new TH1F("h_btag_weight","h_btag_weight",500,0,50);
  fOutput->Add(h_btag_weight);

  h_el_trigg_freq = new TH1D("h_el_trigg_freq", "h_el_trigg_freq",20,0,20);
 fOutput->Add(h_el_trigg_freq);

 h_mu_trigg_freq = new TH1D("h_mu_trigg_freq", "h_mu_trigg_freq",20,0,20);
 fOutput->Add(h_mu_trigg_freq);

  h_cutFlow = new TH1D("h_cutFlow", "Cut flow", 81, 0, 80);
  fOutput->Add(h_cutFlow);
  
  h_cutFlow_w = new TH1D("h_cutFlow_w", "weighted Cut flow", 81, 0, 80);
  fOutput->Add(h_cutFlow_w);  

  h_jet_n = new TH1F("h_jet_n", "h_jet_n",20, 0, 19); 
  fOutput->Add(h_jet_n);
  h_jet_n->Sumw2();

  h_jet_ht = new TH1F("h_jet_ht", "h_jet_ht",200, 0, 800); 
  fOutput->Add(h_jet_ht);
  h_jet_ht->Sumw2();
  
  h_jet_jet_dphi= new TH1F("h_jet_jet_dphi", "h_jet_jet_dphi",100, 0, TMath::Pi()); 
  fOutput->Add(h_jet_jet_dphi); 
  h_jet_jet_dphi->Sumw2();

  h_mu_n = new TH1F("h_mu_n", "h_mu_n",10, 0, 9); 
  fOutput->Add(h_mu_n);
  h_mu_n->Sumw2();
  
  h_el_n = new TH1F("h_el_n", "h_el_n",10, 0, 9); 
  fOutput->Add(h_el_n);
  h_el_n->Sumw2();


  h_elmu_n = new TH1F("h_elmu_n", "h_elmu_n",10, 0, 9);
  fOutput->Add(h_elmu_n);
  h_elmu_n->Sumw2();
  //
  h_jet_or_n = new TH1F("h_jet_or_n", "h_jet_or_n",20, 0, 19); 
  fOutput->Add(h_jet_or_n);
  h_jet_or_n->Sumw2();

  h_jet_b_or_n = new TH1F("h_jet_b_or_n", "h_jet_b_or_n",20, 0, 19); 
  fOutput->Add(h_jet_b_or_n);
  h_jet_b_or_n->Sumw2();  

  h_mu_or_n = new TH1F("h_mu_or_n", "h_mu_or_n",20, 0, 19); 
  fOutput->Add(h_mu_or_n);
  h_mu_or_n->Sumw2();
  
  h_el_or_n = new TH1F("h_el_or_n", "h_el_or_n",20, 0, 19); 
  fOutput->Add(h_el_or_n);
  h_el_or_n->Sumw2();
 
  h_vertex_n = new TH1F("h_vertex_n", "h_vertex_n",20, 0, 19); 
  fOutput->Add(h_vertex_n);
  h_vertex_n->Sumw2();
  
  h_etmiss_mod = new TH1F("h_etmiss_mod", "h_etmiss_mod",50, 0, 500); 
  fOutput->Add(h_etmiss_mod);
  h_etmiss_mod->Sumw2();

  h_etmiss = new TH1F("h_etmiss", "h_etmiss",50, 0, 500); 
  fOutput->Add(h_etmiss); 
  h_etmiss->Sumw2();
  
  h_etmiss_phi = new TH1F("h_etmiss_phi", "h_etmiss_phi",50, -TMath::Pi(), TMath::Pi()); 
  fOutput->Add(h_etmiss_phi);
  h_etmiss_phi->Sumw2();
  
  
  h_etmiss_e = new TH1F("h_etmiss_e", "h_etmiss_e",50, 0, 500);
  fOutput->Add(h_etmiss_e);
  h_etmiss_e->Sumw2();
  
  h_etmiss_emu = new TH1F("h_etmiss_emu", "h_etmiss_emu",50, 0, 500);
  fOutput->Add(h_etmiss_emu);
  h_etmiss_emu->Sumw2();
  
  h_etmiss_mu = new TH1F("h_etmiss_mu", "h_etmiss_mu",50, 0, 500);
  fOutput->Add(h_etmiss_mu);
  h_etmiss_mu->Sumw2();
  
  // dtheta between two leptons

  h_mu_os_dtheta= new TH1F("h_mu_os_dtheta", "h_mu_os_dtheta",50, 0, TMath::Pi());
  fOutput->Add(h_mu_os_dtheta);
  h_mu_os_dtheta->Sumw2();

  h_el_os_dtheta= new TH1F("h_el_os_dtheta", "h_el_os_dtheta",50, 0, TMath::Pi());
  fOutput->Add(h_el_os_dtheta);
  h_el_os_dtheta->Sumw2();

  h_elmu_os_dtheta= new TH1F("h_elmu_os_dtheta", "h_elmu_os_dtheta",50, 0, TMath::Pi());
  fOutput->Add(h_elmu_os_dtheta);
  h_elmu_os_dtheta->Sumw2();
 

  h_mu_os_deta= new TH1F("h_mu_os_deta", "h_mu_os_deta",100, 0,6 );
  fOutput->Add(h_mu_os_deta);
  h_mu_os_deta->Sumw2();

  h_el_os_deta= new TH1F("h_el_os_deta", "h_el_os_deta",100, 0, 6);
  fOutput->Add(h_el_os_deta);
  h_el_os_deta->Sumw2();

  h_elmu_os_deta= new TH1F("h_elmu_os_deta", "h_elmu_os_deta",100, 0, 6);
  fOutput->Add(h_elmu_os_deta);
  h_elmu_os_deta->Sumw2();


  //Pbll only

  h_mu_os_Pbll= new TH1F("h_mu_os_Pbll", "h_mu_os_Pbll",50, 0,500 );
  fOutput->Add(h_mu_os_Pbll);
  h_mu_os_Pbll->Sumw2();

  h_el_os_Pbll= new TH1F("h_el_os_Pbll", "h_el_os_Pbll",50, 0, 500);
  fOutput->Add(h_el_os_Pbll);
  h_el_os_Pbll->Sumw2();

  h_elmu_os_Pbll= new TH1F("h_elmu_os_Pbll", "h_elmu_os_Pbll",50, 0, 500);
  fOutput->Add(h_elmu_os_Pbll);
  h_elmu_os_Pbll->Sumw2();



  h_el_os_DPhiMetPbll= new TH1F("h_el_os_DPhiMetPbll", "h_mu_os_DPhiMetPbll",50, 0,5 );
  fOutput->Add(h_el_os_DPhiMetPbll);
  h_el_os_DPhiMetPbll->Sumw2();
  h_mu_os_DPhiMetPbll= new TH1F("h_mu_os_DPhiMetPbll", "h_mu_os_DPhiMetPbll",50, 0,5 );
  fOutput->Add(h_mu_os_DPhiMetPbll);
  h_mu_os_DPhiMetPbll->Sumw2();
  h_elmu_os_DPhiMetPbll= new TH1F("h_elmu_os_DPhiMetPbll", "h_elmu_os_DPhiMetPbll",50, 0,5 );
  fOutput->Add(h_elmu_os_DPhiMetPbll);
  h_elmu_os_DPhiMetPbll->Sumw2();



  h_el_os_DRll= new TH1F("h_el_os_DRll", "h_el_os_DRll",50, 0,5 );
  fOutput->Add(h_el_os_DRll);
  h_el_os_DRll->Sumw2();
  h_mu_os_DRll= new TH1F("h_mu_os_DRll", "h_mu_os_DRll",50, 0,5 );
  fOutput->Add(h_mu_os_DRll);
  h_mu_os_DRll->Sumw2();
  h_elmu_os_DRll= new TH1F("h_elmu_os_DRll", "h_elmu_os_DRll",50, 0,5 );
  fOutput->Add(h_elmu_os_DRll);
  h_elmu_os_DRll->Sumw2();


  h_el_os_PtLep= new TH1F("h_el_os_PtLep", "h_el_os_PtLep",50, 0,500 );
  fOutput->Add(h_el_os_PtLep);
  h_el_os_PtLep->Sumw2();
  h_mu_os_PtLep= new TH1F("h_mu_os_PtLep", "h_mu_os_PtLep",50, 0,500 );
  fOutput->Add(h_mu_os_PtLep);
  h_mu_os_PtLep->Sumw2();
  h_elmu_os_PtLep= new TH1F("h_elmu_os_PtLep", "h_elmu_os_PtLep",50, 0,500 );
  fOutput->Add(h_elmu_os_PtLep);
  h_elmu_os_PtLep->Sumw2();


  h_el_os_PtRatio= new TH1F("h_el_os_PtRatio", "h_el_os_PtRatio",50, 0,500 );
  fOutput->Add(h_el_os_PtRatio);
  h_el_os_PtRatio->Sumw2();
  h_mu_os_PtRatio= new TH1F("h_mu_os_PtRatio", "h_mu_os_PtRatio",50, 0,500 );
  fOutput->Add(h_mu_os_PtRatio);
  h_mu_os_PtRatio->Sumw2();
  h_elmu_os_PtRatio= new TH1F("h_elmu_os_PtRatio", "h_elmu_os_PtRatio",50, 0,500 );
  fOutput->Add(h_elmu_os_PtRatio);
  h_elmu_os_PtRatio->Sumw2();



  //hadronic mt2



  h_mu_os_had_mt2= new TH1F("h_mu_os_had_mt2", "h_mu_os_had_mt2",50, 0,500 );
  fOutput->Add(h_mu_os_had_mt2);
  h_mu_os_had_mt2->Sumw2();

  h_el_os_had_mt2= new TH1F("h_el_os_had_mt2", "h_el_os_had_mt2",50, 0, 500);
  fOutput->Add(h_el_os_had_mt2);
  h_el_os_had_mt2->Sumw2();

  h_elmu_os_had_mt2= new TH1F("h_elmu_os_had_mt2", "h_elmu_os_had_mt2",50, 0, 500);
  fOutput->Add(h_elmu_os_had_mt2);
  h_elmu_os_had_mt2->Sumw2();





  //********************************************************************************************************
  //*******************************************************************************************************  start C1
  //plots for
 h_mu_mt2_C1 = new TH1F("h_mu_mt2_C1", "h_mu_mt2_C1",60, 0, 600); 
  fOutput->Add(h_mu_mt2_C1); 
  h_mu_mt2_C1->Sumw2();

  h_el_mt2_C1 = new TH1F("h_el_mt2_C1", "h_el_mt2_C1",60, 0, 600); 
  fOutput->Add(h_el_mt2_C1); 
  h_el_mt2_C1->Sumw2();
  
  h_elmu_mt2_C1 = new TH1F("h_elmu_mt2_C1", "h_elmu_mt2_C1",60, 0, 600);  //change name
  fOutput->Add(h_elmu_mt2_C1); 
  h_elmu_mt2_C1->Sumw2();


  h_jet_mt2_C1 = new TH1F("h_jet_mt2_C1", "h_jet_mt2_C1",50, 0, 500); 
  fOutput->Add(h_jet_mt2_C1); 
  h_jet_mt2_C1->Sumw2();

  
  // invariant mass 
 
  
  h_mu_os_minv_C1 = new TH1F("h_mu_os_minv_C1", "h_mu_os_minv_C1",50, 0, 500); 
  fOutput->Add(h_mu_os_minv_C1);
  h_mu_os_minv_C1->Sumw2();

 
  
  h_el_os_minv_C1 = new TH1F("h_el_os_minv_C1", "h_el_os_minv_C1",50, 0, 500); 
  fOutput->Add(h_el_os_minv_C1);
  h_el_os_minv_C1->Sumw2();
  
 
  h_elmu_os_minv_C1 = new TH1F("h_elmu_os_minv_C1", "h_elmu_os_minv_C1",50, 0, 500); //change name 
  fOutput->Add(h_elmu_os_minv_C1);  
  h_elmu_os_minv_C1->Sumw2(); 
  
  // asymmetry 
 
  
  h_mu_os_asym_C1 = new TH1F("h_mu_os_asym_C1", "h_mu_os_asym_C1",50, -3,3); 
  fOutput->Add(h_mu_os_asym_C1);
  h_mu_os_asym_C1->Sumw2();
  

  
  h_el_os_asym_C1 = new TH1F("h_el_os_asym_C1", "h_el_os_asym_C1",50, -3,3); 
  fOutput->Add(h_el_os_asym_C1);
  h_el_os_asym_C1->Sumw2();
  
 
 
  h_elmu_os_asym_C1 = new TH1F("h_elmu_os_asym_C1", "h_elmu_os_asym_C1",50, -3,3); 
  fOutput->Add(h_elmu_os_asym_C1);  
  h_elmu_os_asym_C1->Sumw2();



 // dphi between two leptons 
  h_mu_os_dphi_C1= new TH1F("h_mu_os_dphi_C1", "h_mu_os_dphi_C1",50, 0, TMath::Pi()); 
  fOutput->Add(h_mu_os_dphi_C1); 
  h_mu_os_dphi_C1->Sumw2();

  h_el_os_dphi_C1= new TH1F("h_el_os_dphi_C1", "h_el_os_dphi_C1",50, 0, TMath::Pi()); 
  fOutput->Add(h_el_os_dphi_C1); 
  h_el_os_dphi_C1->Sumw2();
  
  h_elmu_os_dphi_C1= new TH1F("h_elmu_os_dphi_C1", "h_elmu_os_dphi_C1",50, 0, TMath::Pi()); 
  fOutput->Add(h_elmu_os_dphi_C1); 
  h_elmu_os_dphi_C1->Sumw2();

  {
    char hname[30];
    for (int i=0; i<2; i++) {
      sprintf(hname,"h_el_os_pt_C1_%i",i);
      h_el_os_pt_C1[i]= new TH1F(hname,hname,50, 0, 500);
      fOutput->Add(h_el_os_pt_C1[i]);
      h_el_os_pt_C1[i]->Sumw2();


      sprintf(hname,"h_el_os_eta_C1_%i",i);
      h_el_os_eta_C1[i]= new TH1F(hname,hname,100, -3, 3);
      fOutput->Add(h_el_os_eta_C1[i]);
      h_el_os_eta_C1[i]->Sumw2();

      sprintf(hname,"h_mu_os_pt_C1_%i",i);
      h_mu_os_pt_C1[i]= new TH1F(hname,hname,50, 0, 500);
      fOutput->Add(h_mu_os_pt_C1[i]);
      h_mu_os_pt_C1[i]->Sumw2();


      sprintf(hname,"h_mu_os_eta_C1_%i",i);  
      h_mu_os_eta_C1[i]= new TH1F(hname,hname,100, -3, 3);
      fOutput->Add(h_mu_os_eta_C1[i]);
      h_mu_os_eta_C1[i]->Sumw2();

      sprintf(hname,"h_elmu_os_pt_C1_%i",i);
      h_elmu_os_pt_C1[i]= new TH1F(hname,hname,50, 0, 500);
      fOutput->Add(h_elmu_os_pt_C1[i]);
      h_elmu_os_pt_C1[i]->Sumw2();


      sprintf(hname,"h_elmu_os_eta_C1_%i",i);
      h_elmu_os_eta_C1[i]= new TH1F(hname,hname,100, -3, 3);
      fOutput->Add(h_elmu_os_eta_C1[i]);
      h_elmu_os_eta_C1[i]->Sumw2();

    }
  }



  {
    char hname[30];
    for (int i=0; i<2; i++) {
      sprintf(hname,"h_mu_jet_dphi_C1_%i",i);
      h_mu_jet_dphi_C1[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_mu_jet_dphi_C1[i]); 
      h_mu_jet_dphi_C1[i]->Sumw2();
      
      sprintf(hname,"h_el_jet_dphi_C1_%i",i);
      h_el_jet_dphi_C1[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_el_jet_dphi_C1[i]); 
      h_el_jet_dphi_C1[i]->Sumw2();
      
      sprintf(hname,"h_elmu_jet_dphi_C1_%i",i);
      h_elmu_jet_dphi_C1[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_elmu_jet_dphi_C1[i]); 
      h_elmu_jet_dphi_C1[i]->Sumw2();

      sprintf(hname,"h_mu_etm_dphi_C1_%i",i);
      h_mu_etm_dphi_C1[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_mu_etm_dphi_C1[i]); 
      h_mu_etm_dphi_C1[i]->Sumw2();
      
      sprintf(hname,"h_el_etm_dphi_C1_%i",i);
      h_el_etm_dphi_C1[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_el_etm_dphi_C1[i]); 
      h_el_etm_dphi_C1[i]->Sumw2();
      
      sprintf(hname,"h_elmu_etm_dphi_C1_%i",i);
      h_elmu_etm_dphi_C1[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_elmu_etm_dphi_C1[i]); 
      h_elmu_etm_dphi_C1[i]->Sumw2();
    }
  }

  //ETMISS


 h_etmiss_e_C1 = new TH1F("h_etmiss_e_C1", "h_etmiss_e_C1",50, 0, 500);
  fOutput->Add(h_etmiss_e_C1);
  h_etmiss_e_C1->Sumw2();
  
  h_etmiss_emu_C1 = new TH1F("h_etmiss_emu_C1", "h_etmiss_emu_C1",50, 0, 500);//change name
  fOutput->Add(h_etmiss_emu_C1);
  h_etmiss_emu_C1->Sumw2();
  
  h_etmiss_mu_C1 = new TH1F("h_etmiss_mu_C1", "h_etmiss_mu_C1",50, 0, 500);
  fOutput->Add(h_etmiss_mu_C1);
  h_etmiss_mu_C1->Sumw2();



// dtheta between two leptons

  h_mu_os_dtheta_C1= new TH1F("h_mu_os_dtheta_C1", "h_mu_os_dtheta_C1",50, 0, TMath::Pi());
  fOutput->Add(h_mu_os_dtheta_C1);
  h_mu_os_dtheta_C1->Sumw2();

  h_el_os_dtheta_C1= new TH1F("h_el_os_dtheta_C1", "h_el_os_dtheta_C1",50, 0, TMath::Pi());
  fOutput->Add(h_el_os_dtheta_C1);
  h_el_os_dtheta_C1->Sumw2();

  h_elmu_os_dtheta_C1= new TH1F("h_elmu_os_dtheta_C1", "h_elmu_os_dtheta_C1",50, 0, TMath::Pi()); //change name
  fOutput->Add(h_elmu_os_dtheta_C1);
  h_elmu_os_dtheta_C1->Sumw2();


  //deta two leptons



  h_mu_os_deta_C1= new TH1F("h_mu_os_deta_C1", "h_mu_os_deta_C1",100, 0,6 );
  fOutput->Add(h_mu_os_deta_C1);
  h_mu_os_deta_C1->Sumw2();

  h_el_os_deta_C1= new TH1F("h_el_os_deta_C1", "h_el_os_deta_C1",100, 0, 6);
  fOutput->Add(h_el_os_deta_C1);
  h_el_os_deta_C1->Sumw2();

  h_elmu_os_deta_C1= new TH1F("h_elmu_os_deta_C1", "h_elmu_os_deta_C1",100, 0, 6); //change name
  fOutput->Add(h_elmu_os_deta_C1);
  h_elmu_os_deta_C1->Sumw2();




  //*******************************************************************************************************************
  //*****************************************************************************************************************end C1

  //*****************************************************************************************************************start C2
  //******************************************************************************************************************

 h_mu_mt2_C2 = new TH1F("h_mu_mt2_C2", "h_mu_mt2_C2",60, 0, 600); 
  fOutput->Add(h_mu_mt2_C2); 
  h_mu_mt2_C2->Sumw2();

  h_el_mt2_C2 = new TH1F("h_el_mt2_C2", "h_el_mt2_C2",60, 0, 600); 
  fOutput->Add(h_el_mt2_C2); 
  h_el_mt2_C2->Sumw2();
  
  h_elmu_mt2_C2 = new TH1F("h_elmu_mt2_C2", "h_elmu_mt2_C2",60, 0, 600); 
  fOutput->Add(h_elmu_mt2_C2); 
  h_elmu_mt2_C2->Sumw2();


  h_jet_mt2_C2 = new TH1F("h_jet_mt2_C2", "h_jet_mt2_C2",50, 0, 500); 
  fOutput->Add(h_jet_mt2_C2); 
  h_jet_mt2_C2->Sumw2();


  
  // invariant mass 
 
  
  h_mu_os_minv_C2 = new TH1F("h_mu_os_minv_C2", "h_mu_os_minv_C2",50, 0, 500); 
  fOutput->Add(h_mu_os_minv_C2);
  h_mu_os_minv_C2->Sumw2();

 
  
  h_el_os_minv_C2 = new TH1F("h_el_os_minv_C2", "h_el_os_minv_C2",50, 0, 500); 
  fOutput->Add(h_el_os_minv_C2);
  h_el_os_minv_C2->Sumw2();
  
 
  h_elmu_os_minv_C2 = new TH1F("h_elmu_os_minv_C2", "h_elmu_os_minv_C2",50, 0, 500); 
  fOutput->Add(h_elmu_os_minv_C2);  
  h_elmu_os_minv_C2->Sumw2(); 
  
  // asymmetry 
 
  
  h_mu_os_asym_C2 = new TH1F("h_mu_os_asym_C2", "h_mu_os_asym_C2",50, -3,3); 
  fOutput->Add(h_mu_os_asym_C2);
  h_mu_os_asym_C2->Sumw2();
  

  
  h_el_os_asym_C2 = new TH1F("h_el_os_asym_C2", "h_el_os_asym_C2",50, -3,3); 
  fOutput->Add(h_el_os_asym_C2);
  h_el_os_asym_C2->Sumw2();
  
 
 
  h_elmu_os_asym_C2 = new TH1F("h_elmu_os_asym_C2", "h_elmu_os_asym_C2",50, -3,3); 
  fOutput->Add(h_elmu_os_asym_C2);  
  h_elmu_os_asym_C2->Sumw2();



 // dphi between two leptons 
  h_mu_os_dphi_C2= new TH1F("h_mu_os_dphi_C2", "h_mu_os_dphi_C2",50, 0, TMath::Pi()); 
  fOutput->Add(h_mu_os_dphi_C2); 
  h_mu_os_dphi_C2->Sumw2();

  h_el_os_dphi_C2= new TH1F("h_el_os_dphi_C2", "h_el_os_dphi_C2",50, 0, TMath::Pi()); 
  fOutput->Add(h_el_os_dphi_C2); 
  h_el_os_dphi_C2->Sumw2();
  
  h_elmu_os_dphi_C2= new TH1F("h_elmu_os_dphi_C2", "h_elmu_os_dphi_C2",50, 0, TMath::Pi()); 
  fOutput->Add(h_elmu_os_dphi_C2); 
  h_elmu_os_dphi_C2->Sumw2();

  {
    char hname[30];
    for (int i=0; i<2; i++) {
      sprintf(hname,"h_el_os_pt_C2_%i",i);
      h_el_os_pt_C2[i]= new TH1F(hname,hname,50, 0, 500);
      fOutput->Add(h_el_os_pt_C2[i]);
      h_el_os_pt_C2[i]->Sumw2();


      sprintf(hname,"h_el_os_eta_C2_%i",i);
      h_el_os_eta_C2[i]= new TH1F(hname,hname,100, -3, 3);
      fOutput->Add(h_el_os_eta_C2[i]);
      h_el_os_eta_C2[i]->Sumw2();

      sprintf(hname,"h_mu_os_pt_C2_%i",i);
      h_mu_os_pt_C2[i]= new TH1F(hname,hname,50, 0, 500);
      fOutput->Add(h_mu_os_pt_C2[i]);
      h_mu_os_pt_C2[i]->Sumw2();


      sprintf(hname,"h_mu_os_eta_C2_%i",i);
      h_mu_os_eta_C2[i]= new TH1F(hname,hname,100, -3, 3);
      fOutput->Add(h_mu_os_eta_C2[i]);
      h_mu_os_eta_C2[i]->Sumw2();

      sprintf(hname,"h_elmu_os_pt_C2_%i",i);
      h_elmu_os_pt_C2[i]= new TH1F(hname,hname,50, 0, 500);
      fOutput->Add(h_elmu_os_pt_C2[i]);
      h_elmu_os_pt_C2[i]->Sumw2();


      sprintf(hname,"h_elmu_os_eta_C2_%i",i);
      h_elmu_os_eta_C2[i]= new TH1F(hname,hname,100, -3, 3);
      fOutput->Add(h_elmu_os_eta_C2[i]);
      h_elmu_os_eta_C2[i]->Sumw2();

    }
  }



  {
    char hname[30];
    for (int i=0; i<2; i++) {
      sprintf(hname,"h_mu_jet_dphi_C2_%i",i);
      h_mu_jet_dphi_C2[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_mu_jet_dphi_C2[i]); 
      h_mu_jet_dphi_C2[i]->Sumw2();
      
      sprintf(hname,"h_el_jet_dphi_C2_%i",i);
      h_el_jet_dphi_C2[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_el_jet_dphi_C2[i]); 
      h_el_jet_dphi_C2[i]->Sumw2();
      
      sprintf(hname,"h_elmu_jet_dphi_C2_%i",i);
      h_elmu_jet_dphi_C2[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_elmu_jet_dphi_C2[i]); 
      h_elmu_jet_dphi_C2[i]->Sumw2();

      sprintf(hname,"h_mu_etm_dphi_C2_%i",i);
      h_mu_etm_dphi_C2[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_mu_etm_dphi_C2[i]); 
      h_mu_etm_dphi_C2[i]->Sumw2();
      
      sprintf(hname,"h_el_etm_dphi_C2_%i",i);
      h_el_etm_dphi_C2[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_el_etm_dphi_C2[i]); 
      h_el_etm_dphi_C2[i]->Sumw2();
      
      sprintf(hname,"h_elmu_etm_dphi_C2_%i",i);
      h_elmu_etm_dphi_C2[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_elmu_etm_dphi_C2[i]); 
      h_elmu_etm_dphi_C2[i]->Sumw2();
    }
  }

  //ETMISS


 h_etmiss_e_C2 = new TH1F("h_etmiss_e_C2", "h_etmiss_e_C2",50, 0, 500);
  fOutput->Add(h_etmiss_e_C2);
  h_etmiss_e_C2->Sumw2();
  
  h_etmiss_emu_C2 = new TH1F("h_etmiss_emu_C2", "h_etmiss_emu_C2",50, 0, 500);
  fOutput->Add(h_etmiss_emu_C2);
  h_etmiss_emu_C2->Sumw2();
  
  h_etmiss_mu_C2 = new TH1F("h_etmiss_mu_C2", "h_etmiss_mu_C2",50, 0, 500);
  fOutput->Add(h_etmiss_mu_C2);
  h_etmiss_mu_C2->Sumw2();



// dtheta between two leptons

  h_mu_os_dtheta_C2= new TH1F("h_mu_os_dtheta_C2", "h_mu_os_dtheta_C2",50, 0, TMath::Pi());
  fOutput->Add(h_mu_os_dtheta_C2);
  h_mu_os_dtheta_C2->Sumw2();

  h_el_os_dtheta_C2= new TH1F("h_el_os_dtheta_C2", "h_el_os_dtheta_C2",50, 0, TMath::Pi());
  fOutput->Add(h_el_os_dtheta_C2);
  h_el_os_dtheta_C2->Sumw2();

  h_elmu_os_dtheta_C2= new TH1F("h_elmu_os_dtheta_C2", "h_elmu_os_dtheta_C2",50, 0, TMath::Pi());
  fOutput->Add(h_elmu_os_dtheta_C2);
  h_elmu_os_dtheta_C2->Sumw2();



  //deta two leptons                                                                                                                                                                                                                         
  h_mu_os_deta_C2= new TH1F("h_mu_os_deta_C2", "h_mu_os_deta_C2",100, 0,6 );
  fOutput->Add(h_mu_os_deta_C2);
  h_mu_os_deta_C2->Sumw2();

  h_el_os_deta_C2= new TH1F("h_el_os_deta_C2", "h_el_os_deta_C2",100, 0, 6);
  fOutput->Add(h_el_os_deta_C2);
  h_el_os_deta_C2->Sumw2();

  h_elmu_os_deta_C2= new TH1F("h_elmu_os_deta_C2", "h_elmu_os_deta_C2",100, 0, 6);
  fOutput->Add(h_elmu_os_deta_C2);
  h_elmu_os_deta_C2->Sumw2();


                                                                                                                                                                                                               


  //************************************************************************************************************************ENDC2
  //************************************************************************************************************************
  //*********************************************************************************************************************START C3
  //*************************************************************************************************************************

 h_mu_mt2_C3 = new TH1F("h_mu_mt2_C3", "h_mu_mt2_C3",60, 0, 600); 
  fOutput->Add(h_mu_mt2_C3); 
  h_mu_mt2_C3->Sumw2();

  h_el_mt2_C3 = new TH1F("h_el_mt2_C3", "h_el_mt2_C3",60, 0, 600); 
  fOutput->Add(h_el_mt2_C3); 
  h_el_mt2_C3->Sumw2();
  
  h_elmu_mt2_C3 = new TH1F("h_elmu_mt2_C3", "h_elmu_mt2_C3",60, 0, 600); 
  fOutput->Add(h_elmu_mt2_C3); 
  h_elmu_mt2_C3->Sumw2();



  h_jet_mt2_C3 = new TH1F("h_jet_mt2_C3", "h_jet_mt2_C3",50, 0, 500); 
  fOutput->Add(h_jet_mt2_C3); 
  h_jet_mt2_C3->Sumw2();

  
  // invariant mass 
  
  h_mu_os_minv_C3 = new TH1F("h_mu_os_minv_C3", "h_mu_os_minv_C3",50, 0, 500); 
  fOutput->Add(h_mu_os_minv_C3);
  h_mu_os_minv_C3->Sumw2();
  
  h_el_os_minv_C3 = new TH1F("h_el_os_minv_C3", "h_el_os_minv_C3",50, 0, 500); 
  fOutput->Add(h_el_os_minv_C3);
  h_el_os_minv_C3->Sumw2();
  
    
  h_elmu_os_minv_C3 = new TH1F("h_elmu_os_minv_C3", "h_elmu_os_minv_C3",50, 0, 500); 
  fOutput->Add(h_elmu_os_minv_C3);  
  h_elmu_os_minv_C3->Sumw2(); 
  
  // asymmetry  
  
  h_mu_os_asym_C3 = new TH1F("h_mu_os_asym_C3", "h_mu_os_asym_C3",50, -3,3); 
  fOutput->Add(h_mu_os_asym_C3);
  h_mu_os_asym_C3->Sumw2();
  
  h_el_os_asym_C3 = new TH1F("h_el_os_asym_C3", "h_el_os_asym_C3",50, -3,3); 
  fOutput->Add(h_el_os_asym_C3);
  h_el_os_asym_C3->Sumw2();
 
  h_elmu_os_asym_C3 = new TH1F("h_elmu_os_asym_C3", "h_elmu_os_asym_C3",50, -3,3); 
  fOutput->Add(h_elmu_os_asym_C3);  
  h_elmu_os_asym_C3->Sumw2();



 // dphi between two leptons 
  h_mu_os_dphi_C3= new TH1F("h_mu_os_dphi_C3", "h_mu_os_dphi_C3",50, 0, TMath::Pi()); 
  fOutput->Add(h_mu_os_dphi_C3); 
  h_mu_os_dphi_C3->Sumw2();

  h_el_os_dphi_C3= new TH1F("h_el_os_dphi_C3", "h_el_os_dphi_C3",50, 0, TMath::Pi()); 
  fOutput->Add(h_el_os_dphi_C3); 
  h_el_os_dphi_C3->Sumw2();
  
  h_elmu_os_dphi_C3= new TH1F("h_elmu_os_dphi_C3", "h_elmu_os_dphi_C3",50, 0, TMath::Pi()); 
  fOutput->Add(h_elmu_os_dphi_C3); 
  h_elmu_os_dphi_C3->Sumw2();

  {
    char hname[30];
    for (int i=0; i<2; i++) {
      sprintf(hname,"h_el_os_pt_C3_%i",i);
      h_el_os_pt_C3[i]= new TH1F(hname,hname,50, 0, 500);
      fOutput->Add(h_el_os_pt_C3[i]);
      h_el_os_pt_C3[i]->Sumw2();


      sprintf(hname,"h_el_os_eta_C3_%i",i);
      h_el_os_eta_C3[i]= new TH1F(hname,hname,100, -3, 3);
      fOutput->Add(h_el_os_eta_C3[i]);
      h_el_os_eta_C3[i]->Sumw2();

      sprintf(hname,"h_mu_os_pt_C3_%i",i);
      h_mu_os_pt_C3[i]= new TH1F(hname,hname,50, 0, 500);
      fOutput->Add(h_mu_os_pt_C3[i]);
      h_mu_os_pt_C3[i]->Sumw2();


      sprintf(hname,"h_mu_os_eta_C3_%i",i);
      h_mu_os_eta_C3[i]= new TH1F(hname,hname,100, -3, 3);
      fOutput->Add(h_mu_os_eta_C3[i]);
      h_mu_os_eta_C3[i]->Sumw2();

      sprintf(hname,"h_elmu_os_pt_C3_%i",i);
      h_elmu_os_pt_C3[i]= new TH1F(hname,hname,50, 0, 500);
      fOutput->Add(h_elmu_os_pt_C3[i]);
      h_elmu_os_pt_C3[i]->Sumw2();


      sprintf(hname,"h_elmu_os_eta_C3_%i",i);
      h_elmu_os_eta_C3[i]= new TH1F(hname,hname,100, -3, 3);
      fOutput->Add(h_elmu_os_eta_C3[i]);
      h_elmu_os_eta_C3[i]->Sumw2();

    }
  }



  {
    char hname[30];
    for (int i=0; i<2; i++) {
      sprintf(hname,"h_mu_jet_dphi_C3_%i",i);
      h_mu_jet_dphi_C3[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_mu_jet_dphi_C3[i]); 
      h_mu_jet_dphi_C3[i]->Sumw2();
      
      sprintf(hname,"h_el_jet_dphi_C3_%i",i);
      h_el_jet_dphi_C3[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_el_jet_dphi_C3[i]); 
      h_el_jet_dphi_C3[i]->Sumw2();
      
      sprintf(hname,"h_elmu_jet_dphi_C3_%i",i);
      h_elmu_jet_dphi_C3[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_elmu_jet_dphi_C3[i]); 
      h_elmu_jet_dphi_C3[i]->Sumw2();

      sprintf(hname,"h_mu_etm_dphi_C3_%i",i);
      h_mu_etm_dphi_C3[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_mu_etm_dphi_C3[i]); 
      h_mu_etm_dphi_C3[i]->Sumw2();
      
      sprintf(hname,"h_el_etm_dphi_C3_%i",i);
      h_el_etm_dphi_C3[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_el_etm_dphi_C3[i]); 
      h_el_etm_dphi_C3[i]->Sumw2();
      
      sprintf(hname,"h_elmu_etm_dphi_C3_%i",i);
      h_elmu_etm_dphi_C3[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_elmu_etm_dphi_C3[i]); 
      h_elmu_etm_dphi_C3[i]->Sumw2();
    }
  }

  //ETMISS


 h_etmiss_e_C3 = new TH1F("h_etmiss_e_C3", "h_etmiss_e_C3",50, 0, 500);
  fOutput->Add(h_etmiss_e_C3);
  h_etmiss_e_C3->Sumw2();
  
  h_etmiss_emu_C3 = new TH1F("h_etmiss_emu_C3", "h_etmiss_emu_C3",50, 0, 500);
  fOutput->Add(h_etmiss_emu_C3);
  h_etmiss_emu_C3->Sumw2();
  
  h_etmiss_mu_C3 = new TH1F("h_etmiss_mu_C3", "h_etmiss_mu_C3",50, 0, 500);
  fOutput->Add(h_etmiss_mu_C3);
  h_etmiss_mu_C3->Sumw2();



// dtheta between two leptons

  h_mu_os_dtheta_C3= new TH1F("h_mu_os_dtheta_C3", "h_mu_os_dtheta_C3",50, 0, TMath::Pi());
  fOutput->Add(h_mu_os_dtheta_C3);
  h_mu_os_dtheta_C3->Sumw2();

  h_el_os_dtheta_C3= new TH1F("h_el_os_dtheta_C3", "h_el_os_dtheta_C3",50, 0, TMath::Pi());
  fOutput->Add(h_el_os_dtheta_C3);
  h_el_os_dtheta_C3->Sumw2();

  h_elmu_os_dtheta_C3= new TH1F("h_elmu_os_dtheta_C3", "h_elmu_os_dtheta_C3",50, 0, TMath::Pi());
  fOutput->Add(h_elmu_os_dtheta_C3);
  h_elmu_os_dtheta_C3->Sumw2();

  //deta (teo leptons)
  h_mu_os_deta_C3= new TH1F("h_mu_os_deta_C3", "h_mu_os_deta_C3",100, 0,6 );
  fOutput->Add(h_mu_os_deta_C3);
  h_mu_os_deta_C3->Sumw2();

  h_el_os_deta_C3= new TH1F("h_el_os_deta_C3", "h_el_os_deta_C3",100, 0, 6);
  fOutput->Add(h_el_os_deta_C3);
  h_el_os_deta_C3->Sumw2();

  h_elmu_os_deta_C3= new TH1F("h_elmu_os_deta_C3", "h_elmu_os_deta_C3",100, 0, 6);
  fOutput->Add(h_elmu_os_deta_C3);
  h_elmu_os_deta_C3->Sumw2();




  //****************************************************************************************************************************ENDC3
  //*****************************************************************************************************************************
  //**************************************************************************************************************************START C4
  //***************************************************************************************************************************

h_mu_mt2_C4 = new TH1F("h_mu_mt2_C4", "h_mu_mt2_C4",60, 0, 600); 
  fOutput->Add(h_mu_mt2_C4); 
  h_mu_mt2_C4->Sumw2();

  h_el_mt2_C4 = new TH1F("h_el_mt2_C4", "h_el_mt2_C4",60, 0, 600); 
  fOutput->Add(h_el_mt2_C4); 
  h_el_mt2_C4->Sumw2();
  
  h_elmu_mt2_C4 = new TH1F("h_elmu_mt2_C4", "h_elmu_mt2_C4",60, 0, 600); 
  fOutput->Add(h_elmu_mt2_C4); 
  h_elmu_mt2_C4->Sumw2();


  h_jet_mt2_C4 = new TH1F("h_jet_mt2_C4", "h_jet_mt2_C4",50, 0, 500); 
  fOutput->Add(h_jet_mt2_C4); 
  h_jet_mt2_C4->Sumw2();

 
  
  // invariant mass 
  
  
  h_mu_os_minv_C4 = new TH1F("h_mu_os_minv_C4", "h_mu_os_minv_C4",50, 0, 500); 
  fOutput->Add(h_mu_os_minv_C4);
  h_mu_os_minv_C4->Sumw2();
  
  h_el_os_minv_C4 = new TH1F("h_el_os_minv_C4", "h_el_os_minv_C4",50, 0, 500); 
  fOutput->Add(h_el_os_minv_C4);
  h_el_os_minv_C4->Sumw2();
  
  
  h_elmu_os_minv_C4 = new TH1F("h_elmu_os_minv_C4", "h_elmu_os_minv_C4",50, 0, 500); 
  fOutput->Add(h_elmu_os_minv_C4);  
  h_elmu_os_minv_C4->Sumw2(); 
  
  // asymmetry  
  
  h_mu_os_asym_C4 = new TH1F("h_mu_os_asym_C4", "h_mu_os_asym_C4",50, -3,3); 
  fOutput->Add(h_mu_os_asym_C4);
  h_mu_os_asym_C4->Sumw2();
  
  
  h_el_os_asym_C4 = new TH1F("h_el_os_asym_C4", "h_el_os_asym_C4",50, -3,3); 
  fOutput->Add(h_el_os_asym_C4);
  h_el_os_asym_C4->Sumw2();
  
 
  h_elmu_os_asym_C4 = new TH1F("h_elmu_os_asym_C4", "h_elmu_os_asym_C4",50, -3,3); 
  fOutput->Add(h_elmu_os_asym_C4);  
  h_elmu_os_asym_C4->Sumw2();



 // dphi between two leptons 
  h_mu_os_dphi_C4= new TH1F("h_mu_os_dphi_C4", "h_mu_os_dphi_C4",50, 0, TMath::Pi()); 
  fOutput->Add(h_mu_os_dphi_C4); 
  h_mu_os_dphi_C4->Sumw2();

  h_el_os_dphi_C4= new TH1F("h_el_os_dphi_C4", "h_el_os_dphi_C4",50, 0, TMath::Pi()); 
  fOutput->Add(h_el_os_dphi_C4); 
  h_el_os_dphi_C4->Sumw2();
  
  h_elmu_os_dphi_C4= new TH1F("h_elmu_os_dphi_C4", "h_elmu_os_dphi_C4",50, 0, TMath::Pi()); 
  fOutput->Add(h_elmu_os_dphi_C4); 
  h_elmu_os_dphi_C4->Sumw2();

  {
    char hname[30];
    for (int i=0; i<2; i++) {
      sprintf(hname,"h_el_os_pt_C4_%i",i);
      h_el_os_pt_C4[i]= new TH1F(hname,hname,50, 0, 500);
      fOutput->Add(h_el_os_pt_C4[i]);
      h_el_os_pt_C4[i]->Sumw2();


      sprintf(hname,"h_el_os_eta_C4_%i",i);
      h_el_os_eta_C4[i]= new TH1F(hname,hname,100, -3, 3);
      fOutput->Add(h_el_os_eta_C4[i]);
      h_el_os_eta_C4[i]->Sumw2();

      sprintf(hname,"h_mu_os_pt_C4_%i",i);
      h_mu_os_pt_C4[i]= new TH1F(hname,hname,50, 0, 500);
      fOutput->Add(h_mu_os_pt_C4[i]);
      h_mu_os_pt_C4[i]->Sumw2();


      sprintf(hname,"h_mu_os_eta_C4_%i",i);
      h_mu_os_eta_C4[i]= new TH1F(hname,hname,100, -3, 3);
      fOutput->Add(h_mu_os_eta_C4[i]);
      h_mu_os_eta_C4[i]->Sumw2();

      sprintf(hname,"h_elmu_os_pt_C4_%i",i);
      h_elmu_os_pt_C4[i]= new TH1F(hname,hname,50, 0, 500);
      fOutput->Add(h_elmu_os_pt_C4[i]);
      h_elmu_os_pt_C4[i]->Sumw2();


      sprintf(hname,"h_elmu_os_eta_C4_%i",i);
      h_elmu_os_eta_C4[i]= new TH1F(hname,hname,100, -3, 3);
      fOutput->Add(h_elmu_os_eta_C4[i]);
      h_elmu_os_eta_C4[i]->Sumw2();

    }
  }



  {
    char hname[30];
    for (int i=0; i<2; i++) {
      sprintf(hname,"h_mu_jet_dphi_C4_%i",i);
      h_mu_jet_dphi_C4[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_mu_jet_dphi_C4[i]); 
      h_mu_jet_dphi_C4[i]->Sumw2();
      
      sprintf(hname,"h_el_jet_dphi_C4_%i",i);
      h_el_jet_dphi_C4[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_el_jet_dphi_C4[i]); 
      h_el_jet_dphi_C4[i]->Sumw2();
      
      sprintf(hname,"h_elmu_jet_dphi_C4_%i",i);
      h_elmu_jet_dphi_C4[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_elmu_jet_dphi_C4[i]); 
      h_elmu_jet_dphi_C4[i]->Sumw2();

      sprintf(hname,"h_mu_etm_dphi_C4_%i",i);
      h_mu_etm_dphi_C4[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_mu_etm_dphi_C4[i]); 
      h_mu_etm_dphi_C4[i]->Sumw2();
      
      sprintf(hname,"h_el_etm_dphi_C4_%i",i);
      h_el_etm_dphi_C4[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_el_etm_dphi_C4[i]); 
      h_el_etm_dphi_C4[i]->Sumw2();
      
      sprintf(hname,"h_elmu_etm_dphi_C4_%i",i);
      h_elmu_etm_dphi_C4[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_elmu_etm_dphi_C4[i]); 
      h_elmu_etm_dphi_C4[i]->Sumw2();
    }
  }

  //ETMISS


 h_etmiss_e_C4 = new TH1F("h_etmiss_e_C4", "h_etmiss_e_C4",50, 0, 500);
  fOutput->Add(h_etmiss_e_C4);
  h_etmiss_e_C4->Sumw2();
  
  h_etmiss_emu_C4 = new TH1F("h_etmiss_emu_C4", "h_etmiss_emu_C4",50, 0, 500);
  fOutput->Add(h_etmiss_emu_C4);
  h_etmiss_emu_C4->Sumw2();
  
  h_etmiss_mu_C4 = new TH1F("h_etmiss_mu_C4", "h_etmiss_mu_C4",50, 0, 500);
  fOutput->Add(h_etmiss_mu_C4);
  h_etmiss_mu_C4->Sumw2();



// dtheta between two leptons

  h_mu_os_dtheta_C4= new TH1F("h_mu_os_dtheta_C4", "h_mu_os_dtheta_C4",50, 0, TMath::Pi());
  fOutput->Add(h_mu_os_dtheta_C4);
  h_mu_os_dtheta_C4->Sumw2();

  h_el_os_dtheta_C4= new TH1F("h_el_os_dtheta_C4", "h_el_os_dtheta_C4",50, 0, TMath::Pi());
  fOutput->Add(h_el_os_dtheta_C4);
  h_el_os_dtheta_C4->Sumw2();

  h_elmu_os_dtheta_C4= new TH1F("h_elmu_os_dtheta_C4", "h_elmu_os_dtheta_C4",50, 0, TMath::Pi());
  fOutput->Add(h_elmu_os_dtheta_C4);
  h_elmu_os_dtheta_C4->Sumw2();

  //deta two leptons                                                                                                                                                                                                                         



  h_mu_os_deta_C4= new TH1F("h_mu_os_deta_C4", "h_mu_os_deta_C4",100, 0,6 );
  fOutput->Add(h_mu_os_deta_C4);
  h_mu_os_deta_C4->Sumw2();

  h_el_os_deta_C4= new TH1F("h_el_os_deta_C4", "h_el_os_deta_C4",100, 0, 6);
  fOutput->Add(h_el_os_deta_C4);
  h_el_os_deta_C4->Sumw2();

  h_elmu_os_deta_C4= new TH1F("h_elmu_os_deta_C4", "h_elmu_os_deta_C4",100, 0, 6);
  fOutput->Add(h_elmu_os_deta_C4);
  h_elmu_os_deta_C4->Sumw2();




  //**********************************************************************************************************************************
  //************************************************************************************************************************end C4

  // dphi between two leptons 
  h_mu_os_dphi= new TH1F("h_mu_os_dphi", "h_mu_os_dphi",50, 0, TMath::Pi()); 
  fOutput->Add(h_mu_os_dphi); 
  h_mu_os_dphi->Sumw2();

  h_el_os_dphi= new TH1F("h_el_os_dphi", "h_el_os_dphi",50, 0, TMath::Pi()); 
  fOutput->Add(h_el_os_dphi); 
  h_el_os_dphi->Sumw2();
  
  h_elmu_os_dphi= new TH1F("h_elmu_os_dphi", "h_elmu_os_dphi",50, 0, TMath::Pi()); 
  fOutput->Add(h_elmu_os_dphi); 
  h_elmu_os_dphi->Sumw2();



  {
    char hname[30];
    for (int i=0; i<2; i++) {
      sprintf(hname,"h_mu_jet_dphi_%i",i);
      h_mu_jet_dphi[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_mu_jet_dphi[i]); 
      h_mu_jet_dphi[i]->Sumw2();
      
      sprintf(hname,"h_el_jet_dphi_%i",i);
      h_el_jet_dphi[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_el_jet_dphi[i]); 
      h_el_jet_dphi[i]->Sumw2();
      
      sprintf(hname,"h_elmu_jet_dphi_%i",i);
      h_elmu_jet_dphi[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_elmu_jet_dphi[i]); 
      h_elmu_jet_dphi[i]->Sumw2();

      sprintf(hname,"h_mu_etm_dphi_%i",i);
      h_mu_etm_dphi[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_mu_etm_dphi[i]); 
      h_mu_etm_dphi[i]->Sumw2();
      
      sprintf(hname,"h_el_etm_dphi_%i",i);
      h_el_etm_dphi[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_el_etm_dphi[i]); 
      h_el_etm_dphi[i]->Sumw2();
      
      sprintf(hname,"h_elmu_etm_dphi_%i",i);
      h_elmu_etm_dphi[i]= new TH1F(hname,hname,50, 0, TMath::Pi()); 
      fOutput->Add(h_elmu_etm_dphi[i]); 
      h_elmu_etm_dphi[i]->Sumw2();
    }
  }

 

  // mt2 
  h_mu_mt2 = new TH1F("h_mu_mt2", "h_mu_mt2",60, 0, 600); 
  fOutput->Add(h_mu_mt2); 
  h_mu_mt2->Sumw2();

  h_el_mt2 = new TH1F("h_el_mt2", "h_el_mt2",60, 0, 600); 
  fOutput->Add(h_el_mt2); 
  h_el_mt2->Sumw2();
  
  h_elmu_mt2 = new TH1F("h_elmu_mt2", "h_elmu_mt2",60, 0, 600); 
  fOutput->Add(h_elmu_mt2); 
  h_elmu_mt2->Sumw2();

  h_elmu_mt2_cor = new TH1F("h_elmu_mt2_cor", "h_elmu_mt2_cor",60, 0, 600); 
  fOutput->Add(h_elmu_mt2_cor); 
  h_elmu_mt2_cor->Sumw2();

  h_jet_mt2 = new TH1F("h_jet_mt2", "h_jet_mt2",50, 0, 500); 
  fOutput->Add(h_jet_mt2); 
  h_jet_mt2->Sumw2();

  h_jet_elmu_mt2 = new TH1F("h_jet_elmu_mt2", "h_jet_elmu_mt2",60, 0, 600); 
  fOutput->Add(h_jet_elmu_mt2); 
  h_jet_elmu_mt2->Sumw2();
  
  h_el_mt2_bl = new TH1F("h_el_mt2_bl", "h_el_mt2_bl",60, 0, 600);
  fOutput->Add(h_el_mt2_bl);
  h_el_mt2_bl->Sumw2();

  h_elmu_mt2_bl = new TH1F("h_elmu_mt2_bl", "h_elmu_mt2_bl",60, 0, 600);
  fOutput->Add(h_elmu_mt2_bl);
  h_elmu_mt2_bl->Sumw2();

  h_mu_mt2_bl = new TH1F("h_mu_mt2_bl", "h_mu_mt2_bl",60, 0, 600);
  fOutput->Add(h_mu_mt2_bl);
  h_mu_mt2_bl->Sumw2();



  // invariant mass 
  h_mu_ss_minv = new TH1F("h_mu_ss_minv", "h_mu_ss_minv",50, 0, 500); 
  fOutput->Add(h_mu_ss_minv); 
  h_mu_ss_minv->Sumw2();
  
  h_mu_os_minv = new TH1F("h_mu_os_minv", "h_mu_os_minv",50, 0, 500); 
  fOutput->Add(h_mu_os_minv);
  h_mu_os_minv->Sumw2();

  h_el_ss_minv = new TH1F("h_el_ss_minv", "h_el_ss_minv",50, 0, 500); 
  fOutput->Add(h_el_ss_minv);
  h_el_ss_minv->Sumw2();
  
  h_el_os_minv = new TH1F("h_el_os_minv", "h_el_os_minv",50, 0, 500); 
  fOutput->Add(h_el_os_minv);
  h_el_os_minv->Sumw2();
  
  h_elmu1_jet_minv = new TH1F("h_elmu1_jet_minv", "h_elmu1_jet_minv",50, 0, 500); 
  fOutput->Add(h_elmu1_jet_minv);
  h_elmu1_jet_minv->Sumw2();
  
  h_elmu2_jet_minv = new TH1F("h_elmu2_jet_minv", "h_elmu2_jet_minv",50, 0, 500); 
  fOutput->Add(h_elmu2_jet_minv);
  h_elmu2_jet_minv->Sumw2();  

  h_elmu_os_minv = new TH1F("h_elmu_os_minv", "h_elmu_os_minv",50, 0, 500); 
  fOutput->Add(h_elmu_os_minv);  
  h_elmu_os_minv->Sumw2(); 
  
  // asymmetry 
  h_mu_ss_asym = new TH1F("h_mu_ss_asym", "h_mu_ss_asym",50, -3,3); 
  fOutput->Add(h_mu_ss_asym);  
  h_mu_ss_asym->Sumw2(); 
  
  h_mu_os_asym = new TH1F("h_mu_os_asym", "h_mu_os_asym",50, -3,3); 
  fOutput->Add(h_mu_os_asym);
  h_mu_os_asym->Sumw2();
  
  h_el_ss_asym = new TH1F("h_el_ss_asym", "h_el_ss_asym",50, -3,3); 
  fOutput->Add(h_el_ss_asym);
  h_el_ss_asym->Sumw2();
  
  h_el_os_asym = new TH1F("h_el_os_asym", "h_el_os_asym",50, -3,3); 
  fOutput->Add(h_el_os_asym);
  h_el_os_asym->Sumw2();
  
  h_elmu_ss_asym = new TH1F("h_elmu_ss_asym", "h_elmu_ss_asym",50, -3,3); 
  fOutput->Add(h_elmu_ss_asym);
  h_elmu_ss_asym->Sumw2();
 
  h_elmu_os_asym = new TH1F("h_elmu_os_asym", "h_elmu_os_asym",50, -3,3); 
  fOutput->Add(h_elmu_os_asym);  
  h_elmu_os_asym->Sumw2();
  
  char hname[30];
  for (int i=0; i<3; i++) {
    sprintf(hname,"h_el_pt_%i",i);
    h_el_pt[i] = new TH1F(hname, hname, 50, 0, 500);
    fOutput->Add(h_el_pt[i]);
    h_el_pt[i]->Sumw2();
    
    sprintf(hname,"h_el_eta_%i",i);
    h_el_eta[i] = new TH1F(hname,hname, 100, -3, 3);
    fOutput->Add(h_el_eta[i]);
    h_el_eta[i]->Sumw2();
    
    sprintf(hname,"h_el_phi_%i",i);
    h_el_phi[i] = new TH1F(hname, hname, 100, -TMath::Pi(), TMath::Pi());
    fOutput->Add(h_el_phi[i]);
    h_el_phi[i]->Sumw2();
    
    if (i<2) { // os and ss distributions 
      sprintf(hname,"h_el_ss_pt_%i",i);
      h_el_ss_pt[i] = new TH1F(hname, hname, 50, 0, 500);
      fOutput->Add(h_el_ss_pt[i]);
      h_el_ss_pt[i]->Sumw2();
      
      sprintf(hname,"h_el_ss_eta_%i",i);
      h_el_ss_eta[i] = new TH1F(hname,hname, 100, -3, 3);
      fOutput->Add(h_el_ss_eta[i]);
      h_el_ss_eta[i]->Sumw2();
      
      sprintf(hname,"h_el_ss_phi_%i",i);
      h_el_ss_phi[i] = new TH1F(hname, hname, 100, -TMath::Pi(), TMath::Pi());
      fOutput->Add(h_el_ss_phi[i]);
      h_el_ss_phi[i]->Sumw2();
      
      sprintf(hname,"h_el_os_pt_%i",i);
      h_el_os_pt[i] = new TH1F(hname, hname, 50, 0, 500);
      fOutput->Add(h_el_os_pt[i]);
      h_el_os_pt[i]->Sumw2();
      
      sprintf(hname,"h_el_os_eta_%i",i);
      h_el_os_eta[i] = new TH1F(hname,hname, 100, -3, 3);
      fOutput->Add(h_el_os_eta[i]);
      h_el_os_eta[i]->Sumw2();
      
      sprintf(hname,"h_el_os_phi_%i",i);
      h_el_os_phi[i] = new TH1F(hname, hname, 100, -TMath::Pi(), TMath::Pi());
      fOutput->Add(h_el_os_phi[i]);
      h_el_os_phi[i]->Sumw2();

    }

    sprintf(hname,"h_el_trackd0_%i",i);
    h_el_trackd0[i] = new TH1F(hname, hname, 50, 0, 200);
    fOutput->Add(h_el_trackd0[i]);
    h_el_trackd0[i]->Sumw2();    
    
    sprintf(hname,"h_el_etcone20_%i",i);
    h_el_etcone20[i] = new TH1F(hname, hname, 50, 0, 200);
    fOutput->Add(h_el_etcone20[i]);
    h_el_etcone20[i]->Sumw2();

    sprintf(hname,"h_el_etcone30_%i",i);
    h_el_etcone30[i] = new TH1F(hname, hname, 50, 0, 200);
    fOutput->Add(h_el_etcone30[i]);
    h_el_etcone30[i]->Sumw2();
    

    sprintf(hname,"h_mu_pt_%i",i);
    h_mu_pt[i] = new TH1F(hname,hname, 50, 0, 500);
    fOutput->Add(h_mu_pt[i]);
    h_mu_pt[i]->Sumw2();
    
    sprintf(hname,"h_mu_eta_%i",i);
    h_mu_eta[i] = new TH1F(hname,hname, 100, -3,3);
    fOutput->Add(h_mu_eta[i]);
    h_mu_eta[i]->Sumw2();
    
    sprintf(hname,"h_mu_phi_%i",i);
    h_mu_phi[i] = new TH1F(hname,hname,100, -TMath::Pi(), TMath::Pi());
    fOutput->Add(h_mu_phi[i]);
    h_mu_phi[i]->Sumw2();
    

    if (i<2) { // os and ss distributions 
      sprintf(hname,"h_mu_ss_pt_%i",i);
      h_mu_ss_pt[i] = new TH1F(hname, hname, 50, 0, 500);
      fOutput->Add(h_mu_ss_pt[i]);
      h_mu_ss_pt[i]->Sumw2();
      
      sprintf(hname,"h_mu_ss_eta_%i",i);
      h_mu_ss_eta[i] = new TH1F(hname,hname, 100, -3, 3);
      fOutput->Add(h_mu_ss_eta[i]);
      h_mu_ss_eta[i]->Sumw2();
      
      sprintf(hname,"h_mu_ss_phi_%i",i);
      h_mu_ss_phi[i] = new TH1F(hname, hname, 100, -TMath::Pi(), TMath::Pi());
      fOutput->Add(h_mu_ss_phi[i]);
      h_mu_ss_phi[i]->Sumw2();
      
      sprintf(hname,"h_mu_os_pt_%i",i);
      h_mu_os_pt[i] = new TH1F(hname, hname, 50, 0, 500);
      fOutput->Add(h_mu_os_pt[i]);
      h_mu_os_pt[i]->Sumw2();

      sprintf(hname,"h_elmu_os_pt_%i",i);
      h_elmu_os_pt[i] = new TH1F(hname, hname, 50, 0, 500);
      fOutput->Add(h_elmu_os_pt[i]);
      h_elmu_os_pt[i]->Sumw2();

      
      sprintf(hname,"h_mu_os_eta_%i",i);
      h_mu_os_eta[i] = new TH1F(hname,hname, 100, -3, 3);
      fOutput->Add(h_mu_os_eta[i]);
      h_mu_os_eta[i]->Sumw2();
      
      sprintf(hname,"h_mu_os_phi_%i",i);
      h_mu_os_phi[i] = new TH1F(hname, hname, 100, -TMath::Pi(), TMath::Pi());
      fOutput->Add(h_mu_os_phi[i]);
      h_mu_os_phi[i]->Sumw2();

    }


    sprintf(hname,"h_mu_etcone20_%i",i);
    h_mu_etcone20[i] = new TH1F(hname, hname, 50, 0, 200);
    fOutput->Add(h_mu_etcone20[i]);
    h_mu_etcone20[i]->Sumw2();
    
    sprintf(hname,"h_mu_etcone30_%i",i);
    h_mu_etcone30[i] = new TH1F(hname, hname, 50, 0, 200);
    fOutput->Add(h_mu_etcone30[i]);
    h_mu_etcone30[i]->Sumw2();

    sprintf(hname,"h_mu_ptcone20_%i",i);
    h_mu_ptcone20[i] = new TH1F(hname, hname, 50, 0, 200);
    fOutput->Add(h_mu_ptcone20[i]);
    h_mu_ptcone20[i]->Sumw2();
    
    sprintf(hname,"h_mu_ptcone30_%i",i);
    h_mu_ptcone30[i] = new TH1F(hname, hname, 50, 0, 200);
    fOutput->Add(h_mu_ptcone30[i]);
    h_mu_ptcone30[i]->Sumw2();

    sprintf(hname,"h_mu_d0_exPV_%i",i);
    h_mu_d0_exPV[i] = new TH1F(hname, hname, 50, 0, 200);
    fOutput->Add(h_mu_d0_exPV[i]);
    h_mu_d0_exPV[i]->Sumw2();
    
    sprintf(hname,"h_mu_z0_exPV_%i",i);
    h_mu_z0_exPV[i] = new TH1F(hname, hname, 200, 0, 200);
    fOutput->Add(h_mu_z0_exPV[i]);
    h_mu_z0_exPV[i]->Sumw2();
    
    sprintf(hname,"h_jet_pt_%i",i);
    h_jet_pt[i] = new TH1F(hname, hname, 200, 0, 400);
    fOutput->Add(h_jet_pt[i]);
    h_jet_pt[i]->Sumw2();
    
    sprintf(hname,"h_jet_eta_%i",i);
    h_jet_eta[i] = new TH1F(hname,hname, 100, -4, 4);
    fOutput->Add(h_jet_eta[i]);
    h_jet_eta[i]->Sumw2();
    
    sprintf(hname,"h_jet_phi_%i",i);
    h_jet_phi[i] = new TH1F(hname, hname,100, -TMath::Pi(), TMath::Pi());
    fOutput->Add(h_jet_phi[i]);
    h_jet_phi[i]->Sumw2();

    sprintf(hname,"h_jet_m_%i",i);
    h_jet_m[i] = new TH1F(hname, hname,200, 0, 400);
    fOutput->Add(h_jet_m[i]);
    h_jet_m[i]->Sumw2();

    // flavour tagging 
    sprintf(hname,"h_jet_ftag_sv0_%i",i);
    h_jet_ftag_sv0[i] = new TH1F(hname, hname,200, -100, 500);
    fOutput->Add(h_jet_ftag_sv0[i]);
    h_jet_ftag_sv0[i]->Sumw2();

    sprintf(hname,"h_jet_ftag_JProb_%i",i);
    h_jet_ftag_JProb[i] = new TH1F(hname, hname,50, 0, 1);
    fOutput->Add(h_jet_ftag_JProb[i]);
    h_jet_ftag_JProb[i]->Sumw2();

    sprintf(hname,"h_jet_ftag_ntrack_%i",i);
    h_jet_ftag_ntrack[i] = new TH1F(hname, hname,200, -50, 200);
    fOutput->Add(h_jet_ftag_ntrack[i]);
    h_jet_ftag_ntrack[i]->Sumw2();

    // dr 
    sprintf(hname,"h_el_mu_dr_%i",i);
    h_el_mu_dr[i] = new TH1F(hname, hname, 100, 0, 5);
    fOutput->Add(h_el_mu_dr[i]);
    h_el_mu_dr[i]->Sumw2();
    
    sprintf(hname,"h_el_jet_dr_%i",i);
    h_el_jet_dr[i] = new TH1F(hname, hname, 100, 0, 5);
    fOutput->Add(h_el_jet_dr[i]);
    h_el_jet_dr[i]->Sumw2();

    sprintf(hname,"h_el_el_dr_%i",i);
    h_el_el_dr[i] = new TH1F(hname, hname, 100, 0, 5);
    fOutput->Add(h_el_el_dr[i]);
    h_el_el_dr[i]->Sumw2();

    sprintf(hname,"h_el_trig_dr_%i",i);
    h_el_trig_dr[i] = new TH1F(hname, hname, 100, 0, 5);
    fOutput->Add(h_el_trig_dr[i]);
    h_el_trig_dr[i]->Sumw2();

    sprintf(hname,"h_mu_trig_dr_%i",i);
    h_mu_trig_dr[i] = new TH1F(hname, hname, 100, 0, 5);
    fOutput->Add(h_mu_trig_dr[i]);
    h_mu_trig_dr[i]->Sumw2();

    sprintf(hname,"h_mu_etm_dr_%i",i);
    h_mu_etm_dr[i] = new TH1F(hname, hname, 100, 0, 5);
    fOutput->Add(h_mu_etm_dr[i]);
    h_mu_etm_dr[i]->Sumw2();

    sprintf(hname,"h_el_etm_dr_%i",i);
    h_el_etm_dr[i] = new TH1F(hname, hname, 100, 0, 5);
    fOutput->Add(h_el_etm_dr[i]);
    h_el_etm_dr[i]->Sumw2();

    sprintf(hname,"h_jet_etm_dr_%i",i);
    h_jet_etm_dr[i] = new TH1F(hname, hname, 100, 0, 5);
    fOutput->Add(h_jet_etm_dr[i]);
    h_jet_etm_dr[i]->Sumw2();

  }
  
  
}


void read_d4pd::SaveHistos() 
{

 //  TFile *fileOut = new TFile( "histograms.root","RECREATE");
 //  fileOut->cd();


 //  h_cutFlow = dynamic_cast<TH1D *>(fOutput->FindObject(Form("h_cutFlow")));
 //  h_cutFlow->Write();
  
 // h_jet_n = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_jet_n")));
 // h_jet_n->Write();
 
 // h_vertex_n = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_vertex_n")));
 // h_vertex_n->Write();
 
 // h_etmiss_mod = dynamic_cast<TH1F *>(fOutput->FindObject(Form("h_etmiss_mod")));
 // h_etmiss_mod->Write();

 
 // char hname[30]; 
 // for (int i=0; i<3; i++) {
   
 //   // electrons 
 //   sprintf(hname,"h_el_pt_%i",i);
 //   h_el_pt[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_el_pt[i]->Write();
   
 //   sprintf(hname,"h_el_eta_%i",i);   
 //   h_el_eta[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_el_eta[i]->Write();
   
 //   sprintf(hname,"h_el_phi_%i",i);   
 //   h_el_phi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_el_phi[i]->Write();  
   
 //   sprintf(hname,"h_el_trackd0_%i",i);       
 //   h_el_trackd0[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_el_trackd0[i]->Write();
   
 //   sprintf(hname,"h_el_etcone20_%i",i);      
 //   h_el_etcone20[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_el_etcone20[i]->Write();
   
 //   sprintf(hname,"h_el_etcone30_%i",i);      
 //   h_el_etcone30[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_el_etcone30[i]->Write();
   
 //   // muons 
 //   sprintf(hname,"h_mu_pt_%i",i); 
 //   h_mu_pt[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_mu_pt[i]->Write();
  
 //   sprintf(hname,"h_mu_eta_%i",i); 
 //   h_mu_eta[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_mu_eta[i]->Write();
  
 //   sprintf(hname,"h_mu_phi_%i",i); 
 //   h_mu_phi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_mu_phi[i]->Write();

 //   sprintf(hname,"h_mu_etcone20_%i",i);      
 //   h_mu_etcone20[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_mu_etcone20[i]->Write();

 //   sprintf(hname,"h_mu_etcone30_%i",i);      
 //   h_mu_etcone30[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_mu_etcone30[i]->Write();

 //   sprintf(hname,"h_mu_ptcone20_%i",i);      
 //   h_mu_ptcone20[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_mu_ptcone20[i]->Write();

 //   sprintf(hname,"h_mu_ptcone30_%i",i);      
 //   h_mu_ptcone30[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_mu_ptcone30[i]->Write();

 //   sprintf(hname,"h_mu_d0_exPV_%i",i);       
 //   h_mu_d0_exPV[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_mu_d0_exPV[i]->Write();

 //   sprintf(hname,"h_mu_z0_exPV_%i",i);       
 //   h_mu_z0_exPV[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_mu_z0_exPV[i]->Write();

 //   // jets 
 //   sprintf(hname,"h_jet_pt_%i",i);
 //   h_jet_pt[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_jet_pt[i]->Write();
  
 //   sprintf(hname,"h_jet_eta_%i",i);  
 //   h_jet_eta[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_jet_eta[i]->Write();

 //   sprintf(hname,"h_jet_phi_%i",i);  
 //   h_jet_phi[i] = dynamic_cast<TH1F *>(fOutput->FindObject(Form(hname)));
 //   h_jet_phi[i]->Write();  

 // }

 // fileOut->Close();
 
}








bool read_d4pd::isaGoodLumiblock(Long64_t entry){
  b_lbn->GetEntry(entry); 
  b_RunNumber->GetEntry(entry); 
  
  bool answer=false;
  if(m_eventInfo.IsGoodLumiBlock(lbn,RunNumber) ) answer= true;  
  
  //std::cout <<"LBN,RUNNUMBER, passGRL  " <<lbn<<'\t' <<RunNumber <<'\t' <<answear <<'\n'; 
  return answer;  
}

bool read_d4pd::isaMCRun(UInt_t run_number){
  // list of MC runs  
  // for file in $(ls /Volumes/DataA_1/xrd/WW/d4pd/ | grep mc10  | awk -F. ' {print $4}'); do echo -n $file,; done
  // some duplicates here due to different versions:  
  // UInt_t mc_runs[93]={105009,105010,105010,105011,105013,105014,105014,105015,105015,105016,
  // 		     105200,105204,105204,,105986,105987,107650,107650,107650,107651,107652,107653,107654,107655,107660,107660,107660,107661,107661,107662,107663,107664,107665,
  // 		     107670,107670,107670,107671,107672,107673,107674,107675,107680,107680,107681,107682,107683,107684,107685,107690,107691,107692,107693,107694,107695,
  // 		     107700,107701,107702,107703,107704,107705,108340,108340,108341,108341,108342,108343,108344,108344,108345,108346,116250,116251,116252,116253,116254,116255,
  // 		      116260,116261,116262,116262,116263,116264,116265,116270,116270,116271,116271,116272,116272,116272,116273,116274,116275}; 
  
  // m_mc_runs;  contain list of runs 
  bool ismc =false;   
  for (unsigned int i=0; i<m_mc_runs.size(); i++) {
    if (run_number == m_mc_runs.at(i)) {
      ismc = true; 
      break; 
    }
  }   
  return ismc;
}


float read_d4pd::Weight(UInt_t run_number){
  // could handle data e.g. Trigger prescales, here. 
  float w=1; 
  for (unsigned int i=0; i<m_mc_runs.size(); i++) {
    if (run_number == m_mc_runs.at(i)) {
      w=m_mc_scales.at(i);
    }
  }   
  return w;
}

void read_d4pd::ReadWeights(TString mc_info_file) 
{
  ifstream file; 
  file.open(mc_info_file); 
   // mc_info.dat contain MC run number, number of events in the ntuples, and x-section*eff*kfactor(eventually) in pb
  unsigned int run;
  float nevt;
  float xsec; 
  float scale; 
  
  if(file) 
    {	    
      while(!file.eof()) 
	{
	  file >> run >> nevt>> xsec;
	  m_mc_runs.push_back(run);
	  // weight yield for 1fb-1 
	  scale=xsec*1000.0*20.3/(nevt); 
	  m_mc_scales.push_back(scale);   
	  //	  std::cout << "RUN and weight " << run <<'\t' <<scale <<"   "<<"xsec:   "<<xsec<<" "<<"nevt:   "<<nevt<<'\n';  
	}
      file.close();
    }
  else {
    std::cout << "Can not open the file !" << endl;
  }
}

void read_d4pd::fillcutFlow(int id) 
{
  h_cutFlow->Fill(id);
  h_cutFlow_w->Fill(id,w);
}


bool read_d4pd::trigger(Long64_t entry,Int_t m_stream)
{

  
  // el single  
  b_EF_e24vhi_medium1->GetEntry(entry);
  b_EF_e60_medium1->GetEntry(entry); 
  // el pairs
  b_EF_e24vh_medium1_e7_medium1->GetEntry(entry);
  b_EF_2e12Tvh_loose1->GetEntry(entry);
  // mu single 
  b_EF_mu24i_tight->GetEntry(entry); 
  b_EF_mu36_tight->GetEntry(entry); 
  //mu pairs 
  b_EF_mu18_tight_mu8_EFFS->GetEntry(entry);
  b_EF_mu24_tight_mu6_EFFS->GetEntry(entry); 
  b_EF_2mu13->GetEntry(entry); 
  // el mu mixt 
  b_EF_mu18_tight_e7_medium1->GetEntry(entry); 
  b_EF_e12Tvh_medium1_mu8->GetEntry(entry); 
  
  //MET trigger
    if (b_EF_xe80T_tclcw_loose) b_EF_xe80T_tclcw_loose->GetEntry(entry);
  
  if (b_EF_xe80_tclcw_loose)  b_EF_xe80_tclcw_loose->GetEntry(entry);

  //  b_EF_xe80T_tclcw_loose->GetEntry(entry);
  //b_EF_xe80_tclcw_loose->GetEntry(entry);

  bool  Trig_ele(0), Trig_ele_asim(0), Trig_ele_simm(0);
  bool Trig_mu(0), Trig_mu_asim(0), Trig_mu_simm(0);
  bool Trig_emu(0),  Trig_mue(0);
  bool Trig_met(0),  Trig_metT(0);
   if(EF_xe80_tclcw_loose)Trig_met=true;
 
  if(EF_xe80T_tclcw_loose)Trig_metT=true;
  
  if(EF_e24vhi_medium1 || EF_e60_medium1) Trig_ele = true;
  if(EF_e24vh_medium1_e7_medium1)  Trig_ele_asim = true;
  if(EF_2e12Tvh_loose1) Trig_ele_simm = true;
  
  if(EF_mu24i_tight || EF_mu36_tight) Trig_mu = true;
  if(EF_mu18_tight_mu8_EFFS || EF_mu24_tight_mu6_EFFS) Trig_mu_asim = true;
  if(EF_2mu13) Trig_mu_simm = true;
  
  if(EF_e12Tvh_medium1_mu8) Trig_emu = true;
  if(EF_mu18_tight_e7_medium1) Trig_mue = true; 
  
  if(Trig_met==true || Trig_metT==true) m_passTrigger_met =true;
  // if(Trig_ele==true || Trig_ele_simm==true || Trig_ele_asim==true || Trig_emu==true || Trig_mue==true) m_passTrigger_el = true;  //not used
       if(Trig_ele==true || Trig_ele_simm==true || Trig_ele_asim==true || Trig_emu==true) m_passTrigger_el = true;
    if(Trig_mu==true || Trig_mu_simm==true || Trig_mu_asim==true || Trig_mue==true || Trig_emu==true) m_passTrigger_mu = true;
  
  bool TrigPassed=false;
  if(m_stream == 1) TrigPassed =  m_passTrigger_el;  //egamma
  if(m_stream == 2) TrigPassed = m_passTrigger_mu;   //muons
  if(m_stream == 3) TrigPassed = m_passTrigger_mu || m_passTrigger_el;  //combined
  if(m_stream == 4) TrigPassed = m_passTrigger_met; //MET trigger
  if(m_stream == 5) TrigPassed = m_passTrigger_el || m_passTrigger_met ||  m_passTrigger_mu;
  if(m_stream == 6) TrigPassed = m_passTrigger_mu || m_passTrigger_met;
  
  return TrigPassed;
}




void read_d4pd::ReadJobConfiguration(TString conf_file) 
{

  bool dbg;
  bool isdata;
  bool isAtlfast;
  bool ptcut;
  bool ismc;
  Int_t  stream;
  unsigned int  syst_flag;
  bool triglep;
  bool createconfig;


//  namespace SystErr
//  {
//   typedef enum  {
//     NONE,
//     JESDOWN,JESUP, // Return the total JES uncertainty as the quadratic sum of the 7 different uncertainties listed below 
//     BaselineDown,BaselineUp,ForwardDown,ForwardUp,PileupMuDown,PileupMuUp,PileupNPVDown,PileupNPVUp,CloseByDown,CloseByUp,FlavourDown,FlavourUp,BJesDown,BJesUp, //Break down of the JES uncertainties  
//     JER, 
//     EGZEEUP, EGZEEDOWN,EGMATUP,EGMATDOWN,EGPSUP,EGPSDOWN,EGLOWUP,EGLOWDOWN, EGRESDOWN, EGRESUP, EEFFDOWN, EEFFUP, ETRGDOWN, ETRGUP,
//     MMSLOW, MMSUP, MIDLOW, MIDUP, MEFFDOWN, MEFFUP, MSCALE,
//     BJETDOWN, BJETUP, CJETDOWN, CJETUP, BMISTAGDOWN, BMISTAGUP,
//     SCALESTUP,SCALESTDOWN,RESOST,
//     TESUP, TESDOWN /// Tau energy scale and resolution
//     //    RESOSTUP,RESOSTDOWN, // Just have one RESOST, because the UP and DOWN variations are both symmetric smearings.
//     //    SCALEPHUP,SCALEPHDOWN,RESOPHUP,RESOPHDOWN,RESOPHUPDOWN,RESOPHDOWNUP,
//   } Syste;
// }


 


  ifstream file; 
  file.open(conf_file); 
   // config_info.dat contain Job configuration flags 
  
  if(file) 
    {	    
      while(!file.eof()) 
	{
	  //file >> dbg >> isdata >> ismc >> isAtlfast>> stream >> triglep  >> createconfig;
	  file >> dbg >> isdata >> ismc >> isAtlfast>> stream >> ptcut  >> createconfig;
	  
	  //    	  std::cout << "RUNNing with Configuration:   " << dbg <<'\t'<< isdata <<'\t'<< ismc  <<'\t'<< isAtlfast <<'\t'<<   stream <<'\t'<<  triglep  <<'\t' <<  createconfig <<'\n';  
	  std::cout << "RUNNing with Configuration:   " << dbg <<'\t'<< isdata <<'\t'<< ismc  <<'\t'<< isAtlfast <<'\t'<<   stream <<'\t'<<  ptcut  <<'\t' <<  createconfig <<'\n';
	  
	}
      file.close();
    }
  else {
    std::cout << "Can not open the file: " << conf_file << std::endl;
  }
  
  
  debug=dbg;  //false; //true;
  m_isDATA=isdata;  //true; //false; 
  m_isFS=isAtlfast; //false;//true;
  m_isMC=ismc;  //false;//true;
  m_stream=stream;
  //  m_whichsyst = syst_flag; // SystErr::NONE;
  m_isPtcut=ptcut;
  m_istriglep=triglep;
    CreateConfigFile = createconfig; 
  
}




float read_d4pd::MT2(Electron a,Muon b,TLorentzVector Etm) 
{
  // compute mt2... 
  // pa, pb = {mass, px, py}
  // pmiss  = {NULL, pxmiss, pymiss}
  // mn     = invisible particle mass

  //std::cout <<"a.name" <<  a.name <<'\n';

  double pa[3] = { 0.005, 0, 0 };
  double pb[3] = { 0.106, 0, 0 };
  double pmiss[3] = { 0, 0,0 };
  double mn    = 1.;
  pa[0]=a.M();
  pa[1]=a.X();
  pa[2]=a.Y();
  pb[0]=a.M();
  pb[1]=b.X();
  pb[2]=b.Y();
  pmiss[0]=Etm.M();  
  pmiss[1]=Etm.X();
  pmiss[2]=Etm.Y();
  
  mt2_bisect::mt2 mt2_event;
  mt2_event.set_momenta(pa,pb,pmiss);
  mt2_event.set_mn(mn);
  //mt2_event.print();
  //std::cout  << " mt2 = " << mt2_event.get_mt2() << endl;
  return mt2_event.get_mt2();
}


float read_d4pd::MT2(Electron a,Electron b,TLorentzVector Etm) 
{
  // compute mt2... 
  // pa, pb = {mass, px, py}
  // pmiss  = {NULL, pxmiss, pymiss}
  // mn     = invisible particle mass

  //std::cout <<"a.name" <<  a.name <<'\n';

  double pa[3] = { 0.005, 0, 0 };
  double pb[3] = { 0.005, 0, 0 };
  double pmiss[3] = { 0, 0,0 };
  double mn    = 1.;
  pa[0]=a.M();  
  pa[1]=a.X();
  pa[2]=a.Y();
  pb[0]=a.M();
  pb[1]=b.X();
  pb[2]=b.Y();
  pmiss[0]=Etm.M();
  pmiss[1]=Etm.X();
  pmiss[2]=Etm.Y();
  
  mt2_bisect::mt2 mt2_event;
  mt2_event.set_momenta(pa,pb,pmiss);
  mt2_event.set_mn(mn);
  //mt2_event.print();
  //std::cout  << " mt2 = " << mt2_event.get_mt2() << endl;
  return mt2_event.get_mt2();
}


float read_d4pd::MT2(Muon a,Muon b,TLorentzVector Etm) 
{
  // compute mt2... 
  // pa, pb = {mass, px, py}
  // pmiss  = {NULL, pxmiss, pymiss}
  // mn     = invisible particle mass

  //std::cout <<"a.name" <<  a.name <<'\n';

  double pa[3] = { 0.106, 0, 0 };
  double pb[3] = { 0.106, 0, 0 };
  double pmiss[3] = { 0, 0,0 };
  double mn    = 1.;
  pa[0]=a.M();  
  pa[1]=a.X();
  pa[2]=a.Y();
  pb[0]=a.M();
  pb[1]=b.X();
  pb[2]=b.Y();
  pmiss[0]=Etm.M();
  pmiss[1]=Etm.X();
  pmiss[2]=Etm.Y();
  
  mt2_bisect::mt2 mt2_event;
  mt2_event.set_momenta(pa,pb,pmiss);
  mt2_event.set_mn(mn);
  //mt2_event.print();
  //std::cout  << " mt2 = " << mt2_event.get_mt2() << endl;
  return mt2_event.get_mt2();
}

float read_d4pd::MT2_bl( Jet a,Electron b,TLorentzVector Pbll_vec)
{


  double pa[3] = { 0, 0, 0 };
  double pb[3] = { 0.005, 0, 0 };
  double Pbll_1[3] = { 0, 0,0 };
  double mn    = 1.;
  pa[0]=a.M();
  pa[1]=a.X();
  pa[2]=a.Y();
  pb[0]=b.M();
  pb[1]=b.X();
  pb[2]=b.Y();
  Pbll_1[0]=Pbll_vec.M();
  Pbll_1[1]=Pbll_vec.X();
  Pbll_1[2]=Pbll_vec.Y();


  mt2_bisect::mt2 mt2_jetlep;
  mt2_jetlep.set_momenta(pa,pb,Pbll_1);
  mt2_jetlep.set_mn(mn);
  return mt2_jetlep.get_mt2();
}

/*
float read_d4pd::had_MT2( Jet a,Jet b, Muon a,Muon b,TLorentzVector Etm)                                                                                                                                                   
{

  double pa[3] = { 0, 0, 0 };
  double pb[3] = { 0, 0, 0 };
  double pmiss[3] = { 0, 0,0 };
  double mn    = 1.;
  pa[0]=a.M();
  pa[1]=a.X();
  pa[2]=a.Y();
  pb[0]=b.M();
  pb[1]=b.X();
  pb[2]=b.Y();
  pmiss[0]=Etm.M();
  pmiss[1]=Etm.X();
  pmiss[2]=Etm.Y();


  mt2_bisect::mt2 mt2_jet;
  mt2_jet.set_momenta(pa,pb,pmiss);
  mt2_jet.set_mn(mn);
  return mt2_jet.get_mt2();
}


*/


float read_d4pd::had_MT2( Jet a,Jet b,TLorentzVector Pbll_vec)
 //float read_d4pd::had_MT2( Jet a,Jet b,TVector2 Pbll_vec)
					  
{
  
  
  double pa[3] = { 0, 0, 0 };
  double pb[3] = { 0, 0, 0 };
  double Pbll_1[3] = { 0, 0,0 };
  double mn    = 1.;
  pa[0]=a.M();
  pa[1]=a.X();
  pa[2]=a.Y();
  pb[0]=b.M();
  pb[1]=b.X();
  pb[2]=b.Y();
  //  Pbll_1[0]=Pbll_vec.M(); 
  Pbll_1[1]=Pbll_vec.X();
  Pbll_1[2]=Pbll_vec.Y();
  
  mt2_bisect::mt2 mt2_jet;
  mt2_jet.set_momenta(pa,pb,Pbll_1);
  mt2_jet.set_mn(mn);
  return mt2_jet.get_mt2();
}







float read_d4pd::ASYM( Particle a,Particle b,TLorentzVector Etm)  

{
  
  Particle  first=a; 
  Particle second=b; 
  double val= (double) rand() / (double)RAND_MAX;  
  if ( val >0.5) { 
    Particle  first=b; 
    Particle second=a; 
  } 
  float asym=  ( 2*( first* Etm - second * Etm) / (first * Etm + second * Etm) ); 
  
  return asym;
}

