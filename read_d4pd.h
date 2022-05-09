//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul  6 15:31:46 2011 by ROOT version 5.26/00
// from TTree susy/susy 
/////////////////////////////////////////////////////////
#ifndef read_d4pd_h
#define read_d4pd_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include<TH1F.h>
#include<TH1D.h>
#include<TH2F.h>
#include <TGraph.h> 
#include <TSelector.h>
#include <TProofOutputFile.h>
#include <TNtuple.h>
#include <vector>
#include <map>
#include <TRandom.h>
#include "myown_util/EventInfo.h"
#include "myown_util/Particle.h"
#include "SUSYTools/SUSYObjDef.h"
#include "SUSYTools/FakeMetEstimator.h"
#include "SUSYTools/BTagCalib.h"
#include "PileupReweighting/TPileupReweighting.h"
#include "/Users/sdarmora/ST_3_20/SUSYTools/SUSYTools/MuonTriggerSFTool.h"
#include "/Users/sdarmora/ST_3_20/MuonEfficiencyCorrections/MuonEfficiencyCorrections/AnalysisMuonEfficiencyScaleFactors.h"
#include "/Users/sdarmora/ST_3_20/ReweightUtils/ReweightUtils/APWeightEntry.h"  //General weight representation (mandatory)
#include "/Users/sdarmora/ST_3_20/ReweightUtils/ReweightUtils/APReweightND.h"   //3D Reweight class (only for ND parametrization of the weights)
#include "/Users/sdarmora/ST_3_20/ReweightUtils/ReweightUtils/APEvtWeight.h"   //Event Weight calculation class (if single event weights need to be known, e.g. for propagation or usage outside a simple counting/summing)
 

class read_d4pd : public TSelector {
 public :

  APReweightND* m_trigWeighter;
  APEvtWeight* trig_weight_muon;
  MuonTriggerSFTool* MTriggerSF;
  //BTagCalib btag;  
  SUSYObjDef susyobj;
  EventInfo m_eventInfo;
  FakeMetEstimator m_fakemet; 
  Root::TPileupReweighting  pileup;
  BTagCalib *BTagComputation;
  //  Okla_grl *grl; 
  UInt_t Current_RunNumber;
  UInt_t Current_lbn;
  float w; //event weight 
  bool m_bad_lb;

  Int_t m_stream;
  //  Int_t m_whichsyst;
  TString  m_whichsyst;
  
  
  bool doPUreweight;
  bool CreateConfigFile;  
  bool m_isMC;
  bool m_isFS;
  bool m_isDATA;
  bool m_passGRL;    
  bool m_passTrigger_el;
  bool m_passTrigger_mu;
  bool m_passTrigger_met;
  bool m_GoodVertex; 
  bool m_GoodElectron; 
  bool m_GoodMuon; 
  bool m_GoodJet; 
  bool m_BadJetEvent;
  bool m_BadJetEvent1;
  bool debug;
  bool m_JetEvent;
  bool m_ismc12b; 
  bool m_istriglep;
  bool m_isPtcut;
  
  
  TProofOutputFile *fProofFile;
  TNtuple *theNtuple;

  //TFile *f_output;
  //  TDirectory *savedir;
  TDirectory *fFile;

   

  SystErr::Syste whichsyste;
  SUSYBTagger::btagger jetTagger;  
  SUSYMet::met_definition whichmet;

  std::map <int,int> m_totevt; 

  TRandom tr; 
  //  list of mc runs and weigth 
  std::vector<unsigned int> m_mc_runs;
  std::vector<float> m_mc_scales;  

  // count event yield 
  std::vector<unsigned int>  m_runs; 
  std::vector<unsigned int>  m_totevt_runs; 
  std::vector<unsigned int>  m_selevt_runs; 

  // particles: TLorentz + few flags 
  std::vector<Electron> Electrons;
  std::vector<Muon> Muons;
  //std::vector<Photon> Photons; 
  std::vector<Jet> Jets;
  

  void BookHistos();  
  void SaveHistos();
  bool isaGoodLumiblock(Long64_t entry);
  bool isaMCRun(UInt_t run_number);
  void ReadWeights(TString mc_file );
  float Weight(unsigned int runnumber);
  void fillcutFlow(int id);
  bool trigger(Long64_t entry,Int_t m_stream);
  void ReadJobConfiguration(TString conf_file); 
  
  float MT2(Muon a,Muon b,TLorentzVector Etm);
  float MT2(Electron a,Electron b,TLorentzVector Etm);
  float MT2(Electron a,Muon b,TLorentzVector Etm);
   float had_MT2(Jet a,Jet b,TLorentzVector Pbll_vec);

   // float had_MT2(TVector2 a,TVector2 b,TLorentzVector Pbll_vec);
  // float had_MT2(Jet a,Jet b,TVector2 Pbll_vec);
  float MT2_bl(Jet a,Electron b,TLorentzVector Pbll_vec);
  float ASYM(Particle a,Particle b,TLorentzVector Etm);  

  //-------------------------------- Histos ----------------------

 TH1D * h_cutFlow;
  TH1D * h_cutFlow_w;
  //TGraph * g_event_yield; 
  TH1F * h_event_yield; 
  TH1D * h_totevent; 
  TH1F * h_btag_weight;
  TH1D * h_el_trigg_freq;
  TH1D * h_mu_trigg_freq;
  TH1D * h_weight_167803;
  TH1D * h_weight_167804;
  TH1D * h_weight_167805;
  TH1D * h_weight_167815;
  TH1D * h_weight_167816;
  TH1D * h_weight_167817;
  TH1D * h_weight_167821;
  TH1D * h_weight_167822;
  TH1D * h_weight_167823;
  TH1D * h_weight_167824;
  TH1D * h_weight_167825;
  TH1D * h_weight_167826;
  TH1D * h_weight_167827;
  TH1D * h_weight_167828;
  TH1D * h_weight_167829;
  TH1D * h_weight_167833;
  TH1D * h_weight_167834;
  TH1D * h_weight_167835;
  TH1D * h_weight_167839;
  TH1D * h_weight_167840;
  TH1D * h_weight_167841;
 
  TH1D * h_weight_167749;
  TH1D * h_weight_167750;
  TH1D * h_weight_167751;
  TH1D * h_weight_167752;
  TH1D * h_weight_167753;
  TH1D * h_weight_167754;
  TH1D * h_weight_167755;
  TH1D * h_weight_167756;
  TH1D * h_weight_167757;
  TH1D * h_weight_167797;
  TH1D * h_weight_167798;
  TH1D * h_weight_167799;
  TH1D * h_weight_167800;


  TH1D * h_weight_167801;
  TH1D * h_weight_167802;
  TH1D * h_weight_167809;
  TH1D * h_weight_167810;
  TH1D * h_weight_167811;
  TH1D * h_weight_167812;
  TH1D * h_weight_167813;
  TH1D * h_weight_167814;
  TH1D * h_weight_180543;
  TH1D * h_weight_180544;
  TH1D * h_weight_180545;
  TH1D * h_weight_180546;
  TH1D * h_weight_180547;
  TH1D * h_weight_180548;
  TH1D * h_weight_180549;
  TH1D * h_weight_180550;
  TH1D * h_weight_180551;

 


  // electrons 
  TH1F  * h_el_n; 
  TH1F  * h_el_or_n; 
  TH1F  * h_el_pt[3]; 
  TH1F  * h_el_eta[3];
  TH1F  * h_el_phi[3];


  TH1F  * h_el_ss_pt[2]; 
  TH1F  * h_el_ss_eta[2];
  TH1F  * h_el_ss_phi[2];
  TH1F  * h_el_os_pt[2]; 
  TH1F  * h_el_os_eta[2];
  TH1F  * h_el_os_phi[2]; 
  TH1F  * h_el_os_dphi;
  TH1F  * h_el_os_dtheta;
  TH1F  * h_el_jet_dphi[2];
  TH1F  * h_el_etm_dphi[2];


  TH1F * h_el_trackd0[3];
  // z0 missing 
  TH1F * h_el_etcone30[3];
  TH1F * h_el_etcone20[3];


  // muons 
  TH1F  * h_mu_n;  
  TH1F  * h_mu_or_n; 
  TH1F  * h_mu_pt[3]; 
  TH1F  * h_mu_eta[3];
  TH1F  * h_mu_phi[3];
  
  TH1F  * h_mu_ss_pt[2]; 
  TH1F  * h_mu_ss_eta[2];
  TH1F  * h_mu_ss_phi[2];
  TH1F  * h_mu_os_pt[2]; 
  TH1F  * h_mu_os_eta[2];
  TH1F  * h_mu_os_phi[2];
  TH1F  * h_mu_os_dphi;
  TH1F  * h_mu_os_dtheta;
  TH1F  * h_mu_jet_dphi[2];
  TH1F  * h_mu_etm_dphi[2];


  TH1F * h_mu_etcone20[3];
  TH1F * h_mu_etcone30[3];
  TH1F * h_mu_ptcone20[3];
  TH1F * h_mu_ptcone30[3];
  TH1F * h_mu_d0_exPV[3]; 
  TH1F * h_mu_z0_exPV[3]; 
    
  
  // jets 
  TH1F  * h_jet_n;
  TH1F  * h_jet_or_n;
  TH1F  * h_jet_b_or_n;
  TH1F  * h_jet_ht;
  TH1F  * h_jet_jet_dphi;

  TH1F  * h_jet_pt[3]; 
  TH1F  * h_jet_eta[3];
  TH1F  * h_jet_phi[3];
  TH1F  * h_jet_m[3];
  TH1F  * h_jet_ftag_sv0[3];
  TH1F  * h_jet_ftag_JProb[3];
  TH1F  * h_jet_ftag_ntrack[3];
  


    
  // reco lepton to trigger object distance   
  TH1F  * h_el_trig_dr[3]; 
  TH1F  * h_mu_trig_dr[3]; 

  // objects distance  
  TH1F  * h_el_mu_dr[3];   
  TH1F  * h_el_jet_dr[3]; 
  TH1F  * h_el_el_dr[3];
  
  TH1F  * h_el_etm_dr[3];   
  TH1F  * h_mu_etm_dr[3];   
  TH1F  * h_jet_etm_dr[3];   

  TH1F * h_vertex_n; 
  
  // ETMiss 
  TH1F  * h_etmiss_mod;
  TH1F  * h_etmiss_phi; 
  TH1F  * h_etmiss;

  TH1F  * h_etmiss_e;
  TH1F  * h_etmiss_mu;
  TH1F  * h_etmiss_emu;

  TH1F  * h_elmu_n;
  TH1F  * h_mu_mt2;
  TH1F  * h_el_mt2;
  TH1F  * h_elmu_mt2;
  TH1F  * h_elmu_mt2_cor;
  TH1F  * h_jet_mt2;
  TH1F  * h_jet_elmu_mt2; 
  TH1F  * h_el_mt2_bl;
  TH1F  * h_elmu_mt2_bl;
  TH1F  * h_mu_mt2_bl;


  TH1F * h_mu_ss_minv;
  TH1F * h_mu_os_minv;
  TH1F * h_el_ss_minv;
  TH1F * h_el_os_minv;  
  TH1F * h_elmu1_jet_minv;
  TH1F * h_elmu2_jet_minv;
  
  TH1F * h_elmu_os_minv;

  TH1F * h_mu_ss_asym;
  TH1F * h_mu_os_asym;
  TH1F * h_el_ss_asym;
  TH1F * h_el_os_asym;  
  TH1F * h_elmu_ss_asym;
  TH1F * h_elmu_os_asym;

  TH1F  * h_elmu_os_pt[2];
  TH1F  * h_elmu_os_dphi;
  TH1F  * h_elmu_os_dtheta;
  TH1F  * h_elmu_jet_dphi[2];
  TH1F  * h_elmu_etm_dphi[2];


  //deta for leptons

  TH1F  * h_el_os_deta;
  TH1F  * h_mu_os_deta;
  TH1F  * h_elmu_os_deta;
  //Pbll

  TH1F  * h_el_os_Pbll;
  TH1F  * h_mu_os_Pbll;
  TH1F  * h_elmu_os_Pbll;

  TH1F  * h_el_os_DPhiMetPbll;
  TH1F  * h_mu_os_DPhiMetPbll;
  TH1F  * h_elmu_os_DPhiMetPbll;


  //deltaR between leptons
  TH1F  * h_el_os_DRll;
  TH1F  * h_mu_os_DRll;
  TH1F  * h_elmu_os_DRll;

  //sum of lepton pt
  TH1F  * h_el_os_PtLep;
  TH1F  * h_mu_os_PtLep;
  TH1F  * h_elmu_os_PtLep;

  //lep jet pt ratio

  TH1F  * h_el_os_PtRatio;
  TH1F  * h_mu_os_PtRatio;
  TH1F  * h_elmu_os_PtRatio;

  //hadronic mt2

  TH1F  * h_el_os_had_mt2;
  TH1F  * h_mu_os_had_mt2;
  TH1F  * h_elmu_os_had_mt2;

  
  //C1



  TH1F  * h_etmiss_e_C1;
  TH1F  * h_etmiss_mu_C1;
  TH1F  * h_etmiss_emu_C1;
  
  TH1F  * h_mu_mt2_C1;
  TH1F  * h_el_mt2_C1;
  TH1F  * h_elmu_mt2_C1;


  TH1F  * h_jet_mt2_C1;
  
  
  TH1F * h_mu_os_minv_C1;
  TH1F * h_el_os_minv_C1;   
  TH1F * h_elmu_os_minv_C1;
  
  TH1F * h_mu_os_asym_C1;
  TH1F * h_el_os_asym_C1;  
  TH1F * h_elmu_os_asym_C1;
  
  TH1F  * h_elmu_os_dphi_C1;
  TH1F  * h_elmu_os_dtheta_C1;
  TH1F  * h_elmu_jet_dphi_C1[2];
  TH1F  * h_elmu_etm_dphi_C1[2];
  
  TH1F  * h_el_os_dphi_C1;
  TH1F  * h_el_os_dtheta_C1;
  TH1F  * h_el_jet_dphi_C1[2];
  TH1F  * h_el_etm_dphi_C1[2];
  
  
  TH1F  * h_mu_os_dphi_C1;
  TH1F  * h_mu_os_dtheta_C1;
  TH1F  * h_mu_jet_dphi_C1[2];
  TH1F  * h_mu_etm_dphi_C1[2];
  
  
  TH1F  * h_el_os_pt_C1[2];
  TH1F  * h_el_os_eta_C1[2];
  TH1F  * h_mu_os_pt_C1[2];
  TH1F  * h_mu_os_eta_C1[2];
  TH1F  * h_elmu_os_pt_C1[2];
  TH1F  * h_elmu_os_eta_C1[2];

  //deta for leptons
  TH1F  * h_el_os_deta_C1;
  TH1F  * h_mu_os_deta_C1;
  TH1F  * h_elmu_os_deta_C1;
  


//C2


  TH1F  * h_etmiss_e_C2;
  TH1F  * h_etmiss_mu_C2;
  TH1F  * h_etmiss_emu_C2;
  
  TH1F  * h_mu_mt2_C2;
  TH1F  * h_el_mt2_C2;
  TH1F  * h_elmu_mt2_C2;
  TH1F  * h_jet_mt2_C2; 
  
  TH1F * h_mu_os_minv_C2;
  TH1F * h_el_os_minv_C2;   
  TH1F * h_elmu_os_minv_C2;
  
  TH1F * h_mu_os_asym_C2;
  TH1F * h_el_os_asym_C2;  
  TH1F * h_elmu_os_asym_C2;
  
  TH1F  * h_elmu_os_dphi_C2;
  TH1F  * h_elmu_os_dtheta_C2;
  TH1F  * h_elmu_jet_dphi_C2[2];
  TH1F  * h_elmu_etm_dphi_C2[2];
  
  TH1F  * h_el_os_dphi_C2;
  TH1F  * h_el_os_dtheta_C2;
  TH1F  * h_el_jet_dphi_C2[2];
  TH1F  * h_el_etm_dphi_C2[2];
  
  
  TH1F  * h_mu_os_dphi_C2;
  TH1F  * h_mu_os_dtheta_C2;
  TH1F  * h_mu_jet_dphi_C2[2];
  TH1F  * h_mu_etm_dphi_C2[2];
  
  
  TH1F  * h_el_os_pt_C2[2];
  TH1F  * h_el_os_eta_C2[2];
  TH1F  * h_mu_os_pt_C2[2];
  TH1F  * h_mu_os_eta_C2[2];
  TH1F  * h_elmu_os_pt_C2[2];
  TH1F  * h_elmu_os_eta_C2[2];

  //deta for leptons
  TH1F  * h_el_os_deta_C2;
  TH1F  * h_mu_os_deta_C2;
  TH1F  * h_elmu_os_deta_C2;


  //C3


 TH1F  * h_etmiss_e_C3;
  TH1F  * h_etmiss_mu_C3;
  TH1F  * h_etmiss_emu_C3;
  
  TH1F  * h_mu_mt2_C3;
  TH1F  * h_el_mt2_C3;
  TH1F  * h_elmu_mt2_C3;
  TH1F  * h_jet_mt2_C3;
  
  TH1F * h_mu_os_minv_C3;
  TH1F * h_el_os_minv_C3;   
  TH1F * h_elmu_os_minv_C3;
  
  TH1F * h_mu_os_asym_C3;
  TH1F * h_el_os_asym_C3;  
  TH1F * h_elmu_os_asym_C3;
  
  TH1F  * h_elmu_os_dphi_C3;
  TH1F  * h_elmu_os_dtheta_C3;
  TH1F  * h_elmu_jet_dphi_C3[2];
  TH1F  * h_elmu_etm_dphi_C3[2];
  
  TH1F  * h_el_os_dphi_C3;
  TH1F  * h_el_os_dtheta_C3;
  TH1F  * h_el_jet_dphi_C3[2];
  TH1F  * h_el_etm_dphi_C3[2];
  
  
  TH1F  * h_mu_os_dphi_C3;
  TH1F  * h_mu_os_dtheta_C3;
  TH1F  * h_mu_jet_dphi_C3[2];
  TH1F  * h_mu_etm_dphi_C3[2];
  
  
  TH1F  * h_el_os_pt_C3[2];
  TH1F  * h_el_os_eta_C3[2];
  TH1F  * h_mu_os_pt_C3[2];
  TH1F  * h_mu_os_eta_C3[2];
  TH1F  * h_elmu_os_pt_C3[2];
  TH1F  * h_elmu_os_eta_C3[2];

  //deta for leptons                                                                              
  TH1F  * h_el_os_deta_C3;
  TH1F  * h_mu_os_deta_C3;
  TH1F  * h_elmu_os_deta_C3;



  //C4


 TH1F  * h_etmiss_e_C4;
  TH1F  * h_etmiss_mu_C4;
  TH1F  * h_etmiss_emu_C4;
  
  TH1F  * h_mu_mt2_C4;
  TH1F  * h_el_mt2_C4;
  TH1F  * h_elmu_mt2_C4;
  TH1F  * h_jet_mt2_C4;
  
  
  TH1F * h_mu_os_minv_C4;
  TH1F * h_el_os_minv_C4;   
  TH1F * h_elmu_os_minv_C4;
  
  TH1F * h_mu_os_asym_C4;
  TH1F * h_el_os_asym_C4;  
  TH1F * h_elmu_os_asym_C4;
  
  TH1F  * h_elmu_os_dphi_C4;
  TH1F  * h_elmu_os_dtheta_C4;
  TH1F  * h_elmu_jet_dphi_C4[2];
  TH1F  * h_elmu_etm_dphi_C4[2];
  
  TH1F  * h_el_os_dphi_C4;
  TH1F  * h_el_os_dtheta_C4;
  TH1F  * h_el_jet_dphi_C4[2];
  TH1F  * h_el_etm_dphi_C4[2];
  
  
  TH1F  * h_mu_os_dphi_C4;
  TH1F  * h_mu_os_dtheta_C4;
  TH1F  * h_mu_jet_dphi_C4[2];
  TH1F  * h_mu_etm_dphi_C4[2];
  
  
  TH1F  * h_el_os_pt_C4[2];
  TH1F  * h_el_os_eta_C4[2];
  TH1F  * h_mu_os_pt_C4[2];
  TH1F  * h_mu_os_eta_C4[2];
  TH1F  * h_elmu_os_pt_C4[2];
  TH1F  * h_elmu_os_eta_C4[2];

  //deta for leptons                                                                              
  TH1F  * h_el_os_deta_C4;
  TH1F  * h_mu_os_deta_C4;
  TH1F  * h_elmu_os_deta_C4;


 



 

  
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
  // Declaration of leaf types
  Bool_t          EF_2e12Tvh_loose1;
  Bool_t          EF_2mu13;
  Bool_t          EF_e12Tvh_medium1_mu8;
  Bool_t          EF_e24vh_medium1;
  Bool_t          EF_e24vh_medium1_e7_medium1;
  Bool_t          EF_e24vhi_medium1;
  Bool_t          EF_e60_medium1;
  Bool_t          EF_mu18_tight;
  Bool_t          EF_mu18_tight_e7_medium1;
  Bool_t          EF_mu18_tight_mu8_EFFS;
  Bool_t          EF_mu20it_tight;
  Bool_t          EF_mu24_tight;
  Bool_t          EF_mu24_tight_mu6_EFFS;
  Bool_t          EF_mu24i_tight;
  Bool_t          EF_mu36_tight;
  Bool_t          EF_xe80T_tclcw_loose;
  Bool_t          EF_xe80_tclcw_loose;
  
  vector<int>     *trig_EF_el_EF_2e12Tvh_loose1;
  vector<int>     *trig_EF_el_EF_e12Tvh_loose1;
  vector<int>     *trig_EF_el_EF_e12Tvh_medium1;
  vector<int>     *trig_EF_el_EF_e12Tvh_medium1_mu8;
  vector<int>     *trig_EF_el_EF_e24vh_medium1;
  vector<int>     *trig_EF_el_EF_e24vh_medium1_e7_medium1;
  vector<int>     *trig_EF_el_EF_e24vhi_medium1;
  vector<int>     *trig_EF_el_EF_e60_medium1;
  vector<int>     *trig_EF_el_EF_e7T_medium1;
  vector<int>     *trig_EF_trigmuonef_EF_2mu13;
  vector<int>     *trig_EF_trigmuonef_EF_mu13;
  vector<int>     *trig_EF_trigmuonef_EF_mu18_tight;
  vector<int>     *trig_EF_trigmuonef_EF_mu18_tight_e7_medium1;
  vector<int>     *trig_EF_trigmuonef_EF_mu18_tight_mu8_EFFS;
  vector<int>     *trig_EF_trigmuonef_EF_mu20it_tight;
  vector<int>     *trig_EF_trigmuonef_EF_mu24_tight;
  vector<int>     *trig_EF_trigmuonef_EF_mu24_tight_mu6_EFFS;
  vector<int>     *trig_EF_trigmuonef_EF_mu24i_tight;
  vector<int>     *trig_EF_trigmuonef_EF_mu36_tight;
  vector<int>     *trig_EF_trigmuonef_EF_mu8;
  
  UInt_t          RunNumber;
  UInt_t          EventNumber;
  UInt_t          timestamp;
  UInt_t          timestamp_ns;
  UInt_t          lbn;
  UInt_t          bcid;
  Float_t         actualIntPerXing;
  Float_t         averageIntPerXing;
  UInt_t          mc_channel_number;
  UInt_t          coreFlags;
  UInt_t          pixelError;
  UInt_t          sctError;
  UInt_t          trtError;
  UInt_t          larError;
  UInt_t          tileError;
  UInt_t          muonError;
  UInt_t          fwdError;
  UInt_t          coreError;
  Bool_t          isSimulation;
  Int_t           el_n;
  vector<float>   *el_E;
  vector<float>   *el_Et;
  vector<float>   *el_pt;
  vector<float>   *el_eta;
  vector<float>   *el_phi;
  vector<float>   *el_charge;
  vector<int>     *el_author;
  vector<unsigned int> *el_isEM;
  vector<unsigned int> *el_OQ;
  vector<int>     *el_nConv;
  vector<int>     *el_origin;
  vector<float>   *el_truth_E;
  vector<float>   *el_truth_pt;
  vector<float>   *el_truth_eta;
  vector<float>   *el_truth_phi;
  vector<int>     *el_truth_type;
  vector<int>     *el_truth_status;
  vector<int>     *el_truth_barcode;
  vector<int>     *el_truth_mothertype;
  vector<int>     *el_truth_motherbarcode;
  vector<int>     *el_truth_hasHardBrem;
  vector<int>     *el_truth_matched;
  vector<int>     *el_loose;
  vector<int>     *el_medium;
  vector<int>     *el_mediumIso;
  vector<int>     *el_tight;
  vector<int>     *el_tightIso;
  vector<int>     *el_mediumPP;
  vector<int>     *el_tightPP;
  vector<float>   *el_Ethad;
  vector<float>   *el_Ethad1;
  vector<float>   *el_f1;
  vector<float>   *el_Emax2;
  vector<float>   *el_wstot;
  vector<float>   *el_emaxs1;
  vector<float>   *el_E237;
  vector<float>   *el_E277;
  vector<float>   *el_weta2;
  vector<float>   *el_f3;
  vector<float>   *el_Etcone20;
  vector<float>   *el_Etcone30;
  vector<float>   *el_Etcone40;
  vector<float>   *el_ptcone20;
  vector<float>   *el_ptcone30;
  vector<float>   *el_topoEtcone30_corrected;
  vector<float>   *el_deltaeta1;
  vector<float>   *el_deltaphi2;
  vector<float>   *el_expectHitInBLayer;
  vector<float>   *el_trackd0_physics;
  vector<float>   *el_reta;
  vector<float>   *el_etas2;
  vector<float>   *el_cl_E;
  vector<float>   *el_cl_pt;
  vector<float>   *el_cl_eta;
  vector<float>   *el_cl_phi;
  vector<float>   *el_trackd0;
  vector<float>   *el_trackz0;
  vector<float>   *el_trackIPEstimate_d0_unbiasedpvunbiased;
  vector<float>   *el_trackIPEstimate_z0_unbiasedpvunbiased;
  vector<float>   *el_trackIPEstimate_sigd0_unbiasedpvunbiased;
  vector<float>   *el_trackphi;
  vector<float>   *el_trackqoverp;
  vector<float>   *el_trackpt;
  vector<float>   *el_tracketa;
  vector<int>     *el_nBLHits;
  vector<int>     *el_nPixHits;
  vector<int>     *el_nSCTHits;
  vector<int>     *el_nTRTHits;
  vector<int>     *el_nBLayerOutliers;
  vector<int>     *el_nPixelOutliers;
  vector<int>     *el_nSCTOutliers;
  vector<int>     *el_nTRTOutliers;
  vector<int>     *el_nSiHits;
  vector<float>   *el_TRTHighTOutliersRatio;
  vector<float>   *el_trackd0pv;
  vector<float>   *el_trackz0pv;
  vector<int>     *el_hastrack;
  Int_t           mu_staco_n;
  vector<float>   *mu_staco_E;
  vector<float>   *mu_staco_pt;
  vector<float>   *mu_staco_eta;
  vector<float>   *mu_staco_phi;
  vector<float>   *mu_staco_px;
  vector<float>   *mu_staco_py;
  vector<float>   *mu_staco_pz;
  vector<float>   *mu_staco_charge;
  vector<int>     *mu_staco_author;
  vector<float>   *mu_staco_matchchi2;
  vector<float>   *mu_staco_etcone20;
  vector<float>   *mu_staco_etcone30;
  vector<float>   *mu_staco_etcone40;
  vector<float>   *mu_staco_ptcone20;
  vector<float>   *mu_staco_ptcone30;
  vector<float>   *mu_staco_ptcone30_trkelstyle;
  vector<float>   *mu_staco_ptcone40;
  vector<int>     *mu_staco_bestMatch;
  vector<int>     *mu_staco_isCombinedMuon;
  vector<int>     *mu_staco_isLowPtReconstructedMuon;
  vector<int>     *mu_staco_isSegmentTaggedMuon;
  vector<int>     *mu_staco_loose;
  vector<int>     *mu_staco_medium;
  vector<int>     *mu_staco_tight;
  vector<float>   *mu_staco_d0_exPV;
  vector<float>   *mu_staco_z0_exPV;
  vector<float>   *mu_staco_id_phi_exPV;
  vector<float>   *mu_staco_id_theta_exPV;
  vector<float>   *mu_staco_id_qoverp_exPV;
  vector<float>   *mu_staco_me_phi_exPV;
  vector<float>   *mu_staco_me_theta_exPV;
  vector<float>   *mu_staco_me_qoverp_exPV;
  vector<float>   *mu_staco_ms_phi;
  vector<float>   *mu_staco_ms_theta;
  vector<float>   *mu_staco_ms_qoverp;
  vector<float>   *mu_staco_id_phi;
  vector<float>   *mu_staco_id_theta;
  vector<float>   *mu_staco_id_qoverp;
  vector<float>   *mu_staco_me_phi;
  vector<float>   *mu_staco_me_theta;
  vector<float>   *mu_staco_me_qoverp;
  vector<float>   *mu_staco_ie_phi;
  vector<float>   *mu_staco_ie_theta;
  vector<float>   *mu_staco_ie_qoverp;
  vector<int>     *mu_staco_nBLHits;
  vector<int>     *mu_staco_nPixHits;
  vector<int>     *mu_staco_nSCTHits;
  vector<int>     *mu_staco_nTRTHits;
  vector<int>     *mu_staco_nPixHoles;
  vector<int>     *mu_staco_nSCTHoles;
  vector<int>     *mu_staco_nTRTOutliers;
  vector<int>     *mu_staco_nPixelDeadSensors;
  vector<int>     *mu_staco_nSCTDeadSensors;
  vector<int>     *mu_staco_expectBLayerHit;
  vector<int>     *mu_staco_nMDTHits;
  vector<int>     *mu_staco_nCSCEtaHits;
  vector<int>     *mu_staco_nCSCPhiHits;
  vector<int>     *mu_staco_nRPCEtaHits;
  vector<int>     *mu_staco_nRPCPhiHits;
  vector<int>     *mu_staco_nTGCEtaHits;
  vector<int>     *mu_staco_nTGCPhiHits;
  vector<float>   *mu_staco_trackd0;
  vector<float>   *mu_staco_trackz0;
  vector<float>   *mu_staco_trackphi;
  vector<float>   *mu_staco_tracktheta;
  vector<float>   *mu_staco_trackqoverp;
  vector<float>   *mu_staco_trackcov_d0;
  vector<float>   *mu_staco_trackcov_z0;
  vector<float>   *mu_staco_trackIPEstimate_d0_unbiasedpvunbiased;
  vector<float>   *mu_staco_trackIPEstimate_z0_unbiasedpvunbiased;
  vector<float>   *mu_staco_trackIPEstimate_sigd0_unbiasedpvunbiased;
  vector<float>   *mu_staco_truth_dr;
  vector<float>   *mu_staco_truth_E;
  vector<float>   *mu_staco_truth_pt;
  vector<float>   *mu_staco_truth_eta;
  vector<float>   *mu_staco_truth_phi;
  vector<int>     *mu_staco_truth_type;
  vector<int>     *mu_staco_truth_status;
  vector<int>     *mu_staco_truth_barcode;
  vector<int>     *mu_staco_truth_mothertype;
  vector<int>     *mu_staco_truth_motherbarcode;
  vector<int>     *mu_staco_truth_matched;
  Float_t         MET_RefFinal_etx;
  Float_t         MET_RefFinal_ety;
  Float_t         MET_RefFinal_phi;
  Float_t         MET_RefFinal_et;
  Float_t         MET_RefFinal_sumet;
  Float_t         MET_SoftJets_etx;
  Float_t         MET_SoftJets_ety;
  Float_t         MET_SoftJets_sumet;
  Float_t         MET_CellOut_Eflow_etx;
  Float_t         MET_CellOut_Eflow_ety;
  Float_t         MET_CellOut_Eflow_sumet;
  Float_t         MET_Truth_NonInt_etx;
  Float_t         MET_Truth_NonInt_ety;
  Float_t         MET_Truth_NonInt_sumet;
  Float_t         MET_Egamma10NoTau_RefGamma_etx;
  Float_t         MET_Egamma10NoTau_RefGamma_ety;
  Float_t         MET_Egamma10NoTau_RefGamma_sumet;
  Float_t         MET_Egamma10NoTau_RefFinal_sumet;
  Float_t         MET_Egamma10NoTau_CellOut_etx;
  Float_t         MET_Egamma10NoTau_CellOut_ety;
  Float_t         MET_Egamma10NoTau_CellOut_sumet;
  Float_t         MET_Egamma10NoTau_CellOut_Eflow_STVF_etx;
  Float_t         MET_Egamma10NoTau_CellOut_Eflow_STVF_ety;
  Float_t         MET_Egamma10NoTau_CellOut_Eflow_STVF_sumet;
  Float_t         MET_Egamma10NoTau_SoftJets_etx;
  Float_t         MET_Egamma10NoTau_SoftJets_ety;
  Float_t         MET_Egamma10NoTau_SoftJets_sumet;
  vector<vector<float> > *el_MET_Egamma10NoTau_wpx;
  vector<vector<float> > *el_MET_Egamma10NoTau_wpy;
  vector<vector<float> > *el_MET_Egamma10NoTau_wet;
  vector<vector<unsigned int> > *el_MET_Egamma10NoTau_statusWord;
  vector<vector<float> > *jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpx;
  vector<vector<float> > *jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpy;
  vector<vector<float> > *jet_AntiKt4LCTopo_MET_Egamma10NoTau_wet;
  vector<vector<unsigned int> > *jet_AntiKt4LCTopo_MET_Egamma10NoTau_statusWord;
  Int_t           jet_AntiKt4LCTopo_n;
  vector<float>   *jet_AntiKt4LCTopo_E;
  vector<float>   *jet_AntiKt4LCTopo_pt;
  vector<float>   *jet_AntiKt4LCTopo_m;
  vector<float>   *jet_AntiKt4LCTopo_eta;
  vector<float>   *jet_AntiKt4LCTopo_phi;
  vector<float>   *jet_AntiKt4LCTopo_EtaOrigin;
  vector<float>   *jet_AntiKt4LCTopo_PhiOrigin;
  vector<float>   *jet_AntiKt4LCTopo_MOrigin;
  vector<float>   *jet_AntiKt4LCTopo_WIDTH;
  vector<float>   *jet_AntiKt4LCTopo_n90;
  vector<float>   *jet_AntiKt4LCTopo_Timing;
  vector<float>   *jet_AntiKt4LCTopo_LArQuality;
  //  vector<float>   *jet_AntiKt4LCTopo_nTrk;
  vector<float>   *jet_AntiKt4LCTopo_sumPtTrk;
  //   vector<float>   *jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV;
  vector<float>   *jet_AntiKt4LCTopo_HECQuality;
  vector<float>   *jet_AntiKt4LCTopo_NegativeE;
  vector<float>   *jet_AntiKt4LCTopo_AverageLArQF;
  vector<float>   *jet_AntiKt4LCTopo_BCH_CORR_CELL;
  vector<float>   *jet_AntiKt4LCTopo_BCH_CORR_DOTX;
  vector<float>   *jet_AntiKt4LCTopo_BCH_CORR_JET;
  vector<int>     *jet_AntiKt4LCTopo_SamplingMax;
  vector<float>   *jet_AntiKt4LCTopo_fracSamplingMax;
  vector<float>   *jet_AntiKt4LCTopo_hecf;
  vector<float>   *jet_AntiKt4LCTopo_tgap3f;
  vector<float>   *jet_AntiKt4LCTopo_emfrac;
  vector<float>   *jet_AntiKt4LCTopo_EMJES;
  vector<float>   *jet_AntiKt4LCTopo_emscale_E;
  vector<float>   *jet_AntiKt4LCTopo_emscale_pt;
  vector<float>   *jet_AntiKt4LCTopo_emscale_m;
  vector<float>   *jet_AntiKt4LCTopo_emscale_eta;
  vector<float>   *jet_AntiKt4LCTopo_emscale_phi;
  vector<float>   *jet_AntiKt4LCTopo_ActiveAreaPx;
  vector<float>   *jet_AntiKt4LCTopo_ActiveAreaPy;
  vector<float>   *jet_AntiKt4LCTopo_ActiveAreaPz;
  vector<float>   *jet_AntiKt4LCTopo_ActiveAreaE;
  vector<float>   *jet_AntiKt4LCTopo_jvtxf;
  vector<float>   *jet_AntiKt4LCTopo_constscale_E;
  vector<float>   *jet_AntiKt4LCTopo_constscale_m;
  vector<float>   *jet_AntiKt4LCTopo_constscale_eta;
  vector<float>   *jet_AntiKt4LCTopo_constscale_phi;
  vector<float>   *jet_AntiKt4LCTopo_flavor_weight_Comb;
  vector<float>   *jet_AntiKt4LCTopo_flavor_weight_IP3D;
  vector<float>   *jet_AntiKt4LCTopo_flavor_weight_SV0;
  vector<float>   *jet_AntiKt4LCTopo_flavor_weight_SV1;
  vector<float>   *jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN;
  vector<float>   *jet_AntiKt4LCTopo_flavor_weight_MV1;
  vector<int>     *jet_AntiKt4LCTopo_flavor_truth_label;
  Float_t         Eventshape_rhoKt4LC;
  Int_t           vx_n;
  vector<float>   *vx_x;
  vector<float>   *vx_y;
  vector<float>   *vx_z;
  vector<float>   *vx_px;
  vector<float>   *vx_py;
  vector<float>   *vx_pz;
  vector<float>   *vx_E;
  vector<float>   *vx_m;
  vector<int>     *vx_nTracks;
  vector<float>   *vx_sumPt;
  Int_t           top_hfor_type;
  Int_t           jet_AntiKt4TruthJets_n;
  vector<float>   *jet_AntiKt4TruthJets_pt;
  vector<float>   *jet_AntiKt4TruthJets_m;
  vector<float>   *jet_AntiKt4TruthJets_eta;
  vector<float>   *jet_AntiKt4TruthJets_phi;
  vector<int>     *jet_AntiKt4TruthJets_flavor_truth_label;
  Int_t           trig_EF_el_n;
  vector<float>   *trig_EF_el_pt;
  vector<float>   *trig_EF_el_eta;
  vector<float>   *trig_EF_el_phi;
  vector<unsigned int> *trig_L1_TAV;
  Int_t           trig_EF_trigmuonef_n;
  vector<int>     *trig_EF_trigmuonef_track_n;
  vector<vector<float> > *trig_EF_trigmuonef_track_CB_pt;
  vector<vector<float> > *trig_EF_trigmuonef_track_CB_eta;
  vector<vector<float> > *trig_EF_trigmuonef_track_CB_phi;
  vector<vector<int> > *trig_EF_trigmuonef_track_CB_hasCB;
  Int_t           jet_AntiKt4TrackZ_n;
  vector<float>   *jet_AntiKt4TrackZ_pt;
  vector<float>   *jet_AntiKt4TrackZ_eta;
  vector<float>   *jet_AntiKt4TrackZ_phi;
  Int_t           mc_n;
  vector<float>   *mc_pt;
  vector<float>   *mc_m;
  vector<float>   *mc_eta;
  vector<float>   *mc_phi;
  vector<int>     *mc_status;
  vector<int>     *mc_barcode;
  vector<int>     *mc_oldindex;
  vector<int>     *mc_pdgId;
  vector<float>   *mc_charge;
  vector<float>   *mc_vx_x;
  vector<float>   *mc_vx_y;
  vector<float>   *mc_vx_z;
  vector<vector<int> > *mc_children;
  vector<vector<int> > *mc_parents;
  vector<vector<int> > *mc_child_index;
  vector<vector<int> > *mc_parent_index;
  Int_t           mcevt_n;
  vector<int>     *mcevt_signal_process_id;
  vector<int>     *mcevt_event_number;
  vector<double>  *mcevt_event_scale;
  vector<double>  *mcevt_alphaQCD;
  vector<double>  *mcevt_alphaQED;
  vector<int>     *mcevt_pdf_id1;
  vector<int>     *mcevt_pdf_id2;
  vector<double>  *mcevt_pdf_x1;
  vector<double>  *mcevt_pdf_x2;
  vector<double>  *mcevt_pdf_scale;
  vector<double>  *mcevt_pdf1;
  vector<double>  *mcevt_pdf2;
  vector<vector<double> > *mcevt_weight;
  Int_t           idp1;
  Int_t           idp2;

  // List of branches
  TBranch        *b_EF_2e12Tvh_loose1;   //!
  TBranch        *b_EF_2mu13;   //!
  TBranch        *b_EF_e12Tvh_medium1_mu8;   //!
  TBranch        *b_EF_e24vh_medium1;   //!
  TBranch        *b_EF_e24vh_medium1_e7_medium1;   //!
  TBranch        *b_EF_e24vhi_medium1;   //!
  TBranch        *b_EF_e60_medium1;   //!
  TBranch        *b_EF_mu18_tight;   //!
  TBranch        *b_EF_mu18_tight_e7_medium1;   //!
  TBranch        *b_EF_mu18_tight_mu8_EFFS;   //!
  TBranch        *b_EF_mu20it_tight;   //!
  TBranch        *b_EF_mu24_tight;   //!
  TBranch        *b_EF_mu24_tight_mu6_EFFS;   //!
  TBranch        *b_EF_mu24i_tight;   //!
  TBranch        *b_EF_mu36_tight;   //!
  TBranch         *b_EF_xe80_tclcw_loose;
  TBranch         *b_EF_xe80T_tclcw_loose;
  
  TBranch        *b_trig_EF_el_EF_2e12Tvh_loose1;   //!
  TBranch        *b_trig_EF_el_EF_e12Tvh_loose1;   //!
  TBranch        *b_trig_EF_el_EF_e12Tvh_medium1;   //!
  TBranch        *b_trig_EF_el_EF_e12Tvh_medium1_mu8;   //!
  TBranch        *b_trig_EF_el_EF_e24vh_medium1;   //!
  TBranch        *b_trig_EF_el_EF_e24vh_medium1_e7_medium1;   //!
  TBranch        *b_trig_EF_el_EF_e24vhi_medium1;   //!
  TBranch        *b_trig_EF_el_EF_e60_medium1;   //!
  TBranch        *b_trig_EF_el_EF_e7T_medium1;   //!
  TBranch        *b_trig_EF_trigmuonef_EF_2mu13;   //!
  TBranch        *b_trig_EF_trigmuonef_EF_mu13;   //!
  TBranch        *b_trig_EF_trigmuonef_EF_mu18_tight;   //!
  TBranch        *b_trig_EF_trigmuonef_EF_mu18_tight_e7_medium1;   //!
  TBranch        *b_trig_EF_trigmuonef_EF_mu18_tight_mu8_EFFS;   //!
  TBranch        *b_trig_EF_trigmuonef_EF_mu20it_tight;   //!
  TBranch        *b_trig_EF_trigmuonef_EF_mu24_tight;   //!
  TBranch        *b_trig_EF_trigmuonef_EF_mu24_tight_mu6_EFFS;   //!
  TBranch        *b_trig_EF_trigmuonef_EF_mu24i_tight;   //!
  TBranch        *b_trig_EF_trigmuonef_EF_mu36_tight;   //!
  TBranch        *b_trig_EF_trigmuonef_EF_mu8;   //!
  TBranch        *b_RunNumber;   //!
  TBranch        *b_EventNumber;   //!
  TBranch        *b_timestamp;   //!
  TBranch        *b_timestamp_ns;   //!
  TBranch        *b_lbn;   //!
  TBranch        *b_bcid;   //!
  TBranch        *b_actualIntPerXing;   //!
  TBranch        *b_averageIntPerXing;   //!
  TBranch        *b_mc_channel_number;   //!
  TBranch        *b_coreFlags;   //!
  TBranch        *b_pixelError;   //!
  TBranch        *b_sctError;   //!
  TBranch        *b_trtError;   //!
  TBranch        *b_larError;   //!
  TBranch        *b_tileError;   //!
  TBranch        *b_muonError;   //!
  TBranch        *b_fwdError;   //!
  TBranch        *b_coreError;   //!
  TBranch        *b_isSimulation;   //!
  TBranch        *b_el_n;   //!
  TBranch        *b_el_E;   //!
  TBranch        *b_el_Et;   //!
  TBranch        *b_el_pt;   //!
  TBranch        *b_el_eta;   //!
  TBranch        *b_el_phi;   //!
  TBranch        *b_el_charge;   //!
  TBranch        *b_el_author;   //!
  TBranch        *b_el_isEM;   //!
  TBranch        *b_el_OQ;   //!
  TBranch        *b_el_nConv;   //!
  TBranch        *b_el_origin;   //!
  TBranch        *b_el_truth_E;   //!
  TBranch        *b_el_truth_pt;   //!
  TBranch        *b_el_truth_eta;   //!
  TBranch        *b_el_truth_phi;   //!
  TBranch        *b_el_truth_type;   //!
  TBranch        *b_el_truth_status;   //!
  TBranch        *b_el_truth_barcode;   //!
  TBranch        *b_el_truth_mothertype;   //!
  TBranch        *b_el_truth_motherbarcode;   //!
  TBranch        *b_el_truth_hasHardBrem;   //!
  TBranch        *b_el_truth_matched;   //!
  TBranch        *b_el_loose;   //!
  TBranch        *b_el_medium;   //!
  TBranch        *b_el_mediumIso;   //!
  TBranch        *b_el_tight;   //!
  TBranch        *b_el_tightIso;   //!
  TBranch        *b_el_mediumPP;   //!
  TBranch        *b_el_tightPP;   //!
  TBranch        *b_el_Ethad;   //!
  TBranch        *b_el_Ethad1;   //!
  TBranch        *b_el_f1;   //!
  TBranch        *b_el_Emax2;   //!
  TBranch        *b_el_wstot;   //!
  TBranch        *b_el_emaxs1;   //!
  TBranch        *b_el_E237;   //!
  TBranch        *b_el_E277;   //!
  TBranch        *b_el_weta2;   //!
  TBranch        *b_el_f3;   //!
  TBranch        *b_el_Etcone20;   //!
  TBranch        *b_el_Etcone30;   //!
  TBranch        *b_el_Etcone40;   //!
  TBranch        *b_el_ptcone20;   //!
  TBranch        *b_el_ptcone30;   //!
  TBranch        *b_el_deltaeta1;   //!
  TBranch        *b_el_deltaphi2;   //!
  TBranch        *b_el_expectHitInBLayer;   //!
  TBranch        *b_el_trackd0_physics;   //!
  TBranch        *b_el_reta;   //!
  TBranch        *b_el_etas2;   //!
  TBranch        *b_el_cl_E;   //!
  TBranch        *b_el_cl_pt;   //!
  TBranch        *b_el_cl_eta;   //!
  TBranch        *b_el_cl_phi;   //!
  TBranch        *b_el_trackd0;   //!
  TBranch        *b_el_trackz0;   //!
  TBranch        *b_el_trackphi;   //!
  TBranch        *b_el_trackqoverp;   //!
  TBranch        *b_el_trackpt;   //!
  TBranch        *b_el_tracketa;   //!
  TBranch        *b_el_nBLHits;   //!
  TBranch        *b_el_nPixHits;   //!
  TBranch        *b_el_nSCTHits;   //!
  TBranch        *b_el_nTRTHits;   //!
  TBranch        *b_el_nBLayerOutliers;   //!
  TBranch        *b_el_nPixelOutliers;   //!
  TBranch        *b_el_nSCTOutliers;   //!
  TBranch        *b_el_nTRTOutliers;   //!
  TBranch        *b_el_nSiHits;   //!
  TBranch        *b_el_TRTHighTOutliersRatio;   //!
  TBranch        *b_el_trackd0pv;   //!
  TBranch        *b_el_trackz0pv;   //!
  TBranch        *b_el_trackIPEstimate_d0_unbiasedpvunbiased;   //!                                                                                
  TBranch        *b_el_trackIPEstimate_z0_unbiasedpvunbiased;   //!                                                                                
  TBranch        *b_el_trackIPEstimate_sigd0_unbiasedpvunbiased;   //!    
  TBranch        *b_el_hastrack;   //!
  TBranch        *b_el_topoEtcone30_corrected;   //!
  TBranch        *b_mu_staco_n;   //!
  TBranch        *b_mu_staco_E;   //!
  TBranch        *b_mu_staco_pt;   //!
  TBranch        *b_mu_staco_eta;   //!
  TBranch        *b_mu_staco_phi;   //!
  TBranch        *b_mu_staco_px;   //!
  TBranch        *b_mu_staco_py;   //!
  TBranch        *b_mu_staco_pz;   //!
  TBranch        *b_mu_staco_charge;   //!
  TBranch        *b_mu_staco_author;   //!
  TBranch        *b_mu_staco_matchchi2;   //!
  TBranch        *b_mu_staco_etcone20;   //!
  TBranch        *b_mu_staco_etcone30;   //!
  TBranch        *b_mu_staco_etcone40;   //!
  TBranch        *b_mu_staco_ptcone20;   //!
  TBranch        *b_mu_staco_ptcone30;   //!
  TBranch        *b_mu_staco_ptcone30_trkelstyle;   //!
  TBranch        *b_mu_staco_ptcone40;   //!
  TBranch        *b_mu_staco_bestMatch;   //!
  TBranch        *b_mu_staco_isCombinedMuon;   //!
  TBranch        *b_mu_staco_isLowPtReconstructedMuon;   //!
  TBranch        *b_mu_staco_isSegmentTaggedMuon;   //!
  TBranch        *b_mu_staco_loose;   //!
  TBranch        *b_mu_staco_medium;   //!
  TBranch        *b_mu_staco_tight;   //!
  TBranch        *b_mu_staco_d0_exPV;   //!
  TBranch        *b_mu_staco_z0_exPV;   //!
  TBranch        *b_mu_staco_id_phi_exPV;   //!
  TBranch        *b_mu_staco_id_theta_exPV;   //!
  TBranch        *b_mu_staco_id_qoverp_exPV;   //!
  TBranch        *b_mu_staco_me_phi_exPV;   //!
  TBranch        *b_mu_staco_me_theta_exPV;   //!
  TBranch        *b_mu_staco_me_qoverp_exPV;   //!
  TBranch        *b_mu_staco_ms_phi;   //!
  TBranch        *b_mu_staco_ms_theta;   //!
  TBranch        *b_mu_staco_ms_qoverp;   //!
  TBranch        *b_mu_staco_id_phi;   //!
  TBranch        *b_mu_staco_id_theta;   //!
  TBranch        *b_mu_staco_id_qoverp;   //!
  TBranch        *b_mu_staco_me_phi;   //!
  TBranch        *b_mu_staco_me_theta;   //!
  TBranch        *b_mu_staco_me_qoverp;   //!
  TBranch        *b_mu_staco_ie_phi;   //!
  TBranch        *b_mu_staco_ie_theta;   //!
  TBranch        *b_mu_staco_ie_qoverp;   //!
  TBranch        *b_mu_staco_nBLHits;   //!
  TBranch        *b_mu_staco_nPixHits;   //!
  TBranch        *b_mu_staco_nSCTHits;   //!
  TBranch        *b_mu_staco_nTRTHits;   //!
  TBranch        *b_mu_staco_nPixHoles;   //!
  TBranch        *b_mu_staco_nSCTHoles;   //!
  TBranch        *b_mu_staco_nTRTOutliers;   //!
  TBranch        *b_mu_staco_nPixelDeadSensors;   //!
  TBranch        *b_mu_staco_nSCTDeadSensors;   //!
  TBranch        *b_mu_staco_expectBLayerHit;   //!
  TBranch        *b_mu_staco_nMDTHits;   //!
  TBranch        *b_mu_staco_nCSCEtaHits;   //!
  TBranch        *b_mu_staco_nCSCPhiHits;   //!
  TBranch        *b_mu_staco_nRPCEtaHits;   //!
  TBranch        *b_mu_staco_nRPCPhiHits;   //!
  TBranch        *b_mu_staco_nTGCEtaHits;   //!
  TBranch        *b_mu_staco_nTGCPhiHits;   //!
  TBranch        *b_mu_staco_trackd0;   //!
  TBranch        *b_mu_staco_trackz0;   //!
  TBranch        *b_mu_staco_trackphi;   //!
  TBranch        *b_mu_staco_tracktheta;   //!
  TBranch        *b_mu_staco_trackqoverp;   //!
  TBranch        *b_mu_staco_trackcov_d0;   //!
  TBranch        *b_mu_staco_trackcov_z0;   //!
  TBranch        *b_mu_staco_trackIPEstimate_d0_unbiasedpvunbiased;   //!
  TBranch        *b_mu_staco_trackIPEstimate_z0_unbiasedpvunbiased;   //! 
  TBranch        *b_mu_staco_trackIPEstimate_sigd0_unbiasedpvunbiased;   //!
  TBranch        *b_mu_staco_truth_dr;   //!
  TBranch        *b_mu_staco_truth_E;   //!
  TBranch        *b_mu_staco_truth_pt;   //!
  TBranch        *b_mu_staco_truth_eta;   //!
  TBranch        *b_mu_staco_truth_phi;   //!
  TBranch        *b_mu_staco_truth_type;   //!
  TBranch        *b_mu_staco_truth_status;   //!
  TBranch        *b_mu_staco_truth_barcode;   //!
  TBranch        *b_mu_staco_truth_mothertype;   //!
  TBranch        *b_mu_staco_truth_motherbarcode;   //!
  TBranch        *b_mu_staco_truth_matched;   //!
  TBranch        *b_MET_RefFinal_etx;   //!
  TBranch        *b_MET_RefFinal_ety;   //!
  TBranch        *b_MET_RefFinal_phi;   //!
  TBranch        *b_MET_RefFinal_et;   //!
  TBranch        *b_MET_RefFinal_sumet;   //!
  TBranch        *b_MET_SoftJets_etx;   //!
  TBranch        *b_MET_SoftJets_ety;   //!
  TBranch        *b_MET_SoftJets_sumet;   //!
  TBranch        *b_MET_CellOut_Eflow_etx;   //!
  TBranch        *b_MET_CellOut_Eflow_ety;   //!
  TBranch        *b_MET_CellOut_Eflow_sumet;   //!
  TBranch        *b_MET_Truth_NonInt_etx;   //!
  TBranch        *b_MET_Truth_NonInt_ety;   //!
  TBranch        *b_MET_Truth_NonInt_sumet;   //!
  TBranch        *b_MET_Egamma10NoTau_RefGamma_etx;   //!
  TBranch        *b_MET_Egamma10NoTau_RefGamma_ety;   //!
  TBranch        *b_MET_Egamma10NoTau_RefGamma_sumet;   //!
  TBranch        *b_MET_Egamma10NoTau_RefFinal_sumet;   //!
  TBranch        *b_MET_Egamma10NoTau_CellOut_etx;   //!
  TBranch        *b_MET_Egamma10NoTau_CellOut_ety;   //!
  TBranch        *b_MET_Egamma10NoTau_CellOut_sumet;   //!
  TBranch        *b_MET_Egamma10NoTau_CellOut_Eflow_STVF_etx;   //!
  TBranch        *b_MET_Egamma10NoTau_CellOut_Eflow_STVF_ety;   //!
  TBranch        *b_MET_Egamma10NoTau_CellOut_Eflow_STVF_sumet;   //!
  TBranch        *b_MET_Egamma10NoTau_SoftJets_etx;   //!
  TBranch        *b_MET_Egamma10NoTau_SoftJets_ety;   //!
  TBranch        *b_MET_Egamma10NoTau_SoftJets_sumet;   //!
  TBranch        *b_el_MET_Egamma10NoTau_wpx;   //!
  TBranch        *b_el_MET_Egamma10NoTau_wpy;   //!
  TBranch        *b_el_MET_Egamma10NoTau_wet;   //!
  TBranch        *b_el_MET_Egamma10NoTau_statusWord;   //!
  TBranch        *b_jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpx;   //!
  TBranch        *b_jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpy;   //!
  TBranch        *b_jet_AntiKt4LCTopo_MET_Egamma10NoTau_wet;   //!
  TBranch        *b_jet_AntiKt4LCTopo_MET_Egamma10NoTau_statusWord;   //!
  TBranch        *b_jet_AntiKt4LCTopo_n;   //!
  TBranch        *b_jet_AntiKt4LCTopo_E;   //!
  TBranch        *b_jet_AntiKt4LCTopo_pt;   //!
  TBranch        *b_jet_AntiKt4LCTopo_m;   //!
  TBranch        *b_jet_AntiKt4LCTopo_eta;   //!
  TBranch        *b_jet_AntiKt4LCTopo_phi;   //!
  TBranch        *b_jet_AntiKt4LCTopo_EtaOrigin;   //!
  TBranch        *b_jet_AntiKt4LCTopo_PhiOrigin;   //!
  TBranch        *b_jet_AntiKt4LCTopo_MOrigin;   //!
  TBranch        *b_jet_AntiKt4LCTopo_WIDTH;   //!
  TBranch        *b_jet_AntiKt4LCTopo_n90;   //!
  TBranch        *b_jet_AntiKt4LCTopo_Timing;   //!
  TBranch        *b_jet_AntiKt4LCTopo_LArQuality;   //!
  //  TBranch        *b_jet_AntiKt4LCTopo_nTrk;   //!
  TBranch        *b_jet_AntiKt4LCTopo_sumPtTrk;   //!
  // TBranch        *b_jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV;   //!
  TBranch        *b_jet_AntiKt4LCTopo_HECQuality;   //!
  TBranch        *b_jet_AntiKt4LCTopo_NegativeE;   //!
  TBranch        *b_jet_AntiKt4LCTopo_AverageLArQF;   //!
  TBranch        *b_jet_AntiKt4LCTopo_BCH_CORR_CELL;   //!
  TBranch        *b_jet_AntiKt4LCTopo_BCH_CORR_DOTX;   //!
  TBranch        *b_jet_AntiKt4LCTopo_BCH_CORR_JET;   //!
  TBranch        *b_jet_AntiKt4LCTopo_SamplingMax;   //!
  TBranch        *b_jet_AntiKt4LCTopo_fracSamplingMax;   //!
  TBranch        *b_jet_AntiKt4LCTopo_hecf;   //!
  TBranch        *b_jet_AntiKt4LCTopo_tgap3f;   //!
  TBranch        *b_jet_AntiKt4LCTopo_emfrac;   //!
  TBranch        *b_jet_AntiKt4LCTopo_EMJES;   //!
  TBranch        *b_jet_AntiKt4LCTopo_emscale_E;   //!
  TBranch        *b_jet_AntiKt4LCTopo_emscale_pt;   //!
  TBranch        *b_jet_AntiKt4LCTopo_emscale_m;   //!
  TBranch        *b_jet_AntiKt4LCTopo_emscale_eta;   //!
  TBranch        *b_jet_AntiKt4LCTopo_emscale_phi;   //!
  TBranch        *b_jet_AntiKt4LCTopo_ActiveAreaPx;   //!
  TBranch        *b_jet_AntiKt4LCTopo_ActiveAreaPy;   //!
  TBranch        *b_jet_AntiKt4LCTopo_ActiveAreaPz;   //!
  TBranch        *b_jet_AntiKt4LCTopo_ActiveAreaE;   //!
  TBranch        *b_jet_AntiKt4LCTopo_jvtxf;   //!
  TBranch        *b_jet_AntiKt4LCTopo_constscale_E;   //!
  TBranch        *b_jet_AntiKt4LCTopo_constscale_m;   //!
  TBranch        *b_jet_AntiKt4LCTopo_constscale_eta;   //!
  TBranch        *b_jet_AntiKt4LCTopo_constscale_phi;   //!
  TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_Comb;   //!
  TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_IP3D;   //!
  TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_SV0;   //!
  TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_SV1;   //!
  TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN;   //!
  TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_MV1;   //!
  TBranch        *b_jet_AntiKt4LCTopo_flavor_truth_label;   //!
  TBranch        *b_Eventshape_rhoKt4LC;   //!
  TBranch        *b_vx_n;   //!
  TBranch        *b_vx_x;   //!
  TBranch        *b_vx_y;   //!
  TBranch        *b_vx_z;   //!
  TBranch        *b_vx_px;   //!
  TBranch        *b_vx_py;   //!
  TBranch        *b_vx_pz;   //!
  TBranch        *b_vx_E;   //!
  TBranch        *b_vx_m;   //!
  TBranch        *b_vx_nTracks;   //!
  TBranch        *b_vx_sumPt;   //!
  TBranch        *b_top_hfor_type;   //!
  TBranch        *b_jet_AntiKt4TruthJets_n;   //!
  TBranch        *b_jet_AntiKt4TruthJets_pt;   //!
  TBranch        *b_jet_AntiKt4TruthJets_m;   //!
  TBranch        *b_jet_AntiKt4TruthJets_eta;   //!
  TBranch        *b_jet_AntiKt4TruthJets_phi;   //!
  TBranch        *b_jet_AntiKt4TruthJets_flavor_truth_label;   //!
  TBranch        *b_trig_EF_el_n;   //!
  TBranch        *b_trig_EF_el_pt;   //!
  TBranch        *b_trig_EF_el_eta;   //!
  TBranch        *b_trig_EF_el_phi;   //!
  TBranch        *b_trig_L1_TAV;   //!
  TBranch        *b_trig_EF_trigmuonef_n;   //!
  TBranch        *b_trig_EF_trigmuonef_track_n;   //!
  TBranch        *b_trig_EF_trigmuonef_track_CB_pt;   //!
  TBranch        *b_trig_EF_trigmuonef_track_CB_eta;   //!
  TBranch        *b_trig_EF_trigmuonef_track_CB_phi;   //!
  TBranch        *b_trig_EF_trigmuonef_track_CB_hasCB;   //!
  TBranch        *b_jet_AntiKt4TrackZ_n;   //!
  TBranch        *b_jet_AntiKt4TrackZ_pt;   //!
  TBranch        *b_jet_AntiKt4TrackZ_eta;   //!
  TBranch        *b_jet_AntiKt4TrackZ_phi;   //!
  TBranch        *b_mc_n;   //!
  TBranch        *b_mc_pt;   //!
  TBranch        *b_mc_m;   //!
  TBranch        *b_mc_eta;   //!
  TBranch        *b_mc_phi;   //!
  TBranch        *b_mc_status;   //!
  TBranch        *b_mc_barcode;   //!
  TBranch        *b_mc_oldindex;   //!
  TBranch        *b_mc_pdgId;   //!
  TBranch        *b_mc_charge;   //!
  TBranch        *b_mc_vx_x;   //!
  TBranch        *b_mc_vx_y;   //!
  TBranch        *b_mc_vx_z;   //!
  TBranch        *b_mc_children;   //!
  TBranch        *b_mc_parents;   //!
  TBranch        *b_mc_child_index;   //!
  TBranch        *b_mc_parent_index;   //!
  TBranch        *b_mcevt_n;   //!
  TBranch        *b_mcevt_signal_process_id;   //!
  TBranch        *b_mcevt_event_number;   //!
  TBranch        *b_mcevt_event_scale;   //!
  TBranch        *b_mcevt_alphaQCD;   //!
  TBranch        *b_mcevt_alphaQED;   //!
  TBranch        *b_mcevt_pdf_id1;   //!
  TBranch        *b_mcevt_pdf_id2;   //!
  TBranch        *b_mcevt_pdf_x1;   //!
  TBranch        *b_mcevt_pdf_x2;   //!
  TBranch        *b_mcevt_pdf_scale;   //!
  TBranch        *b_mcevt_pdf1;   //!
  TBranch        *b_mcevt_pdf2;   //!
  TBranch        *b_mcevt_weight;   //!
  TBranch        *b_idp1;   //!
  TBranch        *b_idp2;   //!
  
     

  read_d4pd(TTree * /*tree*/ =0) { }
  virtual ~read_d4pd() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();
  
  ClassDef(read_d4pd,0);
};

#endif

#ifdef read_d4pd_cxx
void read_d4pd::Init(TTree *tree)
{
  std::cout<<"Init   "<<std::endl;
  
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set object pointer
 
  trig_EF_el_EF_2e12Tvh_loose1 = 0;
  trig_EF_el_EF_e12Tvh_loose1 = 0;
  trig_EF_el_EF_e12Tvh_medium1 = 0;
  trig_EF_el_EF_e12Tvh_medium1_mu8 = 0;
   trig_EF_el_EF_e24vh_medium1 = 0;
   trig_EF_el_EF_e24vh_medium1_e7_medium1 = 0;
   trig_EF_el_EF_e24vhi_medium1 = 0;
   trig_EF_el_EF_e60_medium1 = 0;
   trig_EF_el_EF_e7T_medium1 = 0;
   trig_EF_trigmuonef_EF_2mu13 = 0;
   trig_EF_trigmuonef_EF_mu13 = 0;
   trig_EF_trigmuonef_EF_mu18_tight = 0;
   trig_EF_trigmuonef_EF_mu18_tight_e7_medium1 = 0;
   trig_EF_trigmuonef_EF_mu18_tight_mu8_EFFS = 0;
   trig_EF_trigmuonef_EF_mu20it_tight = 0;
   trig_EF_trigmuonef_EF_mu24_tight = 0;
   trig_EF_trigmuonef_EF_mu24_tight_mu6_EFFS = 0;
   trig_EF_trigmuonef_EF_mu24i_tight = 0;
   trig_EF_trigmuonef_EF_mu36_tight = 0;
   trig_EF_trigmuonef_EF_mu8 = 0;
   el_E = 0;
   el_Et = 0;
   el_pt = 0;
   el_eta = 0;
   el_phi = 0;
   el_charge = 0;
   el_author = 0;
   el_isEM = 0;
   el_OQ = 0;
   el_nConv = 0;
   el_origin = 0;
   el_truth_E = 0;
   el_truth_pt = 0;
   el_truth_eta = 0;
   el_truth_phi = 0;
   el_truth_type = 0;
   el_truth_status = 0;
   el_truth_barcode = 0;
   el_truth_mothertype = 0;
   el_truth_motherbarcode = 0;
   el_truth_hasHardBrem = 0;
   el_truth_matched = 0;
   el_loose = 0;
   el_medium = 0;
   el_mediumIso = 0;
   el_tight = 0;
   el_tightIso = 0;
   el_mediumPP = 0;
   el_tightPP = 0;
   el_Ethad = 0;
   el_Ethad1 = 0;
   el_f1 = 0;
   el_Emax2 = 0;
   el_wstot = 0;
   el_emaxs1 = 0;
   el_E237 = 0;
   el_E277 = 0;
   el_weta2 = 0;
   el_f3 = 0;
   el_Etcone20 = 0;
   el_Etcone30 = 0;
   el_Etcone40 = 0;
   el_ptcone20 = 0;
   el_ptcone30 = 0;
   el_deltaeta1 = 0;
   el_deltaphi2 = 0;
   el_expectHitInBLayer = 0;
   el_trackd0_physics = 0;
   el_reta = 0;
   el_etas2 = 0;
   el_cl_E = 0;
   el_cl_pt = 0;
   el_cl_eta = 0;
   el_cl_phi = 0;
   el_trackd0 = 0;
   el_trackz0 = 0;
   el_trackphi = 0;
   el_trackqoverp = 0;
   el_trackpt = 0;
   el_tracketa = 0;
   el_nBLHits = 0;
   el_nPixHits = 0;
   el_nSCTHits = 0;
   el_nTRTHits = 0;
   el_nBLayerOutliers = 0;
   el_nPixelOutliers = 0;
   el_nSCTOutliers = 0;
   el_nTRTOutliers = 0;
   el_nSiHits = 0;
   el_TRTHighTOutliersRatio = 0;
   el_trackd0pv = 0;
   el_trackz0pv = 0;
   el_trackIPEstimate_d0_unbiasedpvunbiased = 0;
   el_trackIPEstimate_z0_unbiasedpvunbiased = 0;
   el_trackIPEstimate_sigd0_unbiasedpvunbiased = 0;
   el_hastrack = 0;
   el_topoEtcone30_corrected = 0;
   mu_staco_E = 0;
   mu_staco_pt = 0;
   mu_staco_eta = 0;
   mu_staco_phi = 0;
   mu_staco_px = 0;
   mu_staco_py = 0;
   mu_staco_pz = 0;
   mu_staco_charge = 0;
   mu_staco_author = 0;
   mu_staco_matchchi2 = 0;
   mu_staco_etcone20 = 0;
   mu_staco_etcone30 = 0;
   mu_staco_etcone40 = 0;
   mu_staco_ptcone20 = 0;
   mu_staco_ptcone30 = 0;
   mu_staco_ptcone30_trkelstyle = 0;
   mu_staco_ptcone40 = 0;
   mu_staco_bestMatch = 0;
   mu_staco_isCombinedMuon = 0;
   mu_staco_isLowPtReconstructedMuon = 0;
   mu_staco_isSegmentTaggedMuon = 0;
   mu_staco_loose = 0;
   mu_staco_medium = 0;
   mu_staco_tight = 0;
   mu_staco_d0_exPV = 0;
   mu_staco_z0_exPV = 0;
   mu_staco_id_phi_exPV = 0;
   mu_staco_id_theta_exPV = 0;
   mu_staco_id_qoverp_exPV = 0;
   mu_staco_me_phi_exPV = 0;
   mu_staco_me_theta_exPV = 0;
   mu_staco_me_qoverp_exPV = 0;
   mu_staco_ms_phi = 0;
   mu_staco_ms_theta = 0;
   mu_staco_ms_qoverp = 0;
   mu_staco_id_phi = 0;
   mu_staco_id_theta = 0;
   mu_staco_id_qoverp = 0;
   mu_staco_me_phi = 0;
   mu_staco_me_theta = 0;
   mu_staco_me_qoverp = 0;
   mu_staco_ie_phi = 0;
   mu_staco_ie_theta = 0;
   mu_staco_ie_qoverp = 0;
   mu_staco_nBLHits = 0;
   mu_staco_nPixHits = 0;
   mu_staco_nSCTHits = 0;
   mu_staco_nTRTHits = 0;
   mu_staco_nPixHoles = 0;
   mu_staco_nSCTHoles = 0;
   mu_staco_nTRTOutliers = 0;
   mu_staco_nPixelDeadSensors = 0;
   mu_staco_nSCTDeadSensors = 0;
   mu_staco_expectBLayerHit = 0;
   mu_staco_nMDTHits = 0;
   mu_staco_nCSCEtaHits = 0;
   mu_staco_nCSCPhiHits = 0;
   mu_staco_nRPCEtaHits = 0;
   mu_staco_nRPCPhiHits = 0;
   mu_staco_nTGCEtaHits = 0;
   mu_staco_nTGCPhiHits = 0;
   mu_staco_trackd0 = 0;
   mu_staco_trackz0 = 0;
   mu_staco_trackphi = 0;
   mu_staco_tracktheta = 0;
   mu_staco_trackqoverp = 0;
   mu_staco_trackcov_d0 = 0;
   mu_staco_trackcov_z0 = 0;
   mu_staco_trackIPEstimate_d0_unbiasedpvunbiased = 0;
   mu_staco_trackIPEstimate_z0_unbiasedpvunbiased = 0;
   mu_staco_trackIPEstimate_sigd0_unbiasedpvunbiased = 0;
   mu_staco_truth_dr = 0;
   mu_staco_truth_E = 0;
   mu_staco_truth_pt = 0;
   mu_staco_truth_eta = 0;
   mu_staco_truth_phi = 0;
   mu_staco_truth_type = 0;
   mu_staco_truth_status = 0;
   mu_staco_truth_barcode = 0;
   mu_staco_truth_mothertype = 0;
   mu_staco_truth_motherbarcode = 0;
   mu_staco_truth_matched = 0;
   el_MET_Egamma10NoTau_wpx = 0;
   el_MET_Egamma10NoTau_wpy = 0;
   el_MET_Egamma10NoTau_wet = 0;
   el_MET_Egamma10NoTau_statusWord = 0;
   jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpx = 0;
   jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpy = 0;
   jet_AntiKt4LCTopo_MET_Egamma10NoTau_wet = 0;
   jet_AntiKt4LCTopo_MET_Egamma10NoTau_statusWord = 0;
   jet_AntiKt4LCTopo_E = 0;
   jet_AntiKt4LCTopo_pt = 0;
   jet_AntiKt4LCTopo_m = 0;
   jet_AntiKt4LCTopo_eta = 0;
   jet_AntiKt4LCTopo_phi = 0;
   jet_AntiKt4LCTopo_EtaOrigin = 0;
   jet_AntiKt4LCTopo_PhiOrigin = 0;
   jet_AntiKt4LCTopo_MOrigin = 0;
   jet_AntiKt4LCTopo_WIDTH = 0;
   jet_AntiKt4LCTopo_n90 = 0;
   jet_AntiKt4LCTopo_Timing = 0;
   jet_AntiKt4LCTopo_LArQuality = 0;
   //   jet_AntiKt4LCTopo_nTrk = 0;
   jet_AntiKt4LCTopo_sumPtTrk = 0;
   //jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV = 0;
   jet_AntiKt4LCTopo_HECQuality = 0;
   jet_AntiKt4LCTopo_NegativeE = 0;
   jet_AntiKt4LCTopo_AverageLArQF = 0;
   jet_AntiKt4LCTopo_BCH_CORR_CELL = 0;
   jet_AntiKt4LCTopo_BCH_CORR_DOTX = 0;
   jet_AntiKt4LCTopo_BCH_CORR_JET = 0;
   jet_AntiKt4LCTopo_SamplingMax = 0;
   jet_AntiKt4LCTopo_fracSamplingMax = 0;
   jet_AntiKt4LCTopo_hecf = 0;
   jet_AntiKt4LCTopo_tgap3f = 0;
   jet_AntiKt4LCTopo_emfrac = 0;
   jet_AntiKt4LCTopo_EMJES = 0;
   jet_AntiKt4LCTopo_emscale_E = 0;
   jet_AntiKt4LCTopo_emscale_pt = 0;
   jet_AntiKt4LCTopo_emscale_m = 0;
   jet_AntiKt4LCTopo_emscale_eta = 0;
   jet_AntiKt4LCTopo_emscale_phi = 0;
   jet_AntiKt4LCTopo_ActiveAreaPx = 0;
   jet_AntiKt4LCTopo_ActiveAreaPy = 0;
   jet_AntiKt4LCTopo_ActiveAreaPz = 0;
   jet_AntiKt4LCTopo_ActiveAreaE = 0;
   jet_AntiKt4LCTopo_jvtxf = 0;
   jet_AntiKt4LCTopo_constscale_E = 0;
   jet_AntiKt4LCTopo_constscale_m = 0;
   jet_AntiKt4LCTopo_constscale_eta = 0;
   jet_AntiKt4LCTopo_constscale_phi = 0;
   jet_AntiKt4LCTopo_flavor_weight_Comb = 0;
   jet_AntiKt4LCTopo_flavor_weight_IP3D = 0;
   jet_AntiKt4LCTopo_flavor_weight_SV0 = 0;
   jet_AntiKt4LCTopo_flavor_weight_SV1 = 0;
   jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN = 0;
   jet_AntiKt4LCTopo_flavor_weight_MV1 = 0;
   jet_AntiKt4LCTopo_flavor_truth_label = 0;
   vx_x = 0;
   vx_y = 0;
   vx_z = 0;
   vx_px = 0;
   vx_py = 0;
   vx_pz = 0;
   vx_E = 0;
   vx_m = 0;
   vx_nTracks = 0;
   vx_sumPt = 0;
   jet_AntiKt4TruthJets_pt = 0;
   jet_AntiKt4TruthJets_m = 0;
   jet_AntiKt4TruthJets_eta = 0;
   jet_AntiKt4TruthJets_phi = 0;
   jet_AntiKt4TruthJets_flavor_truth_label = 0;
   trig_EF_el_pt = 0;
   trig_EF_el_eta = 0;
   trig_EF_el_phi = 0;
   trig_L1_TAV = 0;
   trig_EF_trigmuonef_track_n = 0;
   trig_EF_trigmuonef_track_CB_pt = 0;
   trig_EF_trigmuonef_track_CB_eta = 0;
   trig_EF_trigmuonef_track_CB_phi = 0;
   trig_EF_trigmuonef_track_CB_hasCB = 0;
   jet_AntiKt4TrackZ_pt = 0;
   jet_AntiKt4TrackZ_eta = 0;
   jet_AntiKt4TrackZ_phi = 0;
   mc_pt = 0;
   mc_m = 0;
   mc_eta = 0;
   mc_phi = 0;
   mc_status = 0;
   mc_barcode = 0;
   mc_oldindex = 0;
   mc_pdgId = 0;
   mc_charge = 0;
   mc_vx_x = 0;
   mc_vx_y = 0;
   mc_vx_z = 0;
   mc_children = 0;
   mc_parents = 0;
   mc_child_index = 0;
   mc_parent_index = 0;
   mcevt_signal_process_id = 0;
   mcevt_event_number = 0;
   mcevt_event_scale = 0;
   mcevt_alphaQCD = 0;
   mcevt_alphaQED = 0;
   mcevt_pdf_id1 = 0;
   mcevt_pdf_id2 = 0;
   mcevt_pdf_x1 = 0;
   mcevt_pdf_x2 = 0;
   mcevt_pdf_scale = 0;
   mcevt_pdf1 = 0;
   mcevt_pdf2 = 0;
   mcevt_weight = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);
   
   fChain->SetBranchAddress("EF_2e12Tvh_loose1", &EF_2e12Tvh_loose1, &b_EF_2e12Tvh_loose1);
   fChain->SetBranchAddress("EF_2mu13", &EF_2mu13, &b_EF_2mu13);
   fChain->SetBranchAddress("EF_e12Tvh_medium1_mu8", &EF_e12Tvh_medium1_mu8, &b_EF_e12Tvh_medium1_mu8);
   fChain->SetBranchAddress("EF_e24vh_medium1", &EF_e24vh_medium1, &b_EF_e24vh_medium1);
   fChain->SetBranchAddress("EF_e24vh_medium1_e7_medium1", &EF_e24vh_medium1_e7_medium1, &b_EF_e24vh_medium1_e7_medium1);
   fChain->SetBranchAddress("EF_e24vhi_medium1", &EF_e24vhi_medium1, &b_EF_e24vhi_medium1);
   fChain->SetBranchAddress("EF_e60_medium1", &EF_e60_medium1, &b_EF_e60_medium1);
   fChain->SetBranchAddress("EF_mu18_tight", &EF_mu18_tight, &b_EF_mu18_tight);
   fChain->SetBranchAddress("EF_mu18_tight_e7_medium1", &EF_mu18_tight_e7_medium1, &b_EF_mu18_tight_e7_medium1);
   fChain->SetBranchAddress("EF_mu18_tight_mu8_EFFS", &EF_mu18_tight_mu8_EFFS, &b_EF_mu18_tight_mu8_EFFS);
   fChain->SetBranchAddress("EF_mu20it_tight", &EF_mu20it_tight, &b_EF_mu20it_tight);
   fChain->SetBranchAddress("EF_mu24_tight", &EF_mu24_tight, &b_EF_mu24_tight);
   fChain->SetBranchAddress("EF_mu24_tight_mu6_EFFS", &EF_mu24_tight_mu6_EFFS, &b_EF_mu24_tight_mu6_EFFS);
   fChain->SetBranchAddress("EF_mu24i_tight", &EF_mu24i_tight, &b_EF_mu24i_tight);
   fChain->SetBranchAddress("EF_mu36_tight", &EF_mu36_tight, &b_EF_mu36_tight);
   fChain->SetBranchAddress("EF_xe80T_tclcw_loose", &EF_xe80T_tclcw_loose, &b_EF_xe80T_tclcw_loose);
   fChain->SetBranchAddress("EF_xe80_tclcw_loose", &EF_xe80_tclcw_loose, &b_EF_xe80_tclcw_loose);
   fChain->SetBranchAddress("trig_EF_el_EF_2e12Tvh_loose1", &trig_EF_el_EF_2e12Tvh_loose1, &b_trig_EF_el_EF_2e12Tvh_loose1);
   //fChain->SetBranchAddress("trig_EF_el_EF_e12Tvh_loose1", &trig_EF_el_EF_e12Tvh_loose1, &b_trig_EF_el_EF_e12Tvh_loose1);
   //fChain->SetBranchAddress("trig_EF_el_EF_e12Tvh_medium1", &trig_EF_el_EF_e12Tvh_medium1, &b_trig_EF_el_EF_e12Tvh_medium1);
   fChain->SetBranchAddress("trig_EF_el_EF_e12Tvh_medium1_mu8", &trig_EF_el_EF_e12Tvh_medium1_mu8, &b_trig_EF_el_EF_e12Tvh_medium1_mu8);
   fChain->SetBranchAddress("trig_EF_el_EF_e24vh_medium1", &trig_EF_el_EF_e24vh_medium1, &b_trig_EF_el_EF_e24vh_medium1);
   fChain->SetBranchAddress("trig_EF_el_EF_e24vh_medium1_e7_medium1", &trig_EF_el_EF_e24vh_medium1_e7_medium1, &b_trig_EF_el_EF_e24vh_medium1_e7_medium1);
   fChain->SetBranchAddress("trig_EF_el_EF_e24vhi_medium1", &trig_EF_el_EF_e24vhi_medium1, &b_trig_EF_el_EF_e24vhi_medium1);
   fChain->SetBranchAddress("trig_EF_el_EF_e60_medium1", &trig_EF_el_EF_e60_medium1, &b_trig_EF_el_EF_e60_medium1);
   //fChain->SetBranchAddress("trig_EF_el_EF_e7T_medium1", &trig_EF_el_EF_e7T_medium1, &b_trig_EF_el_EF_e7T_medium1);
   fChain->SetBranchAddress("trig_EF_trigmuonef_EF_2mu13", &trig_EF_trigmuonef_EF_2mu13, &b_trig_EF_trigmuonef_EF_2mu13);
  // fChain->SetBranchAddress("trig_EF_trigmuonef_EF_mu13", &trig_EF_trigmuonef_EF_mu13, &b_trig_EF_trigmuonef_EF_mu13);
   fChain->SetBranchAddress("trig_EF_trigmuonef_EF_mu18_tight", &trig_EF_trigmuonef_EF_mu18_tight, &b_trig_EF_trigmuonef_EF_mu18_tight);
   fChain->SetBranchAddress("trig_EF_trigmuonef_EF_mu18_tight_e7_medium1", &trig_EF_trigmuonef_EF_mu18_tight_e7_medium1, &b_trig_EF_trigmuonef_EF_mu18_tight_e7_medium1);
   fChain->SetBranchAddress("trig_EF_trigmuonef_EF_mu18_tight_mu8_EFFS", &trig_EF_trigmuonef_EF_mu18_tight_mu8_EFFS, &b_trig_EF_trigmuonef_EF_mu18_tight_mu8_EFFS);
   fChain->SetBranchAddress("trig_EF_trigmuonef_EF_mu20it_tight", &trig_EF_trigmuonef_EF_mu20it_tight, &b_trig_EF_trigmuonef_EF_mu20it_tight);
   //fChain->SetBranchAddress("trig_EF_trigmuonef_EF_mu24_tight", &trig_EF_trigmuonef_EF_mu24_tight, &b_trig_EF_trigmuonef_EF_mu24_tight);
   //fChain->SetBranchAddress("trig_EF_trigmuonef_EF_mu24_tight_mu6_EFFS", &trig_EF_trigmuonef_EF_mu24_tight_mu6_EFFS, &b_trig_EF_trigmuonef_EF_mu24_tight_mu6_EFFS);
   fChain->SetBranchAddress("trig_EF_trigmuonef_EF_mu24i_tight", &trig_EF_trigmuonef_EF_mu24i_tight, &b_trig_EF_trigmuonef_EF_mu24i_tight);
   fChain->SetBranchAddress("trig_EF_trigmuonef_EF_mu36_tight", &trig_EF_trigmuonef_EF_mu36_tight, &b_trig_EF_trigmuonef_EF_mu36_tight);
   //fChain->SetBranchAddress("trig_EF_trigmuonef_EF_mu8", &trig_EF_trigmuonef_EF_mu8, &b_trig_EF_trigmuonef_EF_mu8);
fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   //fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   //fChain->SetBranchAddress("timestamp_ns", &timestamp_ns, &b_timestamp_ns);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("actualIntPerXing", &actualIntPerXing, &b_actualIntPerXing);
   fChain->SetBranchAddress("averageIntPerXing", &averageIntPerXing, &b_averageIntPerXing);
   fChain->SetBranchAddress("mc_channel_number", &mc_channel_number, &b_mc_channel_number);
   fChain->SetBranchAddress("coreFlags", &coreFlags, &b_coreFlags);
   fChain->SetBranchAddress("pixelError", &pixelError, &b_pixelError);
   fChain->SetBranchAddress("sctError", &sctError, &b_sctError);
   fChain->SetBranchAddress("trtError", &trtError, &b_trtError);
   fChain->SetBranchAddress("larError", &larError, &b_larError);
   fChain->SetBranchAddress("tileError", &tileError, &b_tileError);
   fChain->SetBranchAddress("muonError", &muonError, &b_muonError);
   fChain->SetBranchAddress("fwdError", &fwdError, &b_fwdError);
   fChain->SetBranchAddress("coreError", &coreError, &b_coreError);
   fChain->SetBranchAddress("isSimulation", &isSimulation, &b_isSimulation);
   fChain->SetBranchAddress("el_n", &el_n, &b_el_n);
   fChain->SetBranchAddress("el_E", &el_E, &b_el_E);
   //fChain->SetBranchAddress("el_Et", &el_Et, &b_el_Et);
   fChain->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
   fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
   fChain->SetBranchAddress("el_author", &el_author, &b_el_author);
   fChain->SetBranchAddress("el_isEM", &el_isEM, &b_el_isEM);
   fChain->SetBranchAddress("el_OQ", &el_OQ, &b_el_OQ);
   //fChain->SetBranchAddress("el_nConv", &el_nConv, &b_el_nConv);
   //fChain->SetBranchAddress("el_origin", &el_origin, &b_el_origin);
   //fChain->SetBranchAddress("el_truth_E", &el_truth_E, &b_el_truth_E);
   // fChain->SetBranchAddress("el_truth_pt", &el_truth_pt, &b_el_truth_pt);
   //fChain->SetBranchAddress("el_truth_eta", &el_truth_eta, &b_el_truth_eta);
   // fChain->SetBranchAddress("el_truth_phi", &el_truth_phi, &b_el_truth_phi);
  // fChain->SetBranchAddress("el_truth_type", &el_truth_type, &b_el_truth_type);
  // fChain->SetBranchAddress("el_truth_status", &el_truth_status, &b_el_truth_status);
  // fChain->SetBranchAddress("el_truth_barcode", &el_truth_barcode, &b_el_truth_barcode);
  // fChain->SetBranchAddress("el_truth_mothertype", &el_truth_mothertype, &b_el_truth_mothertype);
 //  fChain->SetBranchAddress("el_truth_motherbarcode", &el_truth_motherbarcode, &b_el_truth_motherbarcode);
  // fChain->SetBranchAddress("el_truth_hasHardBrem", &el_truth_hasHardBrem, &b_el_truth_hasHardBrem);
  // fChain->SetBranchAddress("el_truth_matched", &el_truth_matched, &b_el_truth_matched);
   fChain->SetBranchAddress("el_loose", &el_loose, &b_el_loose);
   fChain->SetBranchAddress("el_medium", &el_medium, &b_el_medium);
  // fChain->SetBranchAddress("el_mediumIso", &el_mediumIso, &b_el_mediumIso);
   fChain->SetBranchAddress("el_tight", &el_tight, &b_el_tight);
  // fChain->SetBranchAddress("el_tightIso", &el_tightIso, &b_el_tightIso);
   fChain->SetBranchAddress("el_mediumPP", &el_mediumPP, &b_el_mediumPP);
   fChain->SetBranchAddress("el_tightPP", &el_tightPP, &b_el_tightPP);
   fChain->SetBranchAddress("el_Ethad", &el_Ethad, &b_el_Ethad);
   fChain->SetBranchAddress("el_Ethad1", &el_Ethad1, &b_el_Ethad1);
   fChain->SetBranchAddress("el_f1", &el_f1, &b_el_f1);
   fChain->SetBranchAddress("el_Emax2", &el_Emax2, &b_el_Emax2);
   fChain->SetBranchAddress("el_wstot", &el_wstot, &b_el_wstot);
   fChain->SetBranchAddress("el_emaxs1", &el_emaxs1, &b_el_emaxs1);
  // fChain->SetBranchAddress("el_E237", &el_E237, &b_el_E237);
  // fChain->SetBranchAddress("el_E277", &el_E277, &b_el_E277);
   fChain->SetBranchAddress("el_weta2", &el_weta2, &b_el_weta2);
   fChain->SetBranchAddress("el_f3", &el_f3, &b_el_f3);
   fChain->SetBranchAddress("el_Etcone20", &el_Etcone20, &b_el_Etcone20);
   fChain->SetBranchAddress("el_Etcone30", &el_Etcone30, &b_el_Etcone30);
   fChain->SetBranchAddress("el_Etcone40", &el_Etcone40, &b_el_Etcone40);
   fChain->SetBranchAddress("el_ptcone20", &el_ptcone20, &b_el_ptcone20);
   fChain->SetBranchAddress("el_ptcone30", &el_ptcone30, &b_el_ptcone30);
   fChain->SetBranchAddress("el_deltaeta1", &el_deltaeta1, &b_el_deltaeta1);
   fChain->SetBranchAddress("el_deltaphi2", &el_deltaphi2, &b_el_deltaphi2);
   fChain->SetBranchAddress("el_expectHitInBLayer", &el_expectHitInBLayer, &b_el_expectHitInBLayer);
   fChain->SetBranchAddress("el_trackd0_physics", &el_trackd0_physics, &b_el_trackd0_physics);
   fChain->SetBranchAddress("el_reta", &el_reta, &b_el_reta);
   fChain->SetBranchAddress("el_etas2", &el_etas2, &b_el_etas2);
   fChain->SetBranchAddress("el_cl_E", &el_cl_E, &b_el_cl_E);
   fChain->SetBranchAddress("el_cl_pt", &el_cl_pt, &b_el_cl_pt);
   fChain->SetBranchAddress("el_cl_eta", &el_cl_eta, &b_el_cl_eta);
   fChain->SetBranchAddress("el_cl_phi", &el_cl_phi, &b_el_cl_phi);
   fChain->SetBranchAddress("el_trackd0", &el_trackd0, &b_el_trackd0);
   fChain->SetBranchAddress("el_trackz0", &el_trackz0, &b_el_trackz0);
   fChain->SetBranchAddress("el_trackphi", &el_trackphi, &b_el_trackphi);
   fChain->SetBranchAddress("el_trackqoverp", &el_trackqoverp, &b_el_trackqoverp);
   fChain->SetBranchAddress("el_trackpt", &el_trackpt, &b_el_trackpt);
   fChain->SetBranchAddress("el_tracketa", &el_tracketa, &b_el_tracketa);
   fChain->SetBranchAddress("el_nBLHits", &el_nBLHits, &b_el_nBLHits);
   fChain->SetBranchAddress("el_nPixHits", &el_nPixHits, &b_el_nPixHits);
   fChain->SetBranchAddress("el_nSCTHits", &el_nSCTHits, &b_el_nSCTHits);
   fChain->SetBranchAddress("el_nTRTHits", &el_nTRTHits, &b_el_nTRTHits);
   fChain->SetBranchAddress("el_nBLayerOutliers", &el_nBLayerOutliers, &b_el_nBLayerOutliers);
   fChain->SetBranchAddress("el_nPixelOutliers", &el_nPixelOutliers, &b_el_nPixelOutliers);
   fChain->SetBranchAddress("el_nSCTOutliers", &el_nSCTOutliers, &b_el_nSCTOutliers);
   fChain->SetBranchAddress("el_nTRTOutliers", &el_nTRTOutliers, &b_el_nTRTOutliers);
   fChain->SetBranchAddress("el_nSiHits", &el_nSiHits, &b_el_nSiHits);
   fChain->SetBranchAddress("el_TRTHighTOutliersRatio", &el_TRTHighTOutliersRatio, &b_el_TRTHighTOutliersRatio);
   fChain->SetBranchAddress("el_trackd0pv", &el_trackd0pv, &b_el_trackd0pv);
   fChain->SetBranchAddress("el_trackz0pv", &el_trackz0pv, &b_el_trackz0pv);
   fChain->SetBranchAddress("el_trackIPEstimate_d0_unbiasedpvunbiased", &el_trackIPEstimate_d0_unbiasedpvunbiased, &b_el_trackIPEstimate_d0_unbiasedpvunbiased);
   fChain->SetBranchAddress("el_trackIPEstimate_z0_unbiasedpvunbiased", &el_trackIPEstimate_z0_unbiasedpvunbiased, &b_el_trackIPEstimate_z0_unbiasedpvunbiased);
   fChain->SetBranchAddress("el_trackIPEstimate_sigd0_unbiasedpvunbiased", &el_trackIPEstimate_sigd0_unbiasedpvunbiased, &b_el_trackIPEstimate_sigd0_unbiasedpvunbiased);
   // fChain->SetBranchAddress("el_hastrack", &el_hastrack, &b_el_hastrack);
   fChain->SetBranchAddress("el_topoEtcone30_corrected", &el_topoEtcone30_corrected, &b_el_topoEtcone30_corrected);
   fChain->SetBranchAddress("mu_staco_n", &mu_staco_n, &b_mu_staco_n);
   fChain->SetBranchAddress("mu_staco_E", &mu_staco_E, &b_mu_staco_E);
   fChain->SetBranchAddress("mu_staco_pt", &mu_staco_pt, &b_mu_staco_pt);
   fChain->SetBranchAddress("mu_staco_eta", &mu_staco_eta, &b_mu_staco_eta);
   fChain->SetBranchAddress("mu_staco_phi", &mu_staco_phi, &b_mu_staco_phi);
   //fChain->SetBranchAddress("mu_staco_px", &mu_staco_px, &b_mu_staco_px);
   // fChain->SetBranchAddress("mu_staco_py", &mu_staco_py, &b_mu_staco_py);
   // fChain->SetBranchAddress("mu_staco_pz", &mu_staco_pz, &b_mu_staco_pz);
   fChain->SetBranchAddress("mu_staco_charge", &mu_staco_charge, &b_mu_staco_charge);
   // fChain->SetBranchAddress("mu_staco_author", &mu_staco_author, &b_mu_staco_author);
   fChain->SetBranchAddress("mu_staco_matchchi2", &mu_staco_matchchi2, &b_mu_staco_matchchi2);
   fChain->SetBranchAddress("mu_staco_etcone20", &mu_staco_etcone20, &b_mu_staco_etcone20);
   fChain->SetBranchAddress("mu_staco_etcone30", &mu_staco_etcone30, &b_mu_staco_etcone30);
   fChain->SetBranchAddress("mu_staco_etcone40", &mu_staco_etcone40, &b_mu_staco_etcone40);
   fChain->SetBranchAddress("mu_staco_ptcone20", &mu_staco_ptcone20, &b_mu_staco_ptcone20);
   fChain->SetBranchAddress("mu_staco_ptcone30", &mu_staco_ptcone30, &b_mu_staco_ptcone30);
   fChain->SetBranchAddress("mu_staco_ptcone30_trkelstyle", &mu_staco_ptcone30_trkelstyle, &b_mu_staco_ptcone30_trkelstyle);
   fChain->SetBranchAddress("mu_staco_ptcone40", &mu_staco_ptcone40, &b_mu_staco_ptcone40);
   // fChain->SetBranchAddress("mu_staco_bestMatch", &mu_staco_bestMatch, &b_mu_staco_bestMatch);
   fChain->SetBranchAddress("mu_staco_isCombinedMuon", &mu_staco_isCombinedMuon, &b_mu_staco_isCombinedMuon);
   fChain->SetBranchAddress("mu_staco_isLowPtReconstructedMuon", &mu_staco_isLowPtReconstructedMuon, &b_mu_staco_isLowPtReconstructedMuon);
   fChain->SetBranchAddress("mu_staco_isSegmentTaggedMuon", &mu_staco_isSegmentTaggedMuon, &b_mu_staco_isSegmentTaggedMuon);
   fChain->SetBranchAddress("mu_staco_loose", &mu_staco_loose, &b_mu_staco_loose);
   fChain->SetBranchAddress("mu_staco_medium", &mu_staco_medium, &b_mu_staco_medium);
   fChain->SetBranchAddress("mu_staco_tight", &mu_staco_tight, &b_mu_staco_tight);
   fChain->SetBranchAddress("mu_staco_d0_exPV", &mu_staco_d0_exPV, &b_mu_staco_d0_exPV);
   fChain->SetBranchAddress("mu_staco_z0_exPV", &mu_staco_z0_exPV, &b_mu_staco_z0_exPV);
   // fChain->SetBranchAddress("mu_staco_id_phi_exPV", &mu_staco_id_phi_exPV, &b_mu_staco_id_phi_exPV);
   fChain->SetBranchAddress("mu_staco_id_theta_exPV", &mu_staco_id_theta_exPV, &b_mu_staco_id_theta_exPV);
   fChain->SetBranchAddress("mu_staco_id_qoverp_exPV", &mu_staco_id_qoverp_exPV, &b_mu_staco_id_qoverp_exPV);
   // fChain->SetBranchAddress("mu_staco_me_phi_exPV", &mu_staco_me_phi_exPV, &b_mu_staco_me_phi_exPV);
   fChain->SetBranchAddress("mu_staco_me_theta_exPV", &mu_staco_me_theta_exPV, &b_mu_staco_me_theta_exPV);
   fChain->SetBranchAddress("mu_staco_me_qoverp_exPV", &mu_staco_me_qoverp_exPV, &b_mu_staco_me_qoverp_exPV);
   fChain->SetBranchAddress("mu_staco_ms_phi", &mu_staco_ms_phi, &b_mu_staco_ms_phi);
   fChain->SetBranchAddress("mu_staco_ms_theta", &mu_staco_ms_theta, &b_mu_staco_ms_theta);
   fChain->SetBranchAddress("mu_staco_ms_qoverp", &mu_staco_ms_qoverp, &b_mu_staco_ms_qoverp);
   fChain->SetBranchAddress("mu_staco_id_phi", &mu_staco_id_phi, &b_mu_staco_id_phi);
   fChain->SetBranchAddress("mu_staco_id_theta", &mu_staco_id_theta, &b_mu_staco_id_theta);
   //fChain->SetBranchAddress("mu_staco_id_qoverp", &mu_staco_id_qoverp, &b_mu_staco_id_qoverp);
  // fChain->SetBranchAddress("mu_staco_me_phi", &mu_staco_me_phi, &b_mu_staco_me_phi);
  // fChain->SetBranchAddress("mu_staco_me_theta", &mu_staco_me_theta, &b_mu_staco_me_theta);
  // fChain->SetBranchAddress("mu_staco_me_qoverp", &mu_staco_me_qoverp, &b_mu_staco_me_qoverp);
  // fChain->SetBranchAddress("mu_staco_ie_phi", &mu_staco_ie_phi, &b_mu_staco_ie_phi);
  // fChain->SetBranchAddress("mu_staco_ie_theta", &mu_staco_ie_theta, &b_mu_staco_ie_theta);
  // fChain->SetBranchAddress("mu_staco_ie_qoverp", &mu_staco_ie_qoverp, &b_mu_staco_ie_qoverp);
   fChain->SetBranchAddress("mu_staco_nBLHits", &mu_staco_nBLHits, &b_mu_staco_nBLHits);
   fChain->SetBranchAddress("mu_staco_nPixHits", &mu_staco_nPixHits, &b_mu_staco_nPixHits);
   fChain->SetBranchAddress("mu_staco_nSCTHits", &mu_staco_nSCTHits, &b_mu_staco_nSCTHits);
   fChain->SetBranchAddress("mu_staco_nTRTHits", &mu_staco_nTRTHits, &b_mu_staco_nTRTHits);
   fChain->SetBranchAddress("mu_staco_nPixHoles", &mu_staco_nPixHoles, &b_mu_staco_nPixHoles);
   fChain->SetBranchAddress("mu_staco_nSCTHoles", &mu_staco_nSCTHoles, &b_mu_staco_nSCTHoles);
   fChain->SetBranchAddress("mu_staco_nTRTOutliers", &mu_staco_nTRTOutliers, &b_mu_staco_nTRTOutliers);
   fChain->SetBranchAddress("mu_staco_nPixelDeadSensors", &mu_staco_nPixelDeadSensors, &b_mu_staco_nPixelDeadSensors);
   fChain->SetBranchAddress("mu_staco_nSCTDeadSensors", &mu_staco_nSCTDeadSensors, &b_mu_staco_nSCTDeadSensors);
   fChain->SetBranchAddress("mu_staco_expectBLayerHit", &mu_staco_expectBLayerHit, &b_mu_staco_expectBLayerHit);
 //  fChain->SetBranchAddress("mu_staco_nMDTHits", &mu_staco_nMDTHits, &b_mu_staco_nMDTHits);
 //  fChain->SetBranchAddress("mu_staco_nCSCEtaHits", &mu_staco_nCSCEtaHits, &b_mu_staco_nCSCEtaHits);
   //fChain->SetBranchAddress("mu_staco_nCSCPhiHits", &mu_staco_nCSCPhiHits, &b_mu_staco_nCSCPhiHits);
   //fChain->SetBranchAddress("mu_staco_nRPCEtaHits", &mu_staco_nRPCEtaHits, &b_mu_staco_nRPCEtaHits);
   //fChain->SetBranchAddress("mu_staco_nRPCPhiHits", &mu_staco_nRPCPhiHits, &b_mu_staco_nRPCPhiHits);
  // fChain->SetBranchAddress("mu_staco_nTGCEtaHits", &mu_staco_nTGCEtaHits, &b_mu_staco_nTGCEtaHits);
  // fChain->SetBranchAddress("mu_staco_nTGCPhiHits", &mu_staco_nTGCPhiHits, &b_mu_staco_nTGCPhiHits);
   fChain->SetBranchAddress("mu_staco_trackd0", &mu_staco_trackd0, &b_mu_staco_trackd0);
   fChain->SetBranchAddress("mu_staco_trackz0", &mu_staco_trackz0, &b_mu_staco_trackz0);
  // fChain->SetBranchAddress("mu_staco_trackphi", &mu_staco_trackphi, &b_mu_staco_trackphi);
  // fChain->SetBranchAddress("mu_staco_tracktheta", &mu_staco_tracktheta, &b_mu_staco_tracktheta);
  // fChain->SetBranchAddress("mu_staco_trackqoverp", &mu_staco_trackqoverp, &b_mu_staco_trackqoverp);
   fChain->SetBranchAddress("mu_staco_trackcov_d0", &mu_staco_trackcov_d0, &b_mu_staco_trackcov_d0);
   fChain->SetBranchAddress("mu_staco_trackcov_z0", &mu_staco_trackcov_z0, &b_mu_staco_trackcov_z0);
   fChain->SetBranchAddress("mu_staco_trackIPEstimate_d0_unbiasedpvunbiased", &mu_staco_trackIPEstimate_d0_unbiasedpvunbiased, &b_mu_staco_trackIPEstimate_d0_unbiasedpvunbiased);
   fChain->SetBranchAddress("mu_staco_trackIPEstimate_z0_unbiasedpvunbiased", &mu_staco_trackIPEstimate_z0_unbiasedpvunbiased, &b_mu_staco_trackIPEstimate_z0_unbiasedpvunbiased);
   fChain->SetBranchAddress("mu_staco_trackIPEstimate_sigd0_unbiasedpvunbiased", &mu_staco_trackIPEstimate_sigd0_unbiasedpvunbiased, &b_mu_staco_trackIPEstimate_sigd0_unbiasedpvunbiased);
   // fChain->SetBranchAddress("mu_staco_truth_dr", &mu_staco_truth_dr, &b_mu_staco_truth_dr);
   // fChain->SetBranchAddress("mu_staco_truth_E", &mu_staco_truth_E, &b_mu_staco_truth_E);
  // fChain->SetBranchAddress("mu_staco_truth_pt", &mu_staco_truth_pt, &b_mu_staco_truth_pt);
  // fChain->SetBranchAddress("mu_staco_truth_eta", &mu_staco_truth_eta, &b_mu_staco_truth_eta);
  // fChain->SetBranchAddress("mu_staco_truth_phi", &mu_staco_truth_phi, &b_mu_staco_truth_phi);
  // fChain->SetBranchAddress("mu_staco_truth_type", &mu_staco_truth_type, &b_mu_staco_truth_type);
  // fChain->SetBranchAddress("mu_staco_truth_status", &mu_staco_truth_status, &b_mu_staco_truth_status);
   //fChain->SetBranchAddress("mu_staco_truth_barcode", &mu_staco_truth_barcode, &b_mu_staco_truth_barcode);
   //fChain->SetBranchAddress("mu_staco_truth_mothertype", &mu_staco_truth_mothertype, &b_mu_staco_truth_mothertype);
  // fChain->SetBranchAddress("mu_staco_truth_motherbarcode", &mu_staco_truth_motherbarcode, &b_mu_staco_truth_motherbarcode);
 //  fChain->SetBranchAddress("mu_staco_truth_matched", &mu_staco_truth_matched, &b_mu_staco_truth_matched);
   fChain->SetBranchAddress("MET_RefFinal_etx", &MET_RefFinal_etx, &b_MET_RefFinal_etx);
   fChain->SetBranchAddress("MET_RefFinal_ety", &MET_RefFinal_ety, &b_MET_RefFinal_ety);
  // fChain->SetBranchAddress("MET_RefFinal_phi", &MET_RefFinal_phi, &b_MET_RefFinal_phi);
  // fChain->SetBranchAddress("MET_RefFinal_et", &MET_RefFinal_et, &b_MET_RefFinal_et);
   fChain->SetBranchAddress("MET_RefFinal_sumet", &MET_RefFinal_sumet, &b_MET_RefFinal_sumet);
   fChain->SetBranchAddress("MET_SoftJets_etx", &MET_SoftJets_etx, &b_MET_SoftJets_etx);
   fChain->SetBranchAddress("MET_SoftJets_ety", &MET_SoftJets_ety, &b_MET_SoftJets_ety);
   fChain->SetBranchAddress("MET_SoftJets_sumet", &MET_SoftJets_sumet, &b_MET_SoftJets_sumet);
   fChain->SetBranchAddress("MET_CellOut_Eflow_etx", &MET_CellOut_Eflow_etx, &b_MET_CellOut_Eflow_etx);
   fChain->SetBranchAddress("MET_CellOut_Eflow_ety", &MET_CellOut_Eflow_ety, &b_MET_CellOut_Eflow_ety);
   fChain->SetBranchAddress("MET_CellOut_Eflow_sumet", &MET_CellOut_Eflow_sumet, &b_MET_CellOut_Eflow_sumet);
   //fChain->SetBranchAddress("MET_Truth_NonInt_etx", &MET_Truth_NonInt_etx, &b_MET_Truth_NonInt_etx);
   //fChain->SetBranchAddress("MET_Truth_NonInt_ety", &MET_Truth_NonInt_ety, &b_MET_Truth_NonInt_ety);
   // fChain->SetBranchAddress("MET_Truth_NonInt_sumet", &MET_Truth_NonInt_sumet, &b_MET_Truth_NonInt_sumet);
   fChain->SetBranchAddress("MET_Egamma10NoTau_RefGamma_etx", &MET_Egamma10NoTau_RefGamma_etx, &b_MET_Egamma10NoTau_RefGamma_etx);
   fChain->SetBranchAddress("MET_Egamma10NoTau_RefGamma_ety", &MET_Egamma10NoTau_RefGamma_ety, &b_MET_Egamma10NoTau_RefGamma_ety);
   fChain->SetBranchAddress("MET_Egamma10NoTau_RefGamma_sumet", &MET_Egamma10NoTau_RefGamma_sumet, &b_MET_Egamma10NoTau_RefGamma_sumet);
   fChain->SetBranchAddress("MET_Egamma10NoTau_RefFinal_sumet", &MET_Egamma10NoTau_RefFinal_sumet, &b_MET_Egamma10NoTau_RefFinal_sumet);
   fChain->SetBranchAddress("MET_Egamma10NoTau_CellOut_etx", &MET_Egamma10NoTau_CellOut_etx, &b_MET_Egamma10NoTau_CellOut_etx);
   fChain->SetBranchAddress("MET_Egamma10NoTau_CellOut_ety", &MET_Egamma10NoTau_CellOut_ety, &b_MET_Egamma10NoTau_CellOut_ety);
   fChain->SetBranchAddress("MET_Egamma10NoTau_CellOut_sumet", &MET_Egamma10NoTau_CellOut_sumet, &b_MET_Egamma10NoTau_CellOut_sumet);
   fChain->SetBranchAddress("MET_Egamma10NoTau_CellOut_Eflow_STVF_etx", &MET_Egamma10NoTau_CellOut_Eflow_STVF_etx, &b_MET_Egamma10NoTau_CellOut_Eflow_STVF_etx);
   fChain->SetBranchAddress("MET_Egamma10NoTau_CellOut_Eflow_STVF_ety", &MET_Egamma10NoTau_CellOut_Eflow_STVF_ety, &b_MET_Egamma10NoTau_CellOut_Eflow_STVF_ety);
   fChain->SetBranchAddress("MET_Egamma10NoTau_CellOut_Eflow_STVF_sumet", &MET_Egamma10NoTau_CellOut_Eflow_STVF_sumet, &b_MET_Egamma10NoTau_CellOut_Eflow_STVF_sumet);
   fChain->SetBranchAddress("MET_Egamma10NoTau_SoftJets_etx", &MET_Egamma10NoTau_SoftJets_etx, &b_MET_Egamma10NoTau_SoftJets_etx);
   fChain->SetBranchAddress("MET_Egamma10NoTau_SoftJets_ety", &MET_Egamma10NoTau_SoftJets_ety, &b_MET_Egamma10NoTau_SoftJets_ety);
   fChain->SetBranchAddress("MET_Egamma10NoTau_SoftJets_sumet", &MET_Egamma10NoTau_SoftJets_sumet, &b_MET_Egamma10NoTau_SoftJets_sumet);
   fChain->SetBranchAddress("el_MET_Egamma10NoTau_wpx", &el_MET_Egamma10NoTau_wpx, &b_el_MET_Egamma10NoTau_wpx);
   fChain->SetBranchAddress("el_MET_Egamma10NoTau_wpy", &el_MET_Egamma10NoTau_wpy, &b_el_MET_Egamma10NoTau_wpy);
   fChain->SetBranchAddress("el_MET_Egamma10NoTau_wet", &el_MET_Egamma10NoTau_wet, &b_el_MET_Egamma10NoTau_wet);
   fChain->SetBranchAddress("el_MET_Egamma10NoTau_statusWord", &el_MET_Egamma10NoTau_statusWord, &b_el_MET_Egamma10NoTau_statusWord);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpx", &jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpx, &b_jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpx);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpy", &jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpy, &b_jet_AntiKt4LCTopo_MET_Egamma10NoTau_wpy);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_MET_Egamma10NoTau_wet", &jet_AntiKt4LCTopo_MET_Egamma10NoTau_wet, &b_jet_AntiKt4LCTopo_MET_Egamma10NoTau_wet);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_MET_Egamma10NoTau_statusWord", &jet_AntiKt4LCTopo_MET_Egamma10NoTau_statusWord, &b_jet_AntiKt4LCTopo_MET_Egamma10NoTau_statusWord);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_n", &jet_AntiKt4LCTopo_n, &b_jet_AntiKt4LCTopo_n);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_E", &jet_AntiKt4LCTopo_E, &b_jet_AntiKt4LCTopo_E);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_pt", &jet_AntiKt4LCTopo_pt, &b_jet_AntiKt4LCTopo_pt);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_m", &jet_AntiKt4LCTopo_m, &b_jet_AntiKt4LCTopo_m);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_eta", &jet_AntiKt4LCTopo_eta, &b_jet_AntiKt4LCTopo_eta);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_phi", &jet_AntiKt4LCTopo_phi, &b_jet_AntiKt4LCTopo_phi);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_EtaOrigin", &jet_AntiKt4LCTopo_EtaOrigin, &b_jet_AntiKt4LCTopo_EtaOrigin);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_PhiOrigin", &jet_AntiKt4LCTopo_PhiOrigin, &b_jet_AntiKt4LCTopo_PhiOrigin);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_MOrigin", &jet_AntiKt4LCTopo_MOrigin, &b_jet_AntiKt4LCTopo_MOrigin);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_WIDTH", &jet_AntiKt4LCTopo_WIDTH, &b_jet_AntiKt4LCTopo_WIDTH);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_n90", &jet_AntiKt4LCTopo_n90, &b_jet_AntiKt4LCTopo_n90);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_Timing", &jet_AntiKt4LCTopo_Timing, &b_jet_AntiKt4LCTopo_Timing);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_LArQuality", &jet_AntiKt4LCTopo_LArQuality, &b_jet_AntiKt4LCTopo_LArQuality);
   //   fChain->SetBranchAddress("jet_AntiKt4LCTopo_nTrk", &jet_AntiKt4LCTopo_nTrk, &b_jet_AntiKt4LCTopo_nTrk);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_sumPtTrk", &jet_AntiKt4LCTopo_sumPtTrk, &b_jet_AntiKt4LCTopo_sumPtTrk);
   //fChain->SetBranchAddress("jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV", &jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV, &b_jet_AntiKt4LCTopo_sumPtTrk_pv0_500MeV);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_HECQuality", &jet_AntiKt4LCTopo_HECQuality, &b_jet_AntiKt4LCTopo_HECQuality);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_NegativeE", &jet_AntiKt4LCTopo_NegativeE, &b_jet_AntiKt4LCTopo_NegativeE);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_AverageLArQF", &jet_AntiKt4LCTopo_AverageLArQF, &b_jet_AntiKt4LCTopo_AverageLArQF);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_BCH_CORR_CELL", &jet_AntiKt4LCTopo_BCH_CORR_CELL, &b_jet_AntiKt4LCTopo_BCH_CORR_CELL);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_BCH_CORR_DOTX", &jet_AntiKt4LCTopo_BCH_CORR_DOTX, &b_jet_AntiKt4LCTopo_BCH_CORR_DOTX);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_BCH_CORR_JET", &jet_AntiKt4LCTopo_BCH_CORR_JET, &b_jet_AntiKt4LCTopo_BCH_CORR_JET);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_SamplingMax", &jet_AntiKt4LCTopo_SamplingMax, &b_jet_AntiKt4LCTopo_SamplingMax);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_fracSamplingMax", &jet_AntiKt4LCTopo_fracSamplingMax, &b_jet_AntiKt4LCTopo_fracSamplingMax);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_hecf", &jet_AntiKt4LCTopo_hecf, &b_jet_AntiKt4LCTopo_hecf);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_tgap3f", &jet_AntiKt4LCTopo_tgap3f, &b_jet_AntiKt4LCTopo_tgap3f);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emfrac", &jet_AntiKt4LCTopo_emfrac, &b_jet_AntiKt4LCTopo_emfrac);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_EMJES", &jet_AntiKt4LCTopo_EMJES, &b_jet_AntiKt4LCTopo_EMJES);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_E", &jet_AntiKt4LCTopo_emscale_E, &b_jet_AntiKt4LCTopo_emscale_E);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_pt", &jet_AntiKt4LCTopo_emscale_pt, &b_jet_AntiKt4LCTopo_emscale_pt);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_m", &jet_AntiKt4LCTopo_emscale_m, &b_jet_AntiKt4LCTopo_emscale_m);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_eta", &jet_AntiKt4LCTopo_emscale_eta, &b_jet_AntiKt4LCTopo_emscale_eta);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_phi", &jet_AntiKt4LCTopo_emscale_phi, &b_jet_AntiKt4LCTopo_emscale_phi);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ActiveAreaPx", &jet_AntiKt4LCTopo_ActiveAreaPx, &b_jet_AntiKt4LCTopo_ActiveAreaPx);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ActiveAreaPy", &jet_AntiKt4LCTopo_ActiveAreaPy, &b_jet_AntiKt4LCTopo_ActiveAreaPy);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ActiveAreaPz", &jet_AntiKt4LCTopo_ActiveAreaPz, &b_jet_AntiKt4LCTopo_ActiveAreaPz);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ActiveAreaE", &jet_AntiKt4LCTopo_ActiveAreaE, &b_jet_AntiKt4LCTopo_ActiveAreaE);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_jvtxf", &jet_AntiKt4LCTopo_jvtxf, &b_jet_AntiKt4LCTopo_jvtxf);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_constscale_E", &jet_AntiKt4LCTopo_constscale_E, &b_jet_AntiKt4LCTopo_constscale_E);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_constscale_m", &jet_AntiKt4LCTopo_constscale_m, &b_jet_AntiKt4LCTopo_constscale_m);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_constscale_eta", &jet_AntiKt4LCTopo_constscale_eta, &b_jet_AntiKt4LCTopo_constscale_eta);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_constscale_phi", &jet_AntiKt4LCTopo_constscale_phi, &b_jet_AntiKt4LCTopo_constscale_phi);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_Comb", &jet_AntiKt4LCTopo_flavor_weight_Comb, &b_jet_AntiKt4LCTopo_flavor_weight_Comb);
  // fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_IP3D", &jet_AntiKt4LCTopo_flavor_weight_IP3D, &b_jet_AntiKt4LCTopo_flavor_weight_IP3D);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_SV0", &jet_AntiKt4LCTopo_flavor_weight_SV0, &b_jet_AntiKt4LCTopo_flavor_weight_SV0);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_SV1", &jet_AntiKt4LCTopo_flavor_weight_SV1, &b_jet_AntiKt4LCTopo_flavor_weight_SV1);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN", &jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN, &b_jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_MV1", &jet_AntiKt4LCTopo_flavor_weight_MV1, &b_jet_AntiKt4LCTopo_flavor_weight_MV1);
  fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_truth_label", &jet_AntiKt4LCTopo_flavor_truth_label, &b_jet_AntiKt4LCTopo_flavor_truth_label);
   fChain->SetBranchAddress("Eventshape_rhoKt4LC", &Eventshape_rhoKt4LC, &b_Eventshape_rhoKt4LC);
   fChain->SetBranchAddress("vx_n", &vx_n, &b_vx_n);
  fChain->SetBranchAddress("vx_x", &vx_x, &b_vx_x);
   fChain->SetBranchAddress("vx_y", &vx_y, &b_vx_y);
   fChain->SetBranchAddress("vx_z", &vx_z, &b_vx_z);
  // fChain->SetBranchAddress("vx_px", &vx_px, &b_vx_px);
  // fChain->SetBranchAddress("vx_py", &vx_py, &b_vx_py);
  // fChain->SetBranchAddress("vx_pz", &vx_pz, &b_vx_pz);
  // fChain->SetBranchAddress("vx_E", &vx_E, &b_vx_E);
  // fChain->SetBranchAddress("vx_m", &vx_m, &b_vx_m);
   fChain->SetBranchAddress("vx_nTracks", &vx_nTracks, &b_vx_nTracks);
   fChain->SetBranchAddress("vx_sumPt", &vx_sumPt, &b_vx_sumPt);
   fChain->SetBranchAddress("top_hfor_type", &top_hfor_type, &b_top_hfor_type);
 //  fChain->SetBranchAddress("jet_AntiKt4TruthJets_n", &jet_AntiKt4TruthJets_n, &b_jet_AntiKt4TruthJets_n);
   // fChain->SetBranchAddress("jet_AntiKt4TruthJets_pt", &jet_AntiKt4TruthJets_pt, &b_jet_AntiKt4TruthJets_pt);
   //fChain->SetBranchAddress("jet_AntiKt4TruthJets_m", &jet_AntiKt4TruthJets_m, &b_jet_AntiKt4TruthJets_m);
   //fChain->SetBranchAddress("jet_AntiKt4TruthJets_eta", &jet_AntiKt4TruthJets_eta, &b_jet_AntiKt4TruthJets_eta);
   //fChain->SetBranchAddress("jet_AntiKt4TruthJets_phi", &jet_AntiKt4TruthJets_phi, &b_jet_AntiKt4TruthJets_phi);
   //fChain->SetBranchAddress("jet_AntiKt4TruthJets_flavor_truth_label", &jet_AntiKt4TruthJets_flavor_truth_label, &b_jet_AntiKt4TruthJets_flavor_truth_label);
  // fChain->SetBranchAddress("trig_EF_el_n", &trig_EF_el_n, &b_trig_EF_el_n);
   fChain->SetBranchAddress("trig_EF_el_pt", &trig_EF_el_pt, &b_trig_EF_el_pt);
   fChain->SetBranchAddress("trig_EF_el_eta", &trig_EF_el_eta, &b_trig_EF_el_eta);
   fChain->SetBranchAddress("trig_EF_el_phi", &trig_EF_el_phi, &b_trig_EF_el_phi);
  // fChain->SetBranchAddress("trig_L1_TAV", &trig_L1_TAV, &b_trig_L1_TAV);
   fChain->SetBranchAddress("trig_EF_trigmuonef_n", &trig_EF_trigmuonef_n, &b_trig_EF_trigmuonef_n);
   fChain->SetBranchAddress("trig_EF_trigmuonef_track_n", &trig_EF_trigmuonef_track_n, &b_trig_EF_trigmuonef_track_n);
   fChain->SetBranchAddress("trig_EF_trigmuonef_track_CB_pt", &trig_EF_trigmuonef_track_CB_pt, &b_trig_EF_trigmuonef_track_CB_pt);
   fChain->SetBranchAddress("trig_EF_trigmuonef_track_CB_eta", &trig_EF_trigmuonef_track_CB_eta, &b_trig_EF_trigmuonef_track_CB_eta);
   fChain->SetBranchAddress("trig_EF_trigmuonef_track_CB_phi", &trig_EF_trigmuonef_track_CB_phi, &b_trig_EF_trigmuonef_track_CB_phi);
   fChain->SetBranchAddress("trig_EF_trigmuonef_track_CB_hasCB", &trig_EF_trigmuonef_track_CB_hasCB, &b_trig_EF_trigmuonef_track_CB_hasCB);
   fChain->SetBranchAddress("jet_AntiKt4TrackZ_n", &jet_AntiKt4TrackZ_n, &b_jet_AntiKt4TrackZ_n);
   fChain->SetBranchAddress("jet_AntiKt4TrackZ_pt", &jet_AntiKt4TrackZ_pt, &b_jet_AntiKt4TrackZ_pt);
   fChain->SetBranchAddress("jet_AntiKt4TrackZ_eta", &jet_AntiKt4TrackZ_eta, &b_jet_AntiKt4TrackZ_eta);
   fChain->SetBranchAddress("jet_AntiKt4TrackZ_phi", &jet_AntiKt4TrackZ_phi, &b_jet_AntiKt4TrackZ_phi);
  // fChain->SetBranchAddress("mc_n", &mc_n, &b_mc_n);
  // fChain->SetBranchAddress("mc_pt", &mc_pt, &b_mc_pt);
  // fChain->SetBranchAddress("mc_m", &mc_m, &b_mc_m);
 //  fChain->SetBranchAddress("mc_eta", &mc_eta, &b_mc_eta);
  // fChain->SetBranchAddress("mc_phi", &mc_phi, &b_mc_phi);
   //fChain->SetBranchAddress("mc_status", &mc_status, &b_mc_status);
   //fChain->SetBranchAddress("mc_barcode", &mc_barcode, &b_mc_barcode);
  // fChain->SetBranchAddress("mc_oldindex", &mc_oldindex, &b_mc_oldindex);
  // fChain->SetBranchAddress("mc_pdgId", &mc_pdgId, &b_mc_pdgId);
  // fChain->SetBranchAddress("mc_charge", &mc_charge, &b_mc_charge);
  // fChain->SetBranchAddress("mc_vx_x", &mc_vx_x, &b_mc_vx_x);
  // fChain->SetBranchAddress("mc_vx_y", &mc_vx_y, &b_mc_vx_y);
  // fChain->SetBranchAddress("mc_vx_z", &mc_vx_z, &b_mc_vx_z);
  // fChain->SetBranchAddress("mc_children", &mc_children, &b_mc_children);
  // fChain->SetBranchAddress("mc_parents", &mc_parents, &b_mc_parents);
  // fChain->SetBranchAddress("mc_child_index", &mc_child_index, &b_mc_child_index);
  // fChain->SetBranchAddress("mc_parent_index", &mc_parent_index, &b_mc_parent_index);
  fChain->SetBranchAddress("mcevt_n", &mcevt_n, &b_mcevt_n);
  // fChain->SetBranchAddress("mcevt_signal_process_id", &mcevt_signal_process_id, &b_mcevt_signal_process_id);
  // fChain->SetBranchAddress("mcevt_event_number", &mcevt_event_number, &b_mcevt_event_number);
  // fChain->SetBranchAddress("mcevt_event_scale", &mcevt_event_scale, &b_mcevt_event_scale);
  // fChain->SetBranchAddress("mcevt_alphaQCD", &mcevt_alphaQCD, &b_mcevt_alphaQCD);
  // fChain->SetBranchAddress("mcevt_alphaQED", &mcevt_alphaQED, &b_mcevt_alphaQED);
   //fChain->SetBranchAddress("mcevt_pdf_id1", &mcevt_pdf_id1, &b_mcevt_pdf_id1);
  // fChain->SetBranchAddress("mcevt_pdf_id2", &mcevt_pdf_id2, &b_mcevt_pdf_id2);
  // fChain->SetBranchAddress("mcevt_pdf_x1", &mcevt_pdf_x1, &b_mcevt_pdf_x1);
  // fChain->SetBranchAddress("mcevt_pdf_x2", &mcevt_pdf_x2, &b_mcevt_pdf_x2);
  // fChain->SetBranchAddress("mcevt_pdf_scale", &mcevt_pdf_scale, &b_mcevt_pdf_scale);
  // fChain->SetBranchAddress("mcevt_pdf1", &mcevt_pdf1, &b_mcevt_pdf1);
  // fChain->SetBranchAddress("mcevt_pdf2", &mcevt_pdf2, &b_mcevt_pdf2);
   fChain->SetBranchAddress("mcevt_weight", &mcevt_weight, &b_mcevt_weight);
   // fChain->SetBranchAddress("idp1", &idp1, &b_idp1);
   //fChain->SetBranchAddress("idp2", &idp2, &b_idp2);   
   
   //   std::cout<<"Entries :" << fChain->GetEntriesFast()  << std::endl; 
   
   std::cout<<"Init terminate   "<<std::endl;
}

Bool_t read_d4pd::Notify()
{
  std::cout<<"Notify " <<std::endl;
  
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
   std::cout<<"Notify " <<std::endl;
 
   return kTRUE;
}

#endif // #ifdef read_d4pd_cxx
