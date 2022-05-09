cut_flow()
{
  //TString hname="(h_cutFlow";
  //            TFile * file_stop= TFile::Open("histograms.signal_177857.root");
        	  TFile * file_stop= TFile::Open("/Users/sdarmora/histograms.B_MET_july_one.root");
  //       	TFile * file_stop= TFile::Open("files_after_final_check_xsection/histograms.ttbar_110001.root");
  //	TFile * file_stop= TFile::Open("histograms.B_MET.root");
  //		  TFile * file_stop= TFile::Open("ntuples_for_TNVA_lep_filter_sept/histograms.signal_166440.root");
  //           	         TFile * file_stop= TFile::Open("histograms.ZJet_single.root");
  //  TFile * file_stop= TFile::Open("histograms.signal_179690.root");
  //                TFile * file_ttbar= TFile::Open("histograms.ttbar_110001.root");
  //  TFile * fileW= TFile::Open("histograms.WW.root");
  //TFile * fileZ= TFile::Open("histograms.WZ.root");
  //TFile * file_ZZ= TFile::Open("histograms.ZZ.root");
  //   TFile * file_ttbar_bos= TFile::Open("histograms.ttbar_bos.root");
  //      TFile * file_Zjet= TFile::Open("/Users/sdarmora/ST_3_20/run/ntuples_for_TMVA_july3/histograms.ZJet.root");
  //TFile * file_DrellYan= TFile::Open("histograms.DrellYan.root");
  //             TFile * file_H= TFile::Open("histograms.H_full.root");
  //    TFile * file_H= TFile::Open("histograms.H_old_susyTool.root");

  //     TFile * filet= TFile::Open("histograms.Wt_mediumIso.root");
  //    TFile * file_dibos= TFile::Open("histograms.dibos_pileup.root");

  //  TFile * file_dibos= TFile::Open("histograms.ttbar_bos.root");

  //   TH1F * h_stop  = (TH1F*)file_stop->Get(hname);
  //TH1F * h_ttbar  = (TH1F*)file_ttbar->Get(hname);
  //TH1F * hW  = (TH1F*)fileW->Get(hname);
  //TH1F * hZ  = (TH1F*)fileZ->Get(hname);
  //TH1F * h_ZZ  = (TH1F*)file_ZZ->Get(hname);
  //  TH1F * h_Zjet  = (TH1F*)file_Zjet->Get(hname);
  //TH1F * h_ttbar_bos  = (TH1F*)file_ttbar_bos->Get(hname);
  //TH1F * h_DrellYan  = (TH1F*)file_DrellYan->Get(hname);
  //  TH1F * ht  = (TH1F*)filet->Get(hname);
  //TH1F * h_H  = (TH1F*)file_H->Get(hname);
  //TH1F * h_dibos  = (TH1F*)file_dibos->Get(hname);

 
  
  cout<<"1.  All events:                                              "<<(h_cutFlow->GetBinContent(1))<<endl;
  cout<<"2.  GRL (only data):                                         "<<(h_cutFlow->GetBinContent(2))<<endl;
  cout<<"3.  coreFlag       :                                         "<<(h_cutFlow->GetBinContent(3))<<endl;
  cout<<"4.  LAr ERROR      :                                         "<<(h_cutFlow->GetBinContent(4))<<endl;
  cout<<"5.  Mu Trigger:                                              "<<(h_cutFlow->GetBinContent(5))<<endl;
  cout<<"6.  EL Trigger:                                              "<<(h_cutFlow->GetBinContent(6))<<endl;
  cout<<"7.  OR Trigger:                                              "<<(h_cutFlow->GetBinContent(7))<<endl;
  cout<<"8.  Primary vertex:                                          "<<(h_cutFlow->GetBinContent(8))<<endl;
  cout<<"9.  cosmic veto:                                             "<<(h_cutFlow->GetBinContent(9))<<endl;
  cout<<"10. good Muon:                                               "<<(h_cutFlow->GetBinContent(10))<<endl;
  cout<<"11. good electron:                                           "<<(h_cutFlow->GetBinContent(11))<<endl;
  cout<<"12. two good lepton:                                         "<<(h_cutFlow->GetBinContent(12))<<endl;
  cout<<"13.  Jet cleaning:                                           "<<(h_cutFlow->GetBinContent(13))<<endl;
  cout<<"14.  MET cut:                                                "<<(h_cutFlow->GetBinContent(14))<<endl;
  cout<<"15.  Smart Veto:                                             "<<(h_cutFlow->GetBinContent(15))<<endl;
  cout<<"16.  jet size<2:                                             "<<(h_cutFlow->GetBinContent(16))<<endl;
    //cout<<"17. Beauty tagging:                                      "<<(h_cutFlow->GetBinContent(17))<<endl;
  cout<<"18.  atleast 2 baselibe leptons:                             "<<(h_cutFlow->GetBinContent(18))<<endl;
  cout<<"19.  2 baseline muons:                                       "<<(h_cutFlow->GetBinContent(19))<<endl;
    cout<<"20.  2 tight muons :                                         "<<(h_cutFlow->GetBinContent(20))<<endl;
  cout<<"21.  2 tight muon+0S+minv>20 :                               "<<(h_cutFlow->GetBinContent(21))<<endl;
  cout<<"22.  2 tight muon+0S:                                        "<<(h_cutFlow->GetBinContent(22))<<endl;
  cout<<"23.  2 baseline electrons :                                  "<<(h_cutFlow->GetBinContent(23))<<endl;
  cout<<"24.  2 tight electron :                                      "<<(h_cutFlow->GetBinContent(24))<<endl;
  cout<<"25.  2 tight el+minv>20+elpt>25+jet>2:                       "<<(h_cutFlow->GetBinContent(25))<<endl;
  cout<<"26.  2 tight elecetron + OS:                                 "<<(h_cutFlow->GetBinContent(26))<<endl;
  cout<<"27.  2 base lepton:                                          "<<(h_cutFlow->GetBinContent(27))<<endl;
  cout<<"28.  2 tight lepton:                                         "<<(h_cutFlow->GetBinContent(28))<<endl;
  cout<<"29.  2 tight lepton + OS:                                    "<<(h_cutFlow->GetBinContent(29))<<endl;
  cout<<"30.  2 baseline lepton + OS:                                 "<<(h_cutFlow->GetBinContent(61))<<endl;
  cout<<"32.  2 tight mu+OS+muPt>20+minv>20:                          "<<(h_cutFlow->GetBinContent(32))<<endl;
  cout<<"33.  2 tight mu+OS+muPt>20+minv>20+jet>2:                    "<<(h_cutFlow->GetBinContent(33))<<endl;
  cout<<"34.  2 tight el+OS+minv>20:                                  "<<(h_cutFlow->GetBinContent(34))<<endl;
  cout<<"35.  2 tight el +OS+minv>20+el_pt>20:                        "<<(h_cutFlow->GetBinContent(35))<<endl;
  cout<<"36.  2 tight lep+OS+minv>20:                                 "<<(h_cutFlow->GetBinContent(36))<<endl;
  cout<<"37.  2 tight lep+OS+minv>20+el_pt>20:                        "<<(h_cutFlow->GetBinContent(37))<<endl;
  cout<<"38.  2 tight lep+OS+minv>20+el_pt>20+jet>2:                  "<<(h_cutFlow->GetBinContent(38))<<endl;
  cout<<"42.  2 tight muon + C1 cut(met>50):                          "<<(h_cutFlow->GetBinContent(42))<<endl;
  cout<<"43.  2 tight muon + C2 cut(met>80):                          "<<(h_cutFlow->GetBinContent(43))<<endl;
  cout<<"44.  2 tight muon + C3 cut(met>50,MuPt>50):                  "<<(h_cutFlow->GetBinContent(44))<<endl;
  cout<<"45.  2 tight muon + C4 cut(met>80,MuPt>50):                  "<<(h_cutFlow->GetBinContent(45))<<endl;
  cout<<"46.  2 tight el + C1 cut(met>50):                            "<<(h_cutFlow->GetBinContent(46))<<endl;
  cout<<"47.  2 tight el + C2 cut(met>80):                            "<<(h_cutFlow->GetBinContent(47))<<endl;
  cout<<"48.  2 tight el + C3 cut(met>50,MuPt>50):                    "<<(h_cutFlow->GetBinContent(48))<<endl;
  cout<<"49.  2 tight el + C4 cut(met>80,MuPt>50):                    "<<(h_cutFlow->GetBinContent(49))<<endl;
  cout<<"50.  2 tight lep + C1 cut(met>50):                           "<<(h_cutFlow->GetBinContent(50))<<endl;
  cout<<"51.  2 tight lep + C2 cut(met>80):                           "<<(h_cutFlow->GetBinContent(51))<<endl;
  cout<<"52.  2 tight lep + C3 cut(met>50,MuPt>50):                   "<<(h_cutFlow->GetBinContent(52))<<endl;
  cout<<"53.  2 tight lep + C4 cut(met>80,MuPt>50):                   "<<(h_cutFlow->GetBinContent(53))<<endl;
  cout<<"54.  2 base mu+opposite sign:                                "<<(h_cutFlow->GetBinContent(54))<<endl;
  cout<<"55.  2 base mu+Same sign:                                    "<<(h_cutFlow->GetBinContent(55))<<endl;
  cout<<"22.  2 tight mu+opposite sign:                               "<<(h_cutFlow->GetBinContent(22))<<endl;
  cout<<"56.  2 tight mu+Same sign:                                   "<<(h_cutFlow->GetBinContent(56))<<endl;
  cout<<"57.  2 base el+opposite sign:                                "<<(h_cutFlow->GetBinContent(57))<<endl;
  cout<<"58.  2 base el+Same sign:                                    "<<(h_cutFlow->GetBinContent(58))<<endl;
  cout<<"26.  2 tight el+opposite sign:                               "<<(h_cutFlow->GetBinContent(26))<<endl;
  cout<<"59.  2 tight el+Same sign:                                   "<<(h_cutFlow->GetBinContent(59))<<endl;
  cout<<"60.  2 base leptons+opposite sign:                           "<<(h_cutFlow->GetBinContent(60))<<endl;
  cout<<"61.  2 base leptons+Same sign:                               "<<(h_cutFlow->GetBinContent(61))<<endl;
  cout<<"29.  2 tight leptons+opposite sign:                          "<<(h_cutFlow->GetBinContent(29))<<endl;
  cout<<"62.  2 tight leptons+Same sign:                              "<<(h_cutFlow->GetBinContent(62))<<endl;





 cout<<(h_cutFlow->GetBinContent(1))<<endl;
  cout<<(h_cutFlow->GetBinContent(2))<<endl;
  cout<<(h_cutFlow->GetBinContent(3))<<endl;
  cout<<(h_cutFlow->GetBinContent(4))<<endl;
  cout<<(h_cutFlow->GetBinContent(5))<<endl;
  cout<<(h_cutFlow->GetBinContent(6))<<endl;
  cout<<(h_cutFlow->GetBinContent(7))<<endl;
  cout<<(h_cutFlow->GetBinContent(8))<<endl;
  cout<<(h_cutFlow->GetBinContent(9))<<endl;
  cout<<(h_cutFlow->GetBinContent(10))<<endl;
  cout<<(h_cutFlow->GetBinContent(11))<<endl;
  cout<<(h_cutFlow->GetBinContent(12))<<endl;
  cout<<(h_cutFlow->GetBinContent(13))<<endl;
  cout<<(h_cutFlow->GetBinContent(14))<<endl;
  cout<<(h_cutFlow->GetBinContent(15))<<endl;
  cout<<(h_cutFlow->GetBinContent(16))<<endl;
    //cout<<"17. Beauty tagging:                                      "<<(h_cutFlow->GetBinContent(17))<<endl;
  cout<<(h_cutFlow->GetBinContent(18))<<endl;
  cout<<(h_cutFlow->GetBinContent(19))<<endl;
  cout<<(h_cutFlow->GetBinContent(20))<<endl;
  cout<<(h_cutFlow->GetBinContent(21))<<endl;
  cout<<(h_cutFlow->GetBinContent(22))<<endl;
  cout<<(h_cutFlow->GetBinContent(23))<<endl;
  cout<<(h_cutFlow->GetBinContent(24))<<endl;
  cout<<(h_cutFlow->GetBinContent(25))<<endl;
  cout<<(h_cutFlow->GetBinContent(26))<<endl;
  cout<<(h_cutFlow->GetBinContent(27))<<endl;
  cout<<(h_cutFlow->GetBinContent(28))<<endl;
  cout<<(h_cutFlow->GetBinContent(29))<<endl;
  cout<<(h_cutFlow->GetBinContent(30))<<endl;
  cout<<(h_cutFlow->GetBinContent(32))<<endl;
  cout<<(h_cutFlow->GetBinContent(33))<<endl;
  cout<<(h_cutFlow->GetBinContent(34))<<endl;
  cout<<(h_cutFlow->GetBinContent(35))<<endl;
  cout<<(h_cutFlow->GetBinContent(36))<<endl;
  cout<<(h_cutFlow->GetBinContent(37))<<endl;
  cout<<(h_cutFlow->GetBinContent(38))<<endl;
  cout<<(h_cutFlow->GetBinContent(42))<<endl;
  cout<<(h_cutFlow->GetBinContent(43))<<endl;
  cout<<(h_cutFlow->GetBinContent(44))<<endl;
  cout<<(h_cutFlow->GetBinContent(45))<<endl;
  cout<<(h_cutFlow->GetBinContent(46))<<endl;
  cout<<(h_cutFlow->GetBinContent(47))<<endl;
  cout<<(h_cutFlow->GetBinContent(48))<<endl;
  cout<<(h_cutFlow->GetBinContent(49))<<endl;
  cout<<(h_cutFlow->GetBinContent(50))<<endl;
  cout<<(h_cutFlow->GetBinContent(51))<<endl;
  cout<<(h_cutFlow->GetBinContent(52))<<endl;
  cout<<(h_cutFlow->GetBinContent(53))<<endl;
  cout<<(h_cutFlow->GetBinContent(54))<<endl;
  cout<<(h_cutFlow->GetBinContent(55))<<endl;
  cout<<(h_cutFlow->GetBinContent(22))<<endl;
  cout<<(h_cutFlow->GetBinContent(56))<<endl;
  cout<<(h_cutFlow->GetBinContent(57))<<endl;
  cout<<(h_cutFlow->GetBinContent(58))<<endl;
  cout<<(h_cutFlow->GetBinContent(26))<<endl;
  cout<<(h_cutFlow->GetBinContent(59))<<endl;
  cout<<(h_cutFlow->GetBinContent(60))<<endl;
  cout<<(h_cutFlow->GetBinContent(61))<<endl;
  cout<<(h_cutFlow->GetBinContent(29))<<endl;
  cout<<(h_cutFlow->GetBinContent(62))<<endl;






    //cout<<"27.  2 baseline leptons OS:  "<<(h_cutFlow->GetBinContent(27))<<endl;
    //cout<<"16.  2 leptons OS:           "<<(h_cutFlow->GetBinContent(16))<<endl;                                                                                                                                                               
    //cout<<"17.  2 leptons SS:           "<<(h_cutFlow->GetBinContent(17)<<endl;                  

    
  /*


  cout<<"1.  All events:              "<<h_cutFlow->GetBinContent(1)<<endl;
  cout<<"2.  GRL (only data):         "<<h_cutFlow->GetBinContent(2)<<endl;
  cout<<"3.  coreFlag       :         "<<h_cutFlow->GetBinContent(3)<<endl;
  cout<<"4.  LAr ERROR      :         "<<h_cutFlow->GetBinContent(4)<<endl;
  cout<<"5.  Mu Trigger:              "<<h_cutFlow->GetBinContent(5)<<endl;
  cout<<"6.  EL Trigger:              "<<h_cutFlow->GetBinContent(6)<<endl;
  cout<<"7.  OR Trigger:              "<<h_cutFlow->GetBinContent(7)<<endl;
  cout<<"8.  Primary vertex:          "<<h_cutFlow->GetBinContent(8)<<endl;
  cout<<"9.  cosmic veto:             "<<h_cutFlow->GetBinContent(9)<<endl;
  cout<<"10. >0 Baseline Mu:          "<<h_cutFlow->GetBinContent(10)<<endl;
  cout<<"11. >0 Baseline EL:          "<<h_cutFlow->GetBinContent(11)<<endl;
  cout<<"12.  Baseline EL or Mu:      "<<h_cutFlow->GetBinContent(12)<<endl;
  cout<<"13.  Jet cleaning:           "<<h_cutFlow->GetBinContent(13)<<endl;
  cout<<"14.  MET cut:                "<<h_cutFlow->GetBinContent(14)<<endl;
  cout<<"15.  Bad jet Veto:           "<<h_cutFlow->GetBinContent(15)<<endl;
  //cout<<"16.  =>2 base leptons:            "<<h_cutFlow->GetBinContent(16)<<endl;
  //cout<<"17. Beauty tagging:          "<<h_cutFlow->GetBinContent(17)<<endl;
  cout<<"18.  2 baseline leptons:     "<<h_cutFlow->GetBinContent(18)<<endl;
  cout<<"19.  2 baseline muons:       "<<h_cutFlow->GetBinContent(19)<<endl;
  cout<<"20.  2 tight muons :        "<<h_cutFlow->GetBinContent(20)<<endl;
  cout<<"21.  2 tight muon+0S+minv :             "<<h_cutFlow->GetBinContent(21)<<endl;
  cout<<"22.  2 tight muon+0S:   "<<h_cutFlow->GetBinContent(22)<<endl;
  cout<<"23.  2 baseline electrons OS:"<<h_cutFlow->GetBinContent(23)<<endl;
  cout<<"24.  2 tight e :             "<<h_cutFlow->GetBinContent(24)<<endl;
  cout<<"25.  2 tight e +minv+elpt+jet>2:           "<<h_cutFlow->GetBinContent(25)<<endl;
  cout<<"26.  2 tight e+OS:        "<<h_cutFlow->GetBinContent(26)<<endl;
  cout<<"32.  2 tight mu_OS_muPt:        "<<h_cutFlow->GetBinContent(32)<<endl;
  cout<<"33.  2 tight mu_OS_muPt+jet>2:        "<<h_cutFlow->GetBinContent(33)<<endl;
  cout<<"34.  2 tight e +OS+minv:        "<<h_cutFlow->GetBinContent(34)<<endl;
  cout<<"35.  2 tight e +OS+minv+E_PT:        "<<h_cutFlow->GetBinContent(35)<<endl;
  cout<<"36.  2 tight e_MU +OS+minv:        "<<h_cutFlow->GetBinContent(36)<<endl;
  cout<<"37.  2 tight e_MU +OS+minv_e_mu_PT:        "<<h_cutFlow->GetBinContent(37)<<endl;
  cout<<"38.  2 tight e_MU +OS+minv_e_mu_PT+JET>2:        "<<h_cutFlow->GetBinContent(38)<<endl;



  */


 cout<<"   All events:                                              "<<endl;
  cout<<"  GRL (only data):                                         "<<endl;
  cout<<"  coreFlag       :                                         "<<endl;
  cout<<"  LAr ERROR      :                                         "<<endl;
  cout<<"  Mu Trigger:                                              "<<endl;
  cout<<"  EL Trigger:                                              "<<endl;
  cout<<"  OR Trigger:                                              "<<endl;
  cout<<"  good Muon:                                               "<<endl;
  cout<<"  good electron:                                           "<<endl;

  cout<<"  Smart Veto:                                             "<<endl;
  cout<<"  Jet cleaning:                                           "<<endl;
  cout<<"  Primary vertex:                                          "<<endl;
  cout<<"  cosmic veto:                                             "<<endl;
  cout<<"  Energy averaged time:                                    "<<endl;
  cout<<"  Bad leading Jet :                                       "<<endl;
  cout<<"  jet size<2:                                             "<<endl;

  cout<<"  atleast 2 baseline leptons:                             "<<endl;
  cout<<"  Exactly 2 baseline leptons:                             "<<endl;
  cout<<"  Exactly 2 tight leptons:                                "<<endl;


  cout<<"  2 baseline muons:                                       "<<endl;
  cout<<"  2 base mu+opposite sign:                                "<<endl;
  cout<<"  2 tight muons :                                         "<<endl;
  cout<<"  2 tight muon+0S:                                        "<<endl;   
  cout<<"  2 tight muon+0S+minv>20 :                               "<<endl;
  cout<<"  2 tight mu+OS+muPt>20+minv>20:                          "<<endl;
  cout<<"  2 tight mu+OS+muPt>20+minv>20+jet>2(Pt>20):             "<<endl;

  cout<<"  2 baseline electrons :                                  "<<endl;
  cout<<"  2 tight electron :                                      "<<endl;
  cout<<"  2 tight elecetron + OS:                                 "<<endl;
  cout<<"  2 tight el+OS+minv>20:                                  "<<endl;
  cout<<"  2 tight el +OS+minv>20+el_pt>20:                        "<<endl;
  cout<<"  2 tight el+minv>20+elpt>25+jet>2(Pt>20):                       "<<endl;

  cout<<"  2 base lepton:                                          "<<endl;
  cout<<"  2 baseline lepton + OS:                                 "<<endl; 
  cout<<"  2 tight lepton:                                         "<<endl;
  cout<<"  2 tight lepton + OS:                                    "<<endl;
  cout<<"  2 tight lep+OS+minv>20:                                 "<<endl;
  cout<<"  2 tight lep+OS+minv>20+el_pt>20:                        "<<endl;
  cout<<"  2 tight lep+OS+minv>20+el_pt>20+jet>2(Pt>20):                  "<<endl;
 
  cout<<"  2 tight muon + C1 cut(met>50):                          "<<endl;
  cout<<"  2 tight muon + C2 cut(met>80):                          "<<endl;
  cout<<"  2 tight muon + C3 cut(met>50,MuPt>50):                  "<<endl;
  cout<<"  2 tight muon + C4 cut(met>80,MuPt>50):                  "<<endl;

  cout<<"  2 tight el + C1 cut(met>50):                            "<<endl;
  cout<<"  2 tight el + C2 cut(met>80):                            "<<endl;
  cout<<"  2 tight el + C3 cut(met>50,MuPt>50):                    "<<endl;
  cout<<"  2 tight el + C4 cut(met>80,MuPt>50):                    "<<endl;

  cout<<"  2 tight lep + C1 cut(met>50):                           "<<endl;
  cout<<"  2 tight lep + C2 cut(met>80):                           "<<endl;
  cout<<"  2 tight lep + C3 cut(met>50,MuPt>50):                   "<<endl;
  cout<<"  2 tight lep + C4 cut(met>80,MuPt>50):                   "<<endl;

  cout<<"  2 base mu+Same sign:                                    "<<endl;
  cout<<"  2 base mu+opposite sign:                                "<<endl;
  cout<<"  2 tight mu+opposite sign:                               "<<endl;
  cout<<"  2 tight mu+Same sign:                                   "<<endl;
  cout<<"  2 base el+opposite sign:                                "<<endl;
  cout<<"  2 base el+Same sign:                                    "<<endl;
  cout<<"  2 tight el+opposite sign:                               "<<endl;
  cout<<"  2 tight el+Same sign:                                   "<<endl;
  cout<<"  2 base leptons+opposite sign:                           "<<endl;
  cout<<"  2 base leptons+Same sign:                               "<<endl;
  cout<<"  2 tight leptons+opposite sign:                          "<<endl;
  cout<<"  2 tight leptons+Same sign:                              "<<endl;





  cout<<"1.  All events:                                              "<<(h_cutFlow->GetBinContent(1))<<endl;
  cout<<"2.  GRL (only data):                                         "<<(h_cutFlow->GetBinContent(2))<<endl;
  cout<<"3.  coreFlag       :                                         "<<(h_cutFlow->GetBinContent(3))<<endl;
  cout<<"4.  LAr ERROR      :                                         "<<(h_cutFlow->GetBinContent(4))<<endl;
  cout<<"5.  Mu Trigger:                                              "<<(h_cutFlow->GetBinContent(5))<<endl;
  cout<<"6.  EL Trigger:                                              "<<(h_cutFlow->GetBinContent(6))<<endl;
  cout<<"7.  OR Trigger:                                              "<<(h_cutFlow->GetBinContent(7))<<endl;
  cout<<"10. good Muon:                                               "<<(h_cutFlow->GetBinContent(10))<<endl;
  cout<<"11. good electron:                                           "<<(h_cutFlow->GetBinContent(11))<<endl;

  cout<<"15.  Smart Veto:                                             "<<(h_cutFlow->GetBinContent(15))<<endl;
  cout<<"13.  Jet cleaning:                                           "<<(h_cutFlow->GetBinContent(13))<<endl;
  cout<<"8.  Primary vertex:                                          "<<(h_cutFlow->GetBinContent(8))<<endl;
  cout<<"9.  cosmic veto:                                             "<<(h_cutFlow->GetBinContent(9))<<endl;
  cout<<"39. Energy averaged time:                                    "<<(h_cutFlow->GetBinContent(39))<<endl;
  cout<<"40. Baad leading Jet :                                       "<<(h_cutFlow->GetBinContent(40))<<endl;

  cout<<"16.  jet size<2:                                             "<<(h_cutFlow->GetBinContent(16))<<endl;

  cout<<"18.  atleast 2 baseline leptons:                             "<<(h_cutFlow->GetBinContent(18))<<endl;
  cout<<"64.  Exactly 2 baseline leptons:                             "<<(h_cutFlow->GetBinContent(64))<<endl;
  cout<<"65.  Exactly 2 tight leptons:                                "<<(h_cutFlow->GetBinContent(65))<<endl;

  cout<<"19.  2 baseline muons:                                       "<<(h_cutFlow->GetBinContent(19))<<endl;
  cout<<"54.  2 base mu+opposite sign:                                "<<(h_cutFlow->GetBinContent(54))<<endl;
  cout<<"20.  2 tight muons :                                         "<<(h_cutFlow->GetBinContent(20))<<endl;
  cout<<"22.  2 tight muon+0S:                                        "<<(h_cutFlow->GetBinContent(22))<<endl;   


  //cout<<"66.  2 tight muons :                                         "<<(h_cutFlow->GetBinContent(66))<<endl;
  //cout<<"67.  2 tight muon+0S:                                        "<<(h_cutFlow->GetBinContent(67))<<endl;
  cout<<"21.  2 tight muon+0S+minv>20 :                               "<<(h_cutFlow->GetBinContent(21))<<endl;
  cout<<"32.  2 tight mu+OS+muPt>20+minv>20:                          "<<(h_cutFlow->GetBinContent(32))<<endl;
  cout<<"33.  2 tight mu+OS+muPt>20+minv>20+jet>2:                    "<<(h_cutFlow->GetBinContent(33))<<endl;

  cout<<"23.  2 baseline electrons :                                  "<<(h_cutFlow->GetBinContent(23))<<endl;
  cout<<"24.  2 tight electron :                                      "<<(h_cutFlow->GetBinContent(24))<<endl;
  cout<<"26.  2 tight elecetron + OS:                                 "<<(h_cutFlow->GetBinContent(26))<<endl;
  //cout<<"68.  2 tight electron :                                      "<<(h_cutFlow->GetBinContent(68))<<endl;
  //cout<<"69.  2 tight elecetron + OS:                                 "<<(h_cutFlow->GetBinContent(69))<<endl;


  cout<<"34.  2 tight el+OS+minv>20:                                  "<<(h_cutFlow->GetBinContent(34))<<endl;
  cout<<"35.  2 tight el +OS+minv>20+el_pt>20:                        "<<(h_cutFlow->GetBinContent(35))<<endl;
  cout<<"25.  2 tight el+minv>20+elpt>25+jet>2:                       "<<(h_cutFlow->GetBinContent(25))<<endl;

  cout<<"27.  2 base lepton:                                          "<<(h_cutFlow->GetBinContent(27))<<endl;
  cout<<"30.  2 baseline lepton + OS:                                 "<<(h_cutFlow->GetBinContent(61))<<endl; 
  cout<<"28.  2 tight lepton:                                         "<<(h_cutFlow->GetBinContent(28))<<endl;
  cout<<"29.  2 tight lepton + OS:                                    "<<(h_cutFlow->GetBinContent(29))<<endl;

  //cout<<"70.  2 tight lepton:                                         "<<(h_cutFlow->GetBinContent(70))<<endl;
  //cout<<"71.  2 tight lepton + OS:                                    "<<(h_cutFlow->GetBinContent(71))<<endl;
  cout<<"36.  2 tight lep+OS+minv>20:                                 "<<(h_cutFlow->GetBinContent(36))<<endl;
  cout<<"37.  2 tight lep+OS+minv>20+el_pt>20:                        "<<(h_cutFlow->GetBinContent(37))<<endl;
  cout<<"38.  2 tight lep+OS+minv>20+el_pt>20+jet>2:                  "<<(h_cutFlow->GetBinContent(38))<<endl;
  /* 
  cout<<"42.  2 tight muon + C1 cut(met>50):                          "<<(h_cutFlow->GetBinContent(42))<<endl;
  cout<<"43.  2 tight muon + C2 cut(met>80):                          "<<(h_cutFlow->GetBinContent(43))<<endl;
  cout<<"44.  2 tight muon + C3 cut(met>50,MuPt>50):                  "<<(h_cutFlow->GetBinContent(44))<<endl;
  cout<<"45.  2 tight muon + C4 cut(met>80,MuPt>50):                  "<<(h_cutFlow->GetBinContent(45))<<endl;

  cout<<"46.  2 tight el + C1 cut(met>50):                            "<<(h_cutFlow->GetBinContent(46))<<endl;
  cout<<"47.  2 tight el + C2 cut(met>80):                            "<<(h_cutFlow->GetBinContent(47))<<endl;
  cout<<"48.  2 tight el + C3 cut(met>50,MuPt>50):                    "<<(h_cutFlow->GetBinContent(48))<<endl;
  cout<<"49.  2 tight el + C4 cut(met>80,MuPt>50):                    "<<(h_cutFlow->GetBinContent(49))<<endl;

  cout<<"50.  2 tight lep + C1 cut(met>50):                           "<<(h_cutFlow->GetBinContent(50))<<endl;
  cout<<"51.  2 tight lep + C2 cut(met>80):                           "<<(h_cutFlow->GetBinContent(51))<<endl;
  cout<<"52.  2 tight lep + C3 cut(met>50,MuPt>50):                   "<<(h_cutFlow->GetBinContent(52))<<endl;
  cout<<"53.  2 tight lep + C4 cut(met>80,MuPt>50):                   "<<(h_cutFlow->GetBinContent(53))<<endl;
  */
  cout<<"55.  2 base mu+Same sign:                                    "<<(h_cutFlow->GetBinContent(55))<<endl;
  cout<<"54.  2 base mu+opposite sign:                                "<<(h_cutFlow->GetBinContent(54))<<endl;
  cout<<"22.  2 tight mu+opposite sign:                               "<<(h_cutFlow->GetBinContent(22))<<endl;
  cout<<"56.  2 tight mu+Same sign:                                   "<<(h_cutFlow->GetBinContent(56))<<endl;
  cout<<"57.  2 base el+opposite sign:                                "<<(h_cutFlow->GetBinContent(57))<<endl;
  cout<<"58.  2 base el+Same sign:                                    "<<(h_cutFlow->GetBinContent(58))<<endl;
  cout<<"26.  2 tight el+opposite sign:                               "<<(h_cutFlow->GetBinContent(26))<<endl;
  cout<<"59.  2 tight el+Same sign:                                   "<<(h_cutFlow->GetBinContent(59))<<endl;
  cout<<"60.  2 base leptons+opposite sign:                           "<<(h_cutFlow->GetBinContent(60))<<endl;
  cout<<"61.  2 base leptons+Same sign:                               "<<(h_cutFlow->GetBinContent(61))<<endl;
  cout<<"29.  2 tight leptons+opposite sign:                          "<<(h_cutFlow->GetBinContent(29))<<endl;
  cout<<"62.  2 tight leptons+Same sign:                              "<<(h_cutFlow->GetBinContent(62))<<endl;



  /*  
    
 cout<<(h_cutFlow->GetBinContent(1))<<endl;
  cout<<(h_cutFlow->GetBinContent(2))<<endl;
  cout<<(h_cutFlow->GetBinContent(3))<<endl;
  cout<<(h_cutFlow->GetBinContent(4))<<endl;
  cout<<(h_cutFlow->GetBinContent(5))<<endl;
  cout<<(h_cutFlow->GetBinContent(6))<<endl;
  cout<<(h_cutFlow->GetBinContent(7))<<endl;
  cout<<(h_cutFlow->GetBinContent(10))<<endl;
  cout<<(h_cutFlow->GetBinContent(11))<<endl;
  cout<<(h_cutFlow->GetBinContent(15))<<endl;
  cout<<(h_cutFlow->GetBinContent(13))<<endl;

  cout<<(h_cutFlow->GetBinContent(8))<<endl;
  cout<<(h_cutFlow->GetBinContent(9))<<endl;
  cout<<(h_cutFlow->GetBinContent(39))<<endl;
  cout<<(h_cutFlow->GetBinContent(40))<<endl;
 
   
  cout<<(h_cutFlow->GetBinContent(16))<<endl;

  cout<<(h_cutFlow->GetBinContent(18))<<endl;

  cout<<(h_cutFlow->GetBinContent(64))<<endl;
  cout<<(h_cutFlow->GetBinContent(65))<<endl;


  cout<<(h_cutFlow->GetBinContent(19))<<endl;
  cout<<(h_cutFlow->GetBinContent(54))<<endl;
  cout<<(h_cutFlow->GetBinContent(22))<<endl;
  cout<<(h_cutFlow->GetBinContent(20))<<endl;   
  cout<<(h_cutFlow->GetBinContent(21))<<endl;
  cout<<(h_cutFlow->GetBinContent(32))<<endl;
  cout<<(h_cutFlow->GetBinContent(33))<<endl;

  cout<<(h_cutFlow->GetBinContent(23))<<endl;
  cout<<(h_cutFlow->GetBinContent(68))<<endl;
  cout<<(h_cutFlow->GetBinContent(69))<<endl;
  cout<<(h_cutFlow->GetBinContent(34))<<endl;
  cout<<(h_cutFlow->GetBinContent(35))<<endl;
  cout<<(h_cutFlow->GetBinContent(25))<<endl;

  cout<<(h_cutFlow->GetBinContent(27))<<endl;
  cout<<(h_cutFlow->GetBinContent(30))<<endl; 
  cout<<(h_cutFlow->GetBinContent(28))<<endl;
  cout<<(h_cutFlow->GetBinContent(29))<<endl;
  cout<<(h_cutFlow->GetBinContent(36))<<endl;
  cout<<(h_cutFlow->GetBinContent(37))<<endl;
  cout<<(h_cutFlow->GetBinContent(38))<<endl;
 
  cout<<(h_cutFlow->GetBinContent(42))<<endl;
  cout<<(h_cutFlow->GetBinContent(43))<<endl;
  cout<<(h_cutFlow->GetBinContent(44))<<endl;
  cout<<(h_cutFlow->GetBinContent(45))<<endl;

  cout<<(h_cutFlow->GetBinContent(46))<<endl;
  cout<<(h_cutFlow->GetBinContent(47))<<endl;
  cout<<(h_cutFlow->GetBinContent(48))<<endl;
  cout<<(h_cutFlow->GetBinContent(49))<<endl;

  cout<<(h_cutFlow->GetBinContent(50))<<endl;
  cout<<(h_cutFlow->GetBinContent(51))<<endl;
  cout<<(h_cutFlow->GetBinContent(52))<<endl;
  cout<<(h_cutFlow->GetBinContent(53))<<endl;

  cout<<(h_cutFlow->GetBinContent(55))<<endl;
  cout<<(h_cutFlow->GetBinContent(54))<<endl;
  cout<<(h_cutFlow->GetBinContent(22))<<endl;
  cout<<(h_cutFlow->GetBinContent(56))<<endl;
  cout<<(h_cutFlow->GetBinContent(57))<<endl;
  cout<<(h_cutFlow->GetBinContent(58))<<endl;
  cout<<(h_cutFlow->GetBinContent(26))<<endl;
  cout<<(h_cutFlow->GetBinContent(59))<<endl;
  cout<<(h_cutFlow->GetBinContent(60))<<endl;
  cout<<(h_cutFlow->GetBinContent(61))<<endl;
  cout<<(h_cutFlow->GetBinContent(29))<<endl;
  cout<<(h_cutFlow->GetBinContent(62))<<endl;

  

  */


  


}

