bool usePROOF = 1; //0; //set to 0 if you don't wan to use PROOF


TString period = "localtest_data";

void run(bool usePROOF_, TString period_ = "localtest_data"){
  usePROOF = usePROOF_;
  period = period_;
  run();
}



void run() {
  //--- Safty : to avoid overwriting the output root file
  
  TChain *c = new TChain("susy"); //create TChain with Tree susy
  
  
  // examples: how to generate list of input files 
  // for file in $(ls /Volumes/Data_2/xrd/ww/ | grep mc10 ); do echo "c->Add(\"/Volumes/Data_2/xrd/ww/"$file/"*root*\");" ; done 
  //for file in $(ls /Volumes/Data_2/xrd/ww/ | grep data | grep Muons); do echo "c->Add(\"/Volumes/Data_2/xrd/ww/"$file/"*root*\");" ; done
  // for file in $(ls /Volumes/Data_1/public/data12/background_skim7/ | grep Ttbar ); do echo "c->Add(\"/Volumes/Data_1/public/data12/background_skim7/"$file/"*root*\");" ; done
 
  if (period =="top") {
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.117800.Sherpa_CT10_TtbarLeptLept.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150511/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.117801.Sherpa_CT10_TtbarLeptTaulept.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150531/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.117802.Sherpa_CT10_TtbarTauleptTaulept.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150558/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.117803.Sherpa_CT10_TtbarLeptHad.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150621/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.117803.Sherpa_CT10_TtbarLeptHad.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150702/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.117803.Sherpa_CT10_TtbarLeptHad.merge.SUSYD4PD.e1434_s1499_s1504/_r3658_r3549_p1181_V0_skim7.121206150731/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.117803.Sherpa_CT10_TtbarLeptHad.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150753/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.117804.Sherpa_CT10_TtbarLeptTauhad.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150812/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.117805.Sherpa_CT10_TtbarTauleptHad.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150830/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.117806.Sherpa_CT10_TtbarTauleptTauhad.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150848/*root*");
  }
  
  if (period=="diboson") {
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.126892.Sherpa_CT10_llnunu_WW.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206145853/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.126893.Sherpa_CT10_lllnu_WZ.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206145912/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.126894.Sherpa_CT10_llll_ZZ.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206145931/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.126895.Sherpa_CT10_llnunu_ZZ.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206145950/*root*");
  }
  
  if (period=="Zjet") {
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.147770.Sherpa_CT10_Zee.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150025/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.147770.Sherpa_CT10_Zee.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150108/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.147771.Sherpa_CT10_Zmumu.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150152/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.147771.Sherpa_CT10_Zmumu.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150228/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.147772.Sherpa_CT10_Ztautau.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150307/*root*");
    c->Add("/Volumes/Data_1/public/data12/background_skim7/user.sdarmora.mc12_8TeV.147772.Sherpa_CT10_Ztautau.merge.SUSYD4PD.e1434_s1499_s1504_r3658_r3549_p1181_V0_skim7.121206150332/*root*");
  }
  
 
  if (period == "newtestdata2") {
    c->Add("/Volumes/Data_2/xrd/stop_data/signal_new/user.Smita.Darmora.mc12_8TeV.118494.Herwigpp_UEEE3_CTEQ6L1_Tt_T300_L100.merge.SUSYD4PD.e1425_a159_a171_r3549_p1181_V0_skim6.121122130113/*.root");  
  }
  
  if (period == "newtestdata1") {
    c->Add("/Volumes/Data_1/public/data12/user.Smita.Darmora.mc12_8TeV.118494.Herwigpp_UEEE3_CTEQ6L1_Tt_T300_L100.merge.SUSYD4PD.e1425_a159_a171_r3549_p1181_V0_skim6.121122130113/*.root");  
  }
 
if (period=="testcutflow") {
c->Add("/Volumes/Data_1/public/data12/user.sdarmora.mc12_8TeV.165884.Herwigpp_UEEE3_CTEQ6L1_pMSSM_2692842_120FT_neut.merge.SUSYD4PD.e1469_a159_a171_r3549_p1181_V0_skim7.121212114426/*.root"); 
}	 
  if (period == "signalall") {
    //laurelin:~ gusai$ for file in $(ls  /Volumes/Data_1/public/data12/signal_skim7/*/*  | grep root ); do  echo "c->Add(\""$file"\");" ; done 
c->Add("/Volumes/Data_1/public/data12/signal_skim7/user.sdarmora.mc12_8TeV.118494.Herwigpp_UEEE3_CTEQ6L1_Tt_T300_L100.merge.SUSYD4PD.e1425_a159_a171_r3549_p1181_V0_skim7.121206135451/user.sdarmora.001948._00001.susy.root");
c->Add("/Volumes/Data_1/public/data12/signal_skim7/user.sdarmora.mc12_8TeV.118495.Herwigpp_UEEE3_CTEQ6L1_Tt_T300_L1.merge.SUSYD4PD.e1425_a159_a171_r3549_p1181_V0_skim7.121206135508/user.sdarmora.001950._00001.susy.root");
c->Add("/Volumes/Data_1/public/data12/signal_skim7/user.sdarmora.mc12_8TeV.164649.Herwigpp_UEEE3_CTEQ6L1_Tt_T600_L1.merge.SUSYD4PD.e1394_a159_a171_r3549_p1181_V0_skim7.121206135240/user.sdarmora.001936._00003.susy.root");
}

   
  gROOT->ProcessLine(".x /Users/gusai/susy12/RootCore/scripts/load_packages.C+");   
 
  if (usePROOF) {
    p = TProof::Open("workers=2");
    c->SetProof();
    p->UploadPackage("root_core.par");
    p->EnablePackage ("root_core.par"); 
    p->ShowEnabledPackages();
  }

  //TSelector *selector = TSelector::GetSelector("/Users/gusai/susy12/Read_d4pd/Root/read_d4pd.cxx+");
  //c->Process(selector);
  c->Process("/Users/gusai/susy12/myown_util/scripts/read_d4pd.cxx++");

  gROOT->ProcessLine(".! mv histograms.root  histograms."+period+".root");
}



void run_singlefile(bool usePROOF, TString name_list) {   // run on  a list of file provided in the inputfiles list.    run on each file separatelly and save individual output histograms 


  //for file in $(ls /Volumes/Data_1/public/data12/signal_skim7/ ); do echo -n "/Volumes/Data_1/public/data12/signal_skim7/"$file"/*root*"; echo $file | awk -F. '{print "\t",$4}'; done >myown_util/scripts/filelist.signal_stopgrid.input

  TString filename="/Users/gusai/susy12/myown_util/scripts/filelist.";
  filename+=+name_list;
  filename+=".input";
  stc::cout <<"reading from: "<< filename<<'\n';

  ifstream *  ifl = new ifstream(filename);
  TString fnm;
  TString IDnm;
  
  while(1)
    {
      *ifl >> fnm >>IDnm ;
      if (ifl->eof()) break;

      TChain *c = new TChain("susy"); //create TChain with Tree susy
      stc::cout <<"adding : "<< fnm<<'\t' <<IDnm <<'\n';
      c->Add(fnm);
      
      
      gROOT->ProcessLine(".x /Users/gusai/susy12/RootCore/scripts/load_packages.C+");   
      
      if (usePROOF) {
	p = TProof::Open("workers=2");
	c->SetProof();
	p->UploadPackage("root_core.par");
	p->EnablePackage ("root_core.par"); 
	p->ShowEnabledPackages();
      }
      
      std::cout<<"All libraries loaded"<<std::endl;
      
      //TSelector *selector = TSelector::GetSelector("/Users/gusai/susy12/Read_d4pd/Root/read_d4pd.cxx+");
      //c->Process(selector);
      c->Process("/Users/gusai/susy12/myown_util/scripts/read_d4pd.cxx++");
      
      
      gROOT->ProcessLine(".! mv histograms.root  histograms."+name_list+"_"+IDnm+".root");
      delete c; 
    }
  
}



count_events(TString name_list)
{

TString filename="/Users/sdarmora/filelist.";
  filename+=+name_list;
  filename+=".input";
  stc::cout <<"reading from: "<< filename<<'\n';

  ifstream *  ifl = new ifstream(filename);
  TString fnm;
  TString IDnm;

  while(1)
    {
      *ifl >> fnm;
      if (ifl->eof()) break;


  
  TChain *c = new TChain("susy");
  TChain *c1 = new TChain("susy");
// stc::cout <<"adding : "<< fnm<<'\t' <<IDnm <<'\n';
c->Add(fnm);
//c->Add(IDnm);
//c->Add("/Volumes/Data_1/public/data12/signal.p1328.15jan_left/user.gusai.mc12.176968.Herwigpp_UEEE3_CTEQ6L1_Tt_T550_L50_left.NTUP.e1721_a188_a171_r3549_p1328.15jan.130130121902/*root*);
//c->Add("/Volumes/Data_1/public/data12/signal.p1328.15jan_left/user.gusai.mc12.176968.Herwigpp_UEEE3_CTEQ6L1_Tt_T550_L50_left.NTUP.e1721_a188_a171_r3549_p1328.15jan.130325135502/*root*");

  TBranch        *b_mc_channel_number;   //!
  UInt_t          mc_channel_number;

  TBranch    *b_RunNumber;
  UInt_t      RunNumber;

   TBranch        *b_mcevt_weight;
  TBranch        *b_mcevt_n;
  Int_t           mcevt_n;
  vector<vector<double> > *mcevt_weight;    

  c->SetBranchAddress("mc_channel_number", &mc_channel_number, &b_mc_channel_number);
  //c->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
  //c->GetEntry(1);
    

  //c->SetBranchAddress("mcevt_weight", &mcevt_weight, &b_mcevt_weight);
  // c->SetBranchAddress("mcevt_n", &mcevt_n, &b_mcevt_n);
        c->GetEntry(1);
    //  c1->GetEntry(1);
    //       b_mcevt_weight->GetEntry(entry);
    double mcevtw=0;  
 
    //mcevt_n++;

    
    //int size_weight=0;
    //size_weight = mcevt_weight->at(0).size();
    //for(int i=0;i<mcevt_n;i++){      
    //if(size_weight!=0){
    //mcevtw=(mcevt_weight->at(0)).at(0); 
    //mcevtw+=mcevtw;
    //}
    //  }
 

      //c->GetEntry(1);
    
        std::cout << mc_channel_number<< "        "<<c->GetEntries() <<'\n';//change
    //    std::cout <<RunNumber<< "        "<<c->GetEntries() <<'\n';
    
    //   std::cout << mc_channel_number<< "        "<<c->mcevtw<<"  "<<mcevt_n<<'\n';
  delete c; 
}
}
