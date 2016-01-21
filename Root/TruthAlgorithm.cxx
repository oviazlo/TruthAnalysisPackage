#include <TruthAnalysis/TruthAlgorithm.h>

#define GEV 0.001

/// this is needed to distribute the algorithm to the workers
ClassImp(TruthAlgorithm)

TruthAlgorithm :: TruthAlgorithm ()
{
}



EL::StatusCode TruthAlgorithm :: setupJob (EL::Job& job)
{
  /// let's initialize the algorithm to use the xAODRootAccess package
  job.useXAOD ();
  xAOD::Init( "TruthAlgorithm" ).ignore(); /// call before opening first file
  
  m_runElectronChannel = false;
  m_debugMode = false;
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthAlgorithm :: histInitialize ()
{
  hMu_pt_off = (TH1D*)WprimeHist::standard("pt","h","","");
  wk()->addOutput(hMu_pt_off); 
  
  hMu_mt_off = (TH1D*)WprimeHist::standard("mt","h","","");
  wk()->addOutput(hMu_mt_off);

  hMu_MET_Muons_off = (TH1D*)WprimeHist::standard("met","h","","");
  wk()->addOutput(hMu_MET_Muons_off); 
  
  hMu_invMass_Muons_off = (TH1D*)WprimeHist::standard("mgen","h","","");
  wk()->addOutput(hMu_invMass_Muons_off); 
  
  hMu_invMass_Muons_alternative = (TH1D*)WprimeHist::standard("mgen","hAlt","",
                                                              "");
  wk()->addOutput(hMu_invMass_Muons_alternative); 
  
  h_event_crossSectionWeight = (TH1D*)WprimeHist::standard("evtwt_xSec","h","",
                                                           "");
  wk()->addOutput(h_event_crossSectionWeight);
  
  h_event_kFactor = (TH1D*)WprimeHist::standard("evtwt_kFactor","h","","");
  wk()->addOutput(h_event_kFactor);
  
  h_event_filterEfficiency = (TH1D*)WprimeHist::standard("evtwt_filterEff","h",
                                                         "","");
  wk()->addOutput(h_event_filterEfficiency);
  
  h_event_totalWeight = (TH1D*)WprimeHist::standard("evtwt","h","","");
  wk()->addOutput(h_event_totalWeight);
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthAlgorithm :: initialize ()
{
  m_pdgIdOfMother = 24; /// W by default
  
  string sampleName = wk()->metaData()->getString ("sample_name");
  
  if (sampleName.find("Pythia8EvtGen_A14NNPDF23LO_Wprime_")!=std::string::npos)
    m_pdgIdOfMother = 34; /// Wprime
  if (sampleName.find("PowhegPythia8EvtGen_AZNLOCTEQ6L1_W")!=
    std::string::npos)
    m_pdgIdOfMother = 24; /// W
    
  m_cut120GeVForInclusiveW = false;
  if ((sampleName.find("361101")!=std::string::npos) || /// incl. Wplusmunu
      (sampleName.find("361100")!=std::string::npos) || /// incl. Wplusenu
      (sampleName.find("361104")!=std::string::npos) || /// incl. Wminmunu
      (sampleName.find("361103")!=std::string::npos))   /// incl. Wminenu
  {
    m_cut120GeVForInclusiveW = true;
    cout << "[WARNING]\tm_cut120GeVForInclusiveW is activated!!!" << endl;
  }
  
  cout << "[JobSetup]\tm_pdgIdOfMother = " << m_pdgIdOfMother << endl;
  
  m_event = wk()->xaodEvent();
  m_store = new xAOD::TStore();
  
  /// as a check, let's see the number of events in our xAOD
  /// print long long int
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); 

  /// Event Info
  const xAOD::EventInfo* m_eventInfo = 0;
  if( ! m_event->retrieve( m_eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", 
          "Failed to retrieve event info collection in initialise. Exiting." );
    return EL::StatusCode::FAILURE;
  }  
    
  /// fill the branches of our trees
  m_isMC = false;
  if(m_eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    m_isMC = true;
  }
  cout << "[JobSetup]\tm_isMC flag = " << m_isMC <<endl;
  
  m_weightkFactor = 1.0;
  m_weighfilterEfficiency = 1.0;
  m_weightCrossSection = 1.0;
  
  if (m_isMC){
    m_LPXKfactorTool = new LPXKfactorTool("LPXKfactorTool");
    EL_RETURN_CHECK("m_LPXKfactorTool_isMC15",
                    m_LPXKfactorTool->setProperty("isMC15", true)); 
    EL_RETURN_CHECK("m_LPXKfactorTool_applyEWCorr",
                    m_LPXKfactorTool->setProperty("applyEWCorr", true)); 
    EL_RETURN_CHECK("m_LPXKfactorTool_applyPICorr",
                    m_LPXKfactorTool->setProperty("applyPICorr", true)); 
    
    EL_RETURN_CHECK( "m_LPXKfactorTool initialize",m_LPXKfactorTool->
    initialize());
  }
  
  if (m_BitsetCutflow)
    m_BitsetCutflow = new BitsetCutflow(wk());
  
  if (m_runElectronChannel){
    cout << "[JobSetup]\trun over ELECTRON channel" <<endl;
    m_leptonPdgId = 11;
  }
  else{
    cout << "[JobSetup]\trun over MUON channel" <<endl;
    m_leptonPdgId = 13;
  }
  
  trueWmass = 0.0;
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthAlgorithm :: execute ()
{
  /// push cutflow bitset to cutflow hist
  m_BitsetCutflow->PushBitSet();
  
  EL_RETURN_CHECK("retrieve EventInfo",
                  m_event->retrieve( m_eventInfo, "EventInfo"));
  m_EventNumber = m_eventInfo->eventNumber();  
    
  const xAOD::TruthVertexContainer* truthVertices = 0;
  EL_RETURN_CHECK("retrieve TruthVertices", 
                  m_event->retrieve( truthVertices, "TruthVertices" ));
 
  bool foundPair = false;
  
  /// Start iterating over truth container
  const xAOD::TruthParticle* truthParticle;
  xAOD::TruthVertexContainer::const_iterator truthV_itr; 
  
//   cout << "[DEBUG]\tm_leptonPdgId = " << m_leptonPdgId << endl;
  
  for (truthV_itr = truthVertices->begin(); 
        truthV_itr != truthVertices->end(); ++truthV_itr )
  {
    for (unsigned int iIn=0; iIn < (*truthV_itr)->nIncomingParticles(); 
          iIn++)
    {
      TVector3 MET(0,0,0);
      TVector3 leptonVec(0,0,0);
      unsigned int counter = 0;
      
      stringstream debugStream;
      double eLepton = 0.0;
      double eNeutrino = 0.0;
      
      const xAOD::TruthParticle* mother;
      const xAOD::TruthParticle* lepton; 
      const xAOD::TruthParticle* neutrino;
      
      if (TMath::Abs((*truthV_itr)->incomingParticle(iIn)->pdgId()) == 
        m_pdgIdOfMother) { 
        
        trueWmass = 
        TMath::Sqrt(TMath::Power((*truthV_itr)->incomingParticle(iIn)->e(),2)
              - TMath::Power((*truthV_itr)->incomingParticle(iIn)->pz(),2)
              - TMath::Power((*truthV_itr)->incomingParticle(iIn)->py(),2)
              - TMath::Power((*truthV_itr)->incomingParticle(iIn)->px(),2)); 
        
        mother = (*truthV_itr)->incomingParticle(iIn);
        
        for (unsigned int iOut=0; iOut < (*truthV_itr)->nOutgoingParticles(); 
              iOut++) {
          truthParticle = (*truthV_itr)->outgoingParticle(iOut);
          unsigned absPdgId = TMath::Abs(truthParticle->pdgId());
          int status = truthParticle->status();
//           if (status!=3)
//             continue;

          if (absPdgId==m_leptonPdgId)
            if (leptonVec.Pt()<(truthParticle->p4()).Vect().Pt()){
              leptonVec = (truthParticle->p4()).Vect();
              eLepton = truthParticle->e();
              lepton = truthParticle;
            }
          if (absPdgId==m_leptonPdgId+1)
            if (MET.Pt()<(truthParticle->p4()).Vect().Pt()){
              MET = (truthParticle->p4()).Vect();
              eNeutrino = truthParticle->e();
              neutrino = truthParticle;
            }
          debugStream << "(" << absPdgId << "," << status << "," << 
          (truthParticle->p4()).Vect().Pt()*GEV << ") "; 
          counter++;
        }
        
      }
//       double invMass = sqrt( 2*leptonVec.Mag()*MET.Mag() * 
//                        (1.0 - TMath::Cos( leptonVec.Phi() - MET.Phi())))*GEV;

      double invMass = TMath::Sqrt((TMath::Power(eLepton+eNeutrino,2) - 
                TMath::Power( (leptonVec+MET).Mag() ,2)))*GEV;
      

                       
      if (counter>=2 && m_debugMode)
        cout << "[DEBUG]\t" << m_EventNumber << ": " << 
        (*truthV_itr)->incomingParticle(iIn)->pdgId() << " --> " << 
        debugStream.str() << " --> " << invMass << endl;
      
      
                       
      if (MET.Pt()>0 && leptonVec.Pt()>0 && 
          ((!m_cut120GeVForInclusiveW) || invMass<120.0) ){
        fillHist(MET,leptonVec,trueWmass);
        fillHist(mother,lepton,neutrino);
        m_BitsetCutflow->FillCutflow("found pair");
        
//           cout << "[DEBUG]\t" << m_EventNumber << ": " << invMass << 
//           " (calc.inv.mass)" << endl;
//           cout << "[DEBUG]\t" << m_EventNumber << ": " << trueWmass*GEV << 
//           " (trueWmass)" << endl;
        
      }
    }
  }
  
  ///***************************************************
  /// TEST mono-W sample
  /*
  for (truthV_itr = truthVertices->begin(); 
      truthV_itr != truthVertices->end(); ++truthV_itr )
  {
    
    stringstream debugStream;
    for (unsigned int iOut=0; iOut < (*truthV_itr)->nOutgoingParticles(); 
            iOut++) {
      const xAOD::TruthParticle* truthOutParticle = 
      (*truthV_itr)->outgoingParticle(iOut);
      if (!truthOutParticle)
        continue;
      int absPdgId = TMath::Abs(truthOutParticle->pdgId());
      int status = truthOutParticle->status();
      if (absPdgId!=m_leptonPdgId)
        continue;
      if ((truthOutParticle->p4()).Vect().Pt()*GEV<100.0)
        continue;
      for (unsigned int iIn=0; iIn < (*truthV_itr)->nIncomingParticles(); 
        iIn++)
      {
        const xAOD::TruthParticle* truthInParticle = 
        (*truthV_itr)->incomingParticle(iIn);
        absPdgId = TMath::Abs(truthInParticle->pdgId());
        debugStream << absPdgId << " "; 
      }
      cout << "[DEBUG]\t" << m_EventNumber << ": " << 
      debugStream.str() << endl;
    }
    
  }*/
  ///***************************************************

  
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TruthAlgorithm :: finalize ()
{
  if (m_LPXKfactorTool){
    delete m_LPXKfactorTool;
    m_LPXKfactorTool = 0;
  }
  
  /// push cutflow for last event
  /// FIXME probably this breaks possibility to use PROOF
  m_BitsetCutflow->PushBitSet();
  
  m_BitsetCutflow->PrintCutflowLocally();
  
  if(m_BitsetCutflow){
    delete m_BitsetCutflow;
    m_BitsetCutflow = 0;
  }
  
  cout << "[INFO]\t\tm_weightCrossSection = " << m_weightCrossSection 
  << " nb" << endl;
  
  return EL::StatusCode::SUCCESS;
}


void TruthAlgorithm :: fillHist (const xAOD::TruthParticle* mother, 
                                 const xAOD::TruthParticle* lepton, 
                                 const xAOD::TruthParticle* neutrino)
{
  
  TVector3 leptonVec = (lepton->p4()).Vect();
  TVector3 neutrinoVec = (neutrino->p4()).Vect();
  
  double trueWmass = 
        TMath::Sqrt(TMath::Power(mother->e(),2)
              - TMath::Power(mother->pz(),2)
              - TMath::Power(mother->py(),2)
              - TMath::Power(mother->px(),2))*GEV;
  
  double leptonEt = TMath::Sqrt(TMath::Power(leptonVec.Pt(),2) + 
                                TMath::Power(lepton->m(),2));
        
  double neutrinoEt = TMath::Sqrt(TMath::Power(neutrinoVec.Pt(),2) + 
                                TMath::Power(neutrino->m(),2));
              
  cout << "lepton->m() = " << lepton->m()*GEV << "; neutrino->m() = " <<
  neutrino->m()*GEV << endl;
  
  cout << "leptonEt = " << leptonEt*GEV << "; leptonVec.Pt() = " << 
  leptonVec.Pt()*GEV << endl;
  
  cout << "neutrinoEt = " << neutrinoEt*GEV << "; neutrinoVec.Pt() = " 
  << neutrinoVec.Pt()*GEV << endl;
    
  double calculatedTrueWmass_v1 = 
        TMath::Sqrt((TMath::Power(lepton->e()+neutrino->e(),2) -
                    TMath::Power( (leptonVec+neutrinoVec).Mag() ,2)))*GEV;
  
  double calculatedTrueWmass_v2 = 
        TMath::Sqrt( 2*leptonVec.Mag()*neutrinoVec.Mag() * 
          (1.0 - TMath::Cos( leptonVec.Angle(neutrinoVec) ) ) )*GEV;
  
  double trueMt = 
        TMath::Sqrt(TMath::Power(leptonEt+neutrinoEt,2) -
                    TMath::Power((leptonVec+neutrinoVec).Pt(),2)
                    )*GEV;
  
                    
  double calculatedTrueMt_v1 = TMath::Sqrt(2*leptonEt*neutrinoEt*
                      (1-TMath::Cos(leptonVec.Phi()-neutrinoVec.Phi())))*GEV;
                      
  double calculatedTrueMt_v2 = TMath::Sqrt(2*leptonVec.Pt()*neutrinoVec.Pt()*
                      (1-TMath::Cos(leptonVec.Phi()-neutrinoVec.Phi())))*GEV;
 
  cout << "trueWmass_v0: " << trueWmass << endl;
  cout << "trueWmass_v1: " << calculatedTrueWmass_v1 << endl;
  cout << "trueWmass_v2: " << calculatedTrueWmass_v2 << endl;
  
  cout << "TrueMt_v0: " << trueMt << endl;
  cout << "TrueMt_v1: " << calculatedTrueMt_v1 << endl;
  cout << "TrueMt_v2: " << calculatedTrueMt_v2 << endl;
  cout << endl;
                      
}

void TruthAlgorithm :: fillHist (TVector3 MET, TVector3 leptonVec, 
                                 double containerInvMass)
{

  double Mt_MET = sqrt( 2*leptonVec.Pt()*MET.Pt() * 
  (1.0 - TMath::Cos( leptonVec.Phi() - MET.Phi())) );      

  double invMass = sqrt( 2*leptonVec.Mag()*MET.Mag() * 
  (1.0 - TMath::Cos( leptonVec.Phi() - MET.Phi())) );

  /// get MC weights
  if (m_isMC){
    m_LPXKfactorTool->execute();
    m_weightkFactor = m_eventInfo->auxdecor<double>("KfactorWeight");
    m_weighfilterEfficiency = m_LPXKfactorTool->getMCFilterEfficiency();
    ///TODO make proper implementation
    m_weightCrossSection = m_LPXKfactorTool->getMCCrossSection(); 
  }
  
  h_event_crossSectionWeight->Fill(m_weightCrossSection);
  h_event_kFactor->Fill(m_weightkFactor);
  h_event_filterEfficiency->Fill(m_weighfilterEfficiency);
  
  double totalWeight = m_weighfilterEfficiency*m_weightkFactor
                                              *m_weightCrossSection;
  
  if (totalWeight==0) totalWeight = 1.0;/// FIXME just for mono-W sample
                                              
  h_event_totalWeight->Fill(totalWeight);
  
  hMu_pt_off->Fill(leptonVec.Pt()*GEV,totalWeight);
  hMu_mt_off->Fill(Mt_MET*GEV,totalWeight);
  hMu_MET_Muons_off->Fill(MET.Pt()*GEV,totalWeight);
  hMu_invMass_Muons_off->Fill(invMass*GEV,totalWeight);
  
  hMu_invMass_Muons_alternative->Fill(containerInvMass*GEV,totalWeight);

}



EL::StatusCode TruthAlgorithm :: fileExecute ()
{
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TruthAlgorithm :: changeInput (bool firstFile)
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TruthAlgorithm :: histFinalize ()
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TruthAlgorithm :: postExecute ()
{
  return EL::StatusCode::SUCCESS;
}