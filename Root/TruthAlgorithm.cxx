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
      
      if (TMath::Abs((*truthV_itr)->incomingParticle(iIn)->pdgId()) == 
        m_pdgIdOfMother) { 
        
        for (unsigned int iOut=0; iOut < (*truthV_itr)->nOutgoingParticles(); 
              iOut++) {
          truthParticle = (*truthV_itr)->outgoingParticle(iOut);
          unsigned absPdgId = TMath::Abs(truthParticle->pdgId());
          int status = truthParticle->status();
//           if (status!=3)
//             continue;

          if (absPdgId==m_leptonPdgId)
            if (leptonVec.Pt()<(truthParticle->p4()).Vect().Pt())
              leptonVec = (truthParticle->p4()).Vect();
          if (absPdgId==m_leptonPdgId+1)
            if (MET.Pt()<(truthParticle->p4()).Vect().Pt())
              MET = (truthParticle->p4()).Vect();
          debugStream << "(" << absPdgId << "," << status << "," << 
          (truthParticle->p4()).Vect().Pt()*GEV << ") "; 
          counter++;
        }
        
      }
      double invMass = sqrt( 2*leptonVec.Mag()*MET.Mag() * 
                       (1.0 - TMath::Cos( leptonVec.Phi() - MET.Phi())))*GEV;
                       
//       if (counter>=2)
//         cout << "[DEBUG]\t" << m_EventNumber << ": " << 
// //         (*truthV_itr)->incomingParticle(iIn)->pdgId() << " --> " << 
//         debugStream.str() << " --> " << invMass << endl;
                       
      if (MET.Pt()>0 && leptonVec.Pt()>0 && 
          ((!m_cut120GeVForInclusiveW) || invMass<120.0) ){
//         cout << "[DEBUG]\tinvMass = " << invMass << endl;
        fillHist(MET,leptonVec);
        m_BitsetCutflow->FillCutflow("found pair");
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
  
  cout << "[INFO]\t\tm_weightCrossSection = " << m_weightCrossSection << endl;
  
  return EL::StatusCode::SUCCESS;
}

void TruthAlgorithm :: fillHist (TVector3 MET, TVector3 leptonVec)
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