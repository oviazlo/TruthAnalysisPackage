#include <TruthAnalysis/TruthAlgorithm.h>

#define GEV 0.001
#define PDGID_DM 1000022

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
  m_sampleName = wk()->metaData()->getString ("sample_name");
    
  string tmpSampleName = "";
  stringstream strStream(m_sampleName);
  getline (strStream, tmpSampleName, '.'); /// return mc15_13TeV
  getline (strStream, tmpSampleName, '.'); /// return DSID
  m_datasetID = atoi(tmpSampleName.c_str());
    
  m_TruthHelper = new TruthHelper();
  
  cout << "[DEBUG]\tm_datasetID = " << m_datasetID << endl;
  
  m_normalizationFactor = 
  m_TruthHelper->getnEventsPerSample(m_datasetID);
  if (m_normalizationFactor==0.0){
    ///FIXME describe this warning
    m_normalizationFactor = 1.0; 
  }
  
  cout << "[DEBUG]\tm_normalizationFactor = " << m_normalizationFactor << endl;
  
  
  
  if (m_sampleName.find("Pythia8EvtGen_A14NNPDF23LO_Wprime_")!=std::string::npos)
    m_pdgIdOfMother = 34; /// Wprime
  if (m_sampleName.find("PowhegPythia8EvtGen_AZNLOCTEQ6L1_W")!=
    std::string::npos)
    m_pdgIdOfMother = 24; /// W
    
  m_cut120GeVForInclusiveW = false;
  if (m_datasetID==361101 || /// incl. Wplusmunu
      m_datasetID==361100 || /// incl. Wplusenu
      m_datasetID==361104 || /// incl. Wminmunu
      m_datasetID==361103)   /// incl. Wminenu
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
                    m_LPXKfactorTool->setProperty("isMC15", false)); /// FIXME make it back to true
    EL_RETURN_CHECK("m_LPXKfactorTool_applyEWCorr",
                    m_LPXKfactorTool->setProperty("applyEWCorr", true)); 
    EL_RETURN_CHECK("m_LPXKfactorTool_applyPICorr",
                    m_LPXKfactorTool->setProperty("applyPICorr", true)); 
    
    EL_RETURN_CHECK( "m_LPXKfactorTool initialize",m_LPXKfactorTool->
    initialize());
    m_weighfilterEfficiency = m_LPXKfactorTool->getMCFilterEfficiency();
    
    m_weightCrossSection = m_LPXKfactorTool->getMCCrossSection()*1000.0; 
    if (m_weightCrossSection==0)
      m_weightCrossSection = m_TruthHelper->getCrossSection(m_datasetID)
      *1000.0;
    
    if (m_weightCrossSection==0){
      m_weightCrossSection = 1.0;
      cout << "[WARNING]\tweightCrossSection = 0; "
      "Assign it equal to 1." << endl;
    }
    if (m_weighfilterEfficiency==0){
      m_weighfilterEfficiency = 1.0;
      cout << "[WARNING]\tweighfilterEfficiency = 0; "
      "Assign it equal to 1." << endl;
    }
  }
    
  
  
  cout << "[INFO]\tm_weighfilterEfficiency = " << m_weighfilterEfficiency <<
  endl;
  cout << "[INFO]\tm_weightCrossSection = " << m_weightCrossSection << " pb" <<
  endl;
  cout << "[INFO]\tN_evnt = " << m_normalizationFactor << endl;
  cout << "[INFO]\tN_evnt/sigma = " << 
  m_normalizationFactor/m_weightCrossSection << endl;
  
  
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
  
  m_trueWmass = 0.0;
  
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

//   EL_RETURN_CHECK("retrieve TruthVertices", 
//                   m_event->retrieve( truthVertices, "TruthVertices" ));
  /// FIXME uncomment above, remove below   
  EL_RETURN_CHECK("retrieve TruthVertices", 
                  m_event->retrieve( truthVertices, "TruthVertex" ));
  
  
  bool foundPair = false;
  
  /// Start iterating over truth container
  const xAOD::TruthParticle* truthParticle;
  xAOD::TruthVertexContainer::const_iterator truthV_itr; 
  const xAOD::TruthParticle* dmPart1 = NULL;
  const xAOD::TruthParticle* dmPart2 = NULL;
  
//   cout << "[DEBUG]\tm_leptonPdgId = " << m_leptonPdgId << endl;
  
//   cout << "m_EventNumber = " << m_EventNumber << endl;
  
  for (truthV_itr = truthVertices->begin(); 
        truthV_itr != truthVertices->end(); ++truthV_itr )
  {
    for (unsigned int iOut=0; iOut < (*truthV_itr)->nOutgoingParticles(); 
              iOut++) {
      truthParticle = (*truthV_itr)->outgoingParticle(iOut);
      if (!truthParticle)
        continue;
      unsigned absPdgId = TMath::Abs(truthParticle->pdgId());
      int status = truthParticle->status();
      if (status==1&&absPdgId==PDGID_DM){
        if (!dmPart1) 
          dmPart1 = truthParticle;
        else if (!dmPart2) 
          dmPart2 = truthParticle;
        else
          cout << "[ERROR]\tFound 3rd DM particle!!! Abort!!!" << endl;
      }
      
//       double pt = sqrt(TMath::Power(truthParticle->py(),2)
//               + TMath::Power(truthParticle->px(),2));
//       double mass = truthParticle->m();
// //       if (mass/1000.0>0.9)
//         cout << "pdg: " << absPdgId << "; status: " << status << "; mass = " <<
//         mass << "; pt = " << pt << endl;
    } 
  }
  if (!dmPart1 || !dmPart2){
    cout << "[ERROR]\tCan't find two DM particles! Skip this event!!!" << endl;
    return EL::StatusCode::SUCCESS;
  }
  
  for (truthV_itr = truthVertices->begin(); 
        truthV_itr != truthVertices->end(); ++truthV_itr )
  {
    
    /// let's find W/W' daughters
    for (unsigned int iIn=0; iIn < (*truthV_itr)->nIncomingParticles(); 
          iIn++)
    {
      TVector3 neutrinoVec(0,0,0);
      TVector3 leptonVec(0,0,0);
      unsigned int counter = 0;
      
      stringstream debugStream;
      
      const xAOD::TruthParticle* mother; 
      const xAOD::TruthParticle* lepton; 
      const xAOD::TruthParticle* neutrino;
      
      if (TMath::Abs((*truthV_itr)->incomingParticle(iIn)->pdgId()) == 
        m_pdgIdOfMother) { 
        
        mother = (*truthV_itr)->incomingParticle(iIn);
        
        m_trueWmass = 
        TMath::Sqrt(TMath::Power((*truthV_itr)->incomingParticle(iIn)->e(),2)
              - TMath::Power((*truthV_itr)->incomingParticle(iIn)->pz(),2)
              - TMath::Power((*truthV_itr)->incomingParticle(iIn)->py(),2)
              - TMath::Power((*truthV_itr)->incomingParticle(iIn)->px(),2)
                   )*GEV; 

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
              lepton = truthParticle;
            }
          if (absPdgId==m_leptonPdgId+1)
            if (neutrinoVec.Pt()<(truthParticle->p4()).Vect().Pt()){
              neutrinoVec = (truthParticle->p4()).Vect();
              neutrino = truthParticle;
            }
          debugStream << "(" << absPdgId << "," << status << "," << 
          (truthParticle->p4()).Vect().Pt()*GEV << ") "; 
          counter++;
        }
        
      }

      if (counter>=2 && m_debugMode)
        cout << "[DEBUG]\t" << m_EventNumber << ": " << 
        (*truthV_itr)->incomingParticle(iIn)->pdgId() << " --> " << 
        debugStream.str() << " --> " << m_trueWmass << endl;
                       
      if (neutrinoVec.Pt()>0 && leptonVec.Pt()>0 && 
          ((!m_cut120GeVForInclusiveW) || m_trueWmass<120.0) ){
        fillHist(mother,lepton,neutrino,dmPart1,dmPart2);
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
  
  return EL::StatusCode::SUCCESS;
}


void TruthAlgorithm :: fillHist (const xAOD::TruthParticle* mother,
                                 const xAOD::TruthParticle* lepton, 
                                 const xAOD::TruthParticle* neutrino,
                                 const xAOD::TruthParticle* dmPart1,
                                 const xAOD::TruthParticle* dmPart2
                                )
{
  
  TVector3 leptonVec = (lepton->p4()).Vect();
  TVector3 neutrinoVec = (neutrino->p4()).Vect();
  if (!dmPart1||!dmPart2)
    cout << "[ERROR]\tCan't find two DM particles!!!" << endl;
  
//   cout << "before\tneutrinoVec.Mag() = " << neutrinoVec.Mag() << endl;
  
  TVector3 dm1Vec = (dmPart1->p4()).Vect();
  TVector3 dm2Vec = (dmPart2->p4()).Vect();
  neutrinoVec += dm1Vec;
  neutrinoVec += dm2Vec;
  
//   cout << "after\tneutrinoVec.Mag() = " << neutrinoVec.Mag() << endl;
  
  double trueWmass = 
        TMath::Sqrt(TMath::Power(mother->e(),2)
              - TMath::Power(mother->pz(),2)
              - TMath::Power(mother->py(),2)
              - TMath::Power(mother->px(),2)
                   )*GEV; 
  
  double calculatedTrueWmass = 
        TMath::Sqrt( 2*leptonVec.Mag()*neutrinoVec.Mag() * 
          (1.0 - TMath::Cos( leptonVec.Angle(neutrinoVec) ) ) );
  
  double calculatedTrueMt = TMath::Sqrt(2*leptonVec.Pt()*neutrinoVec.Pt()*
                      (1-TMath::Cos(leptonVec.Phi()-neutrinoVec.Phi())));
 
  /// get MC weights
  if (m_isMC){
    m_LPXKfactorTool->execute();
    m_weightkFactor = m_eventInfo->auxdecor<double>("KfactorWeight");
  }
  
  h_event_crossSectionWeight->Fill(m_weightCrossSection);
  h_event_kFactor->Fill(m_weightkFactor);
  h_event_filterEfficiency->Fill(m_weighfilterEfficiency);
  
  double totalWeight = m_weighfilterEfficiency*m_weightkFactor
                      *m_weightCrossSection;
   
  totalWeight /= m_normalizationFactor;
  
  if (totalWeight==0) totalWeight = 1.0;/// FIXME just for mono-W sample
                                              
  h_event_totalWeight->Fill(totalWeight);
  
  hMu_pt_off->Fill(leptonVec.Pt()*GEV,totalWeight);
  hMu_mt_off->Fill(calculatedTrueMt*GEV,totalWeight);
  hMu_MET_Muons_off->Fill(neutrinoVec.Pt()*GEV,totalWeight);
  hMu_invMass_Muons_off->Fill(trueWmass,totalWeight);
                      
  hMu_invMass_Muons_alternative->Fill(calculatedTrueWmass*GEV,totalWeight);
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











