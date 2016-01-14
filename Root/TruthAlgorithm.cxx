#include <TruthAnalysis/TruthAlgorithm.h>

// this is needed to distribute the algorithm to the workers
ClassImp(TruthAlgorithm)



TruthAlgorithm :: TruthAlgorithm ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



EL::StatusCode TruthAlgorithm :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  
  /// let's initialize the algorithm to use the xAODRootAccess package
  job.useXAOD ();
  xAOD::Init( "TruthAlgorithm" ).ignore(); /// call before opening first file
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthAlgorithm :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  
  hMu_pt_off = (TH1D*)WprimeHist::standard("pt","h","","");
  wk()->addOutput(hMu_pt_off); 
  
  hMu_mt_off = (TH1D*)WprimeHist::standard("mt","h","","");
  wk()->addOutput(hMu_mt_off);

  hMu_MET_Muons_off = (TH1D*)WprimeHist::standard("met","h","","");
  wk()->addOutput(hMu_MET_Muons_off); 
  
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



EL::StatusCode TruthAlgorithm :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthAlgorithm :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthAlgorithm :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
  
  m_event = wk()->xaodEvent();
  m_store = new xAOD::TStore();
  
  /// as a check, let's see the number of events in our xAOD
  /// print long long int
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); 

  /// Event Info
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", 
          "Failed to retrieve event info collection in initialise. Exiting." );
    return EL::StatusCode::FAILURE;
  }  
  
  /// fill the branches of our trees
  bool isData = true;
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    isData = false;
  }
  cout << "isData flag = " << isData <<endl;
  
  weightkFactor = 1.0;
  weighfilterEfficiency = 1.0;
  weightCrossSection = 1.0;
  
  if (!isData){
    m_LPXKfactorTool = new LPXKfactorTool("LPXKfactorTool");
    EL_RETURN_CHECK("m_LPXKfactorTool_isMC15",
                    m_LPXKfactorTool->setProperty("isMC15", true)); 
    EL_RETURN_CHECK("m_LPXKfactorTool_applyEWCorr",
                    m_LPXKfactorTool->setProperty("applyEWCorr", true)); 
    EL_RETURN_CHECK("m_LPXKfactorTool_applyPICorr",
                    m_LPXKfactorTool->setProperty("applyPICorr", true)); 
    
    EL_RETURN_CHECK( "m_LPXKfactorTool initialize",m_LPXKfactorTool->initialize());
  }
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthAlgorithm :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  
  ///----------------------------
  /// Event information
  ///--------------------------- 

  EL_RETURN_CHECK("retrieve EventInfo",
                  m_event->retrieve( eventInfo, "EventInfo"));
  EventNumber = eventInfo->eventNumber();  

  bool isMC = false;
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    isMC = true;
  }
  
  
  /*
  /// Create truth vertice container
  const xAOD::TruthVertexContainer* truthVertices = 0;
  EL_RETURN_CHECK("retrieve TruthVertices", 
                  m_event->retrieve( truthVertices, "TruthVertices" ));
  
  double highestMt = 0.0;
  double invMass = 0.0;
  double leptonPt = 0.0;
  
  /// Start iterating over truth container
  xAOD::TruthVertexContainer::const_iterator truthV_itr; 
  for (truthV_itr = truthVertices->begin(); 
        truthV_itr != truthVertices->end(); ++truthV_itr )
  {
   
    for (unsigned int iIn=0; iIn < (*truthV_itr)->nIncomingParticles(); 
           iIn++)
    {
      if (TMath::Abs((*truthV_itr)->incomingParticle(iIn)->pdgId()) == 24) { 
         for (unsigned int iOut=0; iOut < (*truthV_itr)->nOutgoingParticles();
               iOut++) {
            if (TMath::Abs((*truthV_itr)->outgoingParticle(iOut)->pdgId()) 
              == 13){
              TVector3 leptonVec = ((*truthV_itr)->outgoingParticle(iOut)->
                                    p4()).Vect();
              leptonPt = leptonVec.Pt() * GeV;
              cout << "leptonPt = " << leptonPt << endl;
            }
         }
            
      }
        
    }
    
  }
  
  */
  const xAOD::TruthVertexContainer* truthVertices = 0;
  /// retrieve arguments: container type, container key
  if ( m_event->retrieve( truthVertices, "truthVertices" ).isSuccess() ){
  EL_RETURN_CHECK("retrieve truthVertices", 
                  m_event->retrieve( truthVertices, "truthVertices" ));
  
  
  unsigned int pdg = 13; /// muons;
  bool foundmuon = false; 

  /// Start iterating over truth container
  xAOD::TruthVertexContainer::const_iterator truthV_itr; 
  const xAOD::TruthParticle* particle;
  for (truthV_itr = truthVertices->begin(); truthV_itr != truthVertices->end();
        ++truthV_itr ) {
    int nout = 0;
    for (unsigned int iOut=0; iOut < (*truthV_itr)->nOutgoingParticles(); 
          iOut++) {
      nout +=1;
      particle = (*truthV_itr)->outgoingParticle(iOut);
      if (particle && TMath::Abs(particle->pdgId())==pdg) {
        particleid = particle->pdgId();
        status = particle->status();
        cout << "found ooutgoing particle: " << nout<< "pdg " << 
         particle->pdgId()<< ", " << particle->pt()<<endl;
        for (unsigned int iIn=0; iIn < (*truthV_itr)->nIncomingParticles(); 
             iIn++) {
          cout << "mother " << 
           TMath::Abs((*truthV_itr)->incomingParticle(iIn)->pdgId())<< endl;
          motherid = (*truthV_itr)->incomingParticle(iIn)->pdgId();
        }
 

      }
      
    }
  }
  
  
  /*
  
  
  TVector3 MET(0,0,0);
  xAOD::TruthParticleContainer::const_iterator truth; 
  for (truth = TruthParticles->begin(); 
         truth != TruthParticles->end(); ++truth ){
        
    Int_t PDG_ID = abs((*truth)->pdgId());
    Int_t Status = (*truth)->status();

//     if (Status != 3) continue; /// hard process particles only
    if (PDG_ID!=12 && PDG_ID!=14 && PDG_ID!=16) continue; /// neutrinos

    MET += ((*truth)->p4()).Vect();
  }
  
  double highestMt = 0.0;
  double invMass = 0.0;
  double leptonPt = 0.0;
  for (truth = TruthParticles->begin(); 
         truth != TruthParticles->end(); ++truth ){

    Int_t PDG_ID = abs((*truth)->pdgId());
    Int_t Status = (*truth)->status();

//     if (Status != 3) continue; /// hard process particles only
    if (PDG_ID!=11 && PDG_ID!=13 && PDG_ID!=15) continue; /// leptons

    TVector3 leptonVec = ((*truth)->p4()).Vect();

    double Mt_MET = sqrt( 2*leptonVec.Pt()*MET.Pt() * 
    (1.0 - TMath::Cos( leptonVec.Phi() - MET.Phi() )) )*GeV;      

    double tmpInvMass = sqrt( 2*leptonVec.Mag()*MET.Mag() * 
    (1.0 - TMath::Cos( leptonVec.Phi() - MET.Phi() )    ) )*GeV;
    
    if (Mt_MET>highestMt){
      leptonPt = leptonVec.Pt() * GeV;
      invMass = tmpInvMass;
      highestMt = Mt_MET;
    }
    
  }*/
  
  /// get MC weights
  if (isMC){
    m_LPXKfactorTool->execute();
    weightkFactor = eventInfo->auxdecor<double>("KfactorWeight");
    weighfilterEfficiency = m_LPXKfactorTool->getMCFilterEfficiency();
    ///TODO make proper implementation
    weightCrossSection = m_LPXKfactorTool->getMCCrossSection(); 
  }
  
  h_event_crossSectionWeight->Fill(weightCrossSection);
  h_event_kFactor->Fill(weightkFactor);
  h_event_filterEfficiency->Fill(weighfilterEfficiency);
  
  double totalWeight = weighfilterEfficiency*weightkFactor*weightCrossSection;
  
  h_event_totalWeight->Fill(totalWeight);
  
  hMu_pt_off->Fill(leptonPt,totalWeight);
  hMu_mt_off->Fill(highestMt,totalWeight);
  hMu_MET_Muons_off->Fill(MET.Pt()*GeV,totalWeight);
  
//   cout << "leptonPt = " << leptonPt << endl;
//   cout << "highestMt = " << highestMt << endl;
//   cout << "MET = " << MET.Pt()*GeV << endl;
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthAlgorithm :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthAlgorithm :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.
  
  if (m_LPXKfactorTool){
    delete m_LPXKfactorTool;
    m_LPXKfactorTool = 0;
  }
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TruthAlgorithm :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}
