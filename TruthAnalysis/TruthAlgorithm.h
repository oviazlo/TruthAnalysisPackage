#ifndef TruthAnalysis_TruthAlgorithm_H
#define TruthAnalysis_TruthAlgorithm_H

/// ROOT
#include "TMath.h"
#include <TLorentzVector.h>
#include <TH1.h>

/// EventLoop
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/Algorithm.h>
#include <EventLoop/Algorithm.h>

/// xAOD
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"

/// std c++
#include <iostream>

/// LPXKfactorTool tool
#include "LPXKfactorTool/LPXKfactorTool.h"

/// private
#include <TruthAnalysis/WprimeHist.h>
#include <TruthAnalysis/WprimeSample.h>

/// Helper macro for checking xAOD::TReturnCode return values
#define EL_RETURN_CHECK( CONTEXT, EXP )                     \
   do {                                                     \
      if( ! EXP.isSuccess() ) {                             \
         Error( CONTEXT,                                    \
                XAOD_MESSAGE( "Failed to execute: %s" ),    \
                #EXP );                                     \
         return EL::StatusCode::FAILURE;                    \
      }                                                     \
   } while( false )

using namespace std;
   
class TruthAlgorithm : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  
  /// Magnar histogramming way
  TH1D* hMu_pt_off; //!
  TH1D* hMu_mt_off; //!
  TH1D* hMu_MET_Muons_off; //!
  
  /// EventInfo
  TH1D* h_event_crossSectionWeight; //!
  TH1D* h_event_kFactor; //!
  TH1D* h_event_filterEfficiency; //!
  TH1D* h_event_totalWeight; //!
  
  xAOD::TEvent *m_event;  //!
  xAOD::TStore *m_store;  //!

  int EventNumber; //!
  
  float weighfilterEfficiency; //!
  float weightkFactor; //!
  float weightCrossSection; //!
  
  Double_t GeV = 0.001; //!
  
  #ifndef __CINT__

  LPXKfactorTool* m_LPXKfactorTool; //!
  const xAOD::EventInfo* eventInfo = 0; //!
    
  #endif /// not __CINT__
  
  // this is a standard constructor
  TruthAlgorithm ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(TruthAlgorithm, 1);
};

#endif