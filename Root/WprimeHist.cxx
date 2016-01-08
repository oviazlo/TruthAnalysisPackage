// WprimeHist.cxx

#include <TruthAnalysis/WprimeHist.h>
#include "TruthAnalysis/WprimeSample.h"
#include <iostream>
#include <cmath>
#include <TH1F.h>
#include <TH1D.h>
#include <TDirectory.h>

std::string WprimeHist_lepton = "l";

//**********************************************************************

TH1* WprimeHist::create(Name name, Name title, int xbin, double xmin, 
                        double xmax) {
  TH1* ph = new TH1D(name.c_str(), title.c_str(), xbin, xmin, xmax);
  ph->Sumw2();
  return ph;
}

TH1* WprimeHist::create(Name name, Name title, int xbin, double* bins) {
  TH1* ph = new TH1D(name.c_str(), title.c_str(), xbin, bins);
  ph->Sumw2();
  return ph;
}

TH1* WprimeHist::createlog(Name name, Name title, int xbin, double xmin, 
                           double xmax) {
  if ( xmin <= 0 ) return 0;
  if ( xmax <= xmin ) return 0;
  if ( xbin <= 0 ) return 0;
  double* bins = new double[xbin+1];
  double fac = pow(xmax/xmin, 1.0/xbin);
  double x = xmin;
  for ( int ibin=0; ibin<=xbin; ++ibin ) {
    bins[ibin] = x;
    x *= fac;
  }
  TH1* ph = new TH1D(name.c_str(), title.c_str(), xbin, bins);
  ph->Sumw2();
  return ph;
}

TH1* WprimeHist::standard(Name name, Name hpname, Name sel, Name hptitle) {
  const int ptbins = 120;
  const double ptmin = 10;
  const double ptmax = 10000;
  const double pi = acos(-1.0);
  Name hname = hpname + "_" + name;
  if ( sel.size() ) hname += "_" + sel;
  TH1* ph = dynamic_cast<TH1*>(gDirectory->Get(hname.c_str()));
  if ( ph != 0 ) return ph;
  Name htitle = hptitle;
  Name ylab = "; Events";
  if ( name == "eta" ) {
    htitle += " #eta; #eta" + ylab;
    return create(hname, htitle, 60, -3.0, 3.0);
  } else if ( name == "phi" ) {
    htitle += " #varphi; #varphi" + ylab;
    return create(hname, htitle, 64, -pi, pi);
  } else if ( name == "pt" ) {
    htitle += " p_{T}; p_{T}^{" + lepton() + "} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, ptmin, ptmax);
  } else if ( name == "wpt" ) {
    htitle += " W p_{T}; W p_{T} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, ptmin, ptmax*1.5);
  } else if ( name == "ptLepMet" ) {
    htitle += " P_{T}^{l,MET}; P_{T}^{l,MET} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, ptmin, ptmax*1.5);
  } else if ( name == "ptLepCaloMet" ) {
    htitle += " P_{T}^{l,Calo MET}; P_{T}^{l,Calo MET} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, ptmin, ptmax*1.5);
  } else if ( name == "ht" ) {
    htitle += " H_{T}; H_{T} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, ptmin, ptmax*1.5);
  } else if ( name == "met" ) {
    htitle += " E_{T}^{miss}; E_{T}^{miss} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, ptmin, ptmax);
  } else if ( name == "sumet" ) {
    htitle += " #Sum E_{T}; #Sum E_{T} [GeV]" + ylab;
    return createlog(hname, htitle, 2*ptbins, ptmin, 2*ptmax);
  } else if ( name == "met_x" ) {
    htitle += " E_{x}^{miss}; E_{x}^{miss} [GeV]" + ylab;
    return create(hname, htitle, 200, -400.0, 400.0);
  } else if ( name == "met_y" ) {
    htitle += " E_{y}^{miss}; E_{y}^{miss} [GeV]" + ylab;
    return create(hname, htitle, 200, -400.0, 400.0);
  } else if ( name == "mt" ) {
    htitle += " m_{T}; m_{T} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, 2*ptmin, 2*ptmax);
  } else if ( name == "mt_cumul" ) {
    htitle += " m_{T} Cumulative; m_{T} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, 2*ptmin, 2*ptmax);
  } else if ( name == "cmet" ) {
    htitle += " calo E_{T,calo}^{miss}; E_{T,calo}^{miss} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, ptmin, ptmax);
  } else if ( name == "cmet_x" ) {
    htitle += " Calo E_{x}^{miss}; Calo E_{x}^{miss} [GeV]" + ylab;
    return create(hname, htitle, 200, -400.0, 400.0);
  } else if ( name == "cmet_y" ) {
    htitle += " Calo E_{y}^{miss}; Calo E_{y}^{miss} [GeV]" + ylab;
    return create(hname, htitle, 200, -400.0, 400.0);
  } else if ( name == "etloss" ) {
    htitle += " E_{Tloss}; E_{Tloss} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, ptmin, ptmax);
  } else if ( name == "evtwt" ) {
    htitle += " event weight; Event weight" + ylab;
    return create(hname, htitle, 200, -10, 10);
  } else if ( name == "evtwt_xSec" ) {
    htitle += " xSec event weight; xSec Event weight" + ylab;
    return create(hname, htitle, 200, -10, 10);
  } else if ( name == "evtwt_filterEff" ) {
    htitle += " filterEff event weight; filterEff Event weight" + ylab;
    return create(hname, htitle, 200, -10, 10);
  } else if ( name == "evtwt_kFactor" ) {
    htitle += " kFactor event weight; kFactor Event weight" + ylab;
    return create(hname, htitle, 200, -10, 10);
  } else if ( name == "mtlimit" ) {
    htitle += " m_{T}; m_{T} [GeV]" + ylab;
    return create(hname, htitle, 40, 0, 5000);
  } else if ( name == "nvtx" ) {
    htitle += " vertex multiplicity; N_{vtx}" + ylab;
    return create(hname, htitle, 40, 0, 40);
  } else if ( name == "njet" ) {
    htitle += " jet multiplicity; N_{jet}" + ylab;
    return create(hname, htitle, 20, 0, 20);
  } else if ( name == "jet1pt" ) {
    htitle += " jet 1 p_{T}; p_{T}^{jet} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, ptmin, ptmax);
  } else if ( name == "lumiblock" ) {
    htitle += " lumiblock; Lumiblock" + ylab;
    return create(hname, htitle, 2000, 0, 2000);
  } else if ( name == "eletisol" ) {
    htitle += " electron isolation; #Sigma E_{T}^{cell};N_{e}/GeV";
    return create(hname, htitle, 120, -20, 100);
  } else if ( name == "muptisol" ) {
    htitle += " muon isolation; #Sigma p_{T}^{trk}/p_{T}^{#mu};N_{#mu}/0.01";
    return create(hname, htitle, 100, 0, 1.0);
  } else if ( name == "elmetoverpt" ) {
    htitle += " MEt/pT; E_{T}^{miss}/p_{T}^{e}" + ylab;
    return createlog(hname, htitle, 100, 0.1, 10.0);
  } else if ( name == "mumetoverpt" ) {
    htitle += " MEt/pT; E_{T}^{miss}/p_{T}^{e}" + ylab;
    return createlog(hname, htitle, 150, 0.1, 15.0);
  } else if ( name == "runnumber" ) {
    htitle += " Run Number; Run Number" + ylab;
    return create(hname, htitle, 9600, 200800, 210400);
  } else if ( name == "cumulative" ) {
    htitle += " Cumulative m_{T}; m_{T} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, 2*ptmin, 2*ptmax);
  // Muon specific.
  } else if ( name == "qpdiff" ) {
    htitle += " MS - ID q/p; (q/p)_{MS}-(q/p)_{ID} [TeV^{-1}]" + ylab;
    return create(hname, htitle, 100, -5.0, 5.0);
  } else if ( name == "qpdifferr" ) {
    htitle += " MS - ID q/p error; #sigma((q/p)_{MS}-(q/p)_{ID}) [TeV^{-1}]" 
              + ylab;
    return create(hname, htitle, 50, 0.0, 5.0);
  } else if ( name == "qpdiffsig" ) {
    htitle += " MS - ID q/p sgnificance; ((q/p)_{MS}-(q/p)_{ID})/#sigma((q/p)"
              "_{MS}-(q/p)_{ID})" + ylab;
    return create(hname, htitle, 200, -10, 10);
  // Filled from generator level information
  } else if ( name == "ptgen" ) {
    htitle += " generated p_{T}; p_{T}^{l} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, ptmin, ptmax);
  } else if ( name == "metgen" ) {
    htitle += " generated E_{T}^{miss}; p_{T}^{#nu} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, ptmin, ptmax);
  } else if ( name == "mtgen" ) {
    htitle += " generated m_{T}; m_{T} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, 2*ptmin, 2*ptmax);
  } else if ( name == "mgen" ) {
    htitle += " generated m_{inv}; m_{inv} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, 2*ptmin, 2*ptmax);
  } else if ( name == "bosonptgen" ) {
    htitle += " generated Boson p_{T}; Boson p_{T}^{l} [GeV]" + ylab;
    return createlog(hname, htitle, ptbins, ptmin, ptmax);
  } else return 0;
}

TH1* WprimeHist::sample(Name name, Name samplename, Name sel) {
  WprimeSample sam(samplename);
  if ( sam.name() == "" ) return 0;
  return standard(name, sam.name(), sel, sam.title());
}
  
WprimeHist::Name WprimeHist::lepton(Name slepton) {
  Name oldval = WprimeHist_lepton;
  if ( slepton.size() ) {
    if ( slepton == "e" ) WprimeHist_lepton = "e";
    if ( slepton == "mu" ) WprimeHist_lepton = "#mu";
    if ( slepton == "#mu" ) WprimeHist_lepton = "#mu";
  }
  return oldval;
}
