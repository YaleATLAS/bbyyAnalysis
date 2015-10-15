/**
 * @file   SignalResonanceMassWindow.cxx
 * @Author James Robinson (james.robinson@cern.ch)
 * @date   11th June 2015
 * @brief  Study effective mass window for H->hh->bbyy
 *
 * Package for Hgamma analysis framework to determine optimal mass window
 */

#include "SignalResonanceMassWindow/SignalResonanceMassWindow.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"

#include <EventLoop/Worker.h>
#include <TFile.h>
#include <TTree.h>
#include <boost/format.hpp>

/// this is needed to distribute the algorithm to the workers
ClassImp(SignalResonanceMassWindow)

SignalResonanceMassWindow::SignalResonanceMassWindow(const char *name)
: HgammaAnalysis(name)
, m_output_tree(0)
, m_m_yyjj_0tag(0)
, m_m_yyjj_1tag(0)
, m_m_yyjj_2tag(0)
, m_m_yyjj_mHconstraint_0tag(0)
, m_m_yyjj_mHconstraint_1tag(0)
, m_m_yyjj_mHconstraint_2tag(0)
{
  /// Here you put any code for the base initialization of variables,
  /// e.g. initialize all pointers to 0.  Note that you should only put
  /// the most basic initialization here, since this method will be
  /// called on both the submission and the worker node.  Most of your
  /// initialization code will go into histInitialize() and
  /// initialize().
}


SignalResonanceMassWindow::~SignalResonanceMassWindow()
{
  /// Here you delete any memory you allocated during your analysis.
}


EL::StatusCode SignalResonanceMassWindow::initialize()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  /// Initialise baseclass
  EL::StatusCode sc = HgammaAnalysis::initialize();

  // Set output tree in output file
  TFile *file = wk()->getOutputFile("MxAOD");
  m_output_tree = new TTree("outputTree","outputTree");
  m_output_tree->SetDirectory(file);
  m_output_tree->Branch("m_yyjj_0tag", &m_m_yyjj_0tag);
  m_output_tree->Branch("m_yyjj_1tag", &m_m_yyjj_1tag);
  m_output_tree->Branch("m_yyjj_2tag", &m_m_yyjj_2tag);
  m_output_tree->Branch("m_yyjj_mHconstraint_0tag", &m_m_yyjj_mHconstraint_0tag);
  m_output_tree->Branch("m_yyjj_mHconstraint_1tag", &m_m_yyjj_mHconstraint_1tag);
  m_output_tree->Branch("m_yyjj_mHconstraint_2tag", &m_m_yyjj_mHconstraint_2tag);

  // Read configuration
  m_m_yy_low  = config()->getNum("SignalResonanceMassWindow.MyyWindow.Low",  120.);
  m_m_yy_high = config()->getNum("SignalResonanceMassWindow.MyyWindow.High", 130.);
  m_debug     = config()->getBool("SignalResonanceMassWindow.DebugMessages", false);
  if( m_debug ) { Info( "initialize()", (boost::format( "Setting m_yy window to %s -> %s" ) % m_m_yy_low % m_m_yy_high ).str().c_str() ); }

  // Get handle to event
  m_event = wk()->xaodEvent();
  return sc;
}

EL::StatusCode SignalResonanceMassWindow::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.
  histoStore()->createTH1F("m_yyjj_0tag", 2501, -1.0, 5001.0,";0 tag, #it{m}_{bbyy} [GeV];N_{events}");
  histoStore()->createTH1F("m_yyjj_1tag", 2501, -1.0, 5001.0,";1 tag, #it{m}_{bbyy} [GeV];N_{events}");
  histoStore()->createTH1F("m_yyjj_2tag", 2501, -1.0, 5001.0,";2 tag, #it{m}_{bbyy} [GeV];N_{events}");
  histoStore()->createTH1F("m_yyjj_mHconstraint_0tag", 2501, -1.0, 5001.0,";0 tag, mHconstraint, #it{m}_{bbyy} [GeV];N_{events}");
  histoStore()->createTH1F("m_yyjj_mHconstraint_1tag", 2501, -1.0, 5001.0,";1 tag, mHconstraint, #it{m}_{bbyy} [GeV];N_{events}");
  histoStore()->createTH1F("m_yyjj_mHconstraint_2tag", 2501, -1.0, 5001.0,";2 tag, mHconstraint, #it{m}_{bbyy} [GeV];N_{events}");
}


EL::StatusCode SignalResonanceMassWindow::execute()
{
  /// Here you do everything that needs to be done on every single
  /// events, e.g. read input variables, apply cuts, and fill
  /// histograms and trees.  This is where most of your actual analysis
  /// code will go.

  /// Important to keep this, so that internal tools / event variables
  /// are filled properly.
  EL::StatusCode sc = HgammaAnalysis::execute();

  // sorry for a late reply! So the mass window is for m_bbyy in the resonant search
  // analysis. Additional correction to m_bb is by applying a mass constraint (scaling
  // it to the Higgs mass) which improves the resolution by quite a bit. Take a look
  // at Fig. 29 in the internal note.  Then the mass window is defined as the smallest
  // window containing 95% of the signal events. This has to be done separately for
  // each mass hypothesis and for 2-tag and 1-tag separately since these will be our
  // signal regions. Also it would be great to compare the mass window before and
  // after applying the mass constraint (Fig. 45). And yes, at reco-level and after
  // all the analysis cuts - so tight photons, pT cuts, b-jets, myy etc.  Another
  // question here is also about the background m_bbyy efficiency —> I think they
  // were derived from a di-jet control sample in Run1.

  // Yes, exactly. We're looking for the m_yybb window for events that are already in
  // the m_yy mass window (for now you could just use 120-130 GeV). And we could do
  // this with and without the m_bb scaling. In terms of the background, it's interesting
  // to see the "efficiency" for background to pass these cuts with and without the m_bb
  // mass constraint. For now of course we have simulation. In run 1 we took this
  // from the < 2 tag data - perhaps we can use 0 tag data this time around when we
  // have it and then compare to our simulation estimate.

  // Get event info
  const xAOD::EventInfo* eventInfo(0);
  HG_CHECK( "execute()", m_event->retrieve(eventInfo, "EventInfo") )

  /// Fetch selected photons
  const xAOD::PhotonContainer *photons(0);
  HG_CHECK( "execute()", m_event->retrieve(photons, "HGamPhotons") );

  // Require at least two photons
  if( photons->size() < 2 ) { return sc; }

  // Require photons to be in the mass window
  TLorentzVector yy = photons->at(0)->p4() + photons->at(1)->p4();
  if( (yy.M() / HG::GeV) < m_m_yy_low || (yy.M() / HG::GeV) > m_m_yy_high ) { return sc; }
  if( m_debug ) { Info( "execute()", (boost::format( "Found m_yy (%s) in range %s -> %s" ) % (yy.M() / HG::GeV) % m_m_yy_low % m_m_yy_high ).str().c_str() ); }

  // Fetch selected jet collections
  const xAOD::JetContainer *b_jets_0(0);
  HG_CHECK( "execute()", m_event->retrieve(b_jets_0, "HGamAntiKt4EMTopoJets_AllSelZeroTag") );
  const xAOD::JetContainer *b_jets_1(0);
  HG_CHECK( "execute()", m_event->retrieve(b_jets_1, "HGamAntiKt4EMTopoJets_AllSelOneTag") );
  const xAOD::JetContainer *b_jets_2(0);
  HG_CHECK( "execute()", m_event->retrieve(b_jets_2, "HGamAntiKt4EMTopoJets_AllSelTwoTag") );

  // Get the uncorrected masses
  m_m_yyjj_0tag = eventInfo->auxdata<double>("HGamAntiKt4EMTopoJets_AllSelZeroTag_m_yyjj");
  m_m_yyjj_1tag = eventInfo->auxdata<double>("HGamAntiKt4EMTopoJets_AllSelOneTag_m_yyjj");
  m_m_yyjj_2tag = eventInfo->auxdata<double>("HGamAntiKt4EMTopoJets_AllSelTwoTag_m_yyjj");

  // Set initial values outside the acceptance
  m_m_yyjj_mHconstraint_0tag = -99;
  m_m_yyjj_mHconstraint_1tag = -99;
  m_m_yyjj_mHconstraint_2tag = -99;

  // Force recalculation with Higgs constraint
  if( b_jets_0->size() >= 2 ) {
    TLorentzVector bb = ApplyHiggsMassScaling( b_jets_0->at(0)->p4() + b_jets_0->at(1)->p4() ) ;
    m_m_yyjj_mHconstraint_0tag = (yy + bb).M() / HG::GeV;
  }
  if( b_jets_1->size() >= 2 ) {
    TLorentzVector bb = ApplyHiggsMassScaling( b_jets_1->at(0)->p4() + b_jets_1->at(1)->p4() ) ;
    m_m_yyjj_mHconstraint_1tag = (yy + bb).M() / HG::GeV;
  }
  if( b_jets_2->size() >= 2 ) {
    TLorentzVector bb = ApplyHiggsMassScaling( b_jets_2->at(0)->p4() + b_jets_2->at(1)->p4() ) ;
    m_m_yyjj_mHconstraint_2tag = (yy + bb).M() / HG::GeV;
  }

  if( m_debug ) { Info( "execute()", (boost::format("m_yyjj 0 tag: %s, with bb constraint %s") % m_m_yyjj_0tag % m_m_yyjj_mHconstraint_0tag).str().c_str() ); }
  if( m_debug ) { Info( "execute()", (boost::format("m_yyjj 1 tag: %s, with bb constraint %s") % m_m_yyjj_1tag % m_m_yyjj_mHconstraint_1tag).str().c_str() ); }
  if( m_debug ) { Info( "execute()", (boost::format("m_yyjj 2 tag: %s, with bb constraint %s") % m_m_yyjj_2tag % m_m_yyjj_mHconstraint_2tag).str().c_str() ); }

  // Fill histograms
  histoStore()->fillTH1F("m_yyjj_0tag",m_m_yyjj_0tag);
  histoStore()->fillTH1F("m_yyjj_1tag",m_m_yyjj_1tag);
  histoStore()->fillTH1F("m_yyjj_2tag",m_m_yyjj_2tag);
  histoStore()->fillTH1F("m_yyjj_mHconstraint_0tag",m_m_yyjj_mHconstraint_0tag);
  histoStore()->fillTH1F("m_yyjj_mHconstraint_1tag",m_m_yyjj_mHconstraint_1tag);
  histoStore()->fillTH1F("m_yyjj_mHconstraint_2tag",m_m_yyjj_mHconstraint_2tag);

  // Fill output tree and return
  m_output_tree->Fill();
  return sc;
}


TLorentzVector SignalResonanceMassWindow::ApplyHiggsMassScaling( TLorentzVector input, double m_H_in_GeV ) {
  TLorentzVector output( input );
  if( output.M() > 0.0 ) {
    output *= ( m_H_in_GeV * HG::GeV / output.M() );
  }
  return output;
}
