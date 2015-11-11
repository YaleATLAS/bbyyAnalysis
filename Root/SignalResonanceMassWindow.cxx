/**
 * @file   SignalResonanceMassWindow.cxx
 * @Author James Robinson (james.robinson@cern.ch)
 * @date   11th June 2015
 * @brief  Study effective mass window for H->hh->bbyy
 *
 * Package for Hgamma analysis framework to determine optimal mass window
 */

#include "bbyyAnalysis/SignalResonanceMassWindow.h"
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
, m_m_yyjj_unscaled_0tag(0)
, m_m_yyjj_unscaled_1tag(0)
, m_m_yyjj_unscaled_2tag(0)
, m_m_yyjj_mHscaled_0tag(0)
, m_m_yyjj_mHscaled_1tag(0)
, m_m_yyjj_mHscaled_2tag(0)
, m_weight_pileup(0)
, m_weight_xslumi(0)
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
  m_output_tree->Branch( "m_yyjj_unscaled_0tag", &m_m_yyjj_unscaled_0tag );
  m_output_tree->Branch( "m_yyjj_unscaled_1tag", &m_m_yyjj_unscaled_1tag );
  m_output_tree->Branch( "m_yyjj_unscaled_2tag", &m_m_yyjj_unscaled_2tag );
  m_output_tree->Branch( "m_yyjj_mHscaled_0tag", &m_m_yyjj_mHscaled_0tag );
  m_output_tree->Branch( "m_yyjj_mHscaled_1tag", &m_m_yyjj_mHscaled_1tag );
  m_output_tree->Branch( "m_yyjj_mHscaled_2tag", &m_m_yyjj_mHscaled_2tag );
  m_output_tree->Branch( "weight_pileup",        &m_weight_pileup        );
  m_output_tree->Branch( "weight_xslumi",        &m_weight_xslumi        );
  Info( "initialize()", "Initialising output tree" );

  // Read configuration
  m_m_yy_low  = config()->getNum("SignalResonanceMassWindow.MyyWindow.Low",  120.);
  m_m_yy_high = config()->getNum("SignalResonanceMassWindow.MyyWindow.High", 130.);
  Info( "initialize()", (boost::format( "Setting m_yy window to %s -> %s GeV" ) % m_m_yy_low % m_m_yy_high ).str().c_str() );

  // Get handle to event
  m_event = wk()->xaodEvent();
  return sc;
}

EL::StatusCode SignalResonanceMassWindow::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.
  for( const auto& constraint : { "unscaled", "mHscaled" } ) {
    for( const auto& tag : { "0tag", "1tag", "2tag" } ) {
      for( const auto& PRW : { "PRW", "noPRW" } ) {
        histoStore()->createTH1F( (boost::format("m_yyjj_%s_%s_%s") % constraint % tag % PRW).str(), 2501, -1.0, 5001.0,";#it{m}_{bbyy} [GeV];N_{events}");
      }
    }
  }
  Info( "createOutput()", "Initialising output histograms" );
  return EL::StatusCode::SUCCESS;
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

  /// Get event info
  const xAOD::EventInfo* HGammaEventInfo(0);
  HG_CHECK( "execute()", m_event->retrieve(HGammaEventInfo, "HH2yybbEventInfo") )

  /// Fetch selected objects
  const xAOD::PhotonContainer *selected_photons(0);
  HG_CHECK( "execute()", m_event->retrieve(selected_photons, "HH2yybbPhotons") );
  const xAOD::JetContainer *selected_jets(0);
  HG_CHECK( "execute()", m_event->retrieve(selected_jets, "HH2yybbAntiKt4EMTopoJets_AllSelAnyTag") );

  // Require at least two photons
  if( selected_photons->size() < 2 ) { return sc; }

  // Require photons to be in the mass window
  TLorentzVector yy = selected_photons->at(0)->p4() + selected_photons->at(1)->p4();
  if( (yy.M() / HG::GeV) < m_m_yy_low || (yy.M() / HG::GeV) > m_m_yy_high ) { return sc; }

  // Retrieve b-tagging category: -1:no jets; 0:two light jets; 1:1 light/1 b-jet; 2:2 b-jets
  int bTagCategory = HGammaEventInfo->auxdata<int>("bTagCategory");
  if( bTagCategory < 0 ) { return sc; }

  // Get the uncorrected masses
  double m_yyjj_unscaled = HGammaEventInfo->auxdata<float>("m_yyjj") / HG::GeV;
  m_m_yyjj_unscaled_0tag = bTagCategory == 0 ? m_yyjj_unscaled : -99;
  m_m_yyjj_unscaled_1tag = bTagCategory == 1 ? m_yyjj_unscaled : -99;
  m_m_yyjj_unscaled_2tag = bTagCategory == 2 ? m_yyjj_unscaled : -99;

  // Construct Higgs-scaled masses
  TLorentzVector jj = ApplyHiggsMassScaling( selected_jets->at(0)->p4() + selected_jets->at(1)->p4() );
  double m_yyjj_mHscaled = (yy + jj).M() / HG::GeV;

  // Set initial values outside the acceptance
  m_m_yyjj_mHscaled_0tag = bTagCategory == 0 ? m_yyjj_mHscaled : -99;
  m_m_yyjj_mHscaled_1tag = bTagCategory == 1 ? m_yyjj_mHscaled : -99;
  m_m_yyjj_mHscaled_2tag = bTagCategory == 2 ? m_yyjj_mHscaled : -99;

  // Overall event weight is pileup weight * xs weight
  m_weight_pileup = HGammaEventInfo->auxdata<float>("weight");
  m_weight_xslumi = HGammaEventInfo->auxdata<float>("weightXsecLumi");

  // Fill histograms - no PRW
  histoStore()->fillTH1F( "m_yyjj_unscaled_0tag_noPRW", m_m_yyjj_unscaled_0tag, m_weight_xslumi );
  histoStore()->fillTH1F( "m_yyjj_unscaled_1tag_noPRW", m_m_yyjj_unscaled_1tag, m_weight_xslumi );
  histoStore()->fillTH1F( "m_yyjj_unscaled_2tag_noPRW", m_m_yyjj_unscaled_2tag, m_weight_xslumi );
  histoStore()->fillTH1F( "m_yyjj_mHscaled_0tag_noPRW", m_m_yyjj_mHscaled_0tag, m_weight_xslumi );
  histoStore()->fillTH1F( "m_yyjj_mHscaled_1tag_noPRW", m_m_yyjj_mHscaled_1tag, m_weight_xslumi );
  histoStore()->fillTH1F( "m_yyjj_mHscaled_2tag_noPRW", m_m_yyjj_mHscaled_2tag, m_weight_xslumi );

  // Fill histograms - with PRW
  histoStore()->fillTH1F( "m_yyjj_unscaled_0tag_PRW", m_m_yyjj_unscaled_0tag, m_weight_xslumi * m_weight_pileup );
  histoStore()->fillTH1F( "m_yyjj_unscaled_1tag_PRW", m_m_yyjj_unscaled_1tag, m_weight_xslumi * m_weight_pileup );
  histoStore()->fillTH1F( "m_yyjj_unscaled_2tag_PRW", m_m_yyjj_unscaled_2tag, m_weight_xslumi * m_weight_pileup );
  histoStore()->fillTH1F( "m_yyjj_mHscaled_0tag_PRW", m_m_yyjj_mHscaled_0tag, m_weight_xslumi * m_weight_pileup );
  histoStore()->fillTH1F( "m_yyjj_mHscaled_1tag_PRW", m_m_yyjj_mHscaled_1tag, m_weight_xslumi * m_weight_pileup );
  histoStore()->fillTH1F( "m_yyjj_mHscaled_2tag_PRW", m_m_yyjj_mHscaled_2tag, m_weight_xslumi * m_weight_pileup );

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
