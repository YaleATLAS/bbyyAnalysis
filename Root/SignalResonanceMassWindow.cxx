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
, m_event_weight(0)
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
  const auto sc = HgammaAnalysis::initialize();
  if( sc != EL::StatusCode::SUCCESS ) { return sc; }

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
  m_output_tree->Branch( "event_weight",         &m_event_weight         );
  Info( "initialize()", "Initialising output tree" );

  // Read configuration
  m_m_yy_low  = config()->getNum("SignalResonanceMassWindow.MyyWindow.Low",  121.9);
  m_m_yy_high = config()->getNum("SignalResonanceMassWindow.MyyWindow.High", 128.1);
  Info( "initialize()", (boost::format( "Setting m_yy window to %s -> %s GeV" ) % m_m_yy_low % m_m_yy_high ).str().c_str() );

  // Get handle to event
  m_event = wk()->xaodEvent();
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode SignalResonanceMassWindow::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.
  for( const auto& constraint : { "unscaled", "mHscaled" } ) {
    for( const auto& tag : { "0tag", "1tag", "2tag" } ) {
      histoStore()->createTH1F( (boost::format("m_yyjj_%s_%s") % constraint % tag).str(), 800, 100, 900,";#it{m}_{yyjj} [GeV];N_{events}");
      histoStore()->createTH1F( (boost::format("m_jj_%s_%s") % constraint % tag).str(), 50, 0, 200,";#it{m}_{jj} [GeV];N_{events}");
      histoStore()->createTH1F( (boost::format("m_yy_%s_%s") % constraint % tag).str(), 55, 105, 160,";#it{m}_{yy} [GeV];N_{events}");
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

  // /// Not strictly needed if all information is correctly in the MxAOD
  // const auto sc = HgammaAnalysis::execute();
  // if( sc != EL::StatusCode::SUCCESS ) { return sc; }

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
  EL_CHECK( "execute()", m_event->retrieve(HGammaEventInfo, "HH2yybbEventInfo") )

  // Require photons to be in the mass window
  double m_yy = HGammaEventInfo->auxdata<float>("m_yy");
  if( (m_yy / HG::GeV) < m_m_yy_low || (m_yy / HG::GeV) > m_m_yy_high ) { return EL::StatusCode::SUCCESS; }

  // Retrieve b-tagging category: -1:no jets; 0:two light jets; 1:1 light/1 b-jet; 2:2 b-jets
  int bTagCategory = HGammaEventInfo->auxdata<int>("bTagCategory");
  if( bTagCategory < 0 ) { return EL::StatusCode::SUCCESS; }

  // Get the uncorrected masses
  double m_yyjj_unscaled = HGammaEventInfo->auxdata<float>("m_yyjj") / HG::GeV;
  m_m_yyjj_unscaled_0tag = bTagCategory == 0 ? m_yyjj_unscaled : -99;
  m_m_yyjj_unscaled_1tag = bTagCategory == 1 ? m_yyjj_unscaled : -99;
  m_m_yyjj_unscaled_2tag = bTagCategory == 2 ? m_yyjj_unscaled : -99;

  // Get Higgs-scaled masses
  double m_yyjj_mHscaled = HGammaEventInfo->auxdata<float>("m_yyjj_constrnd") / HG::GeV;
  m_m_yyjj_mHscaled_0tag = bTagCategory == 0 ? m_yyjj_mHscaled : -99;
  m_m_yyjj_mHscaled_1tag = bTagCategory == 1 ? m_yyjj_mHscaled : -99;
  m_m_yyjj_mHscaled_2tag = bTagCategory == 2 ? m_yyjj_mHscaled : -99;

  // Overall event weight: pileup * PVz * photon scale factors * xsLumi
  m_event_weight = HGammaEventInfo->auxdata<float>("weightFinal");

  // Fill 4-body mass histograms
  histoStore()->fillTH1F( "m_yyjj_unscaled_0tag", m_m_yyjj_unscaled_0tag, m_event_weight );
  histoStore()->fillTH1F( "m_yyjj_unscaled_1tag", m_m_yyjj_unscaled_1tag, m_event_weight );
  histoStore()->fillTH1F( "m_yyjj_unscaled_2tag", m_m_yyjj_unscaled_2tag, m_event_weight );
  histoStore()->fillTH1F( "m_yyjj_mHscaled_0tag", m_m_yyjj_mHscaled_0tag, m_event_weight );
  histoStore()->fillTH1F( "m_yyjj_mHscaled_1tag", m_m_yyjj_mHscaled_1tag, m_event_weight );
  histoStore()->fillTH1F( "m_yyjj_mHscaled_2tag", m_m_yyjj_mHscaled_2tag, m_event_weight );

  // Fill auxiliary histograms
  double m_jj_unscaled = HGammaEventInfo->auxdata<float>("m_jj") / HG::GeV;
  double m_jj_mHscaled = HGammaEventInfo->auxdata<float>("m_yyjj_constrnd") / HG::GeV;
  histoStore()->fillTH1F( (boost::format("m_yy_unscaled_%stag") % bTagCategory).str(), m_yy, m_event_weight );
  histoStore()->fillTH1F( (boost::format("m_yy_mHscaled_%stag") % bTagCategory).str(), m_yy, m_event_weight );
  histoStore()->fillTH1F( (boost::format("m_jj_unscaled_%stag") % bTagCategory).str(), m_jj_unscaled, m_event_weight );
  histoStore()->fillTH1F( (boost::format("m_jj_mHscaled_%stag") % bTagCategory).str(), m_jj_mHscaled, m_event_weight );

  // Fill output tree and return
  m_output_tree->Fill();
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode SignalResonanceMassWindow::finalize() {
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  const auto sc = HgammaAnalysis::finalize();
  if( sc != EL::StatusCode::SUCCESS ) { return sc; }

  Info( "finalise()", (boost::format("Output tree has %s entries")%m_output_tree->GetEntries()).str().c_str() );
  return EL::StatusCode::SUCCESS;
}


TLorentzVector SignalResonanceMassWindow::ApplyHiggsMassScaling( TLorentzVector input, double m_H_in_GeV ) {
  TLorentzVector output( input );
  if( output.M() > 0.0 ) {
    output *= ( m_H_in_GeV * HG::GeV / output.M() );
  }
  return output;
}
