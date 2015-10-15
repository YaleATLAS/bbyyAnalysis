/**
 * @file   ResonanceWidthStudies.cxx
 * @Author James Robinson (james.robinson@cern.ch)
 * @date   22nd May 2015
 * @brief  Study interference effects H->hh->bbyy
 *
 * Package for Hgamma analysis framework to study interference effects
 */

#include "bbyyAnalysis/ResonanceWidthStudies.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/TruthUtils.h"
#include <boost/format.hpp>
#include <TRandom3.h>

// this is needed to distribute the algorithm to the workers
ClassImp(ResonanceWidthStudies)

// bool pTdown ( const xAOD::Jet_v1* i, const xAOD::Jet_v1* j ) { return ( i->pt() > j->pt() ); }
bool pTdown ( const xAOD::IParticle* i, const xAOD::IParticle* j ) { return ( i->pt() > j->pt() ); }

ResonanceWidthStudies::ResonanceWidthStudies(const char *name)
: HgammaAnalysis(name)
, m_categories( { "truth", "reco" } )
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



ResonanceWidthStudies::~ResonanceWidthStudies()
{
  // Here you delete any memory you allocated during your analysis.
}



EL::StatusCode ResonanceWidthStudies::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.
  for( const auto& category : m_categories ) {
    histoStore()->createTH1F(category+"_m_yy", 100, 119.95, 129.95,";"+category+" #it{m}_{#gamma#gamma} [GeV]");
    histoStore()->createTH1F(category+"_pT_yy", 400, -0.5, 399.5,";"+category+" #it{p}_{T#gamma#gamma} [GeV]");
    histoStore()->createTH1F(category+"_y_yy", 100, -5.0, 5.0,";"+category+" #it{y}_{#gamma#gamma} ");
    histoStore()->createTH1F(category+"_m_bb", 100, 119.95, 129.95,";"+category+" #it{m}_{bb} [GeV]");
    histoStore()->createTH1F(category+"_pT_bb", 400, -0.5, 399.5,";"+category+" #it{p}_{Tbb} [GeV]");
    histoStore()->createTH1F(category+"_y_bb", 100, -5.0, 5.0,";"+category+" #it{y}_{bb} ");
    histoStore()->createTH1F(category+"_m_qq", 100, 119.95, 129.95,";"+category+" #it{m}_{qq} [GeV]");
    histoStore()->createTH1F(category+"_pT_qq", 400, -0.5, 399.5,";"+category+" #it{p}_{Tqq} [GeV]");
    histoStore()->createTH1F(category+"_y_qq", 100, -5.0, 5.0,";"+category+" #it{y}_{qq} ");
    histoStore()->createTH1F(category+"_m_bbyy", 5000, -0.5, 4999.5,";"+category+" #it{m}_{bbyy} [GeV]");
    histoStore()->createTH2F(category+"_m_bb_vs_m_yy", 100, 0, 1000, 50, 110, 140,""+category+" #it{m}_{bb} [GeV];truth #it{m}_{#gamma#gamma} [GeV]");
    histoStore()->createTH1F(category+"_m_qqyy", 5000, -0.5, 4999.5,";"+category+" #it{m}_{qqyy} [GeV]");
    histoStore()->createTH2F(category+"_m_qq_vs_m_yy", 100, 0, 1000, 50, 110, 140,""+category+" #it{m}_{qq} [GeV];truth #it{m}_{#gamma#gamma} [GeV]");
  }
  histoStore()->createTH1F("event_weight", 100, -50, 50,";event weight");
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ResonanceWidthStudies::execute()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  HgammaAnalysis::execute();
  static int iteration = 0;
  iteration++;

  float weight = eventHandler()->mcWeight();

  // Retrieve the truth particles
  const xAOD::TruthParticleContainer *truthPtcls = 0;
  if ( event()->retrieve( truthPtcls, "TruthParticles" ).isFailure() ) {
    Error("execute()", "Failed to retrieve TruthParticle container." );
    return EL::StatusCode::FAILURE;
  }

  // Retrieve the truth jets
  const xAOD::JetContainer* ptrTruthJets = 0;
  if ( event()->retrieve( ptrTruthJets,"AntiKt4TruthJets" ).isFailure() ) {
    Error("execute()", "Failed to retrieve AntiKt4TruthJets." );
    return EL::StatusCode::FAILURE;
  }

  // Apply selection - photons
  HG::TruthPtcls photons_selected_truth(SG::VIEW_ELEMENTS);
  for( auto photon : HG::getGoodTruthPhotons(truthPtcls) ) {
    if( photon->pt() > 25*HG::GeV ) { photons_selected_truth.push_back(photon); }
  }
  std::sort( photons_selected_truth.begin(), photons_selected_truth.end(), pTdown );

  // Apply selection - jets
  ConstDataVector<xAOD::JetContainer> b_jets_selected_truth(SG::VIEW_ELEMENTS);
  for( auto truth_jet : *ptrTruthJets ){
    if( truth_jet->pt() <= 20*HG::GeV ) { continue; }
    if( truth_jet->auxdata<int>("HadronConeExclTruthLabelID") == 5 ) {
      b_jets_selected_truth.push_back(truth_jet);
      // if( HG::minDRrap( truth_jet, photons_selected_truth ) < 0.4 ) { Info("execute()", "Too close to photon: %lf", HG::minDRrap(truth_jet, photons_selected_truth) ); }
      // if( HG::minDRrap( truth_jet, electrons_selected_truth ) < 0.4 ) { Info("execute()", "Too close to electron: %lf", HG::minDRrap(truth_jet, electrons_selected_truth) ); }
      // if( HG::minDRrap( truth_jet, b_hadrons_all ) > 0.3 ) { Info("execute()", "Not close enough to b-quark: %lf", HG::minDRrap(truth_jet, b_hadrons_all) ); }
    }
  }
  std::sort( b_jets_selected_truth.begin(), b_jets_selected_truth.end(), pTdown );

  // Apply selection - b-quarks
  HG::TruthPtcls b_quarks_selected_truth(SG::VIEW_ELEMENTS);
  for( auto ptcl : *truthPtcls ) {
    if( ptcl->status() != 23 ) { continue; }
    if( fabs(ptcl->pdgId()) != 5 ) { continue; }
    if( ptcl->pt() > 20*HG::GeV ) { b_quarks_selected_truth.push_back(ptcl); }
  }
  std::sort( b_quarks_selected_truth.begin(), b_quarks_selected_truth.end(), pTdown );

  // Higgs pT, mass and rapidity (truth yy)
  TLorentzVector* yy(0);
  if( photons_selected_truth.size() > 1 ) {
    yy = new TLorentzVector( photons_selected_truth[0]->p4() + photons_selected_truth[1]->p4() );
    histoStore()->fillTH1F("truth_pT_yy", yy->Pt()/HG::GeV, weight);
    histoStore()->fillTH1F("truth_m_yy",  yy->M()/HG::GeV,  weight);
    histoStore()->fillTH1F("truth_y_yy",  yy->Rapidity(),   weight);
  }

  // Higgs pT, mass and rapidity (truth bb)
  TLorentzVector* bb(0);
  if( b_jets_selected_truth.size() > 1 ) {
    bb = new TLorentzVector( b_jets_selected_truth[0]->p4() + b_jets_selected_truth[1]->p4() );
    histoStore()->fillTH1F("truth_m_bb",  bb->M()/HG::GeV,  weight);
    histoStore()->fillTH1F("truth_pT_bb", bb->Pt()/HG::GeV, weight);
    histoStore()->fillTH1F("truth_y_bb",  bb->Rapidity(),   weight);
  }

  // Higgs pT, mass and rapidity (truth qq)
  TLorentzVector* qq(0);
  if( b_quarks_selected_truth.size() > 1 ) {
    qq = new TLorentzVector( b_quarks_selected_truth[0]->p4() + b_quarks_selected_truth[1]->p4() );
    histoStore()->fillTH1F("truth_m_qq",  qq->M()/HG::GeV,  weight);
    histoStore()->fillTH1F("truth_pT_qq", qq->Pt()/HG::GeV, weight);
    histoStore()->fillTH1F("truth_y_qq",  qq->Rapidity(),   weight);
  }

  // Resonance mass
  if( bb != 0 and yy != 0 ) {
    TLorentzVector bbyy = *bb + *yy;
    histoStore()->fillTH1F("truth_m_bbyy",       bbyy.M()/HG::GeV,                weight);
    histoStore()->fillTH2F("truth_m_bb_vs_m_yy", bb->M()/HG::GeV,yy->M()/HG::GeV, weight);
  }
  if( qq != 0 and yy != 0 ) {
    TLorentzVector qqyy = *qq + *yy;
    histoStore()->fillTH1F("truth_m_qqyy",       qqyy.M()/HG::GeV,                weight);
    histoStore()->fillTH2F("truth_m_qq_vs_m_yy", qq->M()/HG::GeV,yy->M()/HG::GeV, weight);
  }

  histoStore()->fillTH1F("event_weight", weight);

  if( bb != 0 ) { delete bb; }
  if( qq != 0 ) { delete qq; }
  if( yy != 0 ) { delete yy; }

  if( iteration % 1000 == 0 ) {
    std::cout << "truth_m_yy: " << histoStore()->getTH1F("truth_m_yy")->GetEntries() << std::endl;
    std::cout << "truth_m_bb: " << histoStore()->getTH1F("truth_m_bb")->GetEntries() << std::endl;
    std::cout << "truth_m_qq: " << histoStore()->getTH1F("truth_m_qq")->GetEntries() << std::endl;
    std::cout << "truth_m_bbyy: " << histoStore()->getTH1F("truth_m_bbyy")->GetEntries() << std::endl;
    std::cout << "truth_m_qqyy: " << histoStore()->getTH1F("truth_m_qqyy")->GetEntries() << std::endl;
    std::cout << "Total events: " << iteration << std::endl;
  }
  return EL::StatusCode::SUCCESS;
}
