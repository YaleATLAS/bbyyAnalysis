/**
 * @file   JetCutStudies.cxx
 * @Author James Robinson (james.robinson@cern.ch)
 * @date   18th July 2016
 * @brief  Study jet pT cuts
 *
 * Package for Hgamma analysis framework to study jet pT cuts
 */

#include "bbyyAnalysis/JetCutStudies.h"
#include "bbyyAnalysis/CommonTools.hpp"
#include "HGamAnalysisFramework/TruthUtils.h"
#include <AsgTools/MsgStream.h>
#include <AsgTools/MsgStreamMacros.h>
#include <boost/format.hpp>

// this is needed to distribute the algorithm to the workers
ClassImp(JetCutStudies)

JetCutStudies::JetCutStudies(const char *name)
  : HgammaAnalysis(name)
  , m_1_tag_WP("")
  , m_2_tag_WP("")
  , m_event_tree(0)
  , m_reader()
  , m_abs_eta_j(0)
  , m_abs_eta_jb(0)
  , m_Delta_eta_jb(0)
  , m_Delta_phi_jb(0)
  , m_idx_by_mH(0)
  , m_idx_by_pT(0)
  , m_m_jb(0)
  , m_pT_j(0)
  , m_pT_jb(0)
  , m_event_weight(0)
  , m_pileup_weight(0)
  , m_sum_mc_weights(0)
  , m_sum_pileup_weights(0)
  , m_cutFlow({{"Events", 0}, {"PassingPreselection", 0}})
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
  this->SetName(name);
}

/**
 * Destructor
 */
JetCutStudies::~JetCutStudies()
{
  // Here you delete any memory you allocated during your analysis.
  // if (m_reader) { delete m_reader; }
}

/**
 * Initialize: inherited from EL::Algorithm
 * @return an EL::StatusCode indicating success/failure
 */
EL::StatusCode JetCutStudies::initialize()
{
  ATH_MSG_INFO("Initialising...");
  const auto sc = HgammaAnalysis::initialize();
  if (sc != EL::StatusCode::SUCCESS) { return sc; }

  // Setup TMVA reader
  ATH_MSG_INFO("Initialising TMVA reader");
  m_reader.AddVariable("abs_eta_j",    &m_abs_eta_j);
  m_reader.AddVariable("abs_eta_jb",   &m_abs_eta_jb);
  m_reader.AddVariable("Delta_eta_jb", &m_Delta_eta_jb);
  m_reader.AddVariable("idx_by_mH",    &m_idx_by_mH);
  m_reader.AddVariable("idx_by_pT",    &m_idx_by_pT);
  m_reader.AddVariable("m_jb",         &m_m_jb);
  m_reader.AddVariable("pT_j",         &m_pT_j);
  m_reader.AddVariable("pT_jb",        &m_pT_jb);
  m_reader.BookMVA("OneTagClassifier", PathResolverFindCalibFile("bbyyAnalysis/skl_BDT_TMVA.weights.xml"));

  // Retrieve b-tagging working point
  ATH_MSG_INFO("Reading in values from config files...");
  m_1_tag_WP = config()->getStr("JetCutStudies.1tag.OperatingPoint", "MV2c10_FixedCutBEff_60");
  m_2_tag_WP = config()->getStr("JetCutStudies.2tag.OperatingPoint", "MV2c10_FixedCutBEff_85");

  return EL::StatusCode::SUCCESS;
}

/**
 * Create histogram/tree output: inherited from EL::Algorithm
 * @return an EL::StatusCode indicating success/failure
 */
EL::StatusCode JetCutStudies::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.
  /// Initialise baseclass
  const auto sc = HgammaAnalysis::createOutput();

  if (sc != EL::StatusCode::SUCCESS) { return sc; }

  // Add correct choice tree to output file
  ATH_MSG_INFO("Creating output trees...");
  TFile *file = wk()->getOutputFile("MxAOD");
  // Add event tree to output file
  m_event_tree = new TTree("events", "events");
  m_event_tree->SetDirectory(file);
  m_event_tree->Branch("photon_n",          &m_photon_n);
  m_event_tree->Branch("photon_pT",         &m_photon_pT);
  m_event_tree->Branch("photon_eta",        &m_photon_eta);
  m_event_tree->Branch("photon_phi",        &m_photon_phi);
  m_event_tree->Branch("photon_E",          &m_photon_E);
  m_event_tree->Branch("photon_isTight",    &m_photon_isTight);
  m_event_tree->Branch("m_yy",              &m_m_yy);
  m_event_tree->Branch("jet_n",             &m_jet_n);
  m_event_tree->Branch("jet_pT",            &m_jet_pT);
  m_event_tree->Branch("jet_eta",           &m_jet_eta);
  m_event_tree->Branch("jet_phi",           &m_jet_phi);
  m_event_tree->Branch("jet_E",             &m_jet_E);
  m_event_tree->Branch("jet_btag_loose",    &m_jet_btag_loose);
  m_event_tree->Branch("jet_btag_tight",    &m_jet_btag_tight);
  m_event_tree->Branch("jet_truth_tag",     &m_jet_truth_tag);
  m_event_tree->Branch("jet_higgs_match",   &m_jet_higgs_match);
  m_event_tree->Branch("jet_JVT",           &m_jet_JVT);
  m_event_tree->Branch("jet_eta_det",       &m_jet_eta_det);
  m_event_tree->Branch("jet_ML_classifier", &m_jet_ML_classifier);
  m_event_tree->Branch("event_weight",      &m_event_weight);
  m_event_tree->Branch("pileup_weight",     &m_pileup_weight);
  return EL::StatusCode::SUCCESS;
}

/**
 * Main event loop: inherited from EL::Algorithm
 * @return an EL::StatusCode indicating success/failure
 */
EL::StatusCode JetCutStudies::execute()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  const auto sc = HgammaAnalysis::execute();

  if (sc != EL::StatusCode::SUCCESS) { return sc; }

  // Get MC weight
  m_sum_mc_weights += eventHandler()->mcWeight();
  m_cutFlow["Events"]++;

  if (!isMC()) {
    m_event_weight = CommonTools::luminosity_invfb() / 12.171;
    m_pileup_weight = 1;
  } else {
    // Get overall event weight, normalised to 1fb-1
    unsigned int mcChannelNumber = eventInfo()->mcChannelNumber();
    m_event_weight = eventHandler()->mcWeight() * CommonTools::luminosity_invfb() *
                     CommonTools::xs_fb(mcChannelNumber, getCrossSection(mcChannelNumber), false) *
                     HgammaAnalysis::getGeneratorEfficiency(mcChannelNumber) *
                     HgammaAnalysis::getKFactor(mcChannelNumber) / CommonTools::sumOfWeights(mcChannelNumber);

    // Get pileup weight
    m_pileup_weight = eventHandler()->pileupWeight();
  }
  m_sum_pileup_weights += m_pileup_weight;

  // ___________________________________________________________________________________________
  // Retrieve truth Higgs bosons
  xAOD::TruthParticleContainer higgsBosons(SG::VIEW_ELEMENTS);
  if (isMC()) { higgsBosons = truthHandler()->getHiggsBosons(); }

  // ___________________________________________________________________________________________
  // Retrieve default jets
  xAOD::JetContainer jets_corrected = jetHandler()->getCorrectedContainer();
  xAOD::JetContainer jets_selected  = jetHandler()->applySelection(jets_corrected);

  // ___________________________________________________________________________________________
  // Retrieve default photons
  xAOD::PhotonContainer photons_corrected = photonHandler()->getCorrectedContainer();
  xAOD::PhotonContainer photons_selected  = photonHandler()->applySelection(photons_corrected);

  // ___________________________________________________________________________________________
  // Retrieve default muons
  xAOD::MuonContainer muons_corrected = muonHandler()->getCorrectedContainer();
  xAOD::MuonContainer muons_selected  = muonHandler()->applySelection(muons_corrected);

  // ___________________________________________________________________________________________
  // Retrieve default electrons
  xAOD::ElectronContainer electrons_corrected = electronHandler()->getCorrectedContainer();
  xAOD::ElectronContainer electrons_selected  = electronHandler()->applySelection(electrons_corrected);

  // ___________________________________________________________________________________________
  // Reject event if it fails the HGamma preselection
  if (!HgammaAnalysis::pass(&photons_selected, &electrons_selected, &muons_selected, &jets_selected)) {
    return StatusCode::SUCCESS;
  }
  m_cutFlow["PassingPreselection"]++;

  // ___________________________________________________________________________________________
  // Decorate muon correction to jets
  CommonTools::decorateMuonCorrection( yybbTool(), jets_selected, muons_corrected );

  // ___________________________________________________________________________________________
  // Perform matching between jet-pairs and Higgs
  CommonTools::decorateHiggsMatching( jets_selected, higgsBosons );

  // Clear vectors
  m_photon_pT.clear(); m_photon_eta.clear(); m_photon_phi.clear(); m_photon_E.clear();
  m_photon_isTight.clear();

  // Fill photon information into tree
  m_photon_n = photons_selected.size();
  for (auto photon : photons_selected) {
    m_photon_pT.push_back( photon->pt() / HG::GeV );
    m_photon_eta.push_back( photon->eta() );
    m_photon_phi.push_back( photon->phi() );
    m_photon_E.push_back( photon->e() / HG::GeV );
    m_photon_isTight.push_back( bool(photonHandler()->passIsoCut(photon, HG::Iso::FixedCutTight)) );
  }

  // Add diphoton mass as cross-check
  if (photons_selected.size() < 2) {
    m_m_yy = -1;
  } else {
    TLorentzVector yy_p4 = photons_selected.at(0)->p4() + photons_selected.at(1)->p4();
    m_m_yy = yy_p4.M() / HG::GeV;
  }

  // Clear vectors
  m_jet_pT.clear(); m_jet_eta.clear(); m_jet_phi.clear(); m_jet_E.clear();
  m_jet_btag_loose.clear(); m_jet_btag_tight.clear(); m_jet_truth_tag.clear();
  m_jet_higgs_match.clear(); m_jet_JVT.clear(); m_jet_eta_det.clear(); m_jet_ML_classifier.clear();

  // ___________________________________________________________________________________________
  // Construct b-jet containers
  xAOD::JetContainer jets_passing_2tag_WP(SG::VIEW_ELEMENTS);
  xAOD::JetContainer jets_passing_1tag_WP(SG::VIEW_ELEMENTS);
  xAOD::JetContainer jets_failing_1tag_WP(SG::VIEW_ELEMENTS);

  // Initialise 1-tag classifier to false and fill b-jet containers
  SG::AuxElement::Accessor<double> accOneTagClassifier("OneTagClassifier");
  for (auto jet : jets_selected) {
    accOneTagClassifier(*jet) = -99;
    if (jet->auxdata<char>(m_2_tag_WP)) { jets_passing_2tag_WP.push_back(jet); }
    if (jet->auxdata<char>(m_1_tag_WP)) { jets_passing_1tag_WP.push_back(jet); }
    else { jets_failing_1tag_WP.push_back(jet); }
  }

  // Decorate only 1-tag events with classifier
  if (jets_passing_2tag_WP.size() < 2 && jets_passing_1tag_WP.size() == 1) {
    decorateWithClassifier(*jets_passing_1tag_WP.at(0), jets_failing_1tag_WP);
  }

  // Fill jet information into tree
  m_jet_n = jets_selected.size();
  for( const auto& jet : jets_selected ) {
    m_jet_pT.push_back( jet->auxdata<double>("muon_pT") / HG::GeV );
    m_jet_eta.push_back( jet->auxdata<double>("muon_eta") );
    m_jet_phi.push_back( jet->auxdata<double>("muon_phi") );
    m_jet_E.push_back( jet->auxdata<double>("muon_E") / HG::GeV );
    m_jet_btag_loose.push_back( jet->auxdata<char>(m_2_tag_WP) );
    m_jet_btag_tight.push_back( jet->auxdata<char>(m_1_tag_WP) );
    m_jet_truth_tag.push_back( jet->auxdata<int>("HadronConeExclTruthLabelID") == 5 );
    m_jet_higgs_match.push_back( jet->auxdata<char>("HiggsMatched") );
    m_jet_JVT.push_back( jet->auxdata<float>("Jvt") );
    m_jet_ML_classifier.push_back( jet->auxdata<double>("OneTagClassifier") );
    if(isMAOD()) {
      m_jet_eta_det.push_back( jet->auxdata<float>("DetectorEta") );
    } else {
      m_jet_eta_det.push_back( jet->getAttribute<xAOD::JetFourMom_t>("JetConstitScaleMomentum").eta() );
    }
    ATH_MSG_DEBUG("... found jet with classifier value: " << jet->auxdata<double>("OneTagClassifier") << " and HiggsMatched: " << (jet->auxdata<char>("HiggsMatched") ? "YES" : "NO"));
  }

  // Fill event-level tree
  m_event_tree->Fill();

  return EL::StatusCode::SUCCESS;
}

/**
 * Add jet pairing to event-level vectors
 * @return nothing
 */
// void JetCutStudies::decorateWithClassifier( const xAOD::Jet& bjet, const xAOD::Jet& otherjet ) {
void JetCutStudies::decorateWithClassifier( const xAOD::Jet& bjet, xAOD::JetContainer& nonbjets ) {
  // Add idx_by_mH and idx_by_pT decorations
  CommonTools::decorateWithIndices(bjet, nonbjets);
  for (auto otherjet : nonbjets) {
    // Construct jb 4-vector
    TLorentzVector b_p4(CommonTools::p4(bjet)), j_p4(CommonTools::p4(*otherjet));
    TLorentzVector jb_p4 = b_p4 + j_p4;
    // Append to vectors for event-level quantities
    m_abs_eta_j = fabs(j_p4.Eta());
    m_abs_eta_jb = fabs(jb_p4.Eta());
    m_Delta_eta_jb = fabs(b_p4.Eta() - j_p4.Eta());
    m_Delta_phi_jb = fabs(b_p4.DeltaPhi(j_p4));
    m_idx_by_mH = otherjet->auxdata<int>("idx_by_mH");
    m_idx_by_pT = otherjet->auxdata<int>("idx_by_pT");
    m_m_jb = jb_p4.M() / HG::GeV;
    m_pT_j = j_p4.Pt() / HG::GeV;
    m_pT_jb = jb_p4.Pt() / HG::GeV;
    otherjet->auxdata<double>("OneTagClassifier") = m_reader.EvaluateMVA("OneTagClassifier");
  }
}

/**
 * Cleanup after last event has been processed: inherited from EL::Algorithm
 * @return an EL::StatusCode indicating success/failure
 */
EL::StatusCode JetCutStudies::finalize() {
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

  if (sc != EL::StatusCode::SUCCESS) { return sc; }

  ATH_MSG_INFO("The sum of MC weights in this job for channel " << eventInfo()->mcChannelNumber() << " was " << m_sum_mc_weights << " from " << m_cutFlow["Events"]++ << " events.");
  ATH_MSG_INFO("The sum of pileup weights in this job was " << m_sum_pileup_weights << " for " << m_cutFlow["Events"]++ << " events.");
  ATH_MSG_INFO(m_cutFlow["PassingPreselection"] << " of " << m_cutFlow["Events"] << " events (" << (m_cutFlow["Events"] > 0 ? 100 * m_cutFlow["PassingPreselection"] / m_cutFlow["Events"] : 0) << "%) passed the HGamma pre-selection.");
  return EL::StatusCode::SUCCESS;
}
