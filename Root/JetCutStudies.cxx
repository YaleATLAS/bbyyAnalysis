/**
 * @file   JetCutStudies.cxx
 * @Author James Robinson (james.robinson@cern.ch)
 * @date   18th July 2016
 * @brief  Study jet pT cuts
 *
 * Package for Hgamma analysis framework to study jet pT cuts
 */

#include "bbyyAnalysis/JetCutStudies.h"
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
  , m_event_weight(0)
  , m_sum_mc_weights(0)
  , m_cutFlow({{"Events", 0}, {"PassingPreselection", 0}, {"ExactlyOneRecoTagged", 0}, {"TruthTagged", 0}, {"TruthMatched", 0}, {"TruthTaggedAndMatched", 0}})
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

  // Retrieve b-tagging working point
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
  m_event_tree->Branch("photon_n",       &m_photon_n);
  m_event_tree->Branch("photon_pT",      &m_photon_pT);
  m_event_tree->Branch("photon_eta",     &m_photon_eta);
  m_event_tree->Branch("photon_phi",     &m_photon_phi);
  m_event_tree->Branch("photon_E",       &m_photon_E);
  m_event_tree->Branch("photon_isTight", &m_photon_isTight);
  m_event_tree->Branch("m_yy",           &m_m_yy);
  m_event_tree->Branch("jet_n",          &m_jet_n);
  m_event_tree->Branch("jet_pT",         &m_jet_pT);
  m_event_tree->Branch("jet_eta",        &m_jet_eta);
  m_event_tree->Branch("jet_phi",        &m_jet_phi);
  m_event_tree->Branch("jet_E",          &m_jet_E);
  m_event_tree->Branch("jet_btag_loose", &m_jet_btag_loose);
  m_event_tree->Branch("jet_btag_tight", &m_jet_btag_tight);
  m_event_tree->Branch("jet_truth_tag",  &m_jet_truth_tag);
  m_event_tree->Branch("event_weight",   &m_event_weight);
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

  // Get overall event weight, normalised to 1fb-1
  unsigned int mcChannelNumber = eventInfo()->mcChannelNumber();
  m_event_weight = eventHandler()->mcWeight() * sampleXS(mcChannelNumber) *
                   HgammaAnalysis::getGeneratorEfficiency(mcChannelNumber) *
                   HgammaAnalysis::getKFactor(mcChannelNumber) / sumOfWeights(mcChannelNumber);

  // ___________________________________________________________________________________________
  // Fetch default jets
  xAOD::JetContainer jets_corrected = jetHandler()->getCorrectedContainer();
  xAOD::JetContainer jets_selected  = jetHandler()->applySelection(jets_corrected);

  // ___________________________________________________________________________________________
  // Fetch default photons
  xAOD::PhotonContainer photons_corrected = photonHandler()->getCorrectedContainer();
  xAOD::PhotonContainer photons_selected  = photonHandler()->applySelection(photons_corrected);

  // ___________________________________________________________________________________________
  // Fetch default muons
  xAOD::MuonContainer muons_corrected = muonHandler()->getCorrectedContainer();
  xAOD::MuonContainer muons_selected  = muonHandler()->applySelection(muons_corrected);

  // ___________________________________________________________________________________________
  // Fetch default electrons
  xAOD::ElectronContainer electrons_corrected = electronHandler()->getCorrectedContainer();
  xAOD::ElectronContainer electrons_selected  = electronHandler()->applySelection(electrons_corrected);

  // ___________________________________________________________________________________________
  // Reject event if it fails the HGamma preselection
  if (!HgammaAnalysis::pass(&photons_selected, &electrons_selected, &muons_selected, &jets_selected)) {
    return StatusCode::SUCCESS;
  }
  m_cutFlow["PassingPreselection"]++;


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


  // ___________________________________________________________________________________________
  // Decorate muon correction to all jets and retrieve corrected jets
  yybbTool()->decorateMuonCorrections(jets_selected, muons_corrected);
  std::vector<TLorentzVector> jets_muon_corrected = yybbTool()->getMuonCorrJet4Vs(jets_selected);
  if ( jets_selected.size() != jets_muon_corrected.size() ) { return StatusCode::FAILURE; }

  // Clear vectors
  m_jet_pT.clear(); m_jet_eta.clear(); m_jet_phi.clear(); m_jet_E.clear();
  m_jet_btag_loose.clear(); m_jet_btag_tight.clear(); m_jet_truth_tag.clear();

  // Fill jet information into tree
  m_jet_n = jets_selected.size();
  for ( unsigned int idx_jet = 0; idx_jet < m_jet_n; ++idx_jet ) {
    m_jet_pT.push_back( jets_muon_corrected.at(idx_jet).Pt() / HG::GeV );
    m_jet_eta.push_back( jets_muon_corrected.at(idx_jet).Eta() );
    m_jet_phi.push_back( jets_muon_corrected.at(idx_jet).Phi() );
    m_jet_E.push_back( jets_muon_corrected.at(idx_jet).E() / HG::GeV );
    m_jet_btag_loose.push_back( jets_selected.at(idx_jet)->auxdata<char>(m_2_tag_WP) );
    m_jet_btag_tight.push_back( jets_selected.at(idx_jet)->auxdata<char>(m_1_tag_WP) );
    m_jet_truth_tag.push_back( jets_selected.at(idx_jet)->auxdata<int>("HadronConeExclTruthLabelID") == 5 );
  }

  // Fill event-level tree
  m_event_tree->Fill();

  return EL::StatusCode::SUCCESS;
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
  ATH_MSG_INFO(m_cutFlow["PassingPreselection"] << " of " << m_cutFlow["Events"] << " events (" << (m_cutFlow["Events"] > 0 ? 100 * m_cutFlow["PassingPreselection"] / m_cutFlow["Events"] : 0) << "%) passed the HGamma pre-selection.");
  return EL::StatusCode::SUCCESS;
}

/**
 * Create histogram output: inherited from EL::Algorithm
 * @return an EL::StatusCode indicating success/failure
 */
double JetCutStudies::sampleXS(int mcID) {
  // Use SM hh cross-section for resonances
  if (mcID == 341173) { return 1e3 * this->getCrossSection(342620) / 5.0; } // X275->hh->yybb
  if (mcID == 341004) { return 1e3 * this->getCrossSection(342620) / 5.0; } // X300->hh->yybb
  if (mcID == 341174) { return 1e3 * this->getCrossSection(342620) / 5.0; } // X325->hh->yybb
  if (mcID == 341175) { return 1e3 * this->getCrossSection(342620) / 5.0; } // X350->hh->yybb
  if (mcID == 341176) { return 1e3 * this->getCrossSection(342620) / 5.0; } // X400->hh->yybb
  return 1e3 * this->getCrossSection(mcID);
}

/**
 * Create histogram output: inherited from EL::Algorithm
 * @return an EL::StatusCode indicating success/failure
 */
double JetCutStudies::sumOfWeights(int mcID) {
  if (mcID == 341061) { return 197000; }     // SM yybb
  if (mcID == 341062) { return 200000; }     // SM yjbb
  if (mcID == 341063) { return 196000; }     // SM yybj
  if (mcID == 341064) { return 200000; }     // SM yjjb
  if (mcID == 341065) { return 180000; }     // SM yyjj
  if (mcID == 341066) { return 198000; }     // SM yjjj
  if (mcID == 341559) { return 100000; }     // SM hh->yybb
  if (mcID == 341173) { return 100000; }     // X275->hh->yybb
  if (mcID == 341004) { return 100000; }     // X300->hh->yybb
  if (mcID == 341174) { return 100000; }     // X325->hh->yybb
  if (mcID == 341175) { return 100000; }     // X350->hh->yybb
  if (mcID == 341176) { return 100000; }     // X400->hh->yybb
  if (mcID == 341939) { return 80604609.6; } // Sherpa photons+jets
  if (mcID == 341559) { return 100000; }     // SM hh->yybb
  if (mcID == 342620) { return 2134.103; }   // SM hh->yybb (from 200000 events)
  return 1.0;
}
