/**
 * @file   OneTagCategorisation.cxx
 * @Author James Robinson (james.robinson@cern.ch)
 * @date   10th March 2016
 * @brief  Study h->bj selection
 *
 * Package for Hgamma analysis framework to select correct non-b-tagged jet
 */

#include "bbyyAnalysis/OneTagCategorisation.h"
#include "bbyyAnalysis/CommonTools.hpp"
#include "HGamAnalysisFramework/TruthUtils.h"
#include <AsgTools/MsgStream.h>
#include <AsgTools/MsgStreamMacros.h>
#include <boost/format.hpp>
// this is needed to distribute the algorithm to the workers
ClassImp(OneTagCategorisation)

OneTagCategorisation::OneTagCategorisation(const char *name)
  : HgammaAnalysis(name)
  , m_1_tag_WP("")
  , m_2_tag_WP("")
  , m_event_tree_1tag(0)
  , m_event_tree_2tag(0)
  , m_v_abs_eta_j(0)
  , m_v_abs_eta_jb(0)
  , m_v_Delta_eta_jb(0)
  , m_v_Delta_phi_jb(0)
  , m_v_idx_by_mH(0)
  , m_v_idx_by_pT(0)
  , m_v_idx_by_pT_jb(0)
  , m_v_m_jb(0)
  , m_v_passes_WP77(0)
  , m_v_passes_WP85(0)
  , m_v_pT_j(0)
  , m_v_pT_jb(0)
  , m_v_isCorrect(0)
  , m_event_weight(0)
  , m_sum_mc_weights(0)
  , m_diphoton_pT(0)
  , m_diphoton_eta(0)
  , m_diphoton_phi(0)
  , m_diphoton_m(0)
  , m_cutFlow({{"Events", 0}, {"PassingPreselection", 0}, {"PassedBTagging", 0}, {"is1tag", 0}, {"is2tag", 0}, {"correctPairs", 0}, {"incorrectPairs", 0}})
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
OneTagCategorisation::~OneTagCategorisation()
{
  // Here you delete any memory you allocated during your analysis.
}

/**
 * Initialize: inherited from EL::Algorithm
 * @return an EL::StatusCode indicating success/failure
 */
EL::StatusCode OneTagCategorisation::initialize()
{
  ATH_MSG_INFO("Initialising...");
  const auto sc = HgammaAnalysis::initialize();

  if (sc != EL::StatusCode::SUCCESS) { return sc; }

  // Retrieve b-tagging working point
  m_1_tag_WP = config()->getStr("OneTagCategorisation.1tag.OperatingPoint", "MV2c10_FixedCutBEff_60");
  m_2_tag_WP = config()->getStr("OneTagCategorisation.2tag.OperatingPoint", "MV2c10_FixedCutBEff_70");

  return EL::StatusCode::SUCCESS;
}

/**
 * Create histogram/tree output: inherited from EL::Algorithm
 * @return an EL::StatusCode indicating success/failure
 */
EL::StatusCode OneTagCategorisation::createOutput()
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

  // Add 1-tag event tree to output file
  m_event_tree_1tag = new TTree("events_1tag", "events_1tag");
  m_event_tree_1tag->SetDirectory(file);
  m_event_tree_1tag->Branch("abs_eta_j",    &m_v_abs_eta_j);
  m_event_tree_1tag->Branch("abs_eta_jb",   &m_v_abs_eta_jb);
  m_event_tree_1tag->Branch("Delta_eta_jb", &m_v_Delta_eta_jb);
  m_event_tree_1tag->Branch("Delta_phi_jb", &m_v_Delta_phi_jb);
  m_event_tree_1tag->Branch("idx_by_mH",    &m_v_idx_by_mH);
  m_event_tree_1tag->Branch("idx_by_pT",    &m_v_idx_by_pT);
  m_event_tree_1tag->Branch("idx_by_pT_jb", &m_v_idx_by_pT_jb);
  m_event_tree_1tag->Branch("m_jb",         &m_v_m_jb);
  m_event_tree_1tag->Branch("passes_WP77",  &m_v_passes_WP77);
  m_event_tree_1tag->Branch("passes_WP85",  &m_v_passes_WP85);
  m_event_tree_1tag->Branch("pT_j",         &m_v_pT_j);
  m_event_tree_1tag->Branch("pT_jb",        &m_v_pT_jb);
  m_event_tree_1tag->Branch("isCorrect",    &m_v_isCorrect);
  m_event_tree_1tag->Branch("event_weight", &m_event_weight);
  m_event_tree_1tag->Branch("diphoton_pT",  &m_diphoton_pT);
  m_event_tree_1tag->Branch("diphoton_eta", &m_diphoton_eta);
  m_event_tree_1tag->Branch("diphoton_phi", &m_diphoton_phi);
  m_event_tree_1tag->Branch("diphoton_m",   &m_diphoton_m);

  // Add 2-tag event tree to output file
  m_event_tree_2tag = new TTree("events_2tag", "events_2tag");
  m_event_tree_2tag->SetDirectory(file);
  m_event_tree_2tag->Branch("abs_eta_j",    &m_v_abs_eta_j);
  m_event_tree_2tag->Branch("abs_eta_jb",   &m_v_abs_eta_jb);
  m_event_tree_2tag->Branch("Delta_eta_jb", &m_v_Delta_eta_jb);
  m_event_tree_2tag->Branch("Delta_phi_jb", &m_v_Delta_phi_jb);
  m_event_tree_2tag->Branch("idx_by_mH",    &m_v_idx_by_mH);
  m_event_tree_2tag->Branch("idx_by_pT",    &m_v_idx_by_pT);
  m_event_tree_2tag->Branch("idx_by_pT_jb", &m_v_idx_by_pT_jb);
  m_event_tree_2tag->Branch("m_jb",         &m_v_m_jb);
  m_event_tree_2tag->Branch("passes_WP77",  &m_v_passes_WP77);
  m_event_tree_2tag->Branch("passes_WP85",  &m_v_passes_WP85);
  m_event_tree_2tag->Branch("pT_j",         &m_v_pT_j);
  m_event_tree_2tag->Branch("pT_jb",        &m_v_pT_jb);
  m_event_tree_2tag->Branch("isCorrect",    &m_v_isCorrect);
  m_event_tree_2tag->Branch("event_weight", &m_event_weight);
  return EL::StatusCode::SUCCESS;
}

/**
 * Main event loop: inherited from EL::Algorithm
 * @return an EL::StatusCode indicating success/failure
 */
EL::StatusCode OneTagCategorisation::execute()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  m_v_abs_eta_j.clear(); m_v_abs_eta_jb.clear();
  m_v_Delta_eta_jb.clear(); m_v_Delta_phi_jb.clear();
  m_v_idx_by_mH.clear(); m_v_idx_by_pT.clear(); m_v_idx_by_pT_jb.clear();
  m_v_m_jb.clear(); m_v_passes_WP77.clear(); m_v_passes_WP85.clear();
  m_v_pT_j.clear(); m_v_pT_jb.clear(); m_v_isCorrect.clear();

  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  const auto sc = HgammaAnalysis::execute();
  if (sc != EL::StatusCode::SUCCESS) { return sc; }

  // Get MC weight
  m_sum_mc_weights += eventHandler()->mcWeight();
  m_cutFlow["Events"]++;

  // Get overall event weight, normalised to 1fb-1
  unsigned int mcChannelNumber = eventInfo()->mcChannelNumber();
  m_event_weight = eventHandler()->mcWeight() * CommonTools::luminosity_invfb() *
                   CommonTools::xs_fb(mcChannelNumber, getCrossSection(mcChannelNumber), true) *
                   HgammaAnalysis::getGeneratorEfficiency(mcChannelNumber) *
                   HgammaAnalysis::getKFactor(mcChannelNumber) / CommonTools::sumOfWeights(mcChannelNumber);

  // Also add pileup and vertex weight
  m_event_weight *= eventHandler()->pileupWeight() * eventHandler()->vertexWeight();

  // ___________________________________________________________________________________________
  // Retrieve truth Higgs bosons
  xAOD::TruthParticleContainer higgsBosons = truthHandler()->getHiggsBosons();

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

  const xAOD::Photon *ph1 = photons_selected[0];
  const xAOD::Photon *ph2 = photons_selected[1];

  TLorentzVector diphoton = ph1->p4() + ph2->p4();
  m_diphoton_pT  = diphoton.Pt();
  m_diphoton_eta = diphoton.Eta();
  m_diphoton_phi = diphoton.Phi();
  m_diphoton_m   = diphoton.M();

  // ___________________________________________________________________________________________
  // Decorate muon correction to jets
  CommonTools::decorateMuonCorrection( yybbTool(), jets_selected, muons_corrected );

  // ___________________________________________________________________________________________
  // Perform matching between jet-pairs and Higgs
  bool isHiggsEvent = CommonTools::decorateHiggsMatching( jets_selected, higgsBosons );

  // ___________________________________________________________________________________________
  // Construct b-jets container
  xAOD::JetContainer jets_passing_2tag_WP(SG::VIEW_ELEMENTS);
  xAOD::JetContainer jets_passing_1tag_WP(SG::VIEW_ELEMENTS);
  xAOD::JetContainer jets_failing_1tag_WP(SG::VIEW_ELEMENTS);

  for (auto jet : jets_selected) {
    if (jet->auxdata<char>(m_2_tag_WP)) { jets_passing_2tag_WP.push_back(jet); }
    if (jet->auxdata<char>(m_1_tag_WP)) { jets_passing_1tag_WP.push_back(jet); }
    else { jets_failing_1tag_WP.push_back(jet); }
  }

  // ___________________________________________________________________________________________
  // Reject 0-tag and 3-tag events
  if (jets_passing_2tag_WP.size() > 2) { return EL::StatusCode::SUCCESS; } // more than two loose b-tags => 3-tag
  if (jets_passing_1tag_WP.size() < 1) { return EL::StatusCode::SUCCESS; } // less than one tight b-tags => 0-tag
  m_cutFlow["PassedBTagging"]++;

  // Consider 1-tag and 2-tag events separately
  if (jets_passing_2tag_WP.size() < 2 && jets_passing_1tag_WP.size() == 1) {
    // ___________________________________________________________________________________________
    // This is a 1-tag event
    m_cutFlow["is1tag"]++;

    // Decorate non-b-jets with their order in pT or distance from mH
    CommonTools::decorateWithIndices(*jets_passing_1tag_WP.at(0), jets_failing_1tag_WP);
    bool passesHiggsMatched(false), passesHadronConeExclTruthLabelID(false);

    // Case (A): this is a Higgs event
    if (isHiggsEvent) {
      // Loop over non b-jets and look for any that are truth-tagged
      for (auto jet : jets_failing_1tag_WP) {
        if (jet->auxdata<char>("HiggsMatched")) { passesHiggsMatched = true; }
        if (jet->auxdata<int>("HadronConeExclTruthLabelID") == 5) { passesHadronConeExclTruthLabelID = true; }
        // Correct pairing
        if (jet->auxdata<char>("HiggsMatched")) {
          appendToOutput( true, *jets_passing_1tag_WP.at(0), *jet );
          m_cutFlow["correctPairs"]++;
        // Incorrect pairing
        } else {
          appendToOutput( false, *jets_passing_1tag_WP.at(0), *jet );
          m_cutFlow["incorrectPairs"]++;
        }
      }

    // Case (B): this is not a Higgs event
    } else {
      // Loop over all non-b jets: all are incorrect pairs
      for (auto jet : jets_failing_1tag_WP) {
        appendToOutput( false, *jets_passing_1tag_WP.at(0), *jet );
        m_cutFlow["incorrectPairs"]++;
      }
    }

    // Fill 1-tag event-level tree
    if (passesHiggsMatched) { m_cutFlow["HiggsMatched"]++; }
    if (passesHadronConeExclTruthLabelID) { m_cutFlow["HadronConeExclTruthLabelID"]++; }
    m_event_tree_1tag->Fill();


  } else if (jets_passing_2tag_WP.size() == 2) {
    // ___________________________________________________________________________________________
    // This is a 2-tag event
    m_cutFlow["is2tag"]++;

    // Set indices to 1 in both cases
    xAOD::JetContainer second_bjet(SG::VIEW_ELEMENTS); second_bjet.push_back(jets_passing_2tag_WP.at(1));
    CommonTools::decorateWithIndices(*jets_passing_2tag_WP.at(0), second_bjet);

    // If both are matched then this is the correct pair
    if( jets_passing_2tag_WP.at(0)->auxdata<char>("HiggsMatched") && jets_passing_2tag_WP.at(1)->auxdata<char>("HiggsMatched") ) {
      appendToOutput( true, *jets_passing_2tag_WP.at(0), *jets_passing_2tag_WP.at(1) );
    } else {
      appendToOutput( false, *jets_passing_2tag_WP.at(0), *jets_passing_2tag_WP.at(1) );
    }

    // Fill 2-tag event-level tree
    m_event_tree_2tag->Fill();
  }

  return EL::StatusCode::SUCCESS;
}


/**
 * Cleanup after last event has been processed: inherited from EL::Algorithm
 * @return an EL::StatusCode indicating success/failure
 */
EL::StatusCode OneTagCategorisation::finalize() {
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
  ATH_MSG_INFO(m_cutFlow["PassedBTagging"] << " of " << m_cutFlow["PassingPreselection"] << " events (" << (m_cutFlow["PassingPreselection"] > 0 ? 100 * m_cutFlow["PassedBTagging"] / m_cutFlow["PassingPreselection"] : 0) << "%) were in the signal categories.");
  ATH_MSG_INFO("...  " << m_cutFlow["is1tag"] << " (" << (m_cutFlow["PassedBTagging"] > 0 ? 100 * m_cutFlow["is1tag"] / m_cutFlow["PassedBTagging"] : 0) << "%) 1-tag events");
  ATH_MSG_INFO("...  " << m_cutFlow["is2tag"] << " (" << (m_cutFlow["PassedBTagging"] > 0 ? 100 * m_cutFlow["is2tag"] / m_cutFlow["PassedBTagging"] : 0) << "%) 2-tag events");
  ATH_MSG_INFO("How many 1-tag events have a matched non-b-tagged jet using: ");
  ATH_MSG_INFO("... Higgs 4-vector matching?                           " << m_cutFlow["HiggsMatched"] << " of " << m_cutFlow["is1tag"] << " events (" << (m_cutFlow["is1tag"] > 0 ? 100 * m_cutFlow["HiggsMatched"] / m_cutFlow["is1tag"] : 0) << "%)");
  ATH_MSG_INFO("... HadronConeExclTruthLabelID truth-tagging?          " << m_cutFlow["HadronConeExclTruthLabelID"] << " of " << m_cutFlow["is1tag"] << " events (" << (m_cutFlow["is1tag"] > 0 ? 100 * m_cutFlow["HadronConeExclTruthLabelID"] / m_cutFlow["is1tag"] : 0) << "%)");
  ATH_MSG_INFO("=> this corresponds to: " << m_cutFlow["correctPairs"] << " correct pairs and " << m_cutFlow["incorrectPairs"] << " incorrect pairs");
  return EL::StatusCode::SUCCESS;
}


/**
 * Add jet pairing information to event-level vectors
 * @return nothing
 */
void OneTagCategorisation::appendToOutput( const bool& isCorrect, const xAOD::Jet& bjet, const xAOD::Jet& otherjet ) {
  // Construct jb 4-vector
  TLorentzVector b_p4(CommonTools::p4(bjet)), j_p4(CommonTools::p4(otherjet));
  TLorentzVector jb_p4 = b_p4 + j_p4;
  // Append to vectors of event-level quantities
  m_v_abs_eta_j.push_back( fabs(j_p4.Eta()) );
  m_v_abs_eta_jb.push_back( fabs(jb_p4.Eta()) );
  m_v_Delta_eta_jb.push_back( fabs(b_p4.Eta() - j_p4.Eta()) );
  m_v_Delta_phi_jb.push_back( fabs(b_p4.DeltaPhi(j_p4)) );
  m_v_idx_by_mH.push_back( otherjet.auxdata<int>("idx_by_mH") );
  m_v_idx_by_pT.push_back( otherjet.auxdata<int>("idx_by_pT") );
  m_v_idx_by_pT_jb.push_back( otherjet.auxdata<int>("idx_by_pT_jb") );
  m_v_m_jb.push_back( jb_p4.M() / HG::GeV );
  m_v_passes_WP77.push_back( (otherjet.auxdata<char>("MV2c10_FixedCutBEff_77") ? 1 : 0) );
  m_v_passes_WP85.push_back( (otherjet.auxdata<char>("MV2c10_FixedCutBEff_85") ? 1 : 0));
  m_v_pT_j.push_back( j_p4.Pt() / HG::GeV );
  m_v_pT_jb.push_back( jb_p4.Pt() / HG::GeV );
  m_v_isCorrect.push_back( isCorrect );
}
