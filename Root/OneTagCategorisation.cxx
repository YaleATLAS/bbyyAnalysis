/**
 * @file   OneTagCategorisation.cxx
 * @Author James Robinson (james.robinson@cern.ch)
 * @date   10th March 2016
 * @brief  Study h->bj selection
 *
 * Package for Hgamma analysis framework to study interference effects
 */

#include "bbyyAnalysis/OneTagCategorisation.h"
#include <AsgTools/MsgStream.h>
#include <AsgTools/MsgStreamMacros.h>
#include <boost/format.hpp>

// this is needed to distribute the algorithm to the workers
ClassImp(OneTagCategorisation)

// // Sorting by pT
// bool pTdown ( const xAOD::IParticle* i, const xAOD::IParticle* j ) { return ( i->pt() > j->pt() ); }

OneTagCategorisation::OneTagCategorisation(const char *name)
: HgammaAnalysis(name)
, m_correct_tree(0)
, m_incorrect_tree(0)
, m_b_tagging_WP("")
, m_sum_mc_weights(0)
, m_nEvents(0)
, m_nPassingPreselection(0)
, m_nPassingReco(0)
, m_nTruthTagged(0)
, m_nTruthMatched(0)
, m_nTruthTaggedAndMatched(0)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
  this->SetName(name);
}



OneTagCategorisation::~OneTagCategorisation()
{
  // Here you delete any memory you allocated during your analysis.
}



EL::StatusCode OneTagCategorisation::initialize ()
{
  const auto sc = HgammaAnalysis::initialize();
  if( sc != EL::StatusCode::SUCCESS ) { return sc; }
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode OneTagCategorisation::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.
  /// Initialise baseclass
  const auto sc = HgammaAnalysis::createOutput();
  if( sc != EL::StatusCode::SUCCESS ) { return sc; }

  // Retrieve b-tagging working point
  m_b_tagging_WP = std::string("MV2c20_") + config()->getStr("OneTagCategorisation.MV2c20.OperatingPoint", "FixedCutBEff_60");

  // Add correct choice tree to output file
  ATH_MSG_INFO( "Initialising output trees..." );
  TFile *file = wk()->getOutputFile("MxAOD");
  m_correct_tree = new TTree("correct","correct");
  m_correct_tree->SetDirectory(file);
  m_correct_tree->Branch( "abs_eta_j",        &m_abs_eta_j );
  m_correct_tree->Branch( "abs_eta_jb",       &m_abs_eta_jb );
  m_correct_tree->Branch( "Delta_eta_jb",     &m_Delta_eta_jb );
  m_correct_tree->Branch( "Delta_phi_jb",     &m_Delta_phi_jb );
  m_correct_tree->Branch( "idx_by_mH",        &m_idx_by_mH );
  m_correct_tree->Branch( "idx_by_pT",        &m_idx_by_pT );
  m_correct_tree->Branch( "m_jb",             &m_m_jb );
  m_correct_tree->Branch( "pT_j",             &m_pT_j );
  m_correct_tree->Branch( "pT_jb",            &m_pT_jb );
  m_correct_tree->Branch( "event_weight",     &m_event_weight );
  if( m_b_tagging_WP != "MV2c20_FixedCutBEff_60" ) { m_correct_tree->Branch( "MV2c20_FCBE_60", &m_MV2c20_FCBE_60 ); }
  if( m_b_tagging_WP != "MV2c20_FixedCutBEff_70" ) { m_correct_tree->Branch( "MV2c20_FCBE_70", &m_MV2c20_FCBE_70 ); }
  if( m_b_tagging_WP != "MV2c20_FixedCutBEff_77" ) { m_correct_tree->Branch( "MV2c20_FCBE_77", &m_MV2c20_FCBE_77 ); }
  if( m_b_tagging_WP != "MV2c20_FixedCutBEff_85" ) { m_correct_tree->Branch( "MV2c20_FCBE_85", &m_MV2c20_FCBE_85 ); }
  // Add incorrect choice tree to output file
  m_incorrect_tree = new TTree("incorrect","incorrect");
  m_incorrect_tree->SetDirectory(file);
  m_incorrect_tree->Branch( "abs_eta_j",        &m_abs_eta_j );
  m_incorrect_tree->Branch( "abs_eta_jb",       &m_abs_eta_jb );
  m_incorrect_tree->Branch( "Delta_eta_jb",     &m_Delta_eta_jb );
  m_incorrect_tree->Branch( "Delta_phi_jb",     &m_Delta_phi_jb );
  m_incorrect_tree->Branch( "idx_by_mH",        &m_idx_by_mH );
  m_incorrect_tree->Branch( "idx_by_pT",        &m_idx_by_pT );
  m_incorrect_tree->Branch( "m_jb",             &m_m_jb );
  m_incorrect_tree->Branch( "pT_j",             &m_pT_j );
  m_incorrect_tree->Branch( "pT_jb",            &m_pT_jb );
  m_incorrect_tree->Branch( "event_weight",     &m_event_weight );
  if( m_b_tagging_WP != "MV2c20_FixedCutBEff_60" ) { m_incorrect_tree->Branch( "MV2c20_FCBE_60", &m_MV2c20_FCBE_60 ); }
  if( m_b_tagging_WP != "MV2c20_FixedCutBEff_70" ) { m_incorrect_tree->Branch( "MV2c20_FCBE_70", &m_MV2c20_FCBE_70 ); }
  if( m_b_tagging_WP != "MV2c20_FixedCutBEff_77" ) { m_incorrect_tree->Branch( "MV2c20_FCBE_77", &m_MV2c20_FCBE_77 ); }
  if( m_b_tagging_WP != "MV2c20_FixedCutBEff_85" ) { m_incorrect_tree->Branch( "MV2c20_FCBE_85", &m_MV2c20_FCBE_85 ); }
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode OneTagCategorisation::execute()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  const auto sc = HgammaAnalysis::execute();
  if( sc != EL::StatusCode::SUCCESS ) { return sc; }

  /// Get MC weight
  m_sum_mc_weights += eventHandler()->mcWeight();
  ++m_nEvents;

  /// Get overall event weight, normalised to 1fb-1
  int mcID = eventInfo()->mcChannelNumber();
  m_event_weight = eventHandler()->mcWeight() * ( 1e3 * HgammaAnalysis::getCrossSection(mcID)) * HgammaAnalysis::getGeneratorEfficiency(mcID) * HgammaAnalysis::getKFactor(mcID) / sumOfWeights(mcID);

  //___________________________________________________________________________________________
  // Fetch default jets
  xAOD::JetContainer jets_corrected = jetHandler()->getCorrectedContainer();
  xAOD::JetContainer jets_selected = jetHandler()->applySelection(jets_corrected);

  //___________________________________________________________________________________________
  // Fetch default photons
  xAOD::PhotonContainer photons_corrected = photonHandler()->getCorrectedContainer();
  xAOD::PhotonContainer photons_selected = photonHandler()->applySelection(photons_corrected);

  //___________________________________________________________________________________________
  // Fetch default muons
  xAOD::MuonContainer muons_corrected = muonHandler()->getCorrectedContainer();
  xAOD::MuonContainer muons_selected = muonHandler()->applySelection(muons_corrected);

  //___________________________________________________________________________________________
  // Fetch default electrons
  xAOD::ElectronContainer electrons_corrected = electronHandler()->getCorrectedContainer();
  xAOD::ElectronContainer electrons_selected = electronHandler()->applySelection(electrons_corrected);

  //___________________________________________________________________________________________
  // Reject event if it fails the HGamma preselection
  if( !HgammaAnalysis::pass( &photons_selected, &electrons_selected, &muons_selected, &jets_selected ) ) {
    return StatusCode::SUCCESS;
  }
  ++m_nPassingPreselection;

  // Retrieve truth-particles and select final-state b-quarks descended from Higgs bosons
  const xAOD::TruthParticleContainer *truthPtcls = 0;
  EL_CHECK( "execute()", event()->retrieve(truthPtcls, "TruthParticles") );
  ConstDataVector<xAOD::TruthParticleContainer> bQuarksSelected = bQuarkHiggsDescendants( truthPtcls );

  // Perform matching to identify jet-quark pairs
  matchQuarksToJets( bQuarksSelected, jets_selected );

  // Construct b-jets container
  xAOD::JetContainer jets_selected_b_tagged( SG::VIEW_ELEMENTS );
  xAOD::JetContainer jets_selected_non_b_tagged( SG::VIEW_ELEMENTS );
  for( auto jet : jets_selected ) {
    bool reco_tagged( jet->auxdata<char>(m_b_tagging_WP) );
    if( reco_tagged ) { jets_selected_b_tagged.push_back( jet ); }
    else { jets_selected_non_b_tagged.push_back( jet ); }
  }

  // We only want events with exactly one b-jet so reject otherwise
  if( jets_selected_b_tagged.size() != 1 ) { return EL::StatusCode::SUCCESS; }
  ++m_nPassingReco;

  // Decorate non-b-jets with their position when ordering by pT or distance from mH
  decorateWithIndices( *jets_selected_b_tagged.at(0), jets_selected_non_b_tagged );

  // Determine whether the b-jet is truth tagged and Higgs matched
  bool b_is_truth_tagged( jets_selected_b_tagged.at(0)->auxdata<int>("HadronConeExclTruthLabelID") == 5 );
  bool b_is_Higgs_matched( jets_selected_b_tagged.at(0)->auxdata<char>("HiggsMatched") );

  // Case (A): matched to a b-quark Higgs descendant
  bool event_has_truth_tagged_jet(false), event_has_Higgs_matched_jet(false);
  if( b_is_Higgs_matched && b_is_truth_tagged ) {
    // Loop over non b-jets and look for any that are also Higgs-matched
    for( auto jet : jets_selected_non_b_tagged ) {
      // Check success rates of truth-tagging and Higgs-matching
      if( jet->auxdata<int>("HadronConeExclTruthLabelID") == 5 ) { event_has_truth_tagged_jet = true; }
      if( jet->auxdata<char>("HiggsMatched") ) { event_has_Higgs_matched_jet = true; }
      // Correct pairing
      if( jet->auxdata<char>("HiggsMatched") ) {
        this->fillOutputTree( m_correct_tree, *jets_selected_b_tagged.at(0), *jet );
        // ->Fill();
      // Incorrect pairing
      } else {
        this->fillOutputTree( m_incorrect_tree, *jets_selected_b_tagged.at(0), *jet );
        // ->Fill();
      }
    }
  // Case (B): not matched to a b-quark Higgs descendant
  } else {
    // Loop over all non-b jets: all are incorrect pairs
    for( auto jet : jets_selected_non_b_tagged ) {
      this->fillOutputTree( m_incorrect_tree, *jets_selected_b_tagged.at(0), *jet );
      // ->Fill();
    }
  }

  // Fill event-level counters
  if( event_has_truth_tagged_jet ) { ++m_nTruthTagged; }
  if( event_has_Higgs_matched_jet ) { ++m_nTruthMatched; }
  if( event_has_truth_tagged_jet && event_has_Higgs_matched_jet ) { ++m_nTruthTaggedAndMatched; }
  return EL::StatusCode::SUCCESS;
}


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
  if( sc != EL::StatusCode::SUCCESS ) { return sc; }

  ATH_MSG_INFO( "The sum of MC weights in this job for channel " << eventInfo()->mcChannelNumber() << " was " << m_sum_mc_weights );
  ATH_MSG_INFO( m_nPassingPreselection << " of " << m_nEvents << " events (" << (m_nEvents > 0 ? 100*m_nPassingPreselection/m_nEvents : 0) << "%) passed the HGamma pre-selection." );
  ATH_MSG_INFO( m_nPassingReco << " of " << m_nPassingPreselection << " events (" << (m_nPassingPreselection > 0 ? 100*m_nPassingReco/m_nPassingPreselection : 0) << "%) had exactly one reconstructed b-jet." );
  ATH_MSG_INFO( "How many events have a matched non-b-tagged jet using..." );
  ATH_MSG_INFO( "... HadronConeExclTruthLabelID truth-tagging?     " << m_nTruthTagged << " (" << (m_nPassingReco > 0 ? 100*m_nTruthTagged/m_nPassingReco : 0) << "%)" );
  ATH_MSG_INFO( "... using b-quarks from Higgs for truth-matching? " << m_nTruthMatched << " (" << (m_nPassingReco > 0 ? 100*m_nTruthMatched/m_nPassingReco : 0) << "%)" );
  ATH_MSG_INFO( "... truth-tagging AND truth-matching?             " << m_nTruthTaggedAndMatched << " (" << (m_nPassingReco > 0 ? 100*m_nTruthTaggedAndMatched/m_nPassingReco : 0) << "%)" );
  ATH_MSG_INFO( "This corresponded to: " << m_correct_tree->GetEntries() << " correct pairs and " << m_incorrect_tree->GetEntries() << " incorrect pairs" );
  return EL::StatusCode::SUCCESS;
}



ConstDataVector<xAOD::TruthParticleContainer> OneTagCategorisation::bQuarkHiggsDescendants( const xAOD::TruthParticleContainer *truthPtcls ) {
  // Select all final-state b-quarks with Higgs ancestors
  ConstDataVector<xAOD::TruthParticleContainer> bQuarkHiggsDescendantsFinalState(SG::VIEW_ELEMENTS);
  for( auto ptcl : *truthPtcls ) {
    if( (ptcl->status() == 23) && (fabs(ptcl->pdgId()) == 5) && HG::isFromHiggs(ptcl) ) { // outgoing ME; b-quarks; from a Higgs decay
      bQuarkHiggsDescendantsFinalState.push_back(ptcl);
    }
  }
  ATH_MSG_DEBUG( "There are " << bQuarkHiggsDescendantsFinalState.size() << " b-quarks outgoing from the hard-process which are descended from Higgs bosons" );
  return bQuarkHiggsDescendantsFinalState;
}



void OneTagCategorisation::matchQuarksToJets( ConstDataVector<xAOD::TruthParticleContainer> bQuarks, xAOD::JetContainer jets ) {
  ATH_MSG_DEBUG( "Matching between " << bQuarks.size() << " quarks and " << jets.size() << " jets" );
  // Initialise 2D-array of dR distances
  std::vector< std::vector<double> > deltaRjq( jets.size()+1, std::vector<double>(bQuarks.size()+1,99) );
  SG::AuxElement::Accessor<char> accHiggsMatched("HiggsMatched");
  // Get deltaR for each pair
  for( unsigned int idx_j = 0; idx_j < jets.size(); ++idx_j ) {
    for( unsigned int idx_q = 0; idx_q < bQuarks.size(); ++idx_q ) {
      deltaRjq.at(idx_j).at(idx_q) = HG::DR( jets.at(idx_j), bQuarks.at(idx_q) );
      deltaRjq.at(idx_j).back() = std::min( deltaRjq.at(idx_j).at(idx_q), deltaRjq.at(idx_j).back() );
      deltaRjq.back().at(idx_q) = std::min( deltaRjq.at(idx_j).at(idx_q), deltaRjq.back().at(idx_q) );
    }
  }
  // Find all jet-quark pairs that are closer to each other than either is to anything else
  for( unsigned int idx_j = 0; idx_j < jets.size(); ++idx_j ) {
    // If there are no quarks then no jets can be matched
    if( bQuarks.size() < 1 ) { accHiggsMatched( *jets.at(idx_j) ) = false; continue; }
    // Otherwise find the closest quark to each jet...
    unsigned int closest_quark_to_jet(-1);
    for( unsigned int idx_q = 0; idx_q < bQuarks.size(); ++idx_q ) {
      if( deltaRjq.at(idx_j).at(idx_q) == deltaRjq.at(idx_j).back() ) { closest_quark_to_jet = idx_q; break; }
    }
    ATH_MSG_DEBUG( "... for jet " << idx_j << " the closest quark is " << closest_quark_to_jet );
    // ... and check that this is also the closest jet to that quark
    if( deltaRjq.at(idx_j).at(closest_quark_to_jet) == deltaRjq.back().at(closest_quark_to_jet) ) {
      ATH_MSG_DEBUG( "... this is also the closest jet to that quark -> match found" );
      accHiggsMatched( *jets.at(idx_j) ) = true;
    } else {
      ATH_MSG_DEBUG( "... but this is not the closest jet to that quark -> no match" );
      accHiggsMatched( *jets.at(idx_j) ) = false;
    }
  }
  return;
}



void OneTagCategorisation::decorateWithIndices( const xAOD::Jet& bjet, xAOD::JetContainer& nonbjets ) {
  // Calculate m_jb
  SG::AuxElement::Accessor<double> accMjb("m_jb");
  for( auto jet : nonbjets ) { accMjb(*jet) = (bjet.p4() + jet->p4()).M() / HG::GeV; }

  // Sort by distance from mH
  std::sort( nonbjets.begin(), nonbjets.end(), [](const xAOD::Jet* i, const xAOD::Jet* j)
    { return fabs( i->auxdata<double>("m_jb") - 125.04 ) < fabs( j->auxdata<double>("m_jb") - 125.04 ); }
  );
  // Add indexing by distance from mH
  SG::AuxElement::Accessor<int> accIdxByMh("idx_by_mH");
  for( unsigned int idx = 0; idx < nonbjets.size(); ++idx ) { accIdxByMh( *nonbjets.at(idx) ) = idx; }

  // Sort by distance from pT
  std::sort( nonbjets.begin(), nonbjets.end(), [](const xAOD::Jet* i, const xAOD::Jet* j)
    { return i->pt() > j->pt(); }
  );
  // Add indexing by pT
  SG::AuxElement::Accessor<int> accIdxByPt("idx_by_pT");
  for( unsigned int idx = 0; idx < nonbjets.size(); ++idx ) { accIdxByPt( *nonbjets.at(idx) ) = idx; }
}



void OneTagCategorisation::fillOutputTree( TTree* outputTree, const xAOD::Jet& bjet, const xAOD::Jet& otherjet ) {
  TLorentzVector jb_p4 = bjet.p4() + otherjet.p4();
  m_m_jb = jb_p4.M() / HG::GeV;
  m_pT_jb = jb_p4.Pt() / HG::GeV;
  m_abs_eta_jb = fabs( jb_p4.Eta() );
  m_Delta_eta_jb = fabs( bjet.eta() - otherjet.eta() );
  m_Delta_phi_jb = fabs( bjet.p4().DeltaPhi( otherjet.p4() ) );
  m_idx_by_mH = otherjet.auxdata<int>("idx_by_mH");
  m_idx_by_pT = otherjet.auxdata<int>("idx_by_pT");
  m_pT_j = otherjet.pt() / HG::GeV;
  m_abs_eta_j = fabs( otherjet.eta() );
  m_MV2c20_FCBE_60 = otherjet.auxdata<char>("MV2c20_FixedCutBEff_60");
  m_MV2c20_FCBE_70 = otherjet.auxdata<char>("MV2c20_FixedCutBEff_70");
  m_MV2c20_FCBE_77 = otherjet.auxdata<char>("MV2c20_FixedCutBEff_77");
  m_MV2c20_FCBE_85 = otherjet.auxdata<char>("MV2c20_FixedCutBEff_85");
  ATH_MSG_DEBUG(  "  Jet with mH: " << m_m_jb << " and pT: " << m_pT_j << " is " << m_idx_by_mH << " by mH and " << m_idx_by_pT << " by pT" );
  outputTree->Fill();
  return;
}



double OneTagCategorisation::sumOfWeights( int mcID ) {
  if( mcID == 341061 ) { return 197000; } // SM yybb
  if( mcID == 341062 ) { return 200000; } // SM yjbb
  if( mcID == 341063 ) { return 196000; } // SM yybj
  if( mcID == 341064 ) { return 200000; } // SM yjjb
  if( mcID == 341065 ) { return 180000; } // SM yyjj
  if( mcID == 341066 ) { return 198000; } // SM yjjj
  if( mcID == 341559 ) { return 100000; } // SM hh->yybb
  if( mcID == 341173 ) { return 100000; } // X275->hh->yybb
  if( mcID == 341004 ) { return 100000; } // X300->hh->yybb
  if( mcID == 341174 ) { return 100000; } // X325->hh->yybb
  if( mcID == 341175 ) { return 100000; } // X350->hh->yybb
  if( mcID == 341176 ) { return 100000; } // X400->hh->yybb
  return 1.0;
}
