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

OneTagCategorisation::OneTagCategorisation(const char *name)
: HgammaAnalysis(name)
, m_correct_tree(0)
, m_incorrect_tree(0)
, m_b_tagging_WP("")
, m_truth_jet_quark_dR(0)
, m_sum_mc_weights(0)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
  // MsgStream::setName(name);
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

  m_truth_jet_quark_dR = config()->getNum("OneTagCategorisation.TruthJetQuarkDR", 0.2);
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
  m_correct_tree->Branch( "m_jb",             &m_m_jb );
  m_correct_tree->Branch( "pT_jb",            &m_pT_jb );
  m_correct_tree->Branch( "eta_jb",           &m_eta_jb );
  m_correct_tree->Branch( "Delta_eta_jb",     &m_Delta_eta_jb );
  m_correct_tree->Branch( "Delta_phi_jb",     &m_Delta_phi_jb );
  m_correct_tree->Branch( "pT_j",             &m_pT_j );
  m_correct_tree->Branch( "eta_j",            &m_eta_j );
  m_correct_tree->Branch( "event_weight",     &m_event_weight );
  if( m_b_tagging_WP != "MV2c20_FixedCutBEff_60" ) { m_correct_tree->Branch( "MV2c20_FCBE_60", &m_MV2c20_FCBE_60 ); }
  if( m_b_tagging_WP != "MV2c20_FixedCutBEff_70" ) { m_correct_tree->Branch( "MV2c20_FCBE_70", &m_MV2c20_FCBE_70 ); }
  if( m_b_tagging_WP != "MV2c20_FixedCutBEff_77" ) { m_correct_tree->Branch( "MV2c20_FCBE_77", &m_MV2c20_FCBE_77 ); }
  if( m_b_tagging_WP != "MV2c20_FixedCutBEff_85" ) { m_correct_tree->Branch( "MV2c20_FCBE_85", &m_MV2c20_FCBE_85 ); }
  // Add incorrect choice tree to output file
  m_incorrect_tree = new TTree("incorrect","incorrect");
  m_incorrect_tree->SetDirectory(file);
  m_incorrect_tree->Branch( "m_jb",             &m_m_jb );
  m_incorrect_tree->Branch( "pT_jb",            &m_pT_jb );
  m_incorrect_tree->Branch( "eta_jb",           &m_eta_jb );
  m_incorrect_tree->Branch( "Delta_eta_jb",     &m_Delta_eta_jb );
  m_incorrect_tree->Branch( "Delta_phi_jb",     &m_Delta_phi_jb );
  m_incorrect_tree->Branch( "pT_j",             &m_pT_j );
  m_incorrect_tree->Branch( "eta_j",            &m_eta_j );
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

  /// Get overall event weight, normalised to 1fb-1
  int mcID = eventInfo()->mcChannelNumber();
  m_event_weight = eventHandler()->mcWeight() * ( 1e3 * HgammaAnalysis::getCrossSection(mcID)) * HgammaAnalysis::getGeneratorEfficiency(mcID) * HgammaAnalysis::getKFactor(mcID) / sumOfWeights(mcID);

  //___________________________________________________________________________________________
  // Fetch default jets
  xAOD::JetContainer jets_corrected = jetHandler()->getCorrectedContainer();
  xAOD::JetContainer jets_selected = jetHandler()->applySelection(jets_corrected);
  xAOD::JetContainer jets_selected_truth_and_tagged( SG::VIEW_ELEMENTS );
  xAOD::JetContainer jets_selected_truth_not_tagged( SG::VIEW_ELEMENTS );
  xAOD::JetContainer jets_selected_not_truth_b( SG::VIEW_ELEMENTS );

  // Retrieve the truth particles
  const xAOD::TruthParticleContainer *truthPtcls = 0;
  EL_CHECK( "execute()", event()->retrieve(truthPtcls, "TruthParticles") );

  // Select final-state b-quarks descended from Higgs bosons
  ConstDataVector<xAOD::TruthParticleContainer> bQuarksSelected = bQuarkHiggsDescendants( truthPtcls );

  //Loop jets and apply muon correction and seperate bjet / untagged jets
  for( auto jet : jets_selected ) {
    bool truth_tagged( jet->auxdata<int>("HadronConeExclTruthLabelID") == 5 );
    bool reco_tagged( jet->auxdata<char>("MV2c20_FixedCutBEff_60") );
    bool truth_matched( HG::minDRrap( jet, bQuarksSelected ) < m_truth_jet_quark_dR );
    if( truth_tagged != truth_matched ) {
      std::cout << "Mismatch: minDRrap: " << HG::minDRrap( jet, bQuarksSelected ) << std::boolalpha << " tagged? " << truth_tagged << ", matched? " << truth_matched << std::endl;
    }
    // Categorise jets
    if( truth_tagged && truth_matched ) {
      if( reco_tagged ) {
        // Truth b-jets from Higgs correctly tagged
        jets_selected_truth_and_tagged.push_back( jet );
      } else {
        // Truth b-jets from Higgs incorrectly tagged
        jets_selected_truth_not_tagged.push_back( jet );
      }
    } else {
      jets_selected_not_truth_b.push_back( jet );
    }
  }
  ATH_MSG_DEBUG( "There are " << jets_selected.size() << " jets, of which "
                              << jets_selected_truth_and_tagged.size() << " are correctly tagged b-jets, "
                              << jets_selected_truth_not_tagged.size() << " are incorrectly tagged b-jets and "
                              << jets_selected_not_truth_b.size() << " are not truth b-jets" );


  // Want events with 2 truth-matched b-jets but only one reco-tagged b-jet
  if( jets_selected_truth_and_tagged.size() == 1 && jets_selected_truth_not_tagged.size() == 1 ) {
    ATH_MSG_DEBUG( "  pT of tagged jet " << jets_selected_truth_and_tagged.at(0)->pt() / HG::GeV );
    ATH_MSG_DEBUG( "  pT of untagged jet " << jets_selected_truth_not_tagged.at(0)->pt() / HG::GeV );

    // Correct pairing
    this->fillVariables( *jets_selected_truth_and_tagged.at(0), *jets_selected_truth_not_tagged.at(0) );
    m_correct_tree->Fill();
    // Incorrect pairing
    for( const auto jet : jets_selected_not_truth_b ) {
      this->fillVariables( *jets_selected_truth_and_tagged.at(0), *jet );
      m_incorrect_tree->Fill();
    }
  }

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

  ATH_MSG_INFO( "Output trees contain: correct " << m_correct_tree->GetEntries() << " incorrect " << m_incorrect_tree->GetEntries() << " pairs" );
  ATH_MSG_INFO( "Sum of MC weights in this job for channel " << eventInfo()->mcChannelNumber() << " is " << m_sum_mc_weights );
  return EL::StatusCode::SUCCESS;
}


ConstDataVector<xAOD::TruthParticleContainer> OneTagCategorisation::bQuarkHiggsDescendants( const xAOD::TruthParticleContainer *truthPtcls ) {
  // Select all b-quarks with Higgs parents
  ConstDataVector<xAOD::TruthParticleContainer> bQuarkHiggsChildren(SG::VIEW_ELEMENTS);
  for( auto ptcl : *truthPtcls ) {
    if( fabs(ptcl->pdgId()) == 5 ) {
      for( unsigned int iParent = 0; iParent < ptcl->nParents(); ++iParent ) {
        if( ptcl->parent(iParent)->isHiggs() ) { bQuarkHiggsChildren.push_back( ptcl ); }
      }
    }
  }
  ATH_MSG_DEBUG( "There are " << bQuarkHiggsChildren.size() << " b-quarks with direct Higgs parents in this event" );

  // Set up two queues of particles to process and a collection of final-state b-quarks descended from Higgses
  ConstDataVector<xAOD::TruthParticleContainer> particlesBeingProcessed( bQuarkHiggsChildren );
  ConstDataVector<xAOD::TruthParticleContainer> particlesInQueue(SG::VIEW_ELEMENTS);
  ConstDataVector<xAOD::TruthParticleContainer> bQuarkHiggsDescendantsFinalState(SG::VIEW_ELEMENTS);

  // Iterate over particles descended from Higgs bosons while any remain in the queue
  while( particlesBeingProcessed.size() > 0 ) {
    for( auto ptcl : particlesBeingProcessed ) {
      // If the particle is a final-state b-quark then save it
      if( fabs(ptcl->pdgId()) == 5 && ptcl->status() == 23 ) {
        bQuarkHiggsDescendantsFinalState.push_back( ptcl );
      // Otherwise add all of its b-quark descendents to the queue
      } else {
        for( unsigned int iChild = 0; iChild < ptcl->nChildren(); ++iChild ) {
          if( fabs(ptcl->child(iChild)->pdgId()) == 5 ) { particlesInQueue.push_back( ptcl->child(iChild) ); }
        }
      }
    }
    // Move the queued particles to be processed next
    // std::cout << "Old sizes: " << particlesBeingProcessed.size() << " done, " << particlesInQueue.size() << " to be done" << std::endl;
    particlesBeingProcessed.swap( particlesInQueue );
    particlesInQueue.clear();
    // std::cout << "New sizes: " << particlesBeingProcessed.size() << " to be done, " << particlesInQueue.size() << " remaining in queue" << std::endl;
  }
  ATH_MSG_DEBUG( "There are " << bQuarkHiggsDescendantsFinalState.size() << " final-state b-quarks which are descended from Higgs bosons" );
  return bQuarkHiggsDescendantsFinalState;
}


double OneTagCategorisation::sumOfWeights( int mcID ) {
  if( mcID == 341063 ) { return 196000; } // SM yybj
  if( mcID == 341559 ) { return 100000; } // SM hh->yybb
  if( mcID == 341173 ) { return 100000; } // X275->hh->yybb
  if( mcID == 341004 ) { return 100000; } // X300->hh->yybb
  if( mcID == 341174 ) { return 100000; } // X325->hh->yybb
  if( mcID == 341175 ) { return 100000; } // X350->hh->yybb
  if( mcID == 341176 ) { return 100000; } // X400->hh->yybb
  return 1.0;
}

void OneTagCategorisation::fillVariables( const xAOD::Jet& bjet, const xAOD::Jet& otherjet ) {
  TLorentzVector tlv_jb = bjet.p4() + otherjet.p4();
  m_m_jb = tlv_jb.M() / HG::GeV;
  m_pT_jb = tlv_jb.Pt() / HG::GeV;
  m_eta_jb = tlv_jb.Eta();
  m_Delta_eta_jb = fabs( bjet.eta() - otherjet.eta() );
  m_Delta_phi_jb = bjet.p4().DeltaPhi( otherjet.p4() );
  m_pT_j = otherjet.pt() / HG::GeV;
  m_eta_j = otherjet.eta();
  m_MV2c20_FCBE_60 = otherjet.auxdata<char>("MV2c20_FixedCutBEff_60");
  m_MV2c20_FCBE_70 = otherjet.auxdata<char>("MV2c20_FixedCutBEff_70");
  m_MV2c20_FCBE_77 = otherjet.auxdata<char>("MV2c20_FixedCutBEff_77");
  m_MV2c20_FCBE_85 = otherjet.auxdata<char>("MV2c20_FixedCutBEff_85");
  return;
}
