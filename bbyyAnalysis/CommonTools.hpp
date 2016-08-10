#pragma once
#include "HGamAnalysisFramework/HHyybbTool.h"
#include <AsgTools/MsgStream.h>
#include <AsgTools/MsgStreamMacros.h>
#include <AsgTools/MessageCheck.h>
#include <cassert>

using namespace asg::msgUserCode;

namespace CommonTools {
  /**
   * Decorate muon-corrected value
   * @return nothing
   */
  inline void correctForMuons(HG::HHyybbTool* yybbTool, xAOD::JetContainer& jets, const xAOD::MuonContainer& muons) {
    yybbTool->decorateMuonCorrections(jets, muons);
    std::vector<TLorentzVector> jets_muon_corrected = yybbTool->getMuonCorrJet4Vs(jets);
    assert( jets_muon_corrected.size() == jets.size() );
    // Define decorations
    SG::AuxElement::Accessor<double> acc_E("muon_E");
    SG::AuxElement::Accessor<double> acc_pT("muon_pT");
    SG::AuxElement::Accessor<double> acc_eta("muon_eta");
    SG::AuxElement::Accessor<double> acc_phi("muon_phi");
    // Add decorations
    for( unsigned int idx = 0; idx < jets.size(); ++idx ) {
      acc_E(*jets.at(idx)) = jets_muon_corrected.at(idx).E();
      acc_pT(*jets.at(idx)) = jets_muon_corrected.at(idx).Pt();
      acc_eta(*jets.at(idx)) = jets_muon_corrected.at(idx).Eta();
      acc_phi(*jets.at(idx)) = jets_muon_corrected.at(idx).Phi();
    }
  }


  /**
   * Construct muon-corrected TLorentzVector
   * @return TLorentzVector of muon-corrected kinematics
   */
  inline TLorentzVector p4(const xAOD::Jet& jet) {
    TLorentzVector output;
    output.SetPtEtaPhiE( jet.auxdata<double>("muon_pT"), jet.auxdata<double>("muon_eta"), jet.auxdata<double>("muon_phi"), jet.auxdata<double>("muon_E") );
    return output;
  }

  /**
   * Return final Higgs before decay to bb
   * @return Higgs if matching was successful, otherwise null pointer
   */
  inline const xAOD::TruthParticle* HbbBeforeDecay( const xAOD::TruthParticleContainer* truthPtcls ) {
    // Look for Higgses among truth particles
    ConstDataVector<xAOD::TruthParticleContainer> higgsBosons = HG::getFinalHiggsBosons(truthPtcls);

    // If at least one is found then look for one decaying to bb
    if( higgsBosons.size() > 0 ) {
      ATH_MSG_DEBUG("Found " << higgsBosons.size() << " Higgses in the last stage before decay");
      for( const auto& higgsBoson : higgsBosons ) {
        for( unsigned int iChild = 0; iChild < higgsBoson->nChildren(); ++iChild ) {
          if(fabs(higgsBoson->child(iChild)->pdgId()) == 5) {
            ATH_MSG_DEBUG("... identified b-quark decay");
            return higgsBoson;
          }
        }
      }
    }
    // Otherwise return null pointer
    ATH_MSG_DEBUG("No Higgses were found decaying to b-quarks!");
    return 0;
  }


  /**
   * Return best matched jet-pair to a given Higgs
   * @return true if matching was successful, otherwise false
   */
  inline bool matchJetsToHiggs( xAOD::JetContainer jets, const xAOD::TruthParticle* higgs ) {
    // Check for null input
    if( higgs == 0 ) { return false; }
    ATH_MSG_DEBUG("Matching " << jets.size() << " jets to a Higgs boson");

    // Initialise decorator and Higgs Lorentz vector
    SG::AuxElement::Accessor<char> accHiggsMatched("HiggsMatched");
    TLorentzVector higgs_p4 = higgs->p4();
    ATH_MSG_DEBUG("Higgs has 4-vector of E = " << higgs_p4.E()/HG::GeV << ", pT = " << higgs_p4.Pt()/HG::GeV << ", eta = " << higgs_p4.Eta() << ", phi = " << higgs_p4.Phi());

    // Initialise counters for current best pairing
    int matched_j1(-1), matched_j2(-1);
    double matched_deltaR(99);

    // Iterate over all jets
    for (unsigned int idx_j1 = 0; idx_j1 < jets.size(); ++idx_j1) {
      // Set initial value of decoration to false
      accHiggsMatched(*jets.at(idx_j1)) = false;
      // Iterate over all possible pairs
      for (unsigned int idx_j2 = idx_j1+1; idx_j2 < jets.size(); ++idx_j2) {
        // Compare DeltaR wrt the Higgs with the best one found so far
        TLorentzVector jj_p4 = p4(*jets.at(idx_j1)) + p4(*jets.at(idx_j2));
        const double deltaR( higgs_p4.DeltaR(jj_p4) );
        ATH_MSG_DEBUG("  ... jet pair has DeltaR wrt Higgs of " << deltaR << " from 4-vector of E = " << jj_p4.E()/HG::GeV << ", pT = " << jj_p4.Pt()/HG::GeV << ", eta = " << jj_p4.Eta() << ", phi = " << jj_p4.Phi() );
        if( deltaR < matched_deltaR ) {
          matched_deltaR = deltaR;
          matched_j1 = idx_j1;
          matched_j2 = idx_j2;
        }
      }
    }

    if( matched_deltaR < 0.4 ) { // 0.6 is 82% efficient; 0.5 is 79% efficient
      ATH_MSG_DEBUG( "=> tagging best matched jets: DeltaR = " << matched_deltaR);
      accHiggsMatched(*jets.at(matched_j1)) = true;
      accHiggsMatched(*jets.at(matched_j2)) = true;
      return true;
    }

    ATH_MSG_DEBUG( "=> no jet-pair candidate had a DeltaR of less than 0.4! Best was " << matched_deltaR);
    return false;
  }
}
