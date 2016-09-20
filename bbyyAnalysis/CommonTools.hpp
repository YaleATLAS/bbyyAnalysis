#pragma once
#include "HGamAnalysisFramework/HHyybbTool.h"
#include <AsgTools/MsgStream.h>
#include <AsgTools/MsgStreamMacros.h>
#include <AsgTools/MessageCheck.h>
#include <cassert>

using namespace asg::msgUserCode;

namespace CommonTools {
  /**
   * Expected luminosity for Moriond 2017
   * @return double
   */
  inline double luminosity_invfb() {
    return 40.0;
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
   * Decorate muon-corrected value
   * @return addition of "muon_*" decoration
   */
  inline void decorateMuonCorrection(HG::HHyybbTool* yybbTool, xAOD::JetContainer& jets, const xAOD::MuonContainer& muons) {
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
    ATH_MSG_DEBUG("Decorated " << jets.size() << " jets with information from " << muons.size() << " muons.");
  }



  /**
   * Decorate jet pair with closest DeltaR match to an identified Higgs to bb decay
   * @return addition of "HiggsMatched" decoration
   */
  inline bool decorateHiggsMatching( xAOD::JetContainer& jets, xAOD::TruthParticleContainer& higgsBosons ) {
    // Initialise decorator to false
    SG::AuxElement::Accessor<char> accHiggsMatched("HiggsMatched");
    for( const auto& jet : jets ) { accHiggsMatched(*jet) = false; }

    // Look for Higgses decaying to bb
    xAOD::TruthParticleContainer higgsBosonsDecayingToBB(SG::VIEW_ELEMENTS);
    for (const auto& higgsBoson : higgsBosons) {
      bool decayingToBB = false;
      for (unsigned int iChild = 0; iChild < higgsBoson->nChildren(); ++iChild) {
        if (fabs(higgsBoson->child(iChild)->pdgId()) == 5) { decayingToBB = true; }
      }
      if (decayingToBB) { higgsBosonsDecayingToBB.push_back(higgsBoson); }
    }
    ATH_MSG_DEBUG("Found " << higgsBosons.size() << " Higgs bosons with " << higgsBosonsDecayingToBB.size() << " decaying to bb");

    // Return false if no Higgs are found decaying to bb
    if( higgsBosonsDecayingToBB.size() == 0 ) { return false; }

    // Construct the Higgs Lorentz vector
    TLorentzVector higgs_p4 = higgsBosonsDecayingToBB.at(0)->p4();
    ATH_MSG_DEBUG("Matching " << jets.size() << " jets to a Higgs with: E = " << higgs_p4.E()/HG::GeV << ", pT = " << higgs_p4.Pt()/HG::GeV << ", eta = " << higgs_p4.Eta() << ", phi = " << higgs_p4.Phi());

    // Initialise counters for current best pairing
    int matched_j1(-1), matched_j2(-1);
    double matched_deltaR(99);

    // Iterate over all possible jet pairs
    for (unsigned int idx_j1 = 0; idx_j1 < jets.size(); ++idx_j1) {
      for (unsigned int idx_j2 = idx_j1+1; idx_j2 < jets.size(); ++idx_j2) {
        // Compare DeltaR wrt the Higgs with the best one found so far
        TLorentzVector jj_p4 = p4(*jets.at(idx_j1)) + p4(*jets.at(idx_j2));
        const double deltaR( higgs_p4.DeltaR(jj_p4) );
        ATH_MSG_DEBUG("  .... jet pair with 4-vector of E = " << jj_p4.E()/HG::GeV << ", pT = " << jj_p4.Pt()/HG::GeV << ", eta = " << jj_p4.Eta() << ", phi = " << jj_p4.Phi() << " has DeltaR wrt Higgs of " << deltaR);
        if( deltaR < matched_deltaR ) {
          matched_deltaR = deltaR;
          matched_j1 = idx_j1;
          matched_j2 = idx_j2;
        }
      }
    }

    if( matched_deltaR < 0.5 ) { // 0.6 is 82% efficient; 0.5 is 79% efficient; 0.4 is 47% efficient
      ATH_MSG_DEBUG( "=> tagging best matched jets: DeltaR = " << matched_deltaR);
      accHiggsMatched(*jets.at(matched_j1)) = true;
      accHiggsMatched(*jets.at(matched_j2)) = true;
      return true;
    }

    ATH_MSG_DEBUG( "=> no jet-pair candidate had a DeltaR of less than 0.4! Best was " << matched_deltaR);
    return false;
  }


  /**
   * Decorate non-b-jets with their order in pT or distance from mH
   * @return addition of "m_jb", "idx_by_mH" and "idx_by_pT" decorations
   */
  inline void decorateWithIndices(const xAOD::Jet& bjet, xAOD::JetContainer& nonbjets) {
    // Calculate m_jb and pT_jb
    SG::AuxElement::Accessor<double> accMjb("m_jb");
    SG::AuxElement::Accessor<double> accpTjb("pT_jb");
    for (auto jet : nonbjets) {
      TLorentzVector jb_p4 = CommonTools::p4(bjet) + CommonTools::p4(*jet);
      accMjb(*jet) = jb_p4.M() / HG::GeV;
      accpTjb(*jet) = jb_p4.Pt() / HG::GeV;
    }

    // Sort by pT_jb and add index
    SG::AuxElement::Accessor<int> accIdxBypTjb("idx_by_pT_jb");
    std::sort(nonbjets.begin(), nonbjets.end(), [](const xAOD::Jet *i, const xAOD::Jet *j) { return i->auxdata<double>("pT_jb") > j->auxdata<double>("pT_jb"); });
    for (unsigned int idx = 0; idx < nonbjets.size(); ++idx) { accIdxBypTjb(*nonbjets.at(idx)) = idx; }

    // Sort by distance from mH and add index
    SG::AuxElement::Accessor<int> accIdxByMh("idx_by_mH");
    std::sort(nonbjets.begin(), nonbjets.end(), [](const xAOD::Jet *i, const xAOD::Jet *j) { return fabs(i->auxdata<double>("m_jb") - 125.09) < fabs(j->auxdata<double>("m_jb") - 125.09); });
    for (unsigned int idx = 0; idx < nonbjets.size(); ++idx) { accIdxByMh(*nonbjets.at(idx)) = idx; }

    // Sort by pT and add index
    SG::AuxElement::Accessor<int> accIdxByPt("idx_by_pT");
    std::sort(nonbjets.begin(), nonbjets.end(), [](const xAOD::Jet *i, const xAOD::Jet *j) { return i->pt() > j->pt(); });
    for (unsigned int idx = 0; idx < nonbjets.size(); ++idx) { accIdxByPt(*nonbjets.at(idx)) = idx; }
  }

  /**
   * Return the sample cross-section
   * @return sample cross-section in fb
   */
  inline double xs_fb(const int& mcID, const double& default_pb, const bool& scaleBSM = false) {
    // Use Moriond 2016 limits
    double xs_fb = 1e3 * default_pb;
    if (mcID == 341173) { xs_fb = 12.89 * (scaleBSM ? 0.2 : 1.0); } // X275->hh->yybb
    if (mcID == 341004) { xs_fb = 12.89 * (scaleBSM ? 0.2 : 1.0); } // X300->hh->yybb
    if (mcID == 341174) { xs_fb = 12.89 * (scaleBSM ? 0.2 : 1.0); } // X325->hh->yybb
    if (mcID == 341175) { xs_fb = 12.89 * (scaleBSM ? 0.2 : 1.0); } // X350->hh->yybb
    if (mcID == 341176) { xs_fb = 12.89 * (scaleBSM ? 0.2 : 1.0); } // X400->hh->yybb
    if (mcID == 342620) { xs_fb = 12.89; } // SM NLO hh->yybb
    return xs_fb; //sherpa default_pb is 4.0127E+001
  }

  /**
   * Return the sum of MC event weights in each sample
   * @return sum of MC weights in the sample
   */
  inline double sumOfWeights(const int& mcID) {
    if (mcID == 341061) { return 197000; }     // SM yybb
    if (mcID == 341062) { return 200000; }     // SM yjbb
    if (mcID == 341063) { return 196000; }     // SM yybj
    if (mcID == 341064) { return 200000; }     // SM yjjb
    if (mcID == 341065) { return 180000; }     // SM yyjj
    if (mcID == 341066) { return 198000; }     // SM yjjj
    if (mcID == 341559) { return 100000; }     // SM LO hh->yybb
    if (mcID == 341173) { return 100000; }     // X275->hh->yybb
    if (mcID == 341004) { return 100000; }     // X300->hh->yybb
    if (mcID == 341174) { return 100000; }     // X325->hh->yybb
    if (mcID == 341175) { return 100000; }     // X350->hh->yybb
    if (mcID == 341176) { return 100000; }     // X400->hh->yybb
    if (mcID == 341939) { return 80604609.6; } // Sherpa photons+jets
    if (mcID == 342620) { return 2134.103; }   // SM NLO hh->yybb (from 200000 events)
    return 1.0;
  }
}
