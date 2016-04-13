#pragma once
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include <vector>

class OneTagCategorisation : public HgammaAnalysis
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)

private:
  /// Output ntuples
  TTree *m_correct_tree; //!
  TTree *m_incorrect_tree; //!

  /// b-tagging working point
  std::string m_b_tagging_WP; //!

  /// Output discriminators
  double m_m_jb, m_pT_jb, m_eta_jb;
  double m_Delta_eta_jb, m_Delta_phi_jb;
  double m_pT_j, m_eta_j;
  int m_MV2c20_FCBE_60, m_MV2c20_FCBE_70, m_MV2c20_FCBE_77, m_MV2c20_FCBE_85;

  /// Event weights
  double m_event_weight, m_sum_mc_weights;

  /// Cutflow counters
  int m_nEvents, m_nPassingReco, m_nTruthTagged, m_nTruthMatched, m_nTruthTaggedAndMatched;

  /// Get final-state b-quarks
  ConstDataVector<xAOD::TruthParticleContainer> bQuarkHiggsDescendants( const xAOD::TruthParticleContainer *truthPtcls );

  /// Get final-state b-quarks
  void matchQuarksToJets( ConstDataVector<xAOD::TruthParticleContainer> bQuarks, xAOD::JetContainer jets );

  /// Calculate quantities for output trees
  void calculateOutputQuantities( const xAOD::Jet& bjet, const xAOD::Jet& otherjet );

  /// Retrieve sum of MC weights
  double sumOfWeights( int mcID );


public:
  // this is a standard constructor
  OneTagCategorisation() { }
  OneTagCategorisation(const char *name);
  virtual ~OneTagCategorisation();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode initialize();
  virtual EL::StatusCode execute();
  virtual EL::StatusCode finalize();

  // this is needed to distribute the algorithm to the workers
  ClassDef(OneTagCategorisation, 1);
};
