#pragma once
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include <map>
#include <vector>

class OneTagCategorisation : public HgammaAnalysis
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)

private:
  /// User-configurables
  /// b-tagging working point
  std::string m_1_tag_WP, m_2_tag_WP;

  /// Not user-configurable
  /// Output ntuples
  TTree *m_event_tree_1tag; //!
  TTree *m_event_tree_2tag; //!

  /// Output discriminator vectors
  std::vector<double> m_v_abs_eta_j; //!
  std::vector<double> m_v_abs_eta_jb; //!
  std::vector<double> m_v_Delta_eta_jb; //!
  std::vector<double> m_v_Delta_phi_jb; //!
  std::vector<int> m_v_idx_by_mH; //!
  std::vector<int> m_v_idx_by_m_jb; //!
  std::vector<int> m_v_idx_by_pT; //!
  std::vector<double> m_v_m_jb; //!
  std::vector<double> m_v_pT_j; //!
  std::vector<double> m_v_pT_jb; //!
  std::vector<double> m_v_isCorrect; //!

  /// Event weights
  double m_event_weight; //!
  double m_sum_mc_weights; //!

  /// Cutflow counters
  std::map< std::string, unsigned int > m_cutFlow; //!

  // /// Decorate with pT and mH indices
  // void decorateWithIndices( const xAOD::Jet& bjet, xAOD::JetContainer& nonbjets );

  /// Calculate quantities for output trees
  void appendToOutput( const bool& isCorrect, const xAOD::Jet& bjet, const xAOD::Jet& otherjet );


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
