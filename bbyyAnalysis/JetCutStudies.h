#pragma once
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include <TMVA/Reader.h>
#include <map>
#include <vector>

class JetCutStudies : public HgammaAnalysis
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
  TTree *m_event_tree; //!

  /// TMVA reader from scikit-learn XML output
  TMVA::Reader m_reader_low_mass, m_reader_high_mass; //!
  TMVA::Reader m_reader_low_mass_without_booleans, m_reader_high_mass_without_booleans; //!
  /// TMVA discriminators
  float m_abs_eta_j; //!
  float m_abs_eta_jb; //!
  float m_Delta_eta_jb; //!
  float m_idx_by_mH; //!
  float m_idx_by_pT; //!
  float m_idx_by_pT_jb; //!
  float m_passes_WP77; //!
  float m_passes_WP85; //!
  float m_m_jb; //!
  float m_pT_j; //!
  float m_pT_jb; //!

  /// Output variables
  double m_m_yy; //!
  unsigned int m_photon_n, m_jet_n; //!
  std::vector<double> m_photon_pT, m_jet_pT; //!
  std::vector<double> m_photon_eta, m_jet_eta, m_jet_eta_det; //!
  std::vector<double> m_photon_phi, m_jet_phi; //!
  std::vector<double> m_photon_E, m_jet_E; //!
  std::vector<int> m_photon_isTight; //!
  std::vector<int> m_jet_btag_1tag; //!
  std::vector<int> m_jet_btag_2tag; //!
  std::vector<int> m_jet_btag_85; //!
  std::vector<double> m_jet_classifier_low_mass; //!
  std::vector<double> m_jet_classifier_high_mass; //!
  std::vector<double> m_jet_classifier_low_mass_without_booleans; //!
  std::vector<double> m_jet_classifier_high_mass_without_booleans; //!
  std::vector<int> m_jet_truth_tag; //!
  std::vector<int> m_jet_higgs_match; //!
  std::vector<double> m_jet_m_jb; //!

  /// Event weights
  double m_event_weight; //!
  double m_sum_mc_weights, m_sum_pileup_weights; //!

  /// Cutflow counters
  std::map< std::string, unsigned int > m_cutFlow; //!


  void decorateWithClassifier( const xAOD::Jet& bjet, xAOD::JetContainer& nonbjets );

public:
  // this is a standard constructor
  JetCutStudies() { }
  JetCutStudies(const char *name);
  virtual ~JetCutStudies();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode initialize();
  virtual EL::StatusCode execute();
  virtual EL::StatusCode finalize();

  // this is needed to distribute the algorithm to the workers
  ClassDef(JetCutStudies, 1);
};
