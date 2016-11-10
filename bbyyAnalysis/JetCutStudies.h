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
  TMVA::Reader m_reader_low_mass_with_booleans, m_reader_high_mass_with_booleans; //!
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
  int m_tag_category; //!
  // 2-tag
  double m_jet_pT1_2tag; //!
  double m_jet_pT2_2tag; //!
  double m_m_jj_2tag; //!
  double m_m_yyjj_2tag; //!
  // 1-tag low mass with booleans
  double m_score_1tag_low_mass_with_booleans; //!
  double m_jet_pT1_1tag_low_mass_with_booleans; //!
  double m_jet_pT2_1tag_low_mass_with_booleans; //!
  double m_m_jj_1tag_low_mass_with_booleans; //!
  double m_m_yyjj_1tag_low_mass_with_booleans; //!
  // 1-tag high mass with booleans
  double m_score_1tag_high_mass_with_booleans; //!
  double m_jet_pT1_1tag_high_mass_with_booleans; //!
  double m_jet_pT2_1tag_high_mass_with_booleans; //!
  double m_m_jj_1tag_high_mass_with_booleans; //!
  double m_m_yyjj_1tag_high_mass_with_booleans; //!
  // 1-tag low mass without booleans
  double m_score_1tag_low_mass_without_booleans; //!
  double m_jet_pT1_1tag_low_mass_without_booleans; //!
  double m_jet_pT2_1tag_low_mass_without_booleans; //!
  double m_m_jj_1tag_low_mass_without_booleans; //!
  double m_m_yyjj_1tag_low_mass_without_booleans; //!
  // 1-tag high mass without booleans
  double m_score_1tag_high_mass_without_booleans; //!
  double m_jet_pT1_1tag_high_mass_without_booleans; //!
  double m_jet_pT2_1tag_high_mass_without_booleans; //!
  double m_m_jj_1tag_high_mass_without_booleans; //!
  double m_m_yyjj_1tag_high_mass_without_booleans; //!
  // 1-tag low mass without booleans with cut on 85% WP
  double m_score_1tag_low_mass_without_booleans_with_cut85; //!
  double m_jet_pT1_1tag_low_mass_without_booleans_with_cut85; //!
  double m_jet_pT2_1tag_low_mass_without_booleans_with_cut85; //!
  double m_m_jj_1tag_low_mass_without_booleans_with_cut85; //!
  double m_m_yyjj_1tag_low_mass_without_booleans_with_cut85; //!
  // 1-tag high mass without booleans with cut on 85% WP
  double m_score_1tag_high_mass_without_booleans_with_cut85; //!
  double m_jet_pT1_1tag_high_mass_without_booleans_with_cut85; //!
  double m_jet_pT2_1tag_high_mass_without_booleans_with_cut85; //!
  double m_m_jj_1tag_high_mass_without_booleans_with_cut85; //!
  double m_m_yyjj_1tag_high_mass_without_booleans_with_cut85; //!

  /// Event weights
  double m_event_weight; //!
  double m_sum_mc_weights, m_sum_pileup_weights; //!

  /// Cutflow counters
  std::map< std::string, unsigned int > m_cutFlow; //!


  void decorateWithClassifiers( const xAOD::Jet& bjet, xAOD::JetContainer& nonbjets );

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
