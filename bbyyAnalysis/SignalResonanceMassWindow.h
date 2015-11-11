/**
 * @file   SignalResonanceMassWindow.h
 * @Author James Robinson (james.robinson@cern.ch)
 * @date   11th June 2015
 * @brief  Study effective mass window for H->hh->bbyy
 *
 * Package for Hgamma analysis framework to determine optimal mass window
 */

#pragma once
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "TLorentzVector.h"

class TTree;

class SignalResonanceMassWindow : public HgammaAnalysis
{
private:
  // Input structures
  xAOD::TEvent *m_event; //!
  TTree *m_output_tree; //!

  // Configurable settings
  double m_m_yy_low, m_m_yy_high;

  // Output masses
  double m_m_yyjj_unscaled_0tag, m_m_yyjj_unscaled_1tag, m_m_yyjj_unscaled_2tag;
  double m_m_yyjj_mHscaled_0tag, m_m_yyjj_mHscaled_1tag, m_m_yyjj_mHscaled_2tag;

  // Event weights
  double m_weight_pileup, m_weight_xslumi;

public:
  /// Constructors/destructors
  SignalResonanceMassWindow() {}
  SignalResonanceMassWindow(const char *name);
  virtual ~SignalResonanceMassWindow();

  /// Functions inherited from HgammaAnalysis
  virtual EL::StatusCode initialize();
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();

  /// this is needed to distribute the algorithm to the workers
  ClassDef(SignalResonanceMassWindow, 1);

private:
  TLorentzVector ApplyHiggsMassScaling( TLorentzVector input, double m_H_in_GeV = 125.0 );
};
