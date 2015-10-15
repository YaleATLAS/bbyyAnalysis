/**
 * @file   SignalResonanceMassWindow.cxx
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
  /// put your configuration variables here as public variables.
  /// that way they can be set directly from CINT and python.
public:
  /// float cutValue;

  /// variables that don't get filled at submission time should be
  /// protected from being send from the submission node to the worker
  /// node (done by the //!)
private:
  // Input structures
  xAOD::TEvent *m_event; //!
  TTree *m_output_tree; //!
  // Configurable settings
  double m_m_yy_low, m_m_yy_high;
  bool m_debug;
  // Output masses
  double m_m_yyjj_0tag, m_m_yyjj_1tag, m_m_yyjj_2tag;
  double m_m_yyjj_mHconstraint_0tag, m_m_yyjj_mHconstraint_1tag, m_m_yyjj_mHconstraint_2tag;

public:
  /// this is a standard constructor
  SignalResonanceMassWindow() { }
  SignalResonanceMassWindow(const char *name);
  virtual ~SignalResonanceMassWindow();

  /// these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode initialize();
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();

  /// this is needed to distribute the algorithm to the workers
  ClassDef(SignalResonanceMassWindow, 1);

private:
  TLorentzVector ApplyHiggsMassScaling( TLorentzVector input, double m_H_in_GeV = 125.0 );
};
