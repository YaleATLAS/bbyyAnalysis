/**
 * @file   runSignalResonanceMassWindow.cxx
 * @Author James Robinson (james.robinson@cern.ch)
 * @date   26th June 2015
 * @brief  Resonant mass window
 *
 * Executable for Hgamma analysis framework to analyse mass windows
 */

#include "SignalResonanceMassWindow/SignalResonanceMassWindow.h"
#include "HGamAnalysisFramework/RunUtils.h"

int main(int argc, char *argv[])
{
  // Set up the job for xAOD access
  xAOD::Init().ignore();

  // Create our algorithm
  SignalResonanceMassWindow *alg = new SignalResonanceMassWindow("SignalResonanceMassWindow");

  // Use helper to start the job
  HG::runJob(alg, argc, argv);

  return 0;
}
