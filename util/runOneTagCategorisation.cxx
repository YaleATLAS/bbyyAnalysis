#include "bbyyAnalysis/OneTagCategorisation.h"
#include "HGamAnalysisFramework/RunUtils.h"

int main(int argc, char *argv[])
{
  // Set up the job for xAOD access
  xAOD::Init().ignore();

  // Create our algorithm
  OneTagCategorisation *alg = new OneTagCategorisation("OneTagCategorisation");

  // Use helper to start the job
  HG::runJob(alg, argc, argv);

  return 0;
}
