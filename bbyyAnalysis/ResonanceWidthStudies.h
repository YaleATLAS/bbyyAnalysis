#pragma once
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include <string>
#include <vector>

class ResonanceWidthStudies : public HgammaAnalysis
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)

private:
  std::vector<std::string> m_categories;

public:
  // this is a standard constructor
  ResonanceWidthStudies() { }
  ResonanceWidthStudies(const char *name);
  virtual ~ResonanceWidthStudies();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();


  // this is needed to distribute the algorithm to the workers
  ClassDef(ResonanceWidthStudies, 1);
};
