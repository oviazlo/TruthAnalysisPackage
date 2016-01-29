#ifndef TruthHelper_H
#define TruthHelper_H  
  
#include <TROOT.h>
#include <map>
#include <iostream>
  
class TruthHelper {

public:

  TruthHelper();
  ~TruthHelper();
  double getnEventsPerSample(unsigned int datasetID);
  double getCrossSection(unsigned int datasetID);
  
private:

  std::map<unsigned int, double> m_nEventsMap; //!
  std::map<unsigned int, double> m_crossSectionMap; //!
  /// this is needed to distribute the algorithm to the workers
  /// WARNING works ONLY if includes contain some ROOT header!!!
  ClassDef(TruthHelper, 1);
  
};
  
#endif