#ifndef defaultConfigurations_H
#define defaultConfigurations_H

#include <boost/assign/list_of.hpp>
#include <map>

std::map<std::string, std::string> configMap = boost::assign::map_list_of
  ("strSamplePattert",          "mc15*Wmintau*") 
  ("pathToExtend_CERN",         "/afs/cern.ch/work/o/oviazlo/Wprime/"
                                "AnalysisFramework/rel20/data/") 
  ("pathToExtend_ALARIK",       "/lunarc/nobackup/users/oviazlo/xAOD/")
  ("pathToExtend_IRIDIUM",      "/nfs/shared/pp/oviazlo/xAOD/")
  ("submitDir",                 "submitDirs/default")
  ;

#endif