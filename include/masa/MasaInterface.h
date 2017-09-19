/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#if defined (NALU_USES_MASA)

#ifndef MasaInterface_h
#define MasaInterface_h

#include <include/masa.h>
#include <FieldTypeDef.h>
#include <EquationSystem.h>
#include <vector>
#include <string>

// yaml for parsing..
#include <yaml-cpp/yaml.h>
#include <NaluParsing.h>

namespace sierra{
namespace nalu{

class Realm;

class MasaInterface
{
public:

  MasaInterface(Realm &realm_);
  ~MasaInterface() {}

  // This was in solutionNorm, should I have it here?
  //Realm &realm_;

  void register_nodal_fields(Realm &realm_, const std::vector<std::string> targetNames);
  void load(const YAML::Node & node);
  void initialize(Realm &realm_);
  void fillExactAndSrcFields(Realm &realm_);
 
  std::string MMStype_;
  std::string masaSoln_;
  std::vector<std::string> masaRequests_;
  
  std::map<std::string, double> masaParamsMap_;
  std::map<std::string, std::string> userFcnsMap_;

  bool masaReq_;
  bool masaStore_;
  bool masaStoreSrc_;

  double masakx_;
  double masaky_;
  double masakz_;

  double masaVisc_;

protected:

};

} // namespace nalu
} // namespace Sierra

#endif

#endif
