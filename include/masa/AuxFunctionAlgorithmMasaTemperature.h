/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#if defined (NALU_USES_MASA)

#ifndef AuxFunctionAlgorithmMasaTemperature_h
#define AuxFunctionAlgorithmMasaTemperature_h

#include <AuxFunctionAlgorithm.h>

#include <vector>
#include <stk_mesh/base/Types.hpp>

namespace stk{
namespace mesh{
class Part;
class FieldBase;
class Selector;

typedef std::vector< Part * > PartVector;
}
}

namespace sierra{
namespace nalu{

class AuxFunction;

class AuxFunctionAlgorithmMasaTemperature : public AuxFunctionAlgorithm
{
public:

  AuxFunctionAlgorithmMasaTemperature(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * field,
    AuxFunction * auxFunction,
    stk::mesh::EntityRank entityRank);

  virtual ~AuxFunctionAlgorithmMasaTemperature();
  virtual void execute();
  bool masaStore_;

private:
  stk::mesh::FieldBase * field_;
  AuxFunction *auxFunction_;
  stk::mesh::EntityRank entityRank_;
  
private:
  // make this non-copyable
  AuxFunctionAlgorithmMasaTemperature(const AuxFunctionAlgorithmMasaTemperature & other);
  AuxFunctionAlgorithmMasaTemperature & operator=(const AuxFunctionAlgorithmMasaTemperature & other);
};

} // namespace nalu
} // namespace Sierra

#endif

#endif
