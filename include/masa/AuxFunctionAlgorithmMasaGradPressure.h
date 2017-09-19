/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#if defined (NALU_USES_MASA)

#ifndef AuxFunctionAlgorithmMasaGradPressure_h
#define AuxFunctionAlgorithmMasaGradPressure_h

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

class AuxFunctionAlgorithmMasaGradPressure : public AuxFunctionAlgorithm
{
public:

  AuxFunctionAlgorithmMasaGradPressure(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * field,
    AuxFunction * auxFunction,
    stk::mesh::EntityRank entityRank);

  virtual ~AuxFunctionAlgorithmMasaGradPressure();
  virtual void execute();
  bool masaStore_;

private:
  stk::mesh::FieldBase * field_;
  AuxFunction *auxFunction_;
  stk::mesh::EntityRank entityRank_;
  
private:
  // make this non-copyable
  AuxFunctionAlgorithmMasaGradPressure(const AuxFunctionAlgorithmMasaGradPressure & other);
  AuxFunctionAlgorithmMasaGradPressure & operator=(const AuxFunctionAlgorithmMasaGradPressure & other);
};

} // namespace nalu
} // namespace Sierra

#endif

#endif
