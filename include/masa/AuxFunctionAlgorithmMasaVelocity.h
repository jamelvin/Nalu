/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef AuxFunctionAlgorithmMasaVelocity_h
#define AuxFunctionAlgorithmMasaVelocity_h

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

class AuxFunctionAlgorithmMasaVelocity : public AuxFunctionAlgorithm
{
public:

  AuxFunctionAlgorithmMasaVelocity(
    Realm & realm,
    stk::mesh::Part * part,
    stk::mesh::FieldBase * field,
    AuxFunction * auxFunction,
    stk::mesh::EntityRank entityRank);

  virtual ~AuxFunctionAlgorithmMasaVelocity();
  virtual void execute();
  bool masaStore_;

private:
  stk::mesh::FieldBase * field_;
  AuxFunction *auxFunction_;
  stk::mesh::EntityRank entityRank_;
  
private:
  // make this non-copyable
  AuxFunctionAlgorithmMasaVelocity(const AuxFunctionAlgorithmMasaVelocity & other);
  AuxFunctionAlgorithmMasaVelocity & operator=(const AuxFunctionAlgorithmMasaVelocity & other);
};

} // namespace nalu
} // namespace Sierra

#endif
