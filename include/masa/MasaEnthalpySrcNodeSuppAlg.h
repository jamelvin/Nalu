/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#if defined (NALU_USES_MASA)

#ifndef MasaEnthalpySrcNodeSuppAlg_h
#define MasaEnthalpySrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MasaEnthalpySrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MasaEnthalpySrcNodeSuppAlg(
    Realm &realm);

  virtual ~MasaEnthalpySrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *coordinates_;
  VectorFieldType *enthalpySrc_;
  ScalarFieldType *dualNodalVolume_;

  const int nDim_;

  bool masaStore_;

};

} // namespace nalu
} // namespace Sierra

#endif

#endif
