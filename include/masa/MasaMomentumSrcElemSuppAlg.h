/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#if defined (NALU_USES_MASA)

#ifndef MasaMomentumSrcElemSuppAlg_h
#define MasaMomentumSrcElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class MasaMomentumSrcElemSuppAlg : public SupplementalAlgorithm
{
public:

  MasaMomentumSrcElemSuppAlg(
    Realm &realm);

  virtual ~MasaMomentumSrcElemSuppAlg() {}

  virtual void setup();

  virtual void elem_resize(
    MasterElement *meSCS,
    MasterElement *meSCV);

  virtual void elem_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    MasterElement *meSCS,
    MasterElement *meSCV);
  
  const stk::mesh::BulkData *bulkData_;

  VectorFieldType *coordinates_;
  VectorFieldType *momentumSrc_;

  double dt_;
  const int nDim_;

  const bool useShifted_;

  // scratch space (at constructor)
  std::vector<double> scvCoords_;
  std::vector<double> scvMomSrc_;
  std::vector<double> srcXi_;
  // at elem_resize
  std::vector<double> ws_shape_function_;
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_momentumSrc_;
  std::vector<double> ws_scv_volume_;

  bool masaStore_;
};

} // namespace nalu
} // namespace Sierra

#endif

#endif
