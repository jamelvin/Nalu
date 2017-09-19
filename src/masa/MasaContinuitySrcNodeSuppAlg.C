/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#if defined (NALU_USES_MASA)

#include <masa/MasaContinuitySrcNodeSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// MasaContinuitySrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MasaContinuitySrcNodeSuppAlg::MasaContinuitySrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    coordinates_(NULL),
    continuitySrc_(NULL),
    dualNodalVolume_(NULL),
    nDim_(realm_.spatialDimension_),
    projTimeScale_(1.0),
    masaStore_(true)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  masaStore_ = realm.masaInterface_.masaStoreSrc_;
  if (masaStore_)
     continuitySrc_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "continuitySrc");  
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MasaContinuitySrcNodeSuppAlg::setup()
{
  projTimeScale_ = realm_.get_time_step()/realm_.get_gamma1();
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MasaContinuitySrcNodeSuppAlg::node_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double *coords = stk::mesh::field_data(*coordinates_, node);
  double *contSrc;
  if (masaStore_) contSrc = stk::mesh::field_data(*continuitySrc_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );

  double src;
  const double x = coords[0];
  const double y = coords[1];
  double z = 0.0;

  if (nDim_ == 3)
      z = coords[2];
  
  if (masaStore_)
     src = contSrc[0];
  else
     src = masa_eval_3d_source_rho(x,y,z);

  rhs[0] += src*dualVolume/projTimeScale_;      
}

} // namespace nalu
} // namespace Sierra

#endif //NALU_USES_MASA

