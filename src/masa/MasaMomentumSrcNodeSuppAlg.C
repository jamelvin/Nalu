/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#if defined (NALU_USES_MASA)

#include <masa/MasaMomentumSrcNodeSuppAlg.h>
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
// MasaMomentumSrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MasaMomentumSrcNodeSuppAlg::MasaMomentumSrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    coordinates_(NULL),
    momentumSrc_(NULL),
    dualNodalVolume_(NULL),
    nDim_(realm_.spatialDimension_),
    masaStore_(true)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  masaStore_ = realm.masaInterface_.masaStoreSrc_;
  if (masaStore_)
     momentumSrc_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "momentumSrc");  
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // internal source
  srcXi_.resize(nDim_);
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MasaMomentumSrcNodeSuppAlg::setup()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MasaMomentumSrcNodeSuppAlg::node_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double *coords = stk::mesh::field_data(*coordinates_, node);
  const double *momSrc;
  if (masaStore_) momSrc  = stk::mesh::field_data(*momentumSrc_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double x = coords[0];
  const double y = coords[1];
  double z = 0.0;

  if (nDim_ == 3)
      z = coords[2];

  if (masaStore_)
  {  
     srcXi_[0] = momSrc[0];
     srcXi_[1] = momSrc[1];

     if (nDim_ == 3) {
         srcXi_[2] = momSrc[2];
     }
  }
  else
  {
     srcXi_[0] = masa_eval_3d_source_u(x,y,z);;
     srcXi_[1] = masa_eval_3d_source_v(x,y,z);;

     if (nDim_ == 3) {
         srcXi_[2] = masa_eval_3d_source_w(x,y,z);;
     }
  }

  for ( int i = 0; i < nDim_; ++i )
    rhs[i] += srcXi_[i]*dualVolume;      
}

} // namespace nalu
} // namespace Sierra

#endif //NALU_USES_MASA

