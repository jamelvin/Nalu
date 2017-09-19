/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#if defined (NALU_USES_MASA)

#include <masa/MasaEnthalpySrcNodeSuppAlg.h>
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
// MasaEnthalpySrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MasaEnthalpySrcNodeSuppAlg::MasaEnthalpySrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    coordinates_(NULL),
    enthalpySrc_(NULL),
    dualNodalVolume_(NULL),
    nDim_(realm_.spatialDimension_),
    masaStore_(true)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  masaStore_ = realm.masaInterface_.masaStoreSrc_;
  if (masaStore_)  
     enthalpySrc_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "enthalpySrc");  
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MasaEnthalpySrcNodeSuppAlg::setup()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MasaEnthalpySrcNodeSuppAlg::node_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double *coords = stk::mesh::field_data(*coordinates_, node);
  const double *entSrc;
  if (masaStore_) entSrc = stk::mesh::field_data(*enthalpySrc_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );

  double src = 0.0;
  const double x = coords[0];
  const double y = coords[1];
  double z = 0.0;
  
  if (nDim_ == 3)
      z = coords[2];
  
  if (masaStore_)  
     src = entSrc[0];
  else
     src = masa_eval_3d_source_e(x,y,z);

  rhs[0] += src*dualVolume;      
}

} // namespace nalu
} // namespace Sierra

#endif //NALU_USES_MASA

