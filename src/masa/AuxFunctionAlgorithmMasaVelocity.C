/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#if defined (NALU_USES_MASA)

#include <masa/AuxFunctionAlgorithmMasaVelocity.h>
#include <AuxFunction.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <Simulation.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <iostream>
#include <fstream>

namespace sierra{
namespace nalu{

AuxFunctionAlgorithmMasaVelocity::AuxFunctionAlgorithmMasaVelocity(
  Realm & realm,
  stk::mesh::Part * part,
  stk::mesh::FieldBase * field,
  AuxFunction * auxFunction,
  stk::mesh::EntityRank entityRank)
  : AuxFunctionAlgorithm(realm, part, field, auxFunction, entityRank),
    field_(field),
    auxFunction_(auxFunction),
    entityRank_(entityRank)
{
  masaStore_ = realm.masaInterface_.masaStore_; 
}

AuxFunctionAlgorithmMasaVelocity::~AuxFunctionAlgorithmMasaVelocity() {
  // delete Aux
  delete auxFunction_;
}

void
AuxFunctionAlgorithmMasaVelocity::execute()
{

  std::ofstream myfile;

  // make sure that partVec_ is size one
  ThrowAssert( partVec_.size() == 1 );

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const unsigned nDim = meta_data.spatial_dimension();
  const double time = realm_.get_current_time();
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  VectorFieldType *velocityExact_;
 
  if (masaStore_) velocityExact_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocityExact");

  stk::mesh::Selector selector = stk::mesh::selectUnion(partVec_) &
    stk::mesh::selectField(*field_);

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( entityRank_, selector );

  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin();
        ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const unsigned fieldSize = field_bytes_per_entity(*field_, b) / sizeof(double);

    const stk::mesh::Bucket::size_type length   = b.size();

    // FIXME: because coordinates are only defined at nodes, this actually
    //        only works for nodal fields. Hrmmm.
    const double * coords = stk::mesh::field_data( *coordinates, *b.begin() );
    double * fieldData = (double*) stk::mesh::field_data( *field_, *b.begin() );
    double * velocityExact;

    if (masaStore_) velocityExact = (double*) stk::mesh::field_data( *velocityExact_, *b.begin() );

    for(unsigned p=0; p < length; ++p) {

      double x = coords[0];
      double y = coords[1];
      double z = 0.0;
 
      if (nDim == 3)
          z = coords[2];

      if (masaStore_)
      {
         fieldData[0] = velocityExact[0]; 
         fieldData[1] = velocityExact[1];
         if (nDim == 3)
            fieldData[2] = velocityExact[2];
         velocityExact += fieldSize;
      }
      else
      {
         fieldData[0] = masa_eval_3d_exact_u(x,y,z);
         fieldData[1] = masa_eval_3d_exact_v(x,y,z);
         if (nDim == 3)
             fieldData[2] = masa_eval_3d_exact_w(x,y,z);
      }

      fieldData += fieldSize;
      coords += nDim;
    }
  }
  
  //myfile.close();
}

} // namespace nalu
} // namespace Sierra

#endif //NALU_USES_MASA

