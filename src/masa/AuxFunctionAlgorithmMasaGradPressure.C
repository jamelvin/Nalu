/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#if defined (NALU_USES_MASA)

#include <masa/AuxFunctionAlgorithmMasaGradPressure.h>
#include <AuxFunction.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <Simulation.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace sierra{
namespace nalu{

AuxFunctionAlgorithmMasaGradPressure::AuxFunctionAlgorithmMasaGradPressure(
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

AuxFunctionAlgorithmMasaGradPressure::~AuxFunctionAlgorithmMasaGradPressure() {
  // delete Aux
  delete auxFunction_;
}

void
AuxFunctionAlgorithmMasaGradPressure::execute()
{

  // make sure that partVec_ is size one
  ThrowAssert( partVec_.size() == 1 );

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const unsigned nDim = meta_data.spatial_dimension();
  const double time = realm_.get_current_time();
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  VectorFieldType *gradPressureExact_;
  
  if (masaStore_) gradPressureExact_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "gradPressureExact");

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
    double * gradPressureExact;
    if (masaStore_) gradPressureExact = (double*) stk::mesh::field_data( *gradPressureExact_, *b.begin() );

    for(unsigned p=0; p < length; ++p) {

      //FIXME: These functions are not a part of MASA yet, using 0.0 as placeholder

      if (masaStore_)
      {
         fieldData[0] = 0.0; //gradPressureExact[0]; 
         fieldData[1] = 0.0; //gradPressureExact[1];
         if (nDim == 3)
            fieldData[2] = 0.0; //gradPressureExact[2];
         gradPressureExact += fieldSize;
      }
      else
      {
         fieldData[0] = 0.0;  
         fieldData[1] = 0.0; 
         if (nDim == 3)
            fieldData[2] = 0.0;
      }

      fieldData += fieldSize;
      coords += nDim;
    }

  }
}

} // namespace nalu
} // namespace Sierra

#endif //NALU_USES_MASA
