/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#if defined (NALU_USES_MASA)

#include <masa/MasaInterface.h>
#include <Realm.h>
#include <Realms.h>
#include <EquationSystems.h>
#include <SolutionOptions.h>

// yaml for parsing..
#include <yaml-cpp/yaml.h>
#include <NaluParsing.h>

// stk_util
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/Comm.hpp>

// stk_io
#include <stk_io/IossBridge.hpp>

// stk_topo
#include <stk_topology/topology.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// Simulation - do some stuff
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MasaInterface::MasaInterface(Realm &realm_) 
  : MMStype_("na"),
    masaSoln_("na"),
    masaReq_(false),
    masaStore_(true),
    masaStoreSrc_(true),
    masaVisc_(0),
    masaParamsMap_(),
    userFcnsMap_(),
    masaRequests_()
{
  // Does Nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
MasaInterface::register_nodal_fields(
  Realm &realm_, const std::vector<std::string> targetNames)
{
  if (masaReq_) {
    stk::mesh::MetaData &meta_data = realm_.meta_data();
  
    const int nDim = meta_data.spatial_dimension();

    // Pull the srcTermMaps created during the Solution Options load to determine what 
    // Masa src term storage fields are required 
    std::map<std::string, std::vector<std::string> > srcTermsMap = realm_.solutionOptions_->srcTermsMap_;
    std::map<std::string, std::vector<std::string> > elemSrcTermsMap = realm_.solutionOptions_->elemSrcTermsMap_;

    for ( size_t itarget = 0; itarget < targetNames.size(); ++itarget ) {
      stk::mesh::Part *targetPart = meta_data.get_part(targetNames[itarget]);
      if ( NULL == targetPart ) {
        NaluEnv::self().naluOutputP0() << "Trouble with part " << targetNames[itarget] << std::endl;
        throw std::runtime_error("Sorry, no part name found by the name " + targetNames[itarget]);
      }
      else {
        // register dof
        if (masaStore_)
        {
          if (masaStoreSrc_)
          {
            std::map<std::string, std::vector<std::string> >::iterator iter;
            for (iter = srcTermsMap.begin(); iter != srcTermsMap.end(); ++iter)
            {
              std::vector<std::string>::iterator it;
              for (it = iter->second.begin(); it != iter->second.end(); ++it)
              {
                std::string srcName = *it;
                // This is a Masa source term, so we want to create storage for it
                if (srcName == "Masa")
                {
                  std::string scaName = iter->first;
                  // Since momentum is the only vector src field we can isolate that one...
                  // FIXME: There is probably a way to combine the Vector and Scalars into a
                  // single statement...maybe following SolutionNormPostProcessing example...
                  if (iter->first == "momentum")
                  {
                    VectorFieldType *momentumSrc =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "momentumSrc"));
                    stk::mesh::put_field(*momentumSrc, *targetPart, nDim);
                  }
                  else 
                  {
                    ScalarFieldType *scaField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, scaName + "Src"));
                    stk::mesh::put_field(*scaField, *targetPart); 
                  }
                }
              }
            }
    
            for (iter = elemSrcTermsMap.begin(); iter != elemSrcTermsMap.end(); ++iter)
            {
              std::vector<std::string>::iterator it;
              for (it = iter->second.begin(); it != iter->second.end(); ++it)
              {
                std::string srcName = *it;
                // This is a Masa source term, so we want to create storage for it
                if (srcName == "Masa")
                {
                  std::string scaName = iter->first;
                  // Since momentum is the only vector src field we can isolate that one...
                  // FIXME: There is probably a way to combine the Vector and Scalars into a
                  // single statement...maybe following SolutionNormPostProcessing example...
                  if (iter->first == "momentum")
                  {
                    VectorFieldType *momentumSrc =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "momentumSrc"));
                    stk::mesh::put_field(*momentumSrc, *targetPart, nDim);
                    masaRequests_.push_back("momentumSrc");
                  }
                  else
                  {
                    ScalarFieldType *scaField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, scaName + "Src"));
                    stk::mesh::put_field(*scaField, *targetPart);
                    masaRequests_.push_back(scaName + "Src");
                  }
                }
              }
            }
          }
 
          std::map<std::string, std::string>::iterator iter2;
          for (iter2 = userFcnsMap_.begin(); iter2 != userFcnsMap_.end(); ++iter2)
          {
            // This is a Masa exact term, so we want to create storage for it
            if (iter2->second == "Masa")
            {
              std::string scaName = iter2->first;
              // FIXME: see above 
              if (iter2->first == "velocity")
              {
                VectorFieldType *velocityExact =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "velocityExact"));
                stk::mesh::put_field(*velocityExact, *targetPart, nDim);
                masaRequests_.push_back("velocityExact");
              }
              else
              {  
                ScalarFieldType *scaField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, scaName + "Exact"));
                stk::mesh::put_field(*scaField, *targetPart);
                masaRequests_.push_back(scaName + "Exact");
              }
            }
          }
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void MasaInterface::load(const YAML::Node & y_node)
{
  const YAML::Node y_MMS = y_node["MMS"];

  std::string thePropType;  
  std::string thePropName;

  if (y_MMS)
  {
    get_if_present(y_MMS, "name", MMStype_, MMStype_);

    if (MMStype_ == "Masa" || MMStype_ == "MASA" || MMStype_ == "masa")
    {
      masaReq_ = true;
      get_if_present(y_MMS, "solution_name", masaSoln_, masaSoln_);
      if (masaSoln_ == "na")
        throw std::runtime_error("ERROR: No MASA solution_name provided");

      get_if_present(y_MMS, "store_masa_terms", masaStore_, masaStore_);
      if (!masaStore_) 
         masaStoreSrc_ = false;
      else
         get_if_present(y_MMS, "store_masa_srcTerms", masaStoreSrc_, masaStoreSrc_);

      const YAML::Node y_masaParams = y_MMS["solution_params"];
      if (y_masaParams)
         masaParamsMap_ = y_masaParams.as<std::map<std::string, double> >();
 
      // This is to pull the viscosity
      const YAML::Node y_matProps = y_node["material_properties"];
      const YAML::Node y_specs = y_matProps["specifications"];
      for (size_t ispec = 0; ispec < y_specs.size(); ++ispec)
      {
        const YAML::Node y_spec = y_specs[ispec];
        get_required(y_spec, "name", thePropName);
        if (thePropName == "viscosity")
        {
          get_required(y_spec, "type", thePropType);
          if (thePropType != "constant")
            throw std::runtime_error("Error: Masa can only handle constant viscosity");
          get_required(y_spec, "value", masaVisc_);
        }
      }

      // This is to pull the initial condition functions required
      const YAML::Node y_IC = expect_sequence(y_node, "initial_conditions");
      for ( size_t iIC = 0; iIC < y_IC.size(); ++iIC)
      {
        const YAML::Node y_userIC = y_IC[iIC];
        if ( y_userIC["user_function"] )
        {
          const YAML::Node y_userFcns = y_userIC["user_function_name"];
          if (y_userFcns)
            userFcnsMap_ = y_userFcns.as<std::map<std::string, std::string> > ();
        }
      }
    }
    // Masa is the only MMS library linked up with Nalu at this time
    else if (MMStype_ != "")
      throw std::runtime_error("ERROR: Only Masa MMS is currently supported!");
  }
}


//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void MasaInterface::initialize(Realm &realm_)
{
   if (masaReq_) {
     const char *masaSoln = masaSoln_.c_str();
     masa_init("default_test", masaSoln);

     std::map<std::string, double>::iterator iter;

     for (iter = masaParamsMap_.begin(); iter != masaParamsMap_.end(); ++iter)
     {
         const char *varName = iter->first.c_str();
         masa_set_param(varName,iter->second);
     }

     masa_set_param("nu",masaVisc_);
  
     masa_sanity_check();
     masa_display_param(); 
     
     if (masaStore_) 
         fillExactAndSrcFields(realm_);
   }
}

//--------------------------------------------------------------------------
//-------- fillExactAndSrcFields -------------------------------------------
//--------------------------------------------------------------------------
void MasaInterface::fillExactAndSrcFields(Realm &realm_)
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const unsigned nDim = meta_data.spatial_dimension();
  const double time = realm_.get_current_time();
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  std::vector<std::string>::iterator it;
  for (it = masaRequests_.begin(); it != masaRequests_.end(); ++it)
  {
    std::string reqName = *it;
    stk::mesh::FieldBase * masaField = meta_data.get_field(stk::topology::NODE_RANK, reqName);

    // JAM: This is what I want to populate right??
    stk::mesh::Selector selector = stk::mesh::selectField(*masaField);
  
    // This defines the quantities on nodes
    stk::mesh::BucketVector const& buckets =
      realm_.get_buckets( stk::topology::NODE_RANK, selector );
  
    for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin();
          ib != buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const unsigned fieldSize = field_bytes_per_entity(*masaField, b) / sizeof(double);
  
      const stk::mesh::Bucket::size_type length   = b.size();
  
      const double * coords = stk::mesh::field_data( *coordinates, *b.begin() );
      double * masaData   = (double*)stk::mesh::field_data( *masaField,   *b.begin() );

      for(unsigned p=0; p < length; ++p) {
        const double x = coords[0];
        const double y = coords[1];
        double z = 0.0;
  
        if (nDim == 3)
            z = coords[2];

        if (reqName == "velocityExact") {
          masaData[0] = masa_eval_3d_exact_u(x,y,z);
          masaData[1] = masa_eval_3d_exact_v(x,y,z);
          if (nDim == 3) masaData[2] = masa_eval_3d_exact_w(x,y,z);
        }
        else if (reqName == "pressureExact")
          masaData[0] = masa_eval_3d_exact_p(x,y,z);
        else if (reqName == "temperatureExact")
          masaData[0] = masa_eval_3d_exact_t(x,y,z);
        else if (reqName == "momentumSrc") {
          masaData[0] = masa_eval_3d_source_u(x,y,z);
          masaData[1] = masa_eval_3d_source_v(x,y,z);
          if (nDim == 3) masaData[2] = masa_eval_3d_source_w(x,y,z);
        }
        else if (reqName == "continuitySrc")
          masaData[0] = masa_eval_3d_source_rho(x,y,z);
        else if (reqName == "enthalpySrc")
          masaData[0] = masa_eval_3d_source_e(x,y,z);
        else
          throw std::runtime_error("Sorry, " + reqName + " has not been setup in MasaInterface.C yet"); 
         
        masaData += fieldSize;
        coords += nDim;
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra

#endif //NALU_USES_MASA

