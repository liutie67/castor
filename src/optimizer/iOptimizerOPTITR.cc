
/*!
  \file
  \ingroup  optimizer
  \brief    Implementation of class iOptimizerOPTITR
*/

#include "iOptimizerOPTITR.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerOPTITR::iOptimizerOPTITR() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // ---------------------------

  // Initial value at 1
  m_initialValue = 1.;
  // Only one backward image for MLEM
  m_nbBackwardImages = 1;
  // OPTITR accepts penalties
  m_requiredPenaltyDerivativesOrder = 1;
  // MLEM is compatible with listmode and histogram data
  m_listmodeCompatibility = true;
  m_histogramCompatibility = true;
  // MLEM is compatible with both emission and log-converted transmission data
  m_emissionCompatibility = true;
  m_transmissionCompatibility = true;

  // --------------------------
  // Specific member parameters
  // --------------------------

  m_dataSpaceDenominatorThreshold = -1.;
  m_minimumImageUpdateFactor = -1.;
  m_maximumImageUpdateFactor = -1.;
  m4p_firstDerivativePenaltyImage = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerOPTITR::~iOptimizerOPTITR()
{
  // Note: there is no need to deallocate the images themselves as they are allocate using the
  //       miscellaneous image function from the image space, which automatically deals with
  //       memory deallocations.
  // Delete the first order derivative penalty image
  if (m4p_firstDerivativePenaltyImage)
  {
    // Loop over time basis functions
    for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
    {
      if (m4p_firstDerivativePenaltyImage[tbf])
      {
        // Loop over respiratory basis functions
        for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
        {
          if (m4p_firstDerivativePenaltyImage[tbf][rbf])
          {
            free(m4p_firstDerivativePenaltyImage[tbf][rbf]);
          }
        }
        free(m4p_firstDerivativePenaltyImage[tbf]);
      }
    }
    free(m4p_firstDerivativePenaltyImage);
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerOPTITR::ShowHelpSpecific()
{
  cout << "This optimizer is the ??????????????" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOPTITR::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  // Read the initial image value option
  key_word = "initial image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_initialValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerOPTITR::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the denominator threshold option
  key_word = "denominator threshold";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_dataSpaceDenominatorThreshold, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerOPTITR::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the minimum image update option
  key_word = "minimum image update";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_minimumImageUpdateFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerOPTITR::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the maximum image update option
  key_word = "maximum image update";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_maximumImageUpdateFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerOPTITR::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOPTITR::ReadOptionsList(const string& a_optionsList)
{
  // There are 4 floating point variables as options
  const int nb_options = 4;
  FLTNB options[nb_options];
  // Read them
  if (ReadStringOption(a_optionsList, options, nb_options, ",", "MLEM configuration"))
  {
    Cerr("***** iOptimizerOPTITR::ReadOptionsList() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  // Affect options
  m_initialValue = options[0];
  m_dataSpaceDenominatorThreshold = options[1];
  m_minimumImageUpdateFactor = options[2];
  m_maximumImageUpdateFactor = options[3];
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOPTITR::CheckSpecificParameters()
{
  // Check that initial image value is strictly positive
  if (m_initialValue<=0.)
  {
    Cerr("***** iOptimizerOPTITR->CheckSpecificParameters() -> Provided initial image value (" << m_initialValue << ") must be strictly positive !" << endl);
    return 1;
  }
  // Check that denominator threshold value is strictly positive
  if (m_dataSpaceDenominatorThreshold<=0.)
  {
    Cerr("***** iOptimizerOPTITR->CheckSpecificParameters() -> Provided data space denominator threshold (" << m_dataSpaceDenominatorThreshold << ") must be strictly positive !" << endl);
    return 1;
  }
  // Check that maximum image update factor is higher than the minimum
  if (m_minimumImageUpdateFactor>0. && m_maximumImageUpdateFactor>0. && m_maximumImageUpdateFactor<m_minimumImageUpdateFactor)
  {
    Cerr("***** iOptimizerOPTITR->CheckSpecificParameters() -> Provided minimum/maximum (" << m_minimumImageUpdateFactor << "/" << m_maximumImageUpdateFactor << " are inconsistent !" << endl);
    return 1;
  }
  // Cannot deal with list-mode transmission data
  if (m_dataSpec==SPEC_TRANSMISSION && m_dataMode==MODE_LIST)
  {
    Cerr("***** iOptimizerOPTITR->CheckSpecificParameters() -> Cannot reconstruct list-mode transmission data !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOPTITR::InitializeSpecific()
{
  // Allocate and create the penalty image
  m4p_firstDerivativePenaltyImage = (FLTNB****)malloc(mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions()*sizeof(FLTNB***));
  // Loop over time basis functions
  for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
  { 
    m4p_firstDerivativePenaltyImage[tbf] = (FLTNB***)malloc(mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions()*sizeof(FLTNB**));
    // Loop over respiratory basis functions
    for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
    {
      m4p_firstDerivativePenaltyImage[tbf][rbf] = (FLTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions()*sizeof(FLTNB*));
      // Loop over cardiac basis functions
      for (int cbf=0; cbf<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cbf++)
      {
        // Get a pointer to a newly allocated image
        m4p_firstDerivativePenaltyImage[tbf][rbf][cbf] = mp_ImageSpace -> AllocateMiscellaneousImage();
        // Initialize to 0, in case the penalty is not used
        for (INTNB v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++) m4p_firstDerivativePenaltyImage[tbf][rbf][cbf][v] = 0.;
      }
    }
  }
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iOptimizerOPTITR::InitializeSpecific() -> Use the OPTITR optimizer" << endl);
    if (m_verbose>=VERBOSE_DETAIL)
    {
      Cout("  --> Initial image value: " << m_initialValue << endl);
      Cout("  --> Data space denominator threshold: " << m_dataSpaceDenominatorThreshold << endl);
      if (m_minimumImageUpdateFactor>0.) Cout("  --> Minimum image update factor: " << m_minimumImageUpdateFactor << endl);
      else Cerr("!!!!! The minimum update value is not set, if using subsets, voxels could be trapped in 0 value causing some negative bias !" << endl);
      if (m_maximumImageUpdateFactor>0.) Cout("  --> Maximum image update factor: " << m_maximumImageUpdateFactor << endl);
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOPTITR::PreImageUpdateSpecificStep()
{
  // ==========================================================================================
  // If no penalty, then exit (the penalty image term has been initialized to 0)
  if (mp_Penalty==NULL) return 0;
  // Set the number of threads
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation());
  #endif
  // Verbose
  if (m_verbose>=1) Cout("iOptimizerPenalizedPreconditionedGradientML::PreImageUpdateSpecificStep() -> Compute penalty term" << endl);
  // ==========================================================================================
  // Global precomputation step if needed by the penalty
  if (mp_Penalty->GlobalPreProcessingStep())
  {
    Cerr("***** iOptimizerPenalizedPreconditionedGradientML::PreImageUpdateSpecificStep() -> A problem occurred while computing the penalty pre-processing step !" << endl);
    return 1;
  }
  // ==========================================================================================
  // Loop over time basis functions
  for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
  {
    // Loop over respiratory basis functions
    for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
    {
      // Loop over cardiac basis functions
      for (int cbf=0; cbf<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cbf++)
      {
        // In order to detect problems in the multi-threaded loop
        bool problem = false;
        // Voxel index
        INTNB v;
        // Multi-threading over voxels
        #pragma omp parallel for private(v) schedule(guided)
        for (v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++)
        { 
          // Get the thread index
          int th = 0;
          #ifdef CASTOR_OMP
          th = omp_get_thread_num();
          #endif
          // Local precomputation step if needed by the penalty
          if (mp_Penalty->LocalPreProcessingStep(tbf,rbf,cbf,v,th))
          {
            Cerr("***** iOptimizerPenalizedPreconditionedGradientML::PreImageUpdateSpecificStep() -> A problem occurred while computing the penalty local pre-processing step for voxel " << v << " !" << endl);
            problem = true;
          }
          // Compute first derivative order penalty terms
          m4p_firstDerivativePenaltyImage[tbf][rbf][cbf][v] = mp_Penalty->ComputeFirstDerivative(tbf,rbf,cbf,v,th);
        }
        // Check for problems
        if (problem)
        {
          Cerr("***** iOptimizerPenalizedPreconditionedGradientML::PreImageUpdateSpecificStep() -> A problem occurred inside the multi-threaded loop, stop now !" << endl);
          return 1;
        }
      }
    }
  }
  // Normal end  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOPTITR::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                                   FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                   FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Line weight here is simply 1
  *ap_weight = 1.;
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOPTITR::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                 FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                 FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Case for emission tomography
  if (m_dataSpec==SPEC_EMISSION)
  {
    // Truncate data to 0 if negative
    if (a_data<0.) a_data = 0.;
    // If the foward model is too close to zero, then ignore this data (set to 1 and not 0 because this line is taken into account in the sensitivity)
    if (a_forwardModel>m_dataSpaceDenominatorThreshold) *ap_backwardValues = a_data / a_forwardModel;
    else *ap_backwardValues = 1.;
  }
  // Case for transmission tomography (equivalent of log-converted MLEM from Nuyts et al 1998)
  else if (m_dataSpec==SPEC_TRANSMISSION)
  {
    // Subtract scatters
    a_data -= a_additiveCorrections;
    a_forwardModel -= a_additiveCorrections;
    // Safely ignore this data (set backward value to 1 and return) in 2 different cases:
    // - if data or model is inferior to 1
    // - if data or model is higher than blank value
    if (a_data<1. || a_forwardModel<1. || a_data>a_blankValue || a_forwardModel>a_blankValue)
    {
      *ap_backwardValues = 1.;
      return 0;
    }
    // Log-convert the data and the model
    a_data = log(a_blankValue/a_data);
    a_forwardModel = log(a_blankValue/a_forwardModel);
    // If the foward model is to close to zero, then ignore this data (set to 1 and not 0 because this line is taken into account in the sensitivity)
    if (a_forwardModel>m_dataSpaceDenominatorThreshold) *ap_backwardValues = a_data / a_forwardModel;
    else *ap_backwardValues = 1.;
  }
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOPTITR::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                  FLTNB a_sensitivity, FLTNB* ap_correctionValues,
                                                  INTNB a_voxel, int a_tbf, int a_rbf, int a_cbf )
{
  // ===== Start by computing the MLEM update
  // Compute image update factor
  FLTNB image_update_factor = *ap_correctionValues / a_sensitivity;
  // Apply minimum image update factor
  if ( m_minimumImageUpdateFactor > 0. && image_update_factor < m_minimumImageUpdateFactor ) image_update_factor = m_minimumImageUpdateFactor;
  // Apply maximum image update factor
  if ( m_maximumImageUpdateFactor > 0. && image_update_factor > m_maximumImageUpdateFactor ) image_update_factor = m_maximumImageUpdateFactor;
  // Update image
  // MLEM update
  HPFLTNB x_EM = (HPFLTNB)(a_currentImageValue) * (HPFLTNB)(image_update_factor);
  // Set the new value to the MLEM value by default
  *ap_newImageValue = x_EM;

  // Optimization transfer
  // ===== If a penalty exists, then include it using De Pierro's derivation for quadratic penalty
  if (mp_Penalty!=NULL)
  {  
    // Scale penalty with respect to the number of subsets to get correct balance between likelihood and penalty
    HPFLTNB penalty = ((HPFLTNB)(m4p_firstDerivativePenaltyImage[a_tbf][a_rbf][a_cbf][a_voxel])) / ((HPFLTNB)(mp_nbSubsets[m_currentIteration]));
    // Get penalty strength and divide it by the number of subsets to get correct balance between likelihood and penalty
    HPFLTNB beta = mp_Penalty->GetPenaltyStrength() / ((HPFLTNB)(mp_nbSubsets[m_currentIteration]));
    // Compute trinome solution
    HPFLTNB b = ((HPFLTNB)a_sensitivity) + penalty;
    HPFLTNB d = sqrt(pow(b,2) + 4 * beta * a_sensitivity * x_EM);
    *ap_newImageValue = 2 * (HPFLTNB)(a_sensitivity) * x_EM / (b + d);
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
