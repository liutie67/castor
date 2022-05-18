
/*!
  \file
  \ingroup  optimizer
  \brief    Implementation of class iPenaltyOPTITR
*/

#include "iPenaltyOPTITR.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iPenaltyOPTITR::iPenaltyOPTITR() : vPenalty()
{
  // ---------------------------
  // Mandatory member parameters
  // (inherited from vPenalty)
  // ---------------------------

  // Specify here the derivative order of the penalty.
  // Most algorithms able to handle penalties require 1 derivative,
  // and some 2 derivatives. If infinite, then set INT_MAX.
  m_penaltyDerivativesOrder = 1;

  // --------------------------
  // Specific member parameters
  // --------------------------
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iPenaltyOPTITR::~iPenaltyOPTITR()
{
  // Note: there is no need to deallocate the images themselves as they are allocate using the
  //       miscellaneous image function from the image space, which automatically deals with
  //       memory deallocations.
  // Delete the first order derivative penalty image
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iPenaltyOPTITR::ShowHelpSpecific()
{
  cout << "This penalty is OPTITR ?????????????" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyOPTITR::ReadConfigurationFile(const string& a_configurationFile)
{
  // This function is designed to read parameters specific to the optimizer through a configuration file.
  // To do that, use the ReadDataASCIIFile() function, and take a look at other optimizers to see how it is done.
  // Do not check the parameters' values, the CheckSpecificParameters() function is designed to do that.
  // Return 1 if any problem. See other penalties' implementation to get guidance.

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
 
int iPenaltyOPTITR::ReadOptionsList(const string& a_optionsList)
{ 
  // This function is designed to read parameters specific to the optimizer through a list of options in a string.
  // To do that, use the ReadStringOption() function, and take a look at other optimizers to see how it is done.
  // Do not check the parameters' values, the CheckSpecificParameters() function is designed to do that.
  // Return 1 if any problem. See other penalties' implementation to get guidance.

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyOPTITR::CheckSpecificParameters()
{
  /*
  // Check potential function type
  if (f_mu==NULL)
  {
    Cerr("***** iPenaltyOPTITR::CheckSpecificParameters() -> Should provide the path for f-mu image !" << endl);
    return 1;
  }
  */

  if (mp_ImageDimensionsAndQuantification->GetNbMultiModalImages()==0)
  {
    Cerr("***** iPenaltyOPTITR::CheckSpecificParameters() -> OPTITR requires a multimodal image (f-mu) !"<< endl);
    return 1;
  }
  else if (mp_ImageDimensionsAndQuantification->GetNbMultiModalImages()>1)
  {
    Cout("***** iPenaltyOPTITR::CheckSpecificParameters() -> Warning : More than one multimodal image, OPTITR will use only the first one !"<< endl);
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyOPTITR::InitializeSpecific()
{
  // Retrieve f-mu image (put in -multimodal argument)
  //for (INTNB v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++) std::cout << std::to_string(mp_ImageSpace-> m2p_multiModalImage[0][v]) << endl;

  //for (INTNB v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++) f_mu[v] = mp_ImageSpace-> m2p_multiModalImage[0][v];
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iOptimizerOPTITR::Initialize() -> Use the OPTITR optimizer" << endl);
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyOPTITR::ComputePenaltyValue(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // This is where you must implement the computation of the penalty value, for the provided parameters.
  // It is used to compute the overall cost function.
  // This function is called inside a parallel loop. The index of the thread is provided.
  // See other penalties' implementation to get guidance.

  // Pointer to the image
  FLTNB* p_image = mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf];
  // Precompute term x - (f-mu), with f-mu the image coming from -multimodal
  FLTNB square = p_image[a_voxel] - mp_ImageSpace-> m2p_multiModalImage[0][a_voxel];
  FLTNB penalty = 0.5 * m_penaltyStrength * square * square;
  // Check for Inf, Nan, etc
  if (fpclassify(penalty) != FP_NORMAL) penalty = 0.;
  // Return result
  return penalty;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyOPTITR::ComputeFirstDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // This is where you must implement the computation of the penalty's first derivative value, for the
  // provided parameters. It is used by the optimizer to compute the next update.
  // This function is called inside a parallel loop. The index of the thread is provided.
  // See other penalties' implementation to get guidance.
  
  // Pointer to the image
  // FLTNB* p_image = mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf];
  // Compute first derivative with f-mu the image coming from -multimodal
  // FLTNB square = p_image[a_voxel] - mp_ImageSpace-> m2p_multiModalImage[0][a_voxel];
  //FLTNB first_derivative = -m_penaltyStrength * (p_image[a_voxel] - mp_ImageSpace-> m2p_multiModalImage[0][a_voxel]); // True derivative
  FLTNB first_derivative = -m_penaltyStrength * mp_ImageSpace-> m2p_multiModalImage[0][a_voxel]; // To be used in OPTITR computation
  // Check for Inf, Nan, etc
  if (fpclassify(first_derivative) != FP_NORMAL) first_derivative = 0.;
  // Return result
  return first_derivative;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyOPTITR::ComputeSecondDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // This is where you must implement the computation of the penalty's second derivative value, for the
  // provided parameters. It may be used by the optimizer to compute the next update.
  // If the penalty does not admit a second derivative, then just let the function as is.
  // This function is called inside a parallel loop. The index of the thread is provided.
  // See other penalties' implementation to get guidance.
  FLTNB second_derivative = 0.;
  return second_derivative;
}
