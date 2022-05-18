
/*!
  \file
  \ingroup  optimizer
  \brief    Implementation of class iPenaltyDIP_ADMM
*/

#include "iPenaltyDIP_ADMM.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iPenaltyDIP_ADMM::iPenaltyDIP_ADMM() : vPenalty()
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

iPenaltyDIP_ADMM::~iPenaltyDIP_ADMM()
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

void iPenaltyDIP_ADMM::ShowHelpSpecific()
{
  cout << "This penalty is coming from the derivation of DIP constraint using ADMM, like in" << endl;
  cout << "K. Gong, C. Catana, J. Qi, and Q. Li., IEEE transactions on medical imaging, 2018." << endl;
  cout << "It corresponds to a L2 distance between the reconstructed image and a reference " << endl;
  cout << "image given in the -multimodal option (computed with DIP output, in a Python script)." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyDIP_ADMM::ReadConfigurationFile(const string& a_configurationFile)
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
 
int iPenaltyDIP_ADMM::ReadOptionsList(const string& a_optionsList)
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

int iPenaltyDIP_ADMM::CheckSpecificParameters()
{
  if (mp_ImageDimensionsAndQuantification->GetNbMultiModalImages()==0)
  {
    Cerr("***** iPenaltyDIP_ADMM::CheckSpecificParameters() -> DIP_ADMM requires a multimodal image (f-mu) !"<< endl);
    return 1;
  }
  else if (mp_ImageDimensionsAndQuantification->GetNbMultiModalImages()>1)
  {
    Cout("***** iPenaltyDIP_ADMM::CheckSpecificParameters() -> Warning : More than one multimodal image, DIP_ADMM will only use the first one !"<< endl);
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyDIP_ADMM::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iOptimizerDIP_ADMM::Initialize() -> Use the DIP_ADMM optimizer" << endl);
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyDIP_ADMM::ComputePenaltyValue(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // Pointer to the image
  FLTNB* p_image = mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf];
  // Precompute term x - (f-mu), with f-mu the image coming from -multimodal (DIP output f - dual variable mu)
  FLTNB difference = p_image[a_voxel] - mp_ImageSpace-> m2p_multiModalImage[0][a_voxel];
  FLTNB penalty = 0.5 * m_penaltyStrength * difference * difference;
  // Check for Inf, Nan, etc
  if (fpclassify(penalty) != FP_NORMAL) penalty = 0.;
  // Return result
  return penalty;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyDIP_ADMM::ComputeFirstDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // Pointer to the image
  FLTNB* p_image = mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf];
  // Compute first derivative with f-mu the image coming from -multimodal (DIP output f - dual variable mu)
  FLTNB difference = p_image[a_voxel] - mp_ImageSpace-> m2p_multiModalImage[0][a_voxel];
  FLTNB first_derivative = -m_penaltyStrength * difference;
  // Check for Inf, Nan, etc
  if (fpclassify(first_derivative) != FP_NORMAL) first_derivative = 0.;
  // Return result
  return first_derivative;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyDIP_ADMM::ComputeSecondDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  FLTNB second_derivative = -m_penaltyStrength;
  return second_derivative;
}
