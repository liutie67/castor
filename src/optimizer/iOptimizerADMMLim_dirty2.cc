/*
This file is part of CASToR.

    CASToR is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    CASToR is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along with
    CASToR (in file GNU_GPL.TXT). If not, see <http://www.gnu.org/licenses/>.

Copyright 2017-2020 all CASToR contributors listed below:

    --> Didier BENOIT, Claude COMTAT, Marina FILIPOVIC, Thibaut MERLIN, Mael MILLARDET, Simon STUTE, Valentin VIELZEUF, Zacharias CHALAMPALAKIS

This is CASToR version 3.1.
*/

/*!
  \file
  \ingroup  optimizer
  \brief    Implementation of class iOptimizerADMMLim_dirty2
*/

#include "iOptimizerADMMLim_dirty2.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================


iOptimizerADMMLim_dirty2::iOptimizerADMMLim_dirty2() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // ---------------------------

  // Initial value at 1
  m_initialValue = 1.;
  // Only one backward image for ADMMLim_dirty2
  m_nbBackwardImages = 1;
  // ADMMLim_dirty2 accepts penalties
  m_requiredPenaltyDerivativesOrder = 1;
  // ADMMLim_dirty2 is compatible with listmode and histogram data (not yet)
  m_listmodeCompatibility = false;
  m_histogramCompatibility = true;
  // ADMMLim_dirty2 is specific to emission data
  m_emissionCompatibility = true;
  m_transmissionCompatibility = false;

  // ADMMLim_dirty2 needs to use the post process loop to compute v and u (one iteration needed for both)
  m_needPreProcessLoop = false;
  m_needPostProcessLoop = true;

  // --------------------------
  // Specific member parameters
  // --------------------------

  m_AxProduct = NULL;
  m_rBackgroundEvents = NULL;
  m_yData = NULL;
  m_alpha = -1;
  m_grad_norm_sum = -1.;
  m_proj_grad_norm_sum = -1.;

  m4p_firstDerivativePenaltyImage = NULL;
  mp_toWrite_vk = NULL;
  mp_toWrite_uk = NULL;
  m_grad_before = NULL;
  m_proj_grad_before = NULL;
  
  m_uk = NULL;
  m_vk = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerADMMLim_dirty2::~iOptimizerADMMLim_dirty2()
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

  // Not clean, but for now, writing files storing v{k+1} and u{k+1} in the destructor
  sOutputManager* p_outputManager = sOutputManager::GetInstance();
  if (mp_toWrite_vk && mp_DataFile->GetNbAdditionalData()!=0) // Do not write v if we initialize v^0 (v^0 is written in u variable...)
  {
    // Write v{k+1} sinogram with provided directory and filename, to be used in next Python iteration
    std::stringstream temp_ss_v;
    temp_ss_v << p_outputManager->GetPathName() << p_outputManager->GetBaseName() << "_v.img";
    std::string a_pathToImg_v = temp_ss_v.str(); 
    IntfWriteImage(a_pathToImg_v, mp_toWrite_vk, mp_DataFile->GetSinogramSize(), m_verbose);
  }
  if (mp_toWrite_uk)
  {
    // Write u{k+1} sinogram with provided directory and filename, to be used in next Python iteration
    std::stringstream temp_ss_u;
    temp_ss_u << p_outputManager->GetPathName() << p_outputManager->GetBaseName() << "_u.img";
    std::string a_pathToImg_u = temp_ss_u.str();
    IntfWriteImage(a_pathToImg_u, mp_toWrite_uk, mp_DataFile->GetSinogramSize(), m_verbose);
  }
  if (mp_toWrite_vk)
  {
    free(mp_toWrite_vk);
  }
  if (mp_toWrite_uk)
  {
    free(mp_toWrite_uk);
  }
  if (m_proj_grad_before)
  {
    free(m_proj_grad_before);
  }
  if (m_grad_before)
  {
    free(m_grad_before);
  }
  if (m_vk)
  {
    free(m_vk);
  }
  if (m_uk)
  {
    free(m_uk);
  }
  if (m_AxProduct)
  {
    free(m_AxProduct);
  }
  if (m_yData)
  {
    free(m_yData);
  }
  if (m_rBackgroundEvents)
  {
    free(m_rBackgroundEvents);
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerADMMLim_dirty2::ShowHelpSpecific()
{
  cout << "This optimizer is the ADMM with non-negativity constraint on projection space." << endl;
  cout << "It was proposed by H. Lim, Y. K. Dewaraja, and J. A. Fessler, Physics " << endl;
  cout << "in Medicine & Biology, 2018. It is for now developed without subsets." << endl; 
  cout << "Stepsize computation is for now only optimized for DIP_ADMM penalty." << endl;
  cout << "It will not converge for other penalties. It is implemented using the same " << endl;
  cout << "equations as in the original paper, except it can be used with as many inner" << endl;
  cout << "iterations as desired. For now, only one outer iteration is computed with CASToR. Thus, " << endl;
  cout << "an external script must be used to call CASToR as many times as there is iterations." << endl;
  cout << "IMPORTANT : the whole code must be compiled with FLTNB = double, otherwise rounding" << endl;
  cout << "errors become too large in and affect the computed images (huge difference between" << endl;
  cout << "different numbers of threads. " << endl;
  cout << "The following option can be used:" << endl;
  cout << "alpha: to set the convergence speed of the ADMM with penalty parameter alpha." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerADMMLim_dirty2::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  // Read the penalty parameter alpha
  key_word = "alpha";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_alpha, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerADMMLim_dirty2::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerADMMLim_dirty2::ReadOptionsList(const string& a_optionsList)
{
  // There are 1 floating point variable as option
  const int nb_options = 1;
  FLTNB options[nb_options];
  // Read them
  if (ReadStringOption(a_optionsList, options, nb_options, ",", "ADMMLim_dirty2 configuration"))
  {
    Cerr("***** iOptimizerADMMLim_dirty2::ReadOptionsList() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  // Affect option
  m_alpha = options[0];
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerADMMLim_dirty2::CheckSpecificParameters()
{
  // Check that penalty parameter alpha is positive
  if (m_alpha<=0.)
  {
    Cerr("***** iOptimizerADMMLim_dirty2->CheckSpecificParameters() -> Provided alpha (" << m_alpha << ") must be strictly positive !" << endl);
    return 1;
  }
  // Cannot deal with list-mode or transmission data
  if (m_dataSpec==SPEC_TRANSMISSION || m_dataMode==MODE_LIST)
  {
    Cerr("***** iOptimizerADMMLim_dirty2->CheckSpecificParameters() -> Cannot reconstruct list-mode or transmission data !" << endl);
    return 1;
  }
  // Check that number of additional data is 2, for v and u. It should only be 0 for v^0 initialization.
  if (mp_DataFile->GetNbAdditionalData()==0)
  {
    Cerr("***** iOptimizerADMMLim_dirty2::CheckSpecificParameters() -> ADMMLim_dirty2 requires 2 additional data for u^k and v^k ! (we will zero them, assuming you want to initialize v^0...)"<< endl);
  }
  else if (mp_DataFile->GetNbAdditionalData()!=2)
  {
    Cerr("***** iOptimizerADMMLim_dirty2::CheckSpecificParameters() -> ADMMLim_dirty2 requires 2 additional data for u^k and v^k !"<< endl);
    return 1;
  }
  // Check that number of multimodal image is 1, for previously computed gradient g^k.
  if (mp_ImageDimensionsAndQuantification->GetNbMultiModalImages()!=2)
  {
    Cerr("***** iOptimizerADMMLim_dirty2::CheckSpecificParameters() -> ADMMLim_dirty2 requires 1 multimodal image for g^k !"<< endl);
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerADMMLim_dirty2::InitializeSpecific()
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
  // Allocate u^k and v^k pointers, depending on current thread
  m_uk = (FLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection()*sizeof(FLTNB));
  m_vk = (FLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection()*sizeof(FLTNB));

  // Allocate and create the whole gradient, which need to be stored to compute forward projection of it
  m_grad_before = (FLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbVoxXYZ()*sizeof(FLTNB));

  // Allocate and create the sinograms v^k and u^k and the copy image
  mp_toWrite_vk = (FLTNB*)malloc(mp_DataFile->GetSinogramSize()*sizeof(FLTNB));
  mp_toWrite_uk = (FLTNB*)malloc(mp_DataFile->GetSinogramSize()*sizeof(FLTNB));
  // Loop over sinogram bins
  for (int lor=0; lor<mp_DataFile->GetSinogramSize(); lor++)
  {
    mp_toWrite_vk[lor] = 0.;
    mp_toWrite_uk[lor] = 0.;
  }

  // Allocate and initialize useful variables for v^k computation
  m_AxProduct = (HPFLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection()*sizeof(HPFLTNB));
  m_yData = (FLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection()*sizeof(FLTNB));
  m_rBackgroundEvents = (FLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection()*sizeof(FLTNB));
  for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection(); th++)
  {
    m_AxProduct[th] = -1.;
    m_yData[th] = -1.;
    m_rBackgroundEvents[th] = -1.;
  }

  // Allocate and initialize gradient projection
  m_proj_grad_before = (HPFLTNB*)malloc(mp_DataFile->GetSinogramSize()*sizeof(HPFLTNB));
  for (int lor=0; lor<mp_DataFile->GetSinogramSize(); lor++)
  {
    m_proj_grad_before[lor] = 0.;
  }

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iOptimizerADMMLim_dirty2::InitializeSpecific() -> Use the ADMMLim_dirty2 optimizer" << endl);
    if (m_verbose>=VERBOSE_DETAIL)
    {
      Cout("  --> Initial image value: " << m_initialValue << endl);
      Cout("  --> Penalty parameter alpha: " << m_alpha << endl);
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerADMMLim_dirty2::DataStep4Optional( oProjectionLine* ap_Line, vEvent* ap_Event,
                                   int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                   int a_th )
{
  // Create u^k and v^k pointers, depending on current thread
  if (mp_DataFile->GetNbAdditionalData()==0) // Initialize u and v to zeros, assuming the user wants to initalize v^0 = Ax^0
  {
    m_uk[a_th] = 0.;
    m_vk[a_th] = 0.;
  }
  else 
  {
    // Retrieve index of current event and get v^k and u^k at this event from -additional-data option
    m_uk[a_th] = mp_DataFile->m2p_additionalData[0][ap_Line->GetEventIndex()];
    m_vk[a_th] = mp_DataFile->m2p_additionalData[1][ap_Line->GetEventIndex()];    
  }

  // Loop over voxels
  for (int v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++)
  {
    m_grad_before[v] = mp_ImageSpace->m2p_multiModalImage[1][v];
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerADMMLim_dirty2::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                                   FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                   FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Line weight here is simply 1
  *ap_weight = 1;
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerADMMLim_dirty2::DataStep5ComputeCorrections( oProjectionLine* ap_Line, vEvent* ap_Event,
                                             int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                             int a_th )
{
  // Loop on TOF bins
  for (int b=0; b<ap_Line->GetNbTOFBins(); b++)
  {
    FLTNB* ap_backwardValues = m3p_backwardValues[a_th][b];                              // the backward values (the result)
    FLTNB a_additiveCorrections = ap_Event->GetAdditiveCorrections(b)*mp_ImageDimensionsAndQuantification->GetFrameDurationInSec(a_bed, a_timeFrame); // the additive corrections
    
    m_AxProduct[a_th] = (HPFLTNB)m2p_forwardValues[a_th][b] - (HPFLTNB)a_additiveCorrections;
    // Backward project (Ax - v^k + u^k), a_forwardModel is Ax here because we have overwritten the Forward Projection computation in this optimizer
    *ap_backwardValues = m_AxProduct[a_th] - (HPFLTNB)m_vk[a_th] + (HPFLTNB)m_uk[a_th];
  
    if (m_isInPostProcessLoop)
    {
      // We need the data and background events without backprojection after, so store them in m_yData and m_rBackgroundEvents
      m_yData[a_th] = (FLTNB)ap_Event->GetEventValue(b);
      m_rBackgroundEvents[a_th] = (FLTNB)ap_Event->GetAdditiveCorrections(b)*mp_ImageDimensionsAndQuantification->GetFrameDurationInSec(a_bed, a_timeFrame);
    }
  }
  
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerADMMLim_dirty2::DataStep6Optional( oProjectionLine* ap_Line, vEvent* ap_Event,
                                   int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                   int a_th )
{
  // Store norm of projection of gradient for next iteration (for stepsize in x computation)
  HPFLTNB proj_grad = (HPFLTNB)ForwardProject(ap_Line,m_grad_before);
  m_proj_grad_before[ap_Line->GetEventIndex()] += proj_grad * proj_grad;

  // Compute v and u as it is only one iteration after all iterations on x
  if (m_isInPostProcessLoop)
  {

    ////////////// v computation //////////////
    m_uk[a_th] =(FLTNB)m_uk[a_th];

    // Update image
    // Define useful variables to compute v_hat solution
    FLTNB beta = 0.5 * (1 / m_alpha + m_rBackgroundEvents[a_th] - m_uk[a_th] - (FLTNB)m_AxProduct[a_th]);
    FLTNB gamma = m_rBackgroundEvents[a_th] * (m_uk[a_th] + (FLTNB)m_AxProduct[a_th]) - (m_rBackgroundEvents[a_th] - m_yData[a_th]) / m_alpha;
    FLTNB v_hat = 0.;

    // Compute v_hat solution given the value of y
    if (m_yData[a_th] == 0.) 
    {
      v_hat = (FLTNB)m_AxProduct[a_th] + m_uk[a_th] - 1 / m_alpha;
    }
    else
    {
      if (beta * beta + gamma >= 0.)
      {

      }
      else // check for numerical instabilities
      {
        cout << "danger sqrt" << endl;
        cout << "beta * beta + gamma = " << beta * beta + gamma << endl;
      }
      
      if (m_yData[a_th] > 0. and beta < 0.)
      {
        v_hat = sqrt(beta * beta + gamma) - beta;
      }
      else if (m_yData[a_th] > 0. and beta >= 0.)
      {
        //v_hat = gamma / (sqrt(beta * beta + gamma) + beta); // formula proposed in Lim et al., but numerically instable
        v_hat = sqrt(beta * beta + gamma) - beta; // same formula written in a different way
      }
    }

    // Correct v^{k+1} if it does not respect the constraint (Ax+r>0) and store image in file
    m_vk[a_th] = ((v_hat + m_rBackgroundEvents[a_th]) > 0.) ? v_hat :  -m_rBackgroundEvents[a_th];
    mp_toWrite_vk[ap_Line->GetEventIndex()] = m_vk[a_th];

    ////////////// u computation //////////////
    if (mp_DataFile->GetNbAdditionalData()==0) // v^0 initialization : v^0 = Ax^0
    {
      mp_toWrite_uk[ap_Line->GetEventIndex()] = (FLTNB)m_AxProduct[a_th];
    }
    else // update image (u^{k+1} := u^k + Ax^{k+1} - v^{k+1}) and store it
    {
      mp_toWrite_uk[ap_Line->GetEventIndex()] = m_uk[a_th] + (FLTNB)m_AxProduct[a_th] - m_vk[a_th];
    }
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerADMMLim_dirty2::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                 FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                 FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Not useful in this optimizer, this function is directly written in DataStep5ComputeCorrections function as m_vk and m_uk depend on the current thread
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerADMMLim_dirty2::PreImageUpdateSpecificStep()
{
  // ==========================================================================================
  // If no penalty, then exit (the penalty image term has been initialized to 0)
  if (mp_Penalty==NULL) return 0;
  // Set the number of threads
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation());
  #endif
  // Verbose
  if (m_verbose>=1) Cout("iOptimizerADMMLim_dirty2::PreImageUpdateSpecificStep() -> Compute penalty term" << endl);
  // ==========================================================================================
  // Global precomputation step if needed by the penalty
  if (mp_Penalty->GlobalPreProcessingStep())
  {
    Cerr("***** iOptimizerADMMLim_dirty2::PreImageUpdateSpecificStep() -> A problem occurred while computing the penalty pre-processing step !" << endl);
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
            Cerr("***** iOptimizerADMMLim_dirty2::PreImageUpdateSpecificStep() -> A problem occurred while computing the penalty local pre-processing step for voxel " << v << " !" << endl);
            problem = true;
          }
          // Compute first derivative order penalty terms
          m4p_firstDerivativePenaltyImage[tbf][rbf][cbf][v] = mp_Penalty->ComputeFirstDerivative(tbf,rbf,cbf,v,th);
        }
        // Check for problems
        if (problem)
        {
          Cerr("***** iOptimizerADMMLim_dirty2::PreImageUpdateSpecificStep() -> A problem occurred inside the multi-threaded loop, stop now !" << endl);
          return 1;
        }
      }
    }
  }

  // Zero norm of gradient projection and norm of gradient projection for this iteration
  m_grad_norm_sum = 0.;
  m_proj_grad_norm_sum = 0.;
  // Compute norms using values from all threads
  for (int lor=0; lor<mp_DataFile->GetSinogramSize(); lor++)
  {
    m_proj_grad_norm_sum += m_proj_grad_before[lor]; // already squared
    // Zero the gradient projection for this LOR
    m_proj_grad_before[lor] = 0.;
  }
  for (int v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++)
  {
    m_grad_norm_sum += (HPFLTNB)m_grad_before[v]*(HPFLTNB)m_grad_before[v];
    // Zero the gradient for this voxel
    m_grad_before[v] = 0.;
  }

  // Normal end  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerADMMLim_dirty2::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                  FLTNB a_sensitivity, FLTNB* ap_correctionValues,
                                                  INTNB a_voxel, int a_tbf, int a_rbf, int a_cbf )
{
  // Store gradient for next iteration
  m_grad_before[a_voxel] = *ap_correctionValues;
  if (!m_isInPostProcessLoop) // Do x computation
  {
    ////////////// Update with preconditioned gradient descent //////////////  
    // Scale penalty with respect to the number of subsets to get correct balance between likelihood and penalty
    HPFLTNB penalty = ((HPFLTNB)(m4p_firstDerivativePenaltyImage[a_tbf][a_rbf][a_cbf][a_voxel])) / ((HPFLTNB)(mp_nbSubsets[m_currentIteration]));
    // Compute conjugate gradient best stepsize after line search using norms from previous iteration
    HPFLTNB stepsize = m_grad_norm_sum / (((HPFLTNB)m_alpha * m_proj_grad_norm_sum) + (mp_Penalty->GetPenaltyStrength()*m_grad_norm_sum));

    // Compute additive image update factor
    HPFLTNB gradient = -(HPFLTNB)m_alpha * (HPFLTNB)*ap_correctionValues + penalty;
    HPFLTNB additive_image_update_factor = stepsize * gradient;

    // Update image value and store it
    *ap_newImageValue = (HPFLTNB)a_currentImageValue + additive_image_update_factor;
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
