
/*!
  \file
  \ingroup algorithm
  \brief Implementation of class iIterativeAlgorithm
*/

#include "gVariables.hh"
#include "iIterativeAlgorithm.hh"
#include "iEventHistoCT.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iIterativeAlgorithm::iIterativeAlgorithm(): vAlgorithm()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iIterativeAlgorithm::~iIterativeAlgorithm()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iIterativeAlgorithm::StepBeforeIterationLoop()
{
  if (vAlgorithm::StepBeforeIterationLoop())
  {
    Cerr("***** iIterativeAlgorithm::StepBeforeIterationLoop() -> A problem occurred while calling StepBeforeIterationLoop() function !" << endl);
    return 1;
  }
  if (m_verbose>=2) Cout("iIterativeAlgorithm::StepBeforeIterationLoop ... " << endl);

  mp_ImageSpace->InstantiateBackwardImageFromDynamicBasis(mp_OptimizerManager->GetNbBackwardImages());
  mp_DeformationManager->InstantiateImageForDeformation(mp_ImageSpace);

  // Main image initialization
  if (mp_ImageSpace->InitImage(m_pathToInitialImg, mp_OptimizerManager->GetInitialValue()) )
  {
    Cerr("***** iIterativeAlgorithm::StepBeforeIterationLoop() -> An error occurred while reading the initialization image !" << endl);
    return 1;
  }

  // double the number of iterations to merge dirty1 and dirty2
  // m_nbIterations = 2*m_nbIterations;

  // Set numbers of iterations and subsets to the optimizer
  mp_OptimizerManager->SetNumbersOfIterationsAndSubsets(m_nbIterations, mp_nbSubsets);
  
  // If the optimizer needs a pre-processing step, do it here
  if (mp_OptimizerManager->NeedPreProcessLoop())
  {
    // Set pre-process flag of the optimizer
    mp_OptimizerManager->EnterPreProcessLoop();
    // Call the independent step function
    if (PreOrPostIterationIndependentStep())
    {
      Cerr("***** iIterativeAlgorithm::StepBeforeIterationLoop() -> A problem occured while calling PreOrPostIterationIndependentStep() function !" << endl);
      return 1;
    }
    // Unset pre-process flag of the optimizer
    mp_OptimizerManager->ExitPreProcessLoop();
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iIterativeAlgorithm::StepBeforeSubsetLoop(int a_iteration)
{
  if (m_verbose>=3) Cout("iIterativeAlgorithm::StepBeforeSubsetLoop ... " << endl);

  // Set the current iteration to the optimizer
  mp_OptimizerManager->SetCurrentIteration(a_iteration);

  // Set the isInDualProcessLoop flag according to the modulo 2 of current iteration
  if (a_iteration%2==0) mp_OptimizerManager->ExitDualProcessLoop();
  if (a_iteration%2==1) mp_OptimizerManager->EnterDualProcessLoop();

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iIterativeAlgorithm::StepPreProcessInsideSubsetLoop(int a_iteration, int a_subset)
{
  if (m_verbose>=3) Cout("iIterativeAlgorithm::StepPreProcessInsideSubsetLoop ... " << endl);

  // Set the current subset to the optimizer
  mp_OptimizerManager->SetCurrentSubset(a_subset);

  // Get the chrono manager singleton pointer
  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();

  // Initialize the correction backward image(s)
  mp_ImageSpace->InitBackwardImage();

  // Copy current image in forward-image buffer
  mp_ImageSpace->PrepareForwardImage() ; 

  // Apply image processing to forward image
  if (mp_ImageProcessingManager->ApplyProcessingForward(mp_ImageSpace))
  {
    Cerr("***** iIterativeAlgorithm::StepPreProcessInsideSubsetLoop() -> A problem occurred while applyin image processing to forward image !" << endl);
    return 1;
  }

  // Apply convolver to forward image
  p_ChronoManager->StartConvolution();
  if (mp_ImageConvolverManager->ConvolveForward(mp_ImageSpace))
  {
    Cerr("***** iIterativeAlgorithm::StepPreProcessInsideSubsetLoop() -> A problem occurred while applying image convolver to forward image !" << endl);
    return 1;
  }
  p_ChronoManager->StopConvolution();
  
  // Call the pre-event loop function from the optimizer manager
  mp_OptimizerManager->PreDataUpdateStep();

  // Initialisation of the backup images for deformation
  // (bu_fward image initialized with current forward image)
  // (bu_bward and bu_sens images set to 0)

  mp_DeformationManager->InitImageForDeformation(mp_ImageSpace);

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iIterativeAlgorithm::StepInnerLoopInsideSubsetLoop(int a_iteration, int a_subset, int a_bed)
{
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    if (m_nbBeds>1) Cout("iIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> Start loop over events for bed " << a_bed+1 << endl << flush);
    else Cout("iIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> Start loop over events" << endl << flush);
  }
    
  // Reinitialize 4D gate indexes
  mp_ID->ResetCurrentDynamicIndices();

  // Apply the bed offset for this bed position
  mp_ProjectorManager->ApplyBedOffset(a_bed);

  // Progression (increments of 2%)
  if (m_verbose>=VERBOSE_NORMAL && mp_ID->GetMPIRank()==0)
  {
    cout << "0   10%  20%  30%  40%  50%  60%  70%  80%  90%  100%" << endl;
    cout << "|----|----|----|----|----|----|----|----|----|----|" << endl;
    cout << "|" << flush;
  }
  int progression_percentage_old = 0;
  int progression_nb_bars = 0;
  uint64_t progression_printing_index = 0;

  // Compute start and stop indices taking MPI into account (the vDataFile does that)
  int64_t index_start = 0;
  int64_t index_stop  = 0;
  m2p_DataFile[a_bed]->GetEventIndexStartAndStop(&index_start, &index_stop, a_subset, mp_nbSubsets[a_iteration]);

  // Synchronize MPI processes
  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // Set the number of threads for projections (right after this loop, we set back the number of threads for image computations)
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ID->GetNbThreadsForProjection());
  #endif

  // This boolean is used to report any problem inside the parallel loop
  bool problem = false;

  // Get the chrono manager singleton pointer
  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();

  // These flags (one per thread) are used to signal when events are beyond the last frame time stop
  bool* p_end_loop_flags = (bool*)malloc(mp_ID->GetNbThreadsForProjection()*sizeof(bool));
  for (int th=0; th<mp_ID->GetNbThreadsForProjection(); th++) p_end_loop_flags[th] = false;

  // Launch the loop with precomputed start and stop and using a step equal to the number of subsets
  int64_t index;
  // Keep the static scheduling with a chunk size at 1, it is important
  #pragma omp parallel for private(index) schedule(static, 1)
  for ( index = index_start  ;  index < index_stop  ;  index += mp_nbSubsets[a_iteration] )
  {              
    // Get the thread index
    int th = 0;
    #ifdef CASTOR_OMP
    th = omp_get_thread_num();
    #endif
    // Print progression (do not log out with Cout here)
    if (m_verbose>=2 && th==0 && mp_ID->GetMPIRank()==0)
    {
      if (progression_printing_index%1000==0)
      {
        int progression_percentage_new = ((int)( (((float)(index-index_start+1))/((float)(index_stop-index_start)) ) * 100.));
        if (progression_percentage_new>=progression_percentage_old+2) // Increments of 2%
        {
          int nb_steps = (progression_percentage_new-progression_percentage_old)/2;
          for (int i=0; i<nb_steps; i++)
          {
            cout << "-" << flush;
            progression_nb_bars++;
          }
          progression_percentage_old += nb_steps*2;
        }
      }
      progression_printing_index++;
    }
    // Skip to the end
    if (p_end_loop_flags[th]) continue;

    // Step 1: Get the current event for that thread index
    p_ChronoManager->StartIterativeDataUpdateStep1(th);
    #ifdef CASTOR_DEBUG
    if (m_verbose>=4) Cout("iIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> Step1: Get current event for that thread index " << endl);
    #endif
    vEvent* event = m2p_DataFile[a_bed]->GetEvent(index, th);
    if (event==NULL)
    {
      Cerr("***** iIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> An error occurred while getting the event from index "
        << index << " (thread " << th << ") !" << endl);
      // Specify that there was a problem
      problem = true;
      // We must continue here because we are inside an OpenMP loop
      continue;
    }
    p_ChronoManager->StopIterativeDataUpdateStep1(th);

    // Step 2: Call the dynamic switch function that updates the current frame and gate numbers, and also detects involuntary patient motion
    p_ChronoManager->StartIterativeDataUpdateStep2(th);
    #ifdef CASTOR_DEBUG
    if (m_verbose>=4)
    {
      Cout("iIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> Step2: Check for Dynamic event (frame/gate switch, image-based deformation " << endl);
    }
    #endif
    int dynamic_switch_value = mp_ID->DynamicSwitch(index, event->GetTimeInMs(), a_bed, th);
    // If the DYNAMIC_SWITCH_CONTINUE is returned, then it means that we are not yet at the first frame
    if ( dynamic_switch_value == DYNAMIC_SWITCH_CONTINUE )
    {
      // Then we just skip this event
      continue;
    }
    // Else, if the DYNAMIC_SWITCH_DEFORMATION is returned, then it means that a change of gate or involuntary motion has occurred
    // and is associated to a deformation
    else if ( dynamic_switch_value == DYNAMIC_SWITCH_DEFORMATION )
    {
      #ifdef CASTOR_OMP
      // set deformation requirement for this thread
      #pragma omp critical (deformation_requirements_first)
      {
        mp_DeformationManager->SetDeformationRequirement(th);
      }
      // block this thread until the deformation is done
      bool wait_for_deformation = true;
      while (wait_for_deformation)
      {
        // the first thread performs the deformation on the forward image
        // as soon as all the threads are ready for deformation
        if (th==0)
        {
          bool all_threads_ready = false;
          //check whether all the threads require deformation
          #pragma omp critical (deformation_requirements_first)
          {
            all_threads_ready = mp_DeformationManager->AllThreadsRequireDeformation();
          }
          if (all_threads_ready)
          {
            // Perform the deformation
            mp_DeformationManager->PerformDeformation(mp_ImageSpace);
            // Free all the threads so that they can continue
            #pragma omp critical (deformation_requirements_second)
            {
              mp_DeformationManager->UnsetDeformationRequirements();
            }
          }
        }
        // check whether the deformation is done and the thread can continue
        #pragma omp critical (deformation_requirements_second)
        {
          wait_for_deformation = mp_DeformationManager->GetDeformationRequirement(th);
        }
      }
      #else
      // Perform here any resp related deformations on the forward image
      mp_DeformationManager->PerformDeformation(mp_ImageSpace);
      #endif
      
      // Restore progression printing
      if (th==0 && m_verbose>=3 && mp_ID->GetMPIRank()==0)
      {
        cout << "0   10%  20%  30%  40%  50%  60%  70%  80%  90%  100%" << endl;
        cout << "|----|----|----|----|----|----|----|----|----|----|" << endl;
        int progression_percentage_new = ((int)( (((float)(index-index_start+1))/((float)(index_stop-index_start)) ) * 100.));
        cout << "|";

        for (int i=0; i< progression_percentage_new / 2 ; i++) cout << "-";
      }
          
    }
    // Else, if the DYNAMIC_SWITCH_END_LOOP is returned, set the local end_loop_flag boolean to true (which allows to quickly end the OMP parallel loop)
    else if ( dynamic_switch_value == DYNAMIC_SWITCH_END_LOOP )
    {
      p_end_loop_flags[th] = true;
      continue;
    }
    p_ChronoManager->StopIterativeDataUpdateStep2(th);

    // Step 3: Compute the projection line
    p_ChronoManager->StartIterativeDataUpdateStep3(th);
    #ifdef CASTOR_DEBUG
    if (m_verbose>=4) Cout("iIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> Step3: Compute the projection line " << endl);
    #endif
    oProjectionLine *line = mp_ProjectorManager->ComputeProjectionLine(event, th);
    if (line==NULL)
    {
      Cerr("***** iIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> A problem occurred while computing the projection line !" << endl);
      // Specify that there was a problem
      problem = true;
      // We must continue here because we are inside an OpenMP loop
      continue;
    }
    p_ChronoManager->StopIterativeDataUpdateStep3(th);

    // Step 4: Optimize in the data space (forward-proj, update, backward-proj)
    p_ChronoManager->StartIterativeDataUpdateStep4(th);
    #ifdef CASTOR_DEBUG
    if (m_verbose>=4) Cout("iIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> Step4: Optimize in the data space " << endl);
    #endif
    if (line->NotEmptyLine()) mp_OptimizerManager->DataUpdateStep( line,
                                                                   event,
                                                                   a_bed,
                                                                   mp_ID->GetCurrentTimeFrame(th),
                                                                   mp_ID->GetCurrentRespImage(th),
                                                                   mp_ID->GetCurrentCardImage(th),
                                                                   th);
    p_ChronoManager->StopIterativeDataUpdateStep4(th);
  } // End of indices loop (OpenMP stops here)

  // Synchronize MPI processes
  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // Set back the number of threads for image computation
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ID->GetNbThreadsForImageComputation());
  #endif

  // End of progression printing (do not log out with Cout here)
  if (m_verbose>=2 && mp_ID->GetMPIRank()==0)
  {
    int progression_total_bars = 49;
    for (int i=0; i<progression_total_bars-progression_nb_bars; i++) cout << "-";
    cout << "|" << endl;
  }

  // If a problem was encountered, then report it here
  if (problem)
  {
    Cerr("***** iIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> A problem occurred inside the parallel loop over events !" << endl);
    return 1;
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iIterativeAlgorithm::StepPostProcessInsideSubsetLoop(int a_iteration, int a_subset)
{
  if (m_verbose>=3) Cout("iIterativeAlgorithm::StepPostProcessInsideSubsetLoop ... " << endl);

  // Get the chrono manager singleton pointer
  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();

  // Merge parallel results
  mp_ImageSpace->Reduce();

  // Perform here any image-based deformations related to motion on the backward image.
  mp_DeformationManager->ApplyDeformationsToBackwardImage(mp_ImageSpace);

  // Apply convolver to backward images and sensitivity-if-needed
  p_ChronoManager->StartConvolution();
  if (mp_ImageConvolverManager->ConvolveBackward(mp_ImageSpace))
  {
    Cerr("***** iIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occurred while applying convolver to backward images !" << endl);
    return 1;
  }
  p_ChronoManager->StopConvolution();
  
  // Call the pre image update step function from the optimizer manager
  if (mp_OptimizerManager->PreImageUpdateStep())
  {
    Cerr("***** iIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occurred while performing the pre-image-update step !" << endl);
    return 1;
  }

  // Optimize in the image space (apply corrections, MAP and sensitivity); pass the number of subsets for list-mode sensitivity scaling
  if (mp_OptimizerManager->ImageUpdateStep())
  {
    Cerr("***** iIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occurred while performing the image update step !" << endl);
    return 1;
  }

  // Apply convolver to current estimed images
  p_ChronoManager->StartConvolution();
  if (mp_ImageConvolverManager->ConvolveIntra(mp_ImageSpace))
  {
    Cerr("***** iIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occurred while applying convolver to current estimate images !" << endl);
    return 1;
  }
  p_ChronoManager->StopConvolution();

  // Post-update Dynamic Modeling step (if enabled). Manage either linear/non linear dynamic model (physiological or not)
  if (mp_DynamicModelManager->ApplyDynamicModel(mp_ImageSpace, a_iteration, a_subset))
  {
    Cerr("***** iIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occurred while applying dynamic model to current estimate images !" << endl);
    return 1;
  }
  
  // Apply image processing to current estime images
  if (mp_ImageProcessingManager->ApplyProcessingIntra(mp_ImageSpace))
  {
    Cerr("***** iIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occurred while applying image processing to current estimate images !" << endl);
    return 1;
  }

  // Save the sensitivity image in histogram mode, if asked for
  if (mp_ID->GetMPIRank()==0 && mp_outputIterations[a_iteration] && m_saveSensitivityHistoFlag && m2p_DataFile[0]->GetDataMode()==MODE_HISTOGRAM)
  {
    // Get output manager to build the file name
    sOutputManager* p_output_manager = sOutputManager::GetInstance();
    // Build the file name
    string temp_sens = p_output_manager->GetPathName() + p_output_manager->GetBaseName();
    stringstream temp_it; temp_it << a_iteration + 1;
    stringstream temp_ss; temp_ss << a_subset + 1;
    temp_sens.append("_it").append(temp_it.str()).append("_ss").append(temp_ss.str()).append("_sensitivity");
    // Save sensitivity
    mp_ImageSpace->SaveSensitivityImage(temp_sens);
  }

  // Save the current image estimate if asked for
  if (mp_ID->GetMPIRank()==0 && mp_outputIterations[a_iteration] && m_saveImageAfterSubsets)
  {
    // Verbose
    if (m_verbose>=1) Cout("iIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> Save image at iteration " << a_iteration+1 << " for subset " << a_subset+1 << endl);
    // Save image
    if (StepImageOutput(a_iteration,a_subset))
    {
      Cerr("***** iIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occurred while saving images at iteration " << a_iteration+1 << " !" << endl);
      return 1;
    }
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iIterativeAlgorithm::StepAfterSubsetLoop(int a_iteration)
{
  if (m_verbose>=3) Cout("iIterativeAlgorithm::StepAfterSubsetLoop() -> Clean never visited voxels and save images if needed" << endl);
  // Clean never visited voxels
  mp_ImageSpace->CleanNeverVisitedVoxels();
  // Save the main image if not already done for each subset
  if (mp_ID->GetMPIRank()==0 && mp_outputIterations[a_iteration] && !m_saveImageAfterSubsets)
  {
    // Verbose
    if (m_verbose>=1) Cout("iIterativeAlgorithm::StepAfterSubsetLoop() -> Save image at iteration " << a_iteration+1 << endl);
    // Save image
    if (StepImageOutput(a_iteration))
    {
      Cerr("***** iIterativeAlgorithm::StepAfterSubsetLoop() -> A problem occurred while saving images at iteration " << a_iteration+1 << " !" << endl);
      return 1;
    }
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iIterativeAlgorithm::StepAfterIterationLoop()
{
  // If the optimizer needs a post-processing step, do it here
  if (mp_OptimizerManager->NeedPostProcessLoop())
  {
    // Set post-process flag of the optimizer
    mp_OptimizerManager->EnterPostProcessLoop();
    // Call the independent step function
    if (PreOrPostIterationIndependentStep())
    {
      Cerr("***** iIterativeAlgorithm::StepAfterIterationLoop() -> A problem occured while calling PreOrPostIterationIndependentStep() function !" << endl);
      return 1;
    }
    // Unset post-process flag of the optimizer
    mp_OptimizerManager->ExitPostProcessLoop();
  }

  if (vAlgorithm::StepAfterIterationLoop())
  {
    Cerr("***** iIterativeAlgorithm::StepAfterIterationLoop() -> A problem occurred while calling StepAfterIterationLoop() function !" << endl);
    return 1;
  }
  if (m_verbose>=2) Cout("iIterativeAlgorithm::StepAfterIterationLoop ... " << endl);
  
  // Deallocate everything
  mp_DeformationManager->DeallocateImageForDeformation(mp_ImageSpace);
  mp_ImageSpace->DeallocateBackwardImageFromDynamicBasis();

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iIterativeAlgorithm::StepImageOutput(int a_iteration, int a_subset)
{
  // Get the chrono manager singleton pointer
  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();
  // =================================================================================================
  // Save dynamic coefficient images from the vDynamicModel
  // =================================================================================================
  if (mp_DynamicModelManager->SaveParametricImages(a_iteration, a_subset))
  {
    Cerr("***** iIterativeAlgorithm::StepImageOutput() -> A problem occurred while saving parametric images related to the dynamic model !" << endl);
    return 1;
  }
  // =================================================================================================
  // Apply pre-processing steps
  // =================================================================================================
  // Copy the current image into the forward image
  mp_ImageSpace->PrepareForwardImage();
  // Apply post-convolution if needed
  p_ChronoManager->StartConvolution();
  if (mp_ImageConvolverManager->ConvolvePost(mp_ImageSpace))
  {
    Cerr("***** iIterativeAlgorithm::StepImageOutput() -> A problem occurred while convolving the output image !" << endl);
    return 1;
  }
  p_ChronoManager->StopConvolution();
  // Apply post-processing if needed
  if (mp_ImageProcessingManager->ApplyProcessingPost(mp_ImageSpace))
  {
    Cerr("***** iIterativeAlgorithm::StepImageOutput() -> A problem occurred while applying image processing the output image !" << endl);
    return 1;
  }
  // Apply output transaxial FOV masking
  if (mp_ImageSpace->ApplyOutputFOVMasking())
  {
    Cerr("***** iIterativeAlgorithm::StepImageOutput() -> A problem occurred while applying output FOV masking !" << endl);
    return 1;
  }
  // Apply output flip
  if (mp_ImageSpace->ApplyOutputFlip())
  {
    Cerr("***** iIterativeAlgorithm::StepImageOutput() -> A problem occurred while applying output flip !" << endl);
    return 1;
  }
  // =================================================================================================
  // If some basis functions are really in use
  // =================================================================================================
  if (!mp_ID->GetTimeStaticFlag() || !mp_ID->GetRespStaticFlag() || !mp_ID->GetCardStaticFlag())
  {
    // Save the basis coefficients images if asked for
    if (m_saveDynamicBasisCoefficients)
    {
      if (mp_ImageSpace->SaveOutputBasisCoefficientImage(a_iteration, a_subset))
      {
        Cerr("***** iIterativeAlgorithm::StepImageOutput() -> A problem occurred while saving output image !" << endl);
        return 1;
      }
    }
    // Compute the output image from the forward image basis functions
    mp_ImageSpace->ComputeOutputImage();
  }
  // =================================================================================================
  // Save frames/gates
  // =================================================================================================

  // Save output image (note that if no basis functions at all are in use, the the output image already points to the forward image)
  if (mp_ImageSpace->SaveOutputImage(a_iteration, a_subset))
  {
    Cerr("***** iIterativeAlgorithm::StepImageOutput() -> A problem occurred while saving output image !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iIterativeAlgorithm::PreOrPostIterationIndependentStep()
{
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("iIterativeAlgorithm::PreOrPostIterationIndependentStep() -> Additional data processing" << endl);

  // Get the chrono manager singleton pointer
  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();

  // Clock start for subset execution time
  clock_t clock_start_step = clock();
  time_t time_start_step = time(NULL);

  // Initialize the correction backward image(s)
  mp_ImageSpace->InitBackwardImage();

  // Copy current image in forward-image buffer
  mp_ImageSpace->PrepareForwardImage() ;

  // Apply image processing to forward image
  if (mp_ImageProcessingManager->ApplyProcessingForward(mp_ImageSpace))
  {
    Cerr("***** iIterativeAlgorithm::PreOrPostIterationIndependentStep() -> A problem occurred while applying image processing to forward image !" << endl);
    return 1;
  }

  // Apply convolver to forward image
  p_ChronoManager->StartConvolution();
  if (mp_ImageConvolverManager->ConvolveForward(mp_ImageSpace))
  {
    Cerr("***** iIterativeAlgorithm::PreOrPostIterationIndependentStep() -> A problem occurred while applying image convolver to forward image !" << endl);
    return 1;
  }
  p_ChronoManager->StopConvolution();

  // Call the pre-event loop function from the optimizer manager
  mp_OptimizerManager->PreDataUpdateStep();

  // Initialisation of the backup images for deformation
  // (bu_fward image initialized with current forward image)
  // (bu_bward and bu_sens images set to 0)
  mp_DeformationManager->InitImageForDeformation(mp_ImageSpace);

  // Synchronize MPI processes
  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // Loop on bed positions
  for (int bed=0 ; bed<m_nbBeds ; bed++)
  {
    // Reinitialize 4D gate indexes
    mp_ID->ResetCurrentDynamicIndices();
    // Apply the bed offset for this bed position
    mp_ProjectorManager->ApplyBedOffset(bed);
    // Progression (increments of 2%)
    if (m_verbose>=VERBOSE_NORMAL && mp_ID->GetMPIRank()==0)
    {
      cout << "0   10%  20%  30%  40%  50%  60%  70%  80%  90%  100%" << endl;
      cout << "|----|----|----|----|----|----|----|----|----|----|" << endl;
      cout << "|" << flush;
    }
    int progression_percentage_old = 0;
    int progression_nb_bars = 0;
    uint64_t progression_printing_index = 0;
    // Compute start and stop indices taking MPI into account (the vDataFile does that)
    int64_t index_start = 0;
    int64_t index_stop  = 0;
    int first_subset = 0;
    int one_subset = 1;
    m2p_DataFile[bed]->GetEventIndexStartAndStop(&index_start, &index_stop, first_subset, one_subset);
    // Synchronize MPI processes
    #ifdef CASTOR_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    // Set the number of threads for projections (right after this loop, we set back the number of threads for image computations)
    #ifdef CASTOR_OMP
    omp_set_num_threads(mp_ID->GetNbThreadsForProjection());
    #endif
    // This boolean is used to report any problem inside the parallel loop
    bool problem = false;
    // These flags (one per thread) are used to signal when events are beyond the last frame time stop
    bool* p_end_loop_flags = (bool*)malloc(mp_ID->GetNbThreadsForProjection()*sizeof(bool));
    for (int th=0; th<mp_ID->GetNbThreadsForProjection(); th++) p_end_loop_flags[th] = false;
    // Launch the loop with precomputed start and stop and using a step equal to the number of subsets
    int64_t index;
    // Keep the static scheduling with a chunk size at 1, it is important
    #pragma omp parallel for private(index) schedule(static, 1)
    for ( index = index_start  ;  index < index_stop  ;  index += one_subset )
    {
      // Get the thread index
      int th = 0;
      #ifdef CASTOR_OMP
      th = omp_get_thread_num();
      #endif
      // Print progression (do not log out with Cout here)
      if (m_verbose>=2 && th==0 && mp_ID->GetMPIRank()==0)
      {
        if (progression_printing_index%1000==0)
        {
          int progression_percentage_new = ((int)( (((float)(index-index_start+1))/((float)(index_stop-index_start)) ) * 100.));
          if (progression_percentage_new>=progression_percentage_old+2) // Increments of 2%
          {
            int nb_steps = (progression_percentage_new-progression_percentage_old)/2;
            for (int i=0; i<nb_steps; i++)
            {
              cout << "-" << flush;
              progression_nb_bars++;
            }
            progression_percentage_old += nb_steps*2;
          }
        }
        progression_printing_index++;
      }
      // Skip to the end
      if (p_end_loop_flags[th]) continue;
      // Step 1: Get the current event for that thread index
      p_ChronoManager->StartIterativeDataUpdateStep1(th);
      #ifdef CASTOR_DEBUG
      if (m_verbose>=4) Cout("iIterativeAlgorithm::PreOrPostIterationIndependentStep() -> Step1: Get current event for that thread index " << endl);
      #endif

      vEvent* event = m2p_DataFile[bed]->GetEvent(index, th);
      if (event==NULL)
      {
        Cerr("***** iIterativeAlgorithm::PreOrPostIterationIndependentStep() -> An error occurred while getting the event from index "
          << index << " (thread " << th << ") !" << endl);
        // Specify that there was a problem
        problem = true;
        // We must continue here because we are inside an OpenMP loop
        continue;
      }
      p_ChronoManager->StopIterativeDataUpdateStep1(th);
      // Step 2: Call the dynamic switch function that updates the current frame and gate numbers, and also detects involuntary patient motion
      p_ChronoManager->StartIterativeDataUpdateStep2(th);
      #ifdef CASTOR_DEBUG
      if (m_verbose>=4)
      {
        Cout("iIterativeAlgorithm::PreOrPostIterationIndependentStep() -> Step2: Check for Dynamic event (frame/gate switch, image-based deformation " << endl);
      }
      #endif
      int dynamic_switch_value = mp_ID->DynamicSwitch(index, event->GetTimeInMs(), bed, th);
      // If the DYNAMIC_SWITCH_CONTINUE is returned, then it means that we are not yet at the first frame
      if ( dynamic_switch_value == DYNAMIC_SWITCH_CONTINUE )
      {
        // Then we just skip this event
        continue;
      }
      // Else, if the DYNAMIC_SWITCH_DEFORMATION is returned, then it means that a change of gate or involuntary motion has occurred
      // and is associated to a deformation
      else if ( dynamic_switch_value == DYNAMIC_SWITCH_DEFORMATION )
      {
        #ifdef CASTOR_OMP
        // set deformation requirement for this thread
        #pragma omp critical (deformation_requirements_first)
        {
          mp_DeformationManager->SetDeformationRequirement(th);
        }
        // block this thread until the deformation is done
        bool wait_for_deformation = true;
        while (wait_for_deformation)
        {
          // the first thread performs the deformation on the forward image
          // as soon as all the threads are ready for deformation
          if (th==0)
          {
            bool all_threads_ready = false;
            //check whether all the threads require deformation
            #pragma omp critical (deformation_requirements_first)
            {
              all_threads_ready = mp_DeformationManager->AllThreadsRequireDeformation();
            }
            if (all_threads_ready)
            {
              // Perform the deformation
              mp_DeformationManager->PerformDeformation(mp_ImageSpace);
              // Free all the threads so that they can continue
              #pragma omp critical (deformation_requirements_second)
              {
                mp_DeformationManager->UnsetDeformationRequirements();
              }
            }
          }
          // check whether the deformation is done and the thread can continue
          #pragma omp critical (deformation_requirements_second)
          {
            wait_for_deformation = mp_DeformationManager->GetDeformationRequirement(th);
          }
        }
        #else
        // Perform here any resp related deformations on the forward image
        mp_DeformationManager->PerformDeformation(mp_ImageSpace);
        #endif
        // Restore progression printing
        if (th==0 && m_verbose>=3 && mp_ID->GetMPIRank()==0)
        {
          cout << "0   10%  20%  30%  40%  50%  60%  70%  80%  90%  100%" << endl;
          cout << "|----|----|----|----|----|----|----|----|----|----|" << endl;
          int progression_percentage_new = ((int)( (((float)(index-index_start+1))/((float)(index_stop-index_start)) ) * 100.));
          cout << "|";
          for (int i=0; i< progression_percentage_new / 2 ; i++) cout << "-";
        }
      }
      // Else, if the DYNAMIC_SWITCH_END_LOOP is returned, set the local end_loop_flag boolean to true (which allows to quickly end the OMP parallel loop)
      else if ( dynamic_switch_value == DYNAMIC_SWITCH_END_LOOP )
      {
        p_end_loop_flags[th] = true;
        continue;
      }
      p_ChronoManager->StopIterativeDataUpdateStep2(th);
      // Step 3: Compute the projection line
      p_ChronoManager->StartIterativeDataUpdateStep3(th);
      #ifdef CASTOR_DEBUG
      if (m_verbose>=4) Cout("iIterativeAlgorithm::PreOrPostIterationIndependentStep() -> Step3: Compute the projection line " << endl);
      #endif
      oProjectionLine *line = mp_ProjectorManager->ComputeProjectionLine(event, th);
      if (line==NULL)
      {
        Cerr("***** iIterativeAlgorithm::PreOrPostIterationIndependentStep() -> A problem occurred while computing the projection line !" << endl);
        // Specify that there was a problem
        problem = true;
        // We must continue here because we are inside an OpenMP loop
        continue;
      }
      p_ChronoManager->StopIterativeDataUpdateStep3(th);
      // Step 4: Optimize in the data space (forward-proj, update, backward-proj)
      p_ChronoManager->StartIterativeDataUpdateStep4(th);
      #ifdef CASTOR_DEBUG
      if (m_verbose>=4) Cout("iIterativeAlgorithm::PreOrPostIterationIndependentStep() -> Step4: Optimize in the data space " << endl);
      #endif
      if (line->NotEmptyLine()) mp_OptimizerManager->DataUpdateStep( line,
                                                                     event,
                                                                     bed,
                                                                     mp_ID->GetCurrentTimeFrame(th),
                                                                     mp_ID->GetCurrentRespImage(th),
                                                                     mp_ID->GetCurrentCardImage(th),
                                                                     th);
      p_ChronoManager->StopIterativeDataUpdateStep4(th);
    } // End of indices loop (OpenMP stops here)

    // Synchronize MPI processes
    #ifdef CASTOR_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    // Set back the number of threads for image computation
    #ifdef CASTOR_OMP
    omp_set_num_threads(mp_ID->GetNbThreadsForImageComputation());
    #endif
    // End of progression printing (do not log out with Cout here)
    if (m_verbose>=2 && mp_ID->GetMPIRank()==0)
    {
      int progression_total_bars = 49;
      for (int i=0; i<progression_total_bars-progression_nb_bars; i++) cout << "-";
      cout << "|" << endl;
    }
    // If a problem was encountered, then report it here
    if (problem)
    {
      Cerr("***** iIterativeAlgorithm::PreOrPostIterationIndependentStep() -> A problem occurred inside the parallel loop over events !" << endl);
      return 1;
    }
  } // End of beds loop

  // Synchronize MPI processes
  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // Merge parallel results
  mp_ImageSpace->Reduce();

  // If mask provided, apply mask to sensitivity, the CleanNeverVisitedVoxels will take care of the images
  mp_ImageSpace->ApplyMaskToSensitivity();

  // Perform here any image-based deformations related to motion on the backward image.
  mp_DeformationManager->ApplyDeformationsToBackwardImage(mp_ImageSpace);

  // Apply convolver to backward images and sensitivity-if-needed
  p_ChronoManager->StartConvolution();
  if (mp_ImageConvolverManager->ConvolveBackward(mp_ImageSpace))
  {
    Cerr("***** iIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occurred while applying convolver to backward images !" << endl);
    return 1;
  }
  p_ChronoManager->StopConvolution();

  // Call the pre image update step function from the optimizer manager
  if (mp_OptimizerManager->PreImageUpdateStep())
  {
    Cerr("***** iIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occurred while performing the pre-image-update step !" << endl);
    return 1;
  }

  // Optimize in the image space (apply corrections, MAP and sensitivity); pass the number of subsets for list-mode sensitivity scaling
  if (mp_OptimizerManager->ImageUpdateStep())
  {
    Cerr("***** iIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occurred while performing the image update step !" << endl);
    return 1;
  }

  // Clock stop for subset execution time
  clock_t clock_stop_step = clock();
  time_t time_stop_step = time(NULL);
  if (m_verbose>=2) Cout ("iIterativeAlgorithm::PreOrPostIterationIndependentStep() -> Time spent for step | User: " << time_stop_step-time_start_step
                       << " sec | CPU: " << (clock_stop_step-clock_start_step)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iIterativeAlgorithm::InitSpecificOptions(string a_specificOptions)
{
  (void)a_specificOptions; // avoid 'unused parameter' compil. warnings
  if (m_verbose>=2) Cout("iIterativeAlgorithm::InitSpecificOptions ... " << endl);
    
  return 0;
}
