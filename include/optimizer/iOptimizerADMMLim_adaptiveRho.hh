
/*!
  \file
  \ingroup  optimizer
  \brief    Declaration of class iOptimizerADMMLim
*/

#ifndef iOptimizerADMMLim_HH
#define iOptimizerADMMLim_HH 1

#include "gVariables.hh"
#include "sAddonManager.hh"
#include "vOptimizer.hh"

/*!
  \class   iOptimizerADMMLim
  \brief   This class implements the ADMM with non-negativity on projection space
  \details This class inherits from vOptimizer and implements the ADMM algorithm
           proposed by Lim et al. to ensure positivity of the projections.
*/
class iOptimizerADMMLim_adaptiveRho : public vOptimizer
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public iOptimizerADMMLim::iOptimizerADMMLim()
      \brief   The constructor of iOptimizerADMMLim
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    iOptimizerADMMLim_adaptiveRho();
    /*!
      \fn      public iOptimizerADMMLim::~iOptimizerADMMLim()
      \brief   The destructor of iOptimizerADMMLim
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
    */
    ~iOptimizerADMMLim_adaptiveRho();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameter and do not add semi-column at the end of the line)
    FUNCTION_OPTIMIZER(iOptimizerADMMLim_adaptiveRho)
    /*!
      \fn      public int iOptimizerADMMLim::ReadConfigurationFile()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to the child optimizer, from
               a configuration file. It is the implementation of the pure virtual function inherited
               from the abstract class vOptimizer. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadConfigurationFile(const string& a_configurationFile);
    /*!
      \fn      public int iOptimizerADMMLim::ReadOptionsList()
      \param   const string& a_optionsList
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to the child optimizer, from
               a list of options. It is the implementation of the pure virtual function inherited
               from the abstract class vOptimizer. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadOptionsList(const string& a_optionsList);
    
    /*!
      \fn      public virtual int vOptimizer::DataStep4Optional()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_bed
      \param   int a_timeFrame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_thread
      \brief   A public function which does nothing but being virtual.
      \details Inside the DataUpdateStep() of the oOptimizerManager, this is the second function called. This function is 
               overloaded ADMMLim optimizer to initialize u^k and v^k variables. 
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int DataStep4Optional( oProjectionLine* ap_Line, vEvent* ap_Event,
                                   int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                   int a_thread );
    /*!
      \fn      public virtual int iOptimizerADMMLim::DataStep5ComputeCorrections()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_bed
      \param   int a_timeFrame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_thread
      \brief   A public function used to compute the correction terms in the data space, for the provided event
      \details Inside the DataUpdateStep() of the oOptimizerManager, this is the fifth function called. Its role is to compute
               the correction terms in the data space, based on the forward model and the data. In order to be specific to each
               optimizer, it calls the pure virtual function DataSpaceSpecificOperations(), where the computation is done. The
               correction terms are put in the m3p_backwardValues (a dimension for threads, one for the number of backward images
               and the last for TOF bins).
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int DataStep5ComputeCorrections( oProjectionLine* ap_Line, vEvent* ap_Event,
                                             int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                             int a_thread );
    /*!
      \fn      public int iOptimizerADMMLim::DataStep6Optional()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_bed
      \param   int a_timeFrame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_thread
      \brief   A public function to compute the analytical solution of ADMM equation on v
      \details Inside the DataUpdateStep() of the oOptimizerManager, this is the sixth function called. 
               For this optimizer, here is computed u^{k+1} and v^{k+1} update thanks to the paper analytical formulas,
               using preivously stored variables and x computation. 
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int DataStep6Optional( oProjectionLine* ap_Line, vEvent* ap_Event,
                                   int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                   int a_thread );
  
  // -------------------------------------------------------------------
  // Private member functions (virtual in vOptimizer)
  private:
  
    /*!
      \fn      private int iOptimizerADMMLim::PreImageUpdateSpecificStep()
      \brief   A private function used to compute the penalty term of the ADMM algorithm
      \details This function implements the virtual eponym function of vOptimizer.
               It computes the penalty term of the ADMM algorithm.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int PreImageUpdateSpecificStep();
    /*!
      \fn      private void iOptimizerADMMLim::ShowHelpSpecific()
      \brief   A function used to show help about the child optimizer
      \details This function must describe what the module does and how to use it. It describes in
               details the different parameters of the optimizer, and how to set them through the use
               of a configuration file or a list of options. It is the implementation of the pure
               virtual function inherited from the abstract class vOptimizer. It is called by the
               public ShowHelp() function.
    */
    void ShowHelpSpecific();
    /*!
      \fn      private int iOptimizerADMMLim::CheckSpecificParameters()
      \brief   A private function used to check the parameters settings specific to the child optimizer
      \details This function is used to check that all parameters specific to the optimizer are correctly set
               within allowed values. It is called by the CheckParameters() function of the mother class.
               It is the implementation of the pure virtual function inherited from the abstract mother
               class vOptimizer.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      private int iOptimizerADMMLim::InitializeSpecific()
      \brief   This function is used to initialize specific stuff to the child optimizer.
      \details It is called by the public Initialize() function from the mother.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int InitializeSpecific();
    /*!
      \fn      private int iOptimizerADMMLim::SensitivitySpecificOperations()
      \param   FLTNB a_data
      \param   FLTNB a_forwardModel
      \param   FLTNB* ap_weight
      \param   FLTNB a_multiplicativeCorrections
      \param   FLTNB a_additiveCorrections
      \param   FLTNB a_blankValue
      \param   FLTNB a_quantificationFactor
      \param   oProjectionLine* ap_Line
      \brief   This function compute the weight associated to the provided event (for sensitivity computation)
      \details It is the implementation of the pure virtual function from vOptimizer. The result is
               put at ap_weight location.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                       FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                       FLTNB a_quantificationFactor, oProjectionLine* ap_Line );
    /*!
      \fn      private int iOptimizerADMMLim::DataSpaceSpecificOperations()
      \param   FLTNB a_data
      \param   FLTNB a_forwardModel
      \param   FLTNB* ap_backwardValues
      \param   FLTNB a_multiplicativeCorrections
      \param   FLTNB a_additiveCorrections
      \param   FLTNB a_blankValue
      \param   FLTNB a_quantificationFactor
      \param   oProjectionLine* ap_Line
      \param   int a_th
      \brief   This function performs the data space operations specific to the optimizer (computes the values
               to be backprojected)
      \details It is the implementation of the pure virtual function from vOptimizer. It is not useful in this optimizer
               this function is directly written in DataStep5ComputeCorrections function as m_vk and m_uk depend on the current thread
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                     FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                     FLTNB a_quantificationFactor, oProjectionLine* ap_Line);
    /*!
      \fn      private int iOptimizerADMMLim::ImageSpaceSpecificOperations()
      \param   FLTNB a_currentImageValue
      \param   FLTNB* ap_newImageValue
      \param   FLTNB a_sensitivity
      \param   FLTNB* ap_correctionValues
      \param   INTNB a_voxel
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \brief   This function perform the image update step specific to the optimizer
      \details It is the implementation of the pure virtual function from vOptimizer. The new image value is
               put at the ap_newImageValue location.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                      FLTNB a_sensitivity, FLTNB* ap_correctionValues,
                                      INTNB a_voxel, int a_tbf = -1, int a_rbf = -1, int a_cbf = -1 );


  // -------------------------------------------------------------------
  // Data members
  private:
    FLTNB**** m4p_firstDerivativePenaltyImage;     /*!< Image containing the penalty terms of the algorithm */
    FLTNB* mp_toWrite_vk;                          /*!< Sinogram v^{k+1} to be written at the end of optimizer computation, to use it in intricated ADMM */
    FLTNB* mp_toWrite_uk;                          /*!< Sinogram u^{k+1} to be written at the end of optimizer computation, to use it in intricated ADMM */
    FLTNB* m_uk;                                   /*!< Pointer containing the Sinogram u^k at index corresponding to current event for each thread, with v=Ax in second ADMM derivation */
    FLTNB* m_vk;                                   /*!< Pointer containing the Sinogram v^k at index corresponding to current event for each thread, with u the scaled dual variable in second ADMM derivation */
    HPFLTNB* m_AxProduct;                          /*!< Number containing the forward projection of x^{k+1} in intricated ADMM*/
    FLTNB* m_rBackgroundEvents;                    /*!< r_bar variable in intricated ADMM, to store background events (scatter + randoms) data to use it in v^{k+1} computation */
    FLTNB* m_yData;                                /*!< y variable in intricated ADMM, to store PET data to use it in v^{k+1} computation */
    FLTNB m_alpha;                                 /*!< alpha variable in intricated ADMM, weighting the strength of the Ax=v constraint */
    HPFLTNB m_grad_norm_sum;                       /*!< Copy of norm of gradient vector, to be used in conjugate gradient in current x computation */
    HPFLTNB m_proj_grad_norm_sum;                  /*!< Copy of norm of projection of gradient vector, to be used in conjugate gradient in current x computation */
    FLTNB* m_grad_before;                          /*!< Gradient vector to compute gradient norm in conjugate gradient in next x computation and needed to compute projection of it */
    HPFLTNB* m_proj_grad_before;                   /*!< Projection of gradient vector in conjugate gradient for each lor, to be stored for next x computation */

    // added variables for adaptive Rho
    FLTNB* mp_previousAx;
    FLTNB* mp_relPrimalResidual;
    FLTNB* mp_relDualResidual;
    // FLTNB* mp_PrimalResidual;
    // FLTNB* mp_DualResidual;
    FLTNB* mp_vectorAx;
    // FLTNB* mp_vectorV;
    // FLTNB* mp_vectorU;

    FLTNB m_square_sum_Ax;
    FLTNB m_square_sum_v;
    FLTNB m_square_sum_u;
    FLTNB m_square_sum_primal;
    FLTNB m_square_sum_dual;

    FLTNB m_adaptiveAlpha;
    FLTNB m_adaptiveTau;

    FLTNB m_xi;
    FLTNB m_mu;
    FLTNB m_tau;


};


// Class for automatic insertion (set here the visible optimizer's name as the first parameter,
// put the class name as the second parameter and do NOT add semi-colon at the end of the line)
CLASS_OPTIMIZER(ADMMLim_adp,iOptimizerADMMLim_adaptiveRho)

#endif

