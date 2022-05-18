
/*!
  \file
  \ingroup  optimizer
  \brief    Declaration of class iOptimizerBSREM
*/

#ifndef IOPTIMIZERBSREM_HH
#define IOPTIMIZERBSREM_HH 1

#include "gVariables.hh"
#include "sAddonManager.hh"
#include "vOptimizer.hh"
#include "oImageSpace.hh"

#define BSREM_RELAXATION_CLASSIC 0
#define BSREM_NOT_DEFINED -1

/*!
  \class   iOptimizerBSREM
  \brief   This class implements the Block Sequential Regularized Expectation Maximization (BSREM) algorithm
  \details This class inherits from vOptimizer and implements the BSREM algorithm. This implementation follows the BSREM II algorithm defined in Ahn and Fessler 2003.
*/
class iOptimizerBSREM : public vOptimizer
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public iOptimizerBSREM::iOptimizerBSREM()
      \brief   The constructor of iOptimizerBSREM
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    iOptimizerBSREM();
    /*!
      \fn      public iOptimizerBSREM::~iOptimizerBSREM()
      \brief   The destructor of iOptimizerBSREM
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
    */
    ~iOptimizerBSREM();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameter and do not add semi-column at the end of the line)
    FUNCTION_OPTIMIZER(iOptimizerBSREM)
    /*!
      \fn      public int iOptimizerBSREM::ReadConfigurationFile()
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
      \fn      public int iOptimizerBSREM::ReadOptionsList()
      \param   const string& a_optionsList
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to the child optimizer, from
               a list of options. It is the implementation of the pure virtual function inherited
               from the abstract class vOptimizer. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadOptionsList(const string& a_optionsList);
    
  // -------------------------------------------------------------------
  // Private member functions (virtual in vOptimizer)
  private:
  
    /*!
      \fn      private int iOptimizerBSREM::PreImageUpdateSpecificStep()
      \brief   A private function used to compute the penalty term of the BSREM algorithm
      \details This function implements the virtual eponym function of vOptimizer.
               It computes the penalty term of the BSREM algorithm.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int PreImageUpdateSpecificStep();


  // -------------------------------------------------------------------
  // Private member functions (pure virtual in vOptimizer)
  private:
    /*!
      \fn      private void iOptimizerBSREM::ShowHelpSpecific()
      \brief   A function used to show help about the child optimizer
      \details This function must describe what the module does and how to use it. It describes in
               details the different parameters of the optimizer, and how to set them through the use
               of a configuration file or a list of options. It is the implementation of the pure
               virtual function inherited from the abstract class vOptimizer. It is called by the
               public ShowHelp() function.
    */
    void ShowHelpSpecific();
    /*!
      \fn      private int iOptimizerBSREM::CheckSpecificParameters()
      \brief   A private function used to check the parameters settings specific to the child optimizer
      \details This function is used to check that all parameters specific to the optimizer are correctly set
               within allowed values. It is called by the CheckParameters() function of the mother class.
               It is the implementation of the pure virtual function inherited from the abstract mother
               class vOptimizer.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      private int iOptimizerBSREM::InitializeSpecific()
      \brief   This function is used to initialize specific stuff to the child optimizer.
      \details It is called by the public Initialize() function from the mother.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int InitializeSpecific();
    /*!
      \fn      private int iOptimizerBSREM::SensitivitySpecificOperations()
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
      \fn      private int iOptimizerBSREM::DataSpaceSpecificOperations()
      \param   FLTNB a_data
      \param   FLTNB a_forwardModel
      \param   FLTNB* ap_backwardValues
      \param   FLTNB a_multiplicativeCorrections
      \param   FLTNB a_additiveCorrections
      \param   FLTNB a_blankValue
      \param   FLTNB a_quantificationFactor
      \param   oProjectionLine* ap_Line
      \brief   This function performs the data space operations specific to the optimizer (computes the values
               to be backprojected)
      \details It is the implementation of the pure virtual function from vOptimizer. The results to be
               backprojected is put at ap_backwardValues location.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                     FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                     FLTNB a_quantificationFactor, oProjectionLine* ap_Line );
    /*!
      \fn      private int iOptimizerBSREM::ImageSpaceSpecificOperations()
      \param   FLTNB a_currentImageValue
      \param   FLTNB* ap_newImageValue
      \param   FLTNB a_sensitivity
      \param   FLTNB* ap_correctionValues
      \param   INTNB a_voxel
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \brief   This function performs the image update step specific to the optimizer
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
    HPFLTNB m_minimumImageValue;               /*!< Threshold applied to the voxels to avoid non-positive values and enforce convergence */
    HPFLTNB m_maximumImageValue;               /*!< Threshold applied to the voxels to avoid too large values values and enforce convergence */
    FLTNB**** m4p_firstDerivativePenaltyImage; /*!< Image containing the penalty terms of the algorithm */
    int m_relaxationFactorType;                /*!< Type of relaxation factors */
    HPFLTNB m_relaxationFactorInitialValue;    /*!< Classic relaxation factor: initial value (alpha_0)*/
    HPFLTNB m_relaxationFactorStepSize;        /*!< Classic relaxation factor: step size (gamma) */
    HPFLTNB m_relaxationFactor;                /*!< Current value of the relaxation factor */

};

// Class for automatic insertion (set here the visible optimizer's name as the first parameter,
// put the class name as the second parameter and do NOT add semi-colon at the end of the line)
CLASS_OPTIMIZER(BSREM,iOptimizerBSREM)

#endif

