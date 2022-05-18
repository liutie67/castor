
/*!
  \file
  \ingroup  optimizer
  \brief    Declaration of class iPenaltyDIP_ADMM
*/

#ifndef IPENALTYDIP_ADMM_HH
#define IPENALTYDIP_ADMM_HH 1

#include "vPenalty.hh"
#include "sAddonManager.hh"
#include "sOutputManager.hh"

/*!
  \class   iPenaltyDIP_ADMM
  \brief   This class is a DIP_ADMM for penalties
  \details This class inherits from vPenalty and provides details on how to implement a penalty.
*/
class iPenaltyDIP_ADMM : public vPenalty
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:

    /*!
      \fn      public iPenaltyDIP_ADMM::iPenaltyDIP_ADMM()
      \brief   The constructor of iPenaltyDIP_ADMM
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    iPenaltyDIP_ADMM();
    /*!
      \fn      public iPenaltyDIP_ADMM::~iPenaltyDIP_ADMM()
      \brief   The destructor of iPenaltyDIP_ADMM
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
    */
    ~iPenaltyDIP_ADMM();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameter and do not add semi-column at the end of the line)
    FUNCTION_PENALTY(iPenaltyDIP_ADMM)
    /*!
      \fn      public int iPenaltyDIP_ADMM::ReadConfigurationFile()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to the child penalty, from
               a configuration file. It is the implementation of the pure virtual function inherited
               from the abstract class vPenalty. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadConfigurationFile(const string& a_configurationFile);
    /*!
      \fn      public int iPenaltyDIP_ADMM::ReadOptionsList()
      \param   const string& a_optionsList
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to the child penalty, from
               a list of options. It is the implementation of the pure virtual function inherited
               from the abstract class vPenalty. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadOptionsList(const string& a_optionsList);
    /*!
      \fn      public int iPenaltyDIP_ADMM::ComputePenaltyValue()
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \param   INTNB a_voxel
      \param   int a_th
      \brief   Implementation of the pure virtual vPenalty::ComputePenaltyValue()
      \details This function computes the value of the penalty function for the provided indices.
               It is the implementation of the pure virtual vPenalty::ComputePenaltyValue().
      \return  The penalty value
    */
    FLTNB ComputePenaltyValue(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th);
    /*!
      \fn      public int iPenaltyDIP_ADMM::ComputeFirstDerivative()
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \param   INTNB a_voxel
      \param   int a_th
      \brief   Implementation of the pure virtual vPenalty::ComputeFirstDerivative()
      \details This function computes the first derivative of the penalty.
               It is the implementation of the pure virtual vPenalty::ComputeFirstDerivative().
      \return  The derivative
    */
    FLTNB ComputeFirstDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th);
    /*!
      \fn      public int iPenaltyDIP_ADMM::ComputeSecondDerivative()
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \param   INTNB a_voxel
      \param   int a_th
      \brief   Implementation of the pure virtual vPenalty::ComputeSecondDerivative()
      \details This function computes the second derivative of the penalty.
               It is the implementation of the pure virtual vPenalty::ComputeSecondDerivative().
      \return  The derivative
    */
    FLTNB ComputeSecondDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th);

  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private void iPenaltyDIP_ADMM::ShowHelpSpecific()
      \brief   A function used to show help about the child penalty
      \details This function must describe what the module does and how to use it. It describes in
               details the different parameters of the penalty, and how to set them through the use
               of a configuration file or a list of options. It is the implementation of the pure
               virtual function inherited from the abstract class vPenalty. It is called by the
               public ShowHelp() function.
    */
    void ShowHelpSpecific();
    /*!
      \fn      private int iPenaltyDIP_ADMM::CheckSpecificParameters()
      \brief   A private function used to check the parameters settings specific to the child penalty
      \details This function is used to check that all parameters specific to the penalty are correctly set
               within allowed values. It is called by the CheckParameters() function of the mother class.
               It is the implementation of the pure virtual function inherited from the abstract mother
               class vPenalty.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      private int iPenaltyDIP_ADMM::InitializeSpecific()
      \brief   This function is used to initialize specific stuff to the child penalty.
      \details It is called by the public Initialize() function from the mother.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int InitializeSpecific();
};

// Class for automatic insertion (set here the visible optimizer's name as the first parameter,
// put the class name as the second parameter and do NOT add semi-colon at the end of the line)
CLASS_PENALTY(DIP_ADMM,iPenaltyDIP_ADMM)

#endif

