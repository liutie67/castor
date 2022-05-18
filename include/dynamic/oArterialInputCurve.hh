
/*!
  \file
  \ingroup  dynamic
  \brief    Declaration of class oArterialInputCurve
*/

#ifndef OARTERIALINPUTCURVE_HH
#define OARTERIALINPUTCURVE_HH 1


#include "gVariables.hh"

/*!
  \class   oArterialInputCurve
  \brief   This class is designed to manage the Arterial Input Curve provided by the user
  \details
*/
class oArterialInputCurve
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
     \fn      oArterialInputCurve::oArterialInputCurve
     \brief   Constructor of oArterialInputCurve. Simply set all data members to default values.
    */
    oArterialInputCurve();

    /*!
     \fn      oArterialInputCurve::~oArterialInputCurve
     \brief   Destructor of oArterialInputCurve. Free memory from all allocated tabs.
    */
    ~oArterialInputCurve();

  // -----------------------------------------------------------------------------------------
  // Public member functions
  public:
    /*!
     \fn      oArterialInputCurve::CheckParameters
     \brief   This function is used to check that the required input parameters are provided
              from the file provided by the user
     \return  0 if success, positive value otherwise.
    */
    int CheckParameters();
    /*!
     \fn      oArterialInputCurve::InitializeInputData
     \brief   This function is used to initialize the input arterial curve samples
              from the file provided by the user
     \return  0 if success, positive value otherwise.
    */
    int InitializeInputData();
    /*!
     \fn      oArterialInputCurve::InterpolateAIC
     \brief   This function performs a linear interpolation within the provided AIC range
              for every discrete unit of millisecond up to the end of the last requested frame
     \return  0 if success, positive value otherwise.
    */
    int InterpolateAIC();
    /*!
     \fn      oArterialInputCurve::Downsample
     \brief   This function downsamples the interpolated arterial input function
              Currently needed to speed up convolution for the SpectralDynamicModel
     \return  0 if success, positive value otherwise.
    */
    int Downsample();


  // -----------------------------------------------------------------------------------------
  // Public Get&Set functions
    /*!
     \fn      public inline void oArterialInputCurve::SetInputFilePath()
     \brief   Set path to file to get the sampled Arterial Input Curve
    */
    inline void SetInputFilePath(const string& a_pathToAICfile) { m_pathToAICfile = a_pathToAICfile; }

    /*!
     \fn      public inline void oArterialInputCurve::SetVerbose()
     \brief   Set the member m_verboseLevel to the provided value
    */
    inline void SetVerbose(int a_verbose) { m_verbose = a_verbose; }

    /*!
     \fn      public void oArterialInputCurve::SetFrames()
     \brief   Set the framing of the reconstruction for the AIC object
    */
    void SetFrames(int a_nbTimeFrames,uint32_t* a_frameTimeStartInMs, uint32_t* a_frameTimeStopInMs);
    /*!
     \fn      public void oArterialInputCurve::GetDownsampledAIC()
     \brief   Set the framing of the reconstruction for the AIC object
    */
    inline HPFLTNB* GetDownsampledAIC () { return mp_AICIntrpY_downsampled; }
    /*!
     \fn      public void oArterialInputCurve::GetInterpolatedAIC()
     \brief   Set the framing of the reconstruction for the AIC object
    */
    inline HPFLTNB* GetInterpolatedAIC () { return mp_AICIntrpY; }


  // -----------------------------------------------------------------------------------------
  // Private member functions

  // -----------------------------------------------------------------------------------------
  // Data members
  private:
    string m_pathToAICfile;                 /*!< The string containing the path to the input AIC file */

    int  m_nbinputDataPoints ;  	         /*!< Number of Arterial Input Function Datapoints */
    HPFLTNB* mp_AICDataPoints ;		         /*!< Pointer to array of Input Function Datapoints */
    uint32_t* mp_AICDataPointsTimes;        /*!< Pointer to array of Input Function Time points */

    HPFLTNB* mp_AICIntrpY;                  /*!< The interpolated points of Y axis (Bq/ml) of the Arterial Input Function
                                                 Each index relates to the time in milliseconds */
    HPFLTNB* mp_AICIntrpY_downsampled;      /*!< The interpolated points of Y axis (Bq/ml) f the Arterial Input Function Downsampled  */

    int m_nbTimeFrames;					     /*!< Number of frames*/
    uint32_t* mp_frameTimeStartInMs;		 /*!< Pointer to the array of frames Start Times */
    uint32_t* mp_frameTimeStopInMs;			 /*!< Pointer to the array of frames End Times */

    // Verbose
    int m_verbose;                          /*!< Verbose level */
};

#endif