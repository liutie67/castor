
/*!
  \file
  \ingroup  datafile
  \brief    Declaration of class iEventHistoSPECT
*/

#ifndef IEVENTHISTOSPECT_HH
#define IEVENTHISTOSPECT_HH 1

#include "iEventSPECT.hh"

/*!
  \class   iEventHistoSPECT
  \brief   Inherit from iEventSPECT. Class for SPECT histogram mode events 
  \details It manages data and functions specific to histo mode SPECT.
*/
class iEventHistoSPECT : public iEventSPECT
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iEventHistoSPECT constructor. 
               Initialize the member variables to their default values.
    */
    iEventHistoSPECT();
    /*!
      \brief   iEventHistoSPECT destructor. 
    */
    ~iEventHistoSPECT();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      inline int iEventHistoSPECT::AllocateSpecificData()
      \brief   Function allowing the allocation of specific data. Return 0 by default for iEventHistoSPECT
      \return  0 is success, positive value otherwise
    */
    inline int AllocateSpecificData()
           {return 0;}
    /*!
      \fn      void iEventHistoSPECT::Describe()
      \brief   This function can be used to get a description of the event printed out
    */
    void Describe();

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      inline FLTNB iEventHistoSPECT::GetEventValue()
      \param   a_bin (0 if noTOF)
      \return  the event value (a_bin is dedicated to PET and ignored for SPECT)
    */
    inline FLTNB GetEventValue(int a_bin)
           {return m_eventValue;}
    /*!
      \fn      inline void iEventHistoSPECT::SetEventValue()
      \param   a_bin (0 if noTOF)
      \param   a_value
      \brief   Cast the FLTNBDATA value passed in parameters in FLTNB, and use it to set the event value (a_bin is dedicated to PET and ignored for SPECT)
    */
    inline void SetEventValue(int a_bin, FLTNBDATA a_value)
           {m_eventValue = (FLTNB)a_value;}
    /*!
      \fn      virtual INTNB vEvent::GetNbValueBins()
      \brief   Get the number of event value bins
      \return  Single value so 1 bin
    */
    inline INTNB GetNbValueBins()
           {return 1;}

  // -------------------------------------------------------------------
  // Private member functions
  private:

  // -------------------------------------------------------------------
  // Data members
  private:
};

#endif
