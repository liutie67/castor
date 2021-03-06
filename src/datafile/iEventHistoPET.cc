
/*!
  \file
  \ingroup datafile

  \brief Implementation of class iEventHistoPET
*/

#include "iEventHistoPET.hh"
#include "vDataFile.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventHistoPET::iEventHistoPET() : iEventPET()
{
  m_dataType = TYPE_PET;
  m_dataSpec = SPEC_EMISSION;
  m_dataMode = MODE_HISTOGRAM;
  m_eventNbTOFBins = 1;
  mp_eventValue = NULL; 
  mp_eventScatRate = NULL;
  m_nbLines = 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventHistoPET::~iEventHistoPET() 
{
  if (mp_eventValue != NULL) free(mp_eventValue);
  if (mp_eventScatRate != NULL) free(mp_eventScatRate);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iEventHistoPET::AllocateSpecificData()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Check that the number of TOF bins is correct
  if (m_eventNbTOFBins<1)
  {
    Cerr("***** iEventHistoPET::AllocateSpecificData() -> Number of TOF bins has not been initialized (<1) !");
    return 1;
  }
  // Allocate values depending on the number of TOF bins
  mp_eventValue = (FLTNB*)malloc(m_eventNbTOFBins*sizeof(FLTNB));
  mp_eventScatRate = (FLTNB*)malloc(m_eventNbTOFBins*sizeof(FLTNB));
  // Deault initialization
  for(int tb=0 ; tb<m_eventNbTOFBins ; tb++)
  {
    mp_eventValue[tb] = 1.;
    mp_eventScatRate[tb] = 0.;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iEventHistoPET::GetAdditiveCorrections(int a_bin)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  return mp_eventScatRate[a_bin] + m_eventRdmRate/((FLTNB)m_eventNbTOFBins);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iEventHistoPET::Describe()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  Cout("iEventHistoPET::Describe() -> Display contents" << endl);
  Cout("Time: " << m_timeInMs << " ms" << endl);
  Cout("Number of lines: " << m_nbLines << endl);
  for (uint16_t l=0; l<m_nbLines; l++) Cout("  --> ID1: " << mp_ID1[l] << " | ID2: " << mp_ID2[l] << endl);
  Cout("Random rate: " << m_eventRdmRate << endl);
  Cout("Normalization factor: " << m_eventNormFactor << endl);
  Cout("ACF: " << m_atnCorrFactor << endl);
  Cout("Number of TOF bins: " << m_eventNbTOFBins << endl);
  for (uint16_t t=0; t<m_eventNbTOFBins; t++) Cout("  --> Event value: " << mp_eventValue[t] << " | Scatter rate: " << mp_eventScatRate[t] << endl);
  Cout(flush);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
