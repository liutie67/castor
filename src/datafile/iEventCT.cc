
/*!
  \file
  \ingroup datafile

  \brief Implementation of class iEventCT
*/

#include "iEventCT.hh"
#include "vDataFile.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventCT::iEventCT() : vEvent()
{
  m_dataType = TYPE_SPECT;
  m_dataSpec = SPEC_TRANSMISSION;
  m_nbLines = 1;
  m_eventScatRate = 0.;
  m_eventBlankValue = 1.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventCT::~iEventCT() {}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iEventCT::Describe()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  Cout("iEventCT::Describe() -> Display contents" << endl);
  Cout("Time: " << m_timeInMs << " ms" << endl);
  Cout("Number of lines: " << m_nbLines << endl);
  for (uint16_t l=0; l<m_nbLines; l++) Cout("  --> ID1: " << mp_ID1[l] << " | ID2: " << mp_ID2[l] << endl);
  Cout("Scatter rate: " << m_eventScatRate << endl);
  Cout("Blank value: " << m_eventBlankValue << endl);
  Cout(flush);
}
