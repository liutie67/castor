
/*!
  \file
  \ingroup datafile

  \brief Implementation of class iEventListCT
*/

#include "iEventListCT.hh"
#include "vDataFile.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventListCT::iEventListCT() : iEventCT()
{
  m_dataType = TYPE_CT;
  m_dataSpec = SPEC_TRANSMISSION;
  m_dataMode = MODE_LIST;
  m_kind = KIND_UNKNOWN;
  m_eventValue = 1.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventListCT::~iEventListCT() {}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iEventListCT::SetEventValue(int a_bin, FLTNBDATA a_value) 
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  Cout("***** iEventListCT::SetEventValue() -> Trying to set the value of a list mode event !");
  Exit(EXIT_FAILURE);
} 

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iEventListCT::Describe()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  Cout("iEventListCT::Describe() -> Display contents" << endl);
  Cout("Time: " << m_timeInMs << " ms" << endl);
  Cout("Number of lines: " << m_nbLines << endl);
  for (uint16_t l=0; l<m_nbLines; l++) Cout("  --> ID1: " << mp_ID1[l] << " | ID2: " << mp_ID2[l] << endl);
  Cout("Scatter rate: " << m_eventScatRate << endl);
  Cout("Blank value: " << m_eventBlankValue << endl);
  Cout("kind: " << m_kind << endl);
  Cout(flush);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
