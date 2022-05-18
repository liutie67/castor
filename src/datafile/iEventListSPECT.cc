
/*!
  \file
  \ingroup datafile

  \brief Implementation of class iEventListSPECT
*/

#include "iEventListSPECT.hh"
#include "vDataFile.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventListSPECT::iEventListSPECT() : iEventSPECT()
{
  m_dataType = TYPE_SPECT;
  m_dataSpec = SPEC_EMISSION;
  m_dataMode = MODE_LIST;
  mp_POI[0] = 0.;
  mp_POI[1] = 0.;
  mp_POI[2] = -1.;
  m_kind = KIND_UNKNOWN;
  m_eventValue = 1.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventListSPECT::~iEventListSPECT() {}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iEventListSPECT::SetEventValue(int a_bin, FLTNBDATA a_value) 
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  Cout("***** iEventListSPECT::SetEventValue() -> Trying to set the value of a list mode event !");
  Exit(EXIT_FAILURE);
} 

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iEventListSPECT::Describe()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  Cout("iEventListSPECT::Describe() -> Display contents" << endl);
  Cout("Time: " << m_timeInMs << " ms" << endl);
  Cout("Number of lines: " << m_nbLines << endl);
  for (uint16_t l=0; l<m_nbLines; l++) Cout("  --> ID1: " << mp_ID1[l] << " | ID2: " << mp_ID2[l] << endl);
  Cout("Scatter rate: " << m_eventScatRate << endl);
  Cout("Normalization factor: " << m_eventNormFactor << endl);
  Cout("kind: " << m_kind << endl);
  Cout("POI (x ; y ; z ) = " << mp_POI[0] <<" ; " << mp_POI[1] <<" ; " << mp_POI[2] <<" ; " << endl);
  Cout(flush);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
