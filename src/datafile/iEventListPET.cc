
/*!
  \file
  \ingroup datafile

  \brief Implementation of class iEventListPET
*/

#include "iEventListPET.hh"
#include "vDataFile.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventListPET::iEventListPET() : iEventPET()
{
  m_dataType = TYPE_PET;
  m_dataSpec = SPEC_EMISSION;
  m_dataMode = MODE_LIST;
  mp_POI1[0] = 0.;
  mp_POI1[1] = 0.;
  mp_POI1[2] = -1.;
  mp_POI2[0] = 0.;
  mp_POI2[1] = 0.;
  mp_POI2[2] = -1.;
  m_kind = KIND_UNKNOWN;
  m_eventScatRate = 0.;
  m_TOFMeasurementInPs = 0.;
  m_eventValue = 1.;
  m_nbLines = 1;
  m_TOFMeasurementRangeInPs = -1.;
  m_hasTOFInfo = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventListPET::~iEventListPET() {}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iEventListPET::SetEventValue(int a_bin, FLTNBDATA a_value) 
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  Cout("***** iEventListPET::SetEventValue() -> Trying to set the value of a list mode event !");
  Exit(EXIT_FAILURE);
} 

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iEventListPET::Describe()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  Cout("iEventListPET::Describe() -> Display contents" << endl);
  Cout("Time: " << m_timeInMs << " ms" << endl);
  Cout("Number of lines: " << m_nbLines << endl);
  for (uint16_t l=0; l<m_nbLines; l++) Cout("  --> ID1: " << mp_ID1[l] << " | ID2: " << mp_ID2[l] << endl);
  Cout("Random rate: " << m_eventRdmRate << endl);
  Cout("Normalization factor: " << m_eventNormFactor << endl);
  Cout("ACF: " << m_atnCorrFactor << endl);
//  Cout("kind: " << m_kind << endl);
  Cout("Scatter rate: " << m_eventScatRate << endl);
  Cout("TOF measurement: " << m_TOFMeasurementInPs << " ps" << endl);
  Cout("TOF measurement range: " << m_TOFMeasurementRangeInPs << " ps" << endl);
  Cout("POI1 (x ; y ; z ) = " << mp_POI1[0] <<" ; " << mp_POI1[1] <<" ; " << mp_POI1[2] <<" ; " << endl);
  Cout("POI2 (x ; y ; z ) = " << mp_POI2[0] <<" ; " << mp_POI2[1] <<" ; " << mp_POI2[2] <<" ; " << endl);
  Cout(flush);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iEventListPET::GetAdditiveCorrections(int a_bin)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  // If TOF, convert the LOR random rate into the random rate matching a single event with continuous TOF measurement
  // (division by the total range of TOF measurements)
  return m_eventScatRate + (m_hasTOFInfo ? (m_eventRdmRate/(m_TOFMeasurementRangeInPs*SPEED_OF_LIGHT_IN_MM_PER_PS*0.5)) : m_eventRdmRate); 
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
