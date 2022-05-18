/*!
  \file
  \ingroup utils
  \brief This program convert a GATE datafile in root format to a list-mode datafile in CASToR format
  \todo Create class specific to GATE geometries (CylindricalPET, SPECT, etc..)
        in order to clean the code, as it currently looks like a huge plate of spaghettis...
*/

#include "gVariables.hh"
#include "gOptions.hh"
#include "iDataFilePET.hh"
#include "iDataFileSPECT.hh"
#include "sOutputManager.hh"
#include "sScannerManager.hh"
#include "gDataConversionUtilities.hh"
#include <iomanip> //std::setprecision

#ifdef CASTOR_ROOT
  #include "TROOT.h"
  #include "TApplication.h"
  #include "TGClient.h"
  #include "TCanvas.h"
  #include "TSystem.h"
  #include "TTree.h"
  #include "TBranch.h"
  #include "TFile.h"
  //#include "RDataFrame.hxx"
#endif


#include "oImageSpace.hh"
#include "oProjectorManager.hh"

/*!
  \fn      ShowHelp()
  \param a_returnCode
  \brief   Display main command line options for castor-GATERootToCastor
*/
void ShowHelp(int a_returnCode)
{
  cout << endl;
  cout << "Usage:  castor-GATERootToCastor      -i   path_to_ifile.root -OR- -il path_to_ifile.txt "<< endl;
  cout << "                                     -s   scanner_alias "<< endl;
  cout << "                                     -o   path_to_out_file "<< endl;
  cout << "                                     -m   path_to_macro.mac " << endl;
  cout << "                                    [-t   only convert true unscattered prompts]" << endl;
  cout << "                                    [-os  only convert true scattered prompts]" << endl;
  cout << "                                    [-or  only convert random prompts]" << endl;
  cout << "                                    [-ots only convert true scatter & unscattered prompts]" << endl;
  cout << "                                    [-otr only convert unscattered prompts]" << endl;
  //cout << "                                    [-sc  compute rates for scatter correction]" << endl;
  //cout << "                                    [-rc  compute rates for random correction]" << endl;
  //cout << "                                    [-src compute rates for scatter & random correction]" << endl;
  cout << "                                    [-cf        set calibration factor]" << endl;
  cout << "                                    [-oh        histogram datafile output]" << endl;
  cout << "                                    [-atn       path_to_atn_image(cm-1)]" << endl;
  cout << "                                    [-k         recover coincidence kind]" << endl;
  cout << "                                    [-ist       isotope alias" << endl;
  cout << "                                    [-TOF_reso  TOF resolution in ps]"  << endl;
  cout << "                                    [-TOF_range maximum range of TOF measurements]"  << endl;
  cout << "                                    [-isrc      image with source XYZ position]"  << endl;
  cout << "                                    [-geo       Generate a CASToR geom file"  << endl;
  cout << "                                    [-sp_bins   nT,nA trans/axial number of bins for SPECT gamma caméra]"  << endl;
  cout << "                                    [-otime     time offsets in s]"  << endl;
  cout << "                                    [-th        number of threads]"  << endl;
  cout << "                                    [-conf      path to config directory]"  << endl;
  cout << "                                    [-vb        verbosity]" << endl;
  cout << endl;
  cout << endl;
  cout << "[Main settings]:" << endl;
  cout << "  -i  path_to_file.root  : give an input root datafile to convert" << endl;
  cout << "  -il path_to_file.txt   : give an input text file containing a list of root files to convert" << endl;
  cout << "                         : they will be concatenated into 1 CASToR datafile" << endl;
  cout << "                         : the path to each datafile to convert must be entered on a newline" << endl;
  cout << "  -m  path_to_macro.mac  : gives the input GATE macro file used for the GATE simulation" << endl;
  cout << "  -o  path_to_out_file   : give the path to the output file which will be created inside this folder (no default)" << endl;
  cout << "  -s  scanner_alias      : provide the name of the scanner used for to acquire the original data" << endl;
  cout << "                         : Must correspond to a .geom or .hscan file in the config/scanner repository." << endl;
  cout << "                         : A geom file can be created from the macro files using the facultative option -geo (see below)" << endl;
  cout << endl;
  cout << "[Optional settings]:" << endl;
  cout << "  -t (or -ot)            : only the true (non-scattered) prompts will be converted" << endl;
  cout << "  -os                    : only the (Compton and Rayleigh) scattered coincidences will be converted" << endl;
  cout << "  -or                    : only the random coincidences will be converted" << endl;
  cout << "  -ots                   : only the true (scattered and unscattered) prompts will be converted" << endl;
  cout << "  -otr                   : only the non scattered prompts (true and random) will be converted" << endl;
  //cout << "  -sc                    : (PET only) scatter correction rates will be computed for each line of response" << endl;
  //cout << "  -rc                    : (PET only) random correction rates will be computed for each line of response" << endl;
  //cout << "  -src                   : (PET only) scatter and random correction rates will be computed for each line of response" << endl;
  cout << "  -cf                    : set a general calibration factor for the data" << endl;
  cout << "  -oh                    : Indicate if the output datafile should be written in histogram format (default : List-mode)" << endl;
  cout << "  -atn path_image:proj   : Provide an attenuation image (cm-1) related to the acquisition to estimate attenuation correction factors (acf)." << endl;
  cout << "                           Analytic projection will be performed during the data conversion in order to estimate PET attenuation correction factors (acf) for each histogram event" << endl;
  cout << "                           path_image : path to the attenuation image" << endl;
  cout << "                           proj       : (optional) projector to use for analytic projection. Defaut projector = Incremental siddon" << endl;
  cout << "                           Notes: " << endl;
  cout << "                           PET histogram output : This is required to perform attenuation correction" << endl;
  cout << "                           PET list-mode output : This is ONLY required if either scatter or random correction rates are estimated (-sc -rc, or -src options) " << endl;
  cout << "                                                  For PET list-mode, the attenuation image MUST be providen to the castor-recon exe as well, in order to perform sensitivity image computation" << endl;
  cout << "                           SPECT : No need to provide the attenuation image for the conversion as the acf will be computed during the reconstruction (atn image must be providen)" << endl;
  cout << "  -k                     : For List-mode output, write kind of coincidence (true/scatter/rdm/...) in the output datafile (disabled by default)" << endl;
  cout << "  -ist isotope_alias     : provide alias of the isotope used during the simulation"<< endl;
  cout << "                           Aliases for supported PET isotopes and their parameters are listed in config/misc/isotopes_pet"<< endl;
  cout << "                           Other isotopes can be added in the same file"<< endl;
  cout << "  -TOF_reso  reso_ps     : For list-mode PET, this option provides a value (in picoseconds) which will be interpreted as the system time-of-flight (ToF) resolution"<< endl;
  cout << "                           This usually corresponds to sqrt(2)*X, when X is the single time resolution defined in GATE."<< endl;
  cout << "                           ->Note that in Gate, there is an additional component S taking part in the time resolution, related to the way the information is stored (last photon interaction in the detector) and the effect of variable depth of interaction. The actual time resolution should be SQRT(2*X^2 + S^2)."<< endl;
  cout << "                             However S value is typically <100ps, and can be negligible depending on the value of X. For the simulation of ultra-fast timing resolution, one might need to consider modifying how the time is set in GATE (e.g the detection of the first photoelectron)"<< endl;
  cout << "                           The ToF Dt for each coincidence is computed from the time information recorded for each detectors"<< endl;
  cout << "  -TOF_range range_ps    : For list-mode PET, this option provides a value (in picoseconds) which will be interpreted as the maximum range of TOF measurements allowed by the scanner"<< endl;
  cout << "                           It is most often equal to the size of the coincidence timing windows. Thus if not provided by the user, the coincidence windown value set in the GATE macro file will be used instead"<< endl;
  cout << "                           If not found, then this value will be estimated from the data as (DtMax - DtMin)"<< endl;
  cout << "  -isrc path_to_img:dims : Provide name and dimensions (separated by a colon) of an image generated with the source (annihilation) XYZ position"<< endl;
  cout << "                           The option must start with the path to the output image which will be generated." << endl;
  cout << "                           Dimensions and voxel sizes of the image must be provided using commas as separators, as in the following template:" << endl;
  cout << "                           path/to/image:dimX,dimY,dimZ,voxSizeX,voxSizeY,voxSizeZ"<< endl;
  cout << "  -geo                   : Generate a CASToR geom file from the provided GATE macro file(s)"<< endl;
  cout << "                           A geom file with the 'scanner_alias' (from the -s option) basename will be generated in the scanner repository (default location : /config/scanner)" << endl;
  cout << "  -sp_bins nbinsT,nbinsA : Option specific to simulation using SPECThead systems, with root output."<< endl;
  cout << "                           Give transaxial and axial number of bins for projections, separated by a comma."<< endl;
  cout << "                           Pixel sizes will be computed from the crystal surface and the transaxial/axial number of bins" << endl;
  cout << "  -otime list            : Provide a serie of time offset in seconds to apply before each input file"<< endl;
  cout << "                           (In the case of converting several datafiles of a dynamic acquisition, timestamps of events may be resetted for each file" << endl;
  cout << "                           This variable allows to manually increment the time between each datafile(s) if required" << endl;
  cout << "                           The number of time offsets must be equal to the number of input files, provided by -i or -if options." << endl;
  cout << "                          'list' is a list of time offsets, separated by ','" << endl;
  #ifdef CASTOR_OMP
  cout << "  -th param              : Set the number of threads for parallel computation (default: 1). If 0 is given, then the maximum number of available threads is automatically determined." << endl;
  cout << "                           Can also give two parameters separated by a comma (e.g. 16,4), to distinguish between the number of threads for projection and image operations respectively." << endl;
  #endif
  cout << endl;
  cout << "[Miscellaneous settings]:" << endl;
  cout << "  -conf                  : Give the path to the config directory (default: located through the CONFIG_DIR environment variable)." << endl;
  cout << endl;
  cout << "  -vb                    : Give the verbosity level, from 0 (no verbose) to above 5 (at the event level) (default: 1)." << endl;
  cout << endl;
  #ifdef CASTOR_VERSION
  cout << "  This program is part of the CASToR release version " << CASTOR_VERSION << "." << endl;
  cout << endl;
  #endif
  Exit(a_returnCode);
}


// --- PET / SPECT structures declarations & definitions  --- //

/*!
  \class   Troot_pet
  \brief   This structure contains variables corresponding to a ROOT
           tree for CylindricalPET and ECAT system.
           Declared in castor-GATERootToCastor.cc
*/
struct Troot_pet
{
  // cylindricalPET specific
  int32_t rsectorID1;
  int32_t rsectorID2;
  int32_t moduleID1;
  int32_t moduleID2;
  int32_t submoduleID1;
  int32_t submoduleID2;
  int32_t layerID1;
  int32_t layerID2;
  // ecat specific
  int32_t blockID1; 
  int32_t blockID2;
  // common
  int32_t  crystalID1;
  int32_t  crystalID2;
  int32_t  eventID1; 
  int32_t  eventID2;
  int32_t  comptonPhantomID1;
  int32_t  comptonPhantomID2;
  int32_t  rayleighPhantomID1;
  int32_t  rayleighPhantomID2;
  double_t time1;
  double_t time2;
  int32_t  sourceID1; 
  int32_t  sourceID2;
  float_t  sourcePosX1;
  float_t  sourcePosY1;
  float_t  sourcePosZ1;
  bool     bool_ecat_mode;
};
/** @} */



/*!
  \class   Troot_spect
  \brief   This structure contains variables corresponding to a ROOT
           tree for SPECTHead system.
           Declared in castor-GATERootToCastor.cc
*/
struct Troot_spect
{
  // SPECThead specific
  float_t rotAngle;
  float_t gPosX;
  float_t gPosY;
  float_t gPosZ;
  int32_t headID;
  int32_t pixelID;
  // common
  int32_t  crystalID;
  int32_t  eventID;
  int32_t  comptonPhantomID;
  int32_t  rayleighPhantomID;
  double_t time;
  int32_t  sourceID;
  float_t  sourcePosX;
  float_t  sourcePosY;
  float_t  sourcePosZ;
};
/** @} */


// --- Function declarations. Definitions are located after main  --- //

/*!
  \fn InitTrootPet
  \param ap_rt : Structure to initialize
  \brief Init the variable of the root tree structure
*/
void InitTrootPet(Troot_pet* ap_rt);


/*!
  \fn InitTrootSpect
  \param ap_rt : Structure to initialize
  \brief Init the variable of the root tree structure
*/
void InitTrootSpect(Troot_spect* ap_rt);


#ifdef CASTOR_ROOT
  /*!
    \fn SetROOTPETBranch
    \param ap_tree : Root tree to use
    \param ap_rt : Structure to use for the branching
    \brief Set the branch of the PET root tree to the structure variables
  */
  void SetROOTPETBranch(TTree* ap_tree, Troot_pet* ap_rt);
  
  /*!
    \fn SetROOTSPECTBranch
    \param ap_tree : Root tree to use
    \param ap_rt : Structure to use for the branching
    \brief Set the branch of the SPECT root tree to the structure variables
  */
  void SetROOTSPECTBranch(TTree* ap_tree, Troot_spect* ap_rt);
#endif


/*!
  \fn RecoverSPECTGeometry
  \param a_inputDataIntf : path to datafile type
  \param a_pathIntf : path to the interfile header
  \param a_pathMac : path to a macro file
  \param distToDetector : distance between center of rotation and detector surface
  \param nHeads : nb of cameras
  \param nPixAxl : nb of transaxial pixels
  \param nPixTrs : nb of axial pixels
  \param crystalSizeAxl : crystal axial dimension
  \param crystalSizeTrs : crystal transaxial dimension
  \param nProjectionsTot : Nb total projections
  \param nProjectionsByHead : Nb total projections by head
  \param head1stAngle : head first angle
  \param headAngPitch : angular pitch between heads
  \param headRotSpeed : angle between projection
  \param headRotDirection : rotation direction of the head (CW or CCW)
  \param start_time_ms : acquisition start time converted in ms
  \param duration_ms : acquisition duration converted in ms
  \param nCrystalsTot : total nb detectors
  \param p_proj_spect_image : array for spect projection image
  \param spect_bin_axl : gamma camera axial (virtual) pixels
  \param spect_bin_trs : gamma camera transaxial (virtual) pixels
  \param nbSimulatedPixels : total nb virtual pixels
  \param a_histoOutFlag : output mode histogram flag
  \param p_bins : prompt histogram
  \param p_scat : scatter histograms
  \param bins_elts : total nb bins elements
  \param a_scatCorrFlag : scatter correction flag
  \param a_nbThreads,
  \param vb : verbosity
  \brief Call dedicated functions to recover SPECT geometry information from the provided datafile
  \return 0 if success, positive value otherwise
*/
int RecoverSPECTGeometry(bool       a_inputDataIntf,
                        string      a_pathIntf,
                        string      a_pathMac,
                        float_t     &distToDetector,
                        uint32_t    &nHeads,
                        uint32_t    &nPixelsAxial,
                        uint32_t    &nPixelsTransaxial,
                        float_t     &crystalSizeAxl,
                        float_t     &crystalSizeTrs,
                        uint32_t    &nProjectionsTot,
                        uint32_t    &nProjectionsByHead,
                        float_t     &head1stAngleDeg,
                        float_t     &headAngPitchDeg,
                        float_t     &headAngStepDeg,
                        int         &headRotDirection,
                        uint32_t    &start_time_ms,
                        uint32_t    &duration_ms,
                        uint32_t    &nCrystalsTot,
                        FLTNB*      &p_proj_spect_image,
                        uint16_t    &spect_bin_axl,
                        uint16_t    &spect_bin_trs,
                        uint32_t    &nbSimulatedPixels,
                        bool        a_histoOutFlag,
                        uint8_t  ***&p_bins,
                        uint8_t  ***&p_scat,
                        uint32_t    &bins_elts,
                        bool        a_scatCorrFlag,
                        int         a_nbThreads,
                        int         vb);


/*!
  \fn ProcessSPECTEvent
  \param Rt
  \param castorID1,
  \param castorID2,
  \param headAngPitchDeg,
  \param headAngStepDeg,
  \param nbSimulatedPixels,
  \param nPixelsTransaxial,
  \param nPixelsAxial,
  \param crystalSizeAxl,
  \param crystalSizeTrs,
  \param nProjectionsByHead,
  \param kind,
  \param nLORs_unknown,
  \param nLORs_trues,
  \param nLORs_scatters,
  \param nLORs_rdms,
  \param vb )
  \brief Compute CASToR IDs for a SPECT event
  \return 0 if success, positive value otherwise
*/
int ProcessSPECTEvent(  Troot_spect *Rt,
                        uint32_t   &castorID1,
                        uint32_t   &castorID2,
                        float_t     headAngPitchDeg,
                        float_t     headAngStepDeg,
                        uint32_t    nbSimulatedPixels,
                        uint32_t    nPixelsTransaxial,
                        uint32_t    nPixelsAxial,
                        float_t     crystalSizeAxl,
                        float_t     crystalSizeTrs,
                        uint32_t    nProjectionsByHead,
                        uint8_t    &kind,
                        uint64_t   &nLORs_unknown,
                        uint64_t   &nLORs_trues,
                        uint64_t   &nLORs_scatters,
                        uint64_t   &nLORs_rdms,
                        int         vb );



/*!
  \fn RecoverPETGeometry
  \param GATE_system_type,
  \param a_pathMac,
  \param nLayers,
  \param nb_crystal_per_layer,
  \param nCrystalsTot,
  \param nCrystalsAxial, 
  \param nCrystalsTransaxial, 
  \param nLayersRptAxial,
  \param nLayersRptTransaxial,
  \param nSubmodulesAxial, 
  \param nSubmodulesTransaxial, 
  \param nModulesAxial, 
  \param nModulesTransaxial, 
  \param nRsectorsAxial,
  \param nRsectorsAngPos,
  \param nBlocksLine, 
  \param nBlocksPerRing,
  \param invert_det_order,
  \param rsector_id_order,
  \param pet_coinc_window,
  \param start_time_ms,
  \param duration_ms,
  \param a_histoOutFlag,
  \param p_bins,
  \param p_scat,
  \param p_rdm,
  \param bins_elts,
  \param a_scatCorrFlag,
  \param a_rdmCorrFlag,
  \param a_nbThreads,
  \param vb : verbosity
  \brief Call dedicated functions to recover PET geometry information from the provided datafile
  \return 0 if success, positive value otherwise
*/
int RecoverPETGeometry( int         GATE_system_type,
                        string      a_pathMac,
                        uint8_t     &nLayers,
                        uint32_t    *nb_crystal_per_layer,
                        uint32_t    &nCrystalsTot,
                        uint32_t    &nCrystalsAxial, 
                        uint32_t    &nCrystalsTransaxial, 
                 vector<uint32_t>   &nLayersRptAxial,
                 vector<uint32_t>   &nLayersRptTransaxial,
                        uint32_t    &nSubmodulesAxial, 
                        uint32_t    &nSubmodulesTransaxial, 
                        uint32_t    &nModulesAxial, 
                        uint32_t    &nModulesTransaxial, 
                        uint32_t    &nRsectorsAxial,
                        uint32_t    &nRsectorsAngPos,
                        uint32_t    &nBlocksLine, 
                        uint32_t    &nBlocksPerRing,
                        bool        &invert_det_order,
                        int         &rsector_id_order,
                        FLTNB       &pet_coinc_window,
                        uint32_t    &start_time_ms,
                        uint32_t    &duration_ms,
                        bool        a_histoOutFlag,
                        uint8_t  ***&p_bins,
                        uint8_t  ***&p_scat,
                        uint8_t  ***&p_rdm,
                        uint32_t    &bins_elts,
                        bool        a_scatCorrFlag,
                        bool        a_rdmCorrFlag,
                        int         a_nbThreads,
                        int         vb
                        );
                        
                        
                        
                        
/*!
  \fn ProcessPETEvent
  \param GATE_system_type,
  \param Rt,
  \param castorID1,
  \param castorID2,
  \param nLayers,
  \param nb_crystal_per_layer,
  \param nCrystalsTot,
  \param nCrystalsAxial, 
  \param nCrystalsTransaxial, 
  \param nLayersRptAxial,
  \param nLayersRptTransaxial,
  \param nSubmodulesAxial, 
  \param nSubmodulesTransaxial, 
  \param nModulesAxial, 
  \param nModulesTransaxial, 
  \param nRsectorsAxial,
  \param nRsectorsAngPos,
  \param nBlocksLine, 
  \param nBlocksPerRing,
  \param invert_det_order,
  \param rsector_id_order,
  \param kind,
  \param nLORs_unknown,
  \param nLORs_trues,
  \param nLORs_scatters,
  \param nLORs_rdms, 
  \param vb 
  \brief Compute CASToR IDs for a PET event
  \return 0 if success, positive value otherwise
*/
int ProcessPETEvent(    int         GATE_system_type,
                        Troot_pet   *Rt,
                        uint32_t   &castorID1,
                        uint32_t   &castorID2,
                        uint8_t     nLayers,
                        uint32_t    *nb_crystal_per_layer,
                        uint32_t    nCrystalsTot,
                        uint32_t    nCrystalsAxial, 
                        uint32_t    nCrystalsTransaxial, 
                 vector<uint32_t>   nLayersRptAxial,
                 vector<uint32_t>   nLayersRptTransaxial,
                        uint32_t    nSubmodulesAxial, 
                        uint32_t    nSubmodulesTransaxial, 
                        uint32_t    nModulesAxial, 
                        uint32_t    nModulesTransaxial, 
                        uint32_t    nRsectorsAxial,
                        uint32_t    nRsectorsAngPos,
                        uint32_t    nBlocksLine, 
                        uint32_t    nBlocksPerRing,
                        bool        invert_det_order,
                        int         rsector_id_order,
                        uint8_t    &kind,
                        uint64_t   &nLORs_unknown,
                        uint64_t   &nLORs_trues,
                        uint64_t   &nLORs_scatters,
                        uint64_t   &nLORs_rdms, 
                        int         vb );
      








/*
  Main program
*/

int main(int argc, char** argv)
{
  // ============================================================================================================
  // MPI stuff (we make all instances but the first one returning 0 directly)
  // ============================================================================================================
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  int mpi_size = 1;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_rank!=0) return 0;
  #endif

  // No argument, then show help
  if (argc==1) ShowHelp(0);

  // ============================================================================================================
  // Parameterized variables with their default values
  // ============================================================================================================

  // String which gathers the path to the input data filename provided by the user. no default 
  string input_file = "";
  vector<string> path_to_input_file;
  
  // Path to user provided data and output files/images
  string path_to_out_file = "";
  string path_to_data_filename = "";
  string path_to_header_filename = "";
  string path_to_mac_file = "";
  string path_to_atn_image = "";
  string path_to_src_image = "";
          
  // Scanner name provided by the user
  string scanner_name = "";

  // GATE system type
  int GATE_system_type = GATE_SYS_UNKNOWN;
  
  // Specific coincidence recovery flag 
  bool true_only_flag = false;
  bool true_scat_flag = false;
  bool true_rdm_flag = false;
  bool scat_only_flag = false;
  bool rdm_only_flag = false;
  
  // Scat/Rdm correction flag 
  bool scat_corr_flag = false;
  bool rdm_corr_flag = false;
  
  // Verbosity
  int vb = 0;

  // PET ToF variables
  HPFLTNB calibration_factor = 1.;
  
  // Histogram output
  bool histo_out_flag = false;

  // Estimate acf (histogram output) 
  bool estimate_acf_flag = false;
  // Projector for analytic projection
  string options_projector  = "incrementalSiddon";
  // Image and voxel dimensions for optional analytic projections
  uint32_t proj_img_dim[3]  = {1,1,1};
  FLTNB    proj_vox_size[3] = {1.,1.,1.};

  // Recover kind of coincidence (list-mode output)
  bool kind_flag = false;

  // Isotope alias
  string istp_alias = "unknown";

  // Variables related to the source image
  bool src_img_flag = false;
  INTNB dim_src_img[3];
  FLTNB vox_src_img[3];
  FLTNB** p_src_img = NULL;

  // Generate geom file
  bool geom_flag = false;

  // Input is interfile (for SPECT simulation) -> not supported (def = false)
  bool input_is_intf = false;
  
  // SPECT bins
  uint16_t spect_bin_axl = 0,
           spect_bin_trs = 0;

  // PET ToF variables
  FLTNB tof_reso  = -1.;
  FLTNB tof_range = -1.;
  FLTNB pet_coinc_window = -1.;
  
  // Time offsets
  string offset_time_str = "";
  uint32_t* offset_time_list = NULL;

  // Path to config directory
  string path_to_config_dir = "";

  // Number of threads
  string nb_threads = "1";
  
  // ============================================================================================================
  // Read command-line parameters
  // ============================================================================================================
  for (int i=1; i<argc; i++)
  {
    string option = (string)argv[i];
    
    if (option=="-h" || option=="--help" || option=="-help") ShowHelp(0);
    
    // Just one file is provided
    if (option=="-i")
    {
      if (path_to_input_file.size() > 0)
      {
        Cerr("***** castor-GATERootToCastor :: the following file names have already been provided (-i/-il options must be used ONCE): " << endl);    
        for (size_t i=0 ; i<path_to_input_file.size() ; i++) 
          Cout(path_to_input_file[i] << endl); 
        Exit(EXIT_FAILURE);
      }
      else
      {
        if (argv[i+1] == NULL)
        {
          Cerr("***** castor-GATERootToCastor :: argument missing for option: " << option << endl);
          Exit(EXIT_FAILURE);
        }
        else
        {
          input_file = argv[i+1];
          path_to_input_file.push_back(input_file);
        }        
        i++;
      }
    }

    // Read a text file containing a list of root datafiles to read
    else if (option=="-il")
    {
      if (path_to_input_file.size() > 0)  
      {
        Cerr("***** castor-GATERootToCastor :: the following file names have already been provided (-i/-il options must be used ONCE) " << endl);   
        for (size_t i=0 ; i<path_to_input_file.size() ; i++) 
          Cout(path_to_input_file[i] << endl);  
        Exit(EXIT_FAILURE); 
      }
      else
      {
        if (argv[i+1] == NULL)
        {
          Cerr("***** castor-GATERootToCastor :: argument missing for option: " << option << endl);
          Exit(EXIT_FAILURE);
        }
        else
        {
          ifstream ifile(argv[i+1], ios::in);
          string line;
        
          // Read list of rootfgiles
          if(ifile)
          {
            // Get the position of the list_file, then append the name of the content datafiles to this position.
            input_file = argv[i+1];
          
            // Recover the files
            while (getline(ifile, line))
            {
              string path_to_datafile = GetPathOfFile(input_file);
              
              path_to_datafile += path_to_datafile.empty() ?  
                                  line:
                                  OS_SEP + line;
                                  
              path_to_input_file.push_back(path_to_datafile);
            }
          }
          else
          {
            Cerr("***** castor-GATERootToCastor :: Error, can't read txt file: " << argv[i+1] << endl);
            Exit(EXIT_FAILURE);
          }
     
          ifile.close();
        }
        i++;
      }
    }
    // Mac file
    else if (option=="-m")
    {
      path_to_mac_file = (string)argv[i+1];
      i++;
    }
    // Scanner alias
    else if (option=="-s")
    {
      if (argv[i+1] == NULL)
      {
        Cerr("***** castor-GATERootToCastor :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        scanner_name = argv[i+1];
      i++;
    }
    // Output CASToR datafile
    else if (option=="-o")
    {
      if (argv[i+1] == NULL)
      {
        Cerr("***** castor-GATERootToCastor :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        path_to_out_file = argv[i+1];
      i++;
    }
    // Recover only trues
    else if (option=="-t" || option=="-ot")
    {
      #ifdef CASTOR_ROOT
        true_only_flag = true;
      #else
        Cerr("***** castor-GATERootToCastor :: Option: " << option << " is only available for dataset generated with Gate in a Root format." << endl);
        Cerr("***** castor-GATERootToCastor :: Root support is currently not enabled (CASTOR_ROOT environnement variable should be set before compilation)" << endl);
        Exit(EXIT_FAILURE);
      #endif
    }
    // Recover trues and scatter
    else if (option=="-ots")
    {
      #ifdef CASTOR_ROOT
        true_scat_flag = true;
      #else
        Cerr("***** castor-GATERootToCastor :: Option: " << option << " is only available for dataset generated with Gate in a Root format." << endl);
        Cerr("***** castor-GATERootToCastor :: Root support is currently not enabled (CASTOR_ROOT environnement variable should be set before compilation)" << endl);
        Exit(EXIT_FAILURE);
      #endif
    }
    // Recover trues and randoms, unscattered
    else if (option=="-otr")
    {
      #ifdef CASTOR_ROOT
        true_rdm_flag = true;
      #else
        Cerr("***** castor-GATERootToCastor :: Option: " << option << " is only available for dataset generated with Gate in a Root format." << endl);
        Cerr("***** castor-GATERootToCastor :: Root support is currently not enabled (CASTOR_ROOT environnement variable should be set before compilation)" << endl);
        Exit(EXIT_FAILURE);
      #endif
    }
    // Recover scatter only
    else if (option=="-os")
    {
      #ifdef CASTOR_ROOT
        scat_only_flag = true;
      #else
        Cerr("***** castor-GATERootToCastor :: Option: " << option << " is only available for dataset generated with Gate in a Root format." << endl);
        Cerr("***** castor-GATERootToCastor :: Root support is currently not enabled (CASTOR_ROOT environnement variable should be set before compilation)" << endl);
        Exit(EXIT_FAILURE);
      #endif
    }
    // Recover random only
    else if (option=="-or")
    {
      #ifdef CASTOR_ROOT
        rdm_only_flag = true;
      #else
        Cerr("***** castor-GATERootToCastor :: Option: " << option << " is only available for dataset generated with Gate in a Root format." << endl);
        Cerr("***** castor-GATERootToCastor :: Root support is currently not enabled (CASTOR_ROOT environnement variable should be set before compilation)" << endl);
        Exit(EXIT_FAILURE);
      #endif
    }
    
    // THESE OPTIONS ARE DISABLED
    
    // Compute scatter correction
    /*
    else if (option=="-sc")
    {
      #ifdef CASTOR_ROOT
        scat_corr_flag = true;
      #else
        Cerr("***** castor-GATERootToCastor :: Option: " << option << " is only available for dataset generated with Gate in a Root format." << endl);
        Cerr("***** castor-GATERootToCastor :: Root support is currently not enabled (CASTOR_ROOT environnement variable should be set before compilation)" << endl);
        Exit(EXIT_FAILURE);
      #endif
    }
    // Recover random correction
    else if (option=="-rc")
    {
      #ifdef CASTOR_ROOT
        rdm_corr_flag = true;
      #else
        Cerr("***** castor-GATERootToCastor :: Option: " << option << " is only available for dataset generated with Gate in a Root format." << endl);
        Cerr("***** castor-GATERootToCastor :: Root support is currently not enabled (CASTOR_ROOT environnement variable should be set before compilation)" << endl);
        Exit(EXIT_FAILURE);
      #endif
    }
    // Recover random correction
    else if (option=="-src" || option=="-rsc" )
    {
      #ifdef CASTOR_ROOT
        scat_corr_flag = true;
        rdm_corr_flag = true;
      #else
        Cerr("***** castor-GATERootToCastor :: Option: " << option << " is only available for dataset generated with Gate in a Root format." << endl);
        Cerr("***** castor-GATERootToCastor :: Root support is currently not enabled (CASTOR_ROOT environnement variable should be set before compilation)" << endl);
        Exit(EXIT_FAILURE);
      #endif
    }
    */

    // PET ToF resolution 
    else if (option=="-TOF_reso")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-GATERootToCastor :: Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &tof_reso))
      {
        Cerr("***** castor-GATERootToCastor :: Exception when trying to value '" << argv[i+1] << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }    
    
    else if (option=="-TOF_range")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-GATERootToCastor :: Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &tof_range))
      {
        Cerr("***** castor-GATERootToCastor :: Exception when trying to value '" << argv[i+1] << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }    
    
    // Time offsets
    else if (option=="-otime")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-GATERootToCastor :: Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      offset_time_str = (string)argv[i+1];
      i++;
    }    
    
    // Provide a calibration factor
    else if (option=="-cf")
    {
      if (argv[i+1] == NULL)
      {
        Cerr("***** castor-GATERootToCastor :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &calibration_factor))
      {
        Cerr("***** castor-GATERootToCastor :: Exception when trying to value '" << argv[i+1] << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
        
      i++;
    }
    
    // Output datafile in histogram mode
    else if (option=="-oh")
    {
      histo_out_flag = true;
    }

    else if (option=="-atn")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-GATERootToCastor :: Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      
      path_to_atn_image = argv[i+1];
      
      // check first if a projector has been provided (':' included)
      size_t pos = path_to_atn_image.find_first_of(":");
              
      if (pos != string::npos)
      {
        options_projector = path_to_atn_image.substr(pos+1);
        path_to_atn_image = path_to_atn_image.substr(0,pos);
      }

      // Recover image dimensions from the provided atn image
      Intf_fields IF;
      IntfKeyInitFields(&IF);
      if(IntfReadHeader(path_to_atn_image, &IF, vb) )
      {
        Cerr("***** castor-GATERootToCastor :: An error occurred while trying to read the interfile header of attenuation file " << path_to_atn_image << " !" << endl);  
        Exit(1);
      }

      for(int d=0 ; d<3 ; d++)
      {
        proj_img_dim[ d ] = IF.mtx_size[ d ];
        proj_vox_size[ d ] = IF.vox_size[ d ];
      }

      estimate_acf_flag = true;
      i++;
    }
    
    // Output datafile in histogram mode
    else if (option=="-k")
    {
      kind_flag = true;
    }

    // Provide an isotope alias
    else if (option=="-ist")
    {
      if (argv[i+1] == NULL)
      {
        Cerr("***** castor-GATERootToCastor :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
        istp_alias = argv[i+1];
        
      i++;
    }

    // Generate image from the sourceID
    else if (option=="-isrc")
    {
      if (argv[i+1] == NULL)
      {
        Cerr("***** castor-GATERootToCastor :: argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
      {
        // Recover the path to the output image
        string input = argv[i+1];
        int pos_colon = input.find_first_of(":");
        path_to_src_image = input.substr(0,pos_colon);
        input = input.substr(pos_colon + 1);
        
        // Get string section related to dimensions
        string p_dims_str[6];
        if (ReadStringOption(input.c_str(), p_dims_str, 6, ",", option))
        {
          Cerr("***** castor-GATERootToCastor :: Invalid argument " << argv[i+1] << " for option " << option << " !" << endl);
          Exit(EXIT_FAILURE);
        }
        
        // Recover dimensions
        for(int i=0 ; i<3 ; i++)
          if (ConvertFromString(p_dims_str[i], &dim_src_img[i]))
          {
            Cerr("***** castor-GATERootToCastor :: Conversion error for elt " << p_dims_str[i] << " for option " << option << " !" << endl);
            Exit(EXIT_FAILURE);
          }

        // Recover dimensions
        for(int i=0 ; i<3 ; i++)
          if (ConvertFromString(p_dims_str[i+3], &vox_src_img[i]))
          {
            Cerr("***** castor-GATERootToCastor :: Conversion error for elt " << p_dims_str[i+3] << " for option " << option << " !" << endl);
            Exit(EXIT_FAILURE);
          }
        
        // Initilize source image
        src_img_flag = true;
      }
      i++;
    }

    // Generate geom file
    else if (option=="-geo")
    {
      geom_flag = true;
    }
    
    // SPECT bins
    else if (option=="-sp_bins")
    {
      string input = argv[i+1];
      uint16_t s_bins[2];
      if (ReadStringOption(input.c_str(), s_bins, 2, ",", option))
      {
        Cerr("***** castor-GATERootToCastor :: Invalid argument " << argv[i+1] << " for option " << option << " !" << endl);
        Exit(EXIT_FAILURE);
      }
      
      spect_bin_trs = s_bins[0];
      spect_bin_axl = s_bins[1];
      
      i++;
    }
   
    // Path to config directory
    else if (option=="-conf")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_to_config_dir = (string)argv[i+1];
      i++;
    }
    
    // Number of threads
    #ifdef CASTOR_OMP
    else if (option=="-th")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      
      nb_threads = (string)argv[i+1];
      i++;
    }
    #else
    else if (option=="-th")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      Cerr("!!!!! castor-recon() -> Option -th is available only if the code is compiled using the CASTOR_OMP environment variable set to 1 !" << endl);
      Cerr("!!!!!                   The execution will continue BUT WITH ONLY ONE THREAD !" << endl);
      Cerr("!!!!!                   We strongly advice to compile CASToR with OpenMP !" << endl);
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      i++;
    }
    #endif
    
    // Verbosity level
    else if (option=="-vb")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-GATERootToCastor :: Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      vb = atoi(argv[i+1]);
      i++;
    }

    else
    {
      Cerr("***** castor-GATERootToCastor :: Unknown option '" << option << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
  }


  // ============================================================================================================
  // Mandatory checks:
  // ============================================================================================================

  // Basic initialization checks (minimal initializations mandatory for the next steps)

  // data files
  if (path_to_input_file.empty() )
  {
    Cerr("***** castor-GATERootToCastor :: Please provide at least one data filename (-i or -if)" << endl);
    ShowHelp(0);
    Exit(EXIT_FAILURE);
  }
  else
  {
    /*
    // Check if we have interfile
    if(path_to_input_file.size() == 1)
    {
      string check;
      int rvalue = 0;
      rvalue = IntfKeyGetValueFromFile(path_to_input_file[0], "interfile", &check, 1, KEYWORD_OPTIONAL);
      
      if(rvalue == 1)
      {
        // error 
        Cerr("***** castor-GATERootToCastor :: Error when trying to read file: " << path_to_input_file[0] << "!" << endl);
        Exit(EXIT_FAILURE);
      }
      else if(rvalue == 0)
        // key found, we have an interfile as input
        input_is_intf = true;
    }
    */
    if(vb >= 2)
    { 
      Cout(" Selected root data-file(s) to convert: " << endl);
      for (size_t i=0 ; i<path_to_input_file.size() ; i++) 
        Cout(path_to_input_file[i] << endl);
    }
  }



  if(!input_is_intf)
  {
    // Check ROOT is enabled if the input file is not interfile (SPECT)
    #ifndef CASTOR_ROOT
    Cerr("***** castor-GATERootToCastor :: CASToR must be compiled with ROOT to read input root files (check installation instructions)" << endl);
    Exit(EXIT_FAILURE);
    #endif
  }
  else if(src_img_flag) // intf input + image of the source -> error
  {
    Cerr("***** castor-GATERootToCastor :: Can't use -isrc with interfile dataset ");
    Cerr(" (image of the source can only be generated from root file) !" << endl);
    Exit(EXIT_FAILURE);
  }
  
  // Check file(s) existence
  for (size_t i=0 ; i<path_to_input_file.size() ; i++) 
  {
    ifstream input_file(path_to_input_file[i], ios::in);
    if (!input_file)
    {
      Cerr("***** castor-GATERootToCastor :: Couldn't find or read data-file '"<< path_to_input_file[i] << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
  }
  
  // output directory
  if (path_to_out_file.empty() )
  {
    Cerr("***** castor-GATERootToCastor :: Please provide the output file name (-o)" << endl);
    ShowHelp(0);
    Exit(EXIT_FAILURE);
  }
  else
    if(vb >= 2) Cout(" selected output file:" << path_to_out_file << endl);

  // macro
  if (path_to_mac_file.empty())
  {
    Cerr("***** castor-GATERootToCastor :: Please provide the macro file associated to the GATE root datafile (-m) :" << endl);
    ShowHelp(0);
    Exit(EXIT_FAILURE);
  }
  else
    if(vb >= 2) Cout(" selected macro file: " << path_to_mac_file << endl);
    
  // scanner
  if (scanner_name.empty())
  {
    Cerr("***** castor-GATERootToCastor :: Please provide a scanner alias (-s) :" << endl);
    ShowHelp(0);
    Exit(EXIT_FAILURE);
  }
  else
    if(vb >= 2) Cout(" selected scanner alias: " << scanner_name << endl);


  #ifdef CASTOR_ROOT
  
  // No tof with scat correction rates
  if (tof_reso>0 && scat_corr_flag)
  {
    Cerr("***** castor-GATERootToCastor :: ToF is not currently compatible with the computation of rates for scatter correction" << endl);
    Cerr("*****                            ( -TOF_reso can't be used with -sc )" << endl);
    Exit(EXIT_FAILURE);
  }
  
  if (  (true_only_flag + true_scat_flag + true_rdm_flag + scat_only_flag + rdm_only_flag) > 1 )
  {
    Cerr("***** castor-GATERootToCastor :: Only use one of this option in the command line according to the type of coincidence to convert: -t -ts -tr -so -ro" << endl);
    Exit(EXIT_FAILURE);
  }
  #endif
  
  // No tof with histogram output
  if (tof_reso>0 &&  histo_out_flag)
  {
    Cerr("***** castor-GATERootToCastor :: ToF is not currently compatible with the histogram datafile output, only with list-mode!" << endl);
    Cerr("*****                            ( -TOF_reso can't be used with -oh )" << endl);
    Exit(EXIT_FAILURE);
  }
  
  // (Required for options using interfile I/O)
  GetUserEndianness();
    
    
  // ============================================================================================================
  // SOutputManager object initialisation:
  // ============================================================================================================
    
  sOutputManager* p_outputManager = sOutputManager::GetInstance();  
  // Set verbose level
  p_outputManager->SetVerbose(vb);
  // Set MPI rank
  p_outputManager->SetMPIRank(0);

  // Set path to the config directory
  if (p_outputManager->CheckConfigDir(path_to_config_dir))
  {
    Cerr("***** castor-GATERootToCastor :: A problem occurred while checking for the config directory path !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize output directory and base name
  if (p_outputManager->InitOutputDirectory(path_to_out_file, ""))
  {
    Cerr("*****castor-GATERootToCastor ::  A problem occurred while initializing output directory !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Log command line
  if (p_outputManager->LogCommandLine(argc,argv))
  {
    Cerr("***** castor-GATERootToCastor :: A problem occurred while logging command line arguments !" << endl);
    Exit(EXIT_FAILURE);
  }
  
  // Output progression options
  cout << std::fixed;
  cout << std::setprecision(2);
  

  // ============================================================================================================
  // Check system type from the macro file
  // ============================================================================================================
  GATE_system_type = GetGATESystemType(path_to_mac_file);

  if(GATE_system_type<0) 
  {
    // Unknown system
    Cerr("***** castor-GATERootToCastor :: GATE system type not found : " << endl);
    Cerr("                           This script only supports conversion for cylindricalPET, ecat and SPECThead systems"  << endl);
    Cerr("                           The system type is recovered from the lines '/gate/systems/...'"  << endl);
    Exit(EXIT_FAILURE);  
  }


  // ============================================================================================================
  // Generate the geom file from the mac file(s) is required
  // ============================================================================================================
  
  if(geom_flag) 
  {
    string scanner_repository = sOutputManager::GetInstance()->GetPathToConfigDir() + "scanner" + OS_SEP;
    string path_to_geom = scanner_repository + scanner_name + ".geom";
         
    // Check system type 
    switch ( GATE_system_type ) 
    {
      case GATE_SYS_CYLINDRICAL:
        if( vb>=2 )Cout(endl << " --- CylindricalPET system detected. Proceeding to conversion... --- " << endl << endl);
        
        if(CreateGeomWithCylindrical(path_to_mac_file , path_to_geom) )
        {
          Cerr("***** castor-GATERootToCastor :: An error occurred while trying to process mac file for cylindrical system: " << path_to_mac_file << endl);
          Exit(EXIT_FAILURE);
        }
        break;
        
      case GATE_SYS_ECAT:
        if( vb>=2 )Cout(endl << " --- ECAT system detected. Proceeding to conversion... --- " << endl << endl);
        
        if(CreateGeomWithECAT(path_to_mac_file , path_to_geom) )
        {
          Cerr("***** castor-GATEMacToGeom :: An error occurred while trying to process mac file for ecat system: " << path_to_mac_file  << endl);
          Exit(EXIT_FAILURE);
        }
        break;
      // TODO
      case GATE_SYS_SPECT:
        if( vb>=2 )Cout(endl << " --- SPECThead system detected. Proceeding to conversion... --- " << endl << endl);
        
        if(CreateGeomWithSPECT(path_to_mac_file , path_to_geom) )
        {
          Cerr("***** castor-GATEMacToGeom :: An error occurred while trying to process mac file for SPECT system: " << path_to_mac_file  << endl);
          Exit(EXIT_FAILURE);
        }
        break;
      
      default: // Unknown system
        Cerr("***** castor-GATERootToCastor :: System type not found : " << endl);
        Cerr("                   This script only supports conversion for cylindricalPET ecat and SPECThead systems"  << endl);
        Cerr("                   The system type is recovered from the lines '/gate/systems/...'"  << endl);
        Exit(EXIT_FAILURE);  
        break;
    }
      
    if( vb>=2 )Cout(endl << " --- Conversion completed --- " << endl << endl);
  }


  // ============================================================================================================
  // ScannerManager object initialisation:
  // ============================================================================================================
  
  sScannerManager* p_scannermanager = sScannerManager::GetInstance();  
  p_scannermanager->SetVerbose(vb);
    
  // Check if the scanner exists and get the name from ScannerManager
  scanner_name = (scanner_name.find(OS_SEP)) ? 
                  scanner_name.substr(scanner_name.find_last_of(OS_SEP)+1) :
                  scanner_name;

  if(p_scannermanager->FindScannerSystem(scanner_name) )
  {
    Cerr("**** castor-GATERootToCastor :: A problem occurred while searching for scanner system !" << endl);
    Exit(EXIT_FAILURE);
  } 

  // Build output file names
  path_to_data_filename = path_to_out_file + ".Cdf";
  path_to_header_filename = path_to_out_file + ".Cdh";
  


  // ============================================================================================================
  // Declaration and Allocation of oImageDimensions object 
  // ============================================================================================================
  
  oImageDimensionsAndQuantification* p_ID = new oImageDimensionsAndQuantification(); 
  
  // --- oImageDimensionsAndQuantification initialization ---
  p_ID->SetNbVoxX(proj_img_dim[ 0 ]);
  p_ID->SetNbVoxY(proj_img_dim[ 1 ]);
  p_ID->SetNbVoxZ(proj_img_dim[ 2 ]);
  p_ID->SetNbBeds(1);
  p_ID->SetVoxSizeX(proj_vox_size[ 0 ]);
  p_ID->SetVoxSizeY(proj_vox_size[ 1 ]);
  p_ID->SetVoxSizeZ(proj_vox_size[ 2 ]);
  p_ID->SetFOVOutMasking(0., 0);
  p_ID->SetFOVSizeX(-1.);
  p_ID->SetFOVSizeY(-1.);
  p_ID->SetFOVSizeZ(-1.);
  p_ID->SetOffsetX(0);
  p_ID->SetOffsetY(0);
  p_ID->SetOffsetZ(0);
  p_ID->SetVerbose(vb);
  p_ID->SetNbRespGates(1);
  p_ID->SetNbCardGates(1);
  p_ID->SetFrames("");  
  p_ID->SetNbThreads(nb_threads);
  
  if (p_ID->CheckParameters())
  {
    Cerr("***** castor-GATERootToCastor :: A problem occurred while checking image dimensions parameters !" << endl);
    Exit(1);
  }
  if (p_ID->Initialize())
  {
    Cerr("***** castor-GATERootToCastor :: A problem occurred while initializing image dimensions !" << endl);
    Exit(1);
  }
  // Initialization of DynamicDataManager class, related 4D data splitting management 
  if (p_ID->InitDynamicData( "", 0, 0, 0, 1, 1 ) )
  {
    Cerr("***** castor-GATERootToCastor :: A problem occurred while initializing Dynamic data manager's class !" << endl);
    Exit(EXIT_FAILURE);
  }




  // ============================================================================================================
  // Instanciate/Initialize CASToR DataFile
  // ============================================================================================================
    
  // Instanciate & Initialize iDataFilePET and Event objects
  vDataFile* Out_data_file = NULL;
  //vEvent** Event = NULL;

  uint16_t max_nb_lines_per_event = 1; // No compression for root files
  
  if(GATE_system_type == GATE_SYS_SPECT)
  {
    Out_data_file = new iDataFileSPECT();
    iDataFileSPECT* p_datafile = (dynamic_cast<iDataFileSPECT*>(Out_data_file));
    p_datafile->SetDataType(TYPE_SPECT);
    p_datafile->SetIsotope(istp_alias);
    histo_out_flag = true; // force histogram output for SPECT
    
    if(scat_corr_flag)
      p_datafile->SetScatterCorrectionFlagOn();
  }
  else
  {
    Out_data_file = new iDataFilePET();
    iDataFilePET* p_datafile = (dynamic_cast<iDataFilePET*>(Out_data_file));
    p_datafile->SetDataType(TYPE_PET);
    p_datafile->SetIsotope(istp_alias);
    
    // Scatter and Random correction estimation
    if(scat_corr_flag)
      p_datafile->SetScatterCorrectionFlagOn();
    if(rdm_corr_flag)
      p_datafile->SetRandomCorrectionFlagOn();
            
    // ACF computed
    if(estimate_acf_flag)
      p_datafile->SetAtnCorrectionFlagOn();
    
    p_datafile->SetMaxNumberOfLinesPerEvent(max_nb_lines_per_event);
  }
    
  Out_data_file->SetImageDimensionsAndQuantification(p_ID);
  Out_data_file->SetHeaderDataFileName(path_to_header_filename);
  Out_data_file->SetVerbose(0);

  // Histogram
  if(histo_out_flag)
    Out_data_file->SetDataMode(MODE_HISTOGRAM);
  // List-mode 
  else 
  {
    Out_data_file->SetDataMode(MODE_LIST);
    
    // Specific list-mode metadata
    if(GATE_system_type == GATE_SYS_SPECT )
    {
      if(kind_flag)
        ((iDataFileSPECT*)Out_data_file)->SetEventKindFlagOn();
    }
    else
    {
      if(kind_flag)
        ((iDataFilePET*)Out_data_file)->SetEventKindFlagOn();
      // Set ToF
      if(tof_reso>=0) 
      {
        ((iDataFilePET*)Out_data_file)->SetTOFInfoFlag();
        ((iDataFilePET*)Out_data_file)->SetTOFResolutionInPs(tof_reso);
        // Set it to an acceptable value (>0) in order to perform projector initialization
        // (Projection with TOF kernel will be disabled anyway)
        // The actual value will be set after conversion during datafile header creation
        ((iDataFilePET*)Out_data_file)->SetTOFMeasurementRangeInPs(1000);
      }
    }
  }
  
  Out_data_file->PROJ_InitFile();
  Out_data_file->SetCalibrationFactor(calibration_factor);
  Out_data_file->ComputeSizeEvent();
  Out_data_file->PrepareDataFile();
  
  
  
  
  // ============================================================================================================
  // If acf estimation is enabled for histogram output, initialize all required objects for analytic projection
  // (Image space, projector and scanner geometry)
  // ============================================================================================================
  oProjectorManager* p_ProjectorManager = new oProjectorManager();
  oImageSpace* p_ImageSpace = new oImageSpace();
  
    
  if(estimate_acf_flag)
  {
    // Check if system is SPECT, throw error in this case 
    if(GATE_system_type == GATE_SYS_SPECT)
    {
      Cerr("***** castor-GATERootToCastor :: Estimation of acf from an attenuation image (-atn option) only available for PET systems ! (detected system is SPECT)" << endl);
      Exit(1);
    }
    
    // --- Image space initialization ---
    p_ImageSpace->SetImageDimensionsAndQuantification(p_ID);
    p_ImageSpace->SetVerbose(vb);
  
    // Allocate memory for main image
    p_ImageSpace->LMS_InstantiateImage();

    // Read attenuation image
    if(p_ImageSpace->PROJ_LoadInitialImage(path_to_atn_image) )
    {
      Cerr("***** castor-GATERootToCastor :: Error during image initialization !" << endl);  
      Exit(EXIT_FAILURE);
    }


    // --- Build Scanner geometry ---
    if(p_scannermanager->BuildScannerObject() )
    {
      Cerr("**** castor-GATERootToCastor :: A problem occurred during scanner object construction !" << endl);
      Exit(EXIT_FAILURE);
    }
    
    if(p_scannermanager->InstantiateScanner() )
    {
      Cerr("**** castor-GATERootToCastor :: A problem occurred while creating Scanner object !" << endl);
      Exit(EXIT_FAILURE);
    } 
    
    if(p_scannermanager->BuildLUT() )
    {
      Cerr("***** castor-GATERootToCastor :: A problem occurred while generating/reading the LUT !" << endl);
      Exit(EXIT_FAILURE);
    }
    
    // Check the scanner manager parameters and initialize the scanner
    if (p_scannermanager->CheckParameters())
    {
      Cerr("***** castor-GATERootToCastor :: A problem occurred while checking scanner manager parameters !" << endl);
      Exit(1);
    }

    if (p_scannermanager->Initialize())
    {
      Cerr("***** castor-GATERootToCastor :: A problem occurred while initializing scanner !" << endl);
      Exit(1);
    }
    
    
    // --- Projector Manager initialization---  
    p_ProjectorManager->SetScanner(p_scannermanager->GetScannerObject());
    p_ProjectorManager->SetImageDimensionsAndQuantification(p_ID);
    p_ProjectorManager->SetDataFile(Out_data_file);
    p_ProjectorManager->SetComputationStrategy(FIXED_LIST_COMPUTATION_STRATEGY);
    p_ProjectorManager->SetOptionsForward(options_projector);
    p_ProjectorManager->SetOptionsBackward(options_projector);
    p_ProjectorManager->SetOptionsCommon("");
    p_ProjectorManager->SetVerbose(vb);

    // Check parameters
    if (p_ProjectorManager->CheckParameters())
    {
      Cerr("***** castor-GATERootToCastor :: A problem occurred while checking projector manager's parameters !" << endl);
      Exit(EXIT_FAILURE);
    }
  
    // Initialize projector manager
    if (p_ProjectorManager->Initialize())
    {
      Cerr("***** castor-GATERootToCastor :: A problem occurred while initializing projector manager !" << endl);
      Exit(EXIT_FAILURE);
    }
    p_ProjectorManager->SetSensitivityModeOn();
  }
  
  
  // Check isotope existence
  
  if(!istp_alias.empty())
  {
    if(p_ID->SetPETIsotope(0, istp_alias) )
    {
      Cerr("**** castor-GATERootToCastor :: A problem occurred while checking isotope name !" << endl);
      Exit(EXIT_FAILURE);  
    }
  }
  
  
  
  // ============================================================================================================
  // Parse Time offsets
  // ============================================================================================================
  offset_time_list = new uint32_t[path_to_input_file.size()];
  for(size_t f=0 ; f<path_to_input_file.size() ; f++)
    offset_time_list[f] = 0;
  
  // Parse offset_time_str, if it contains any data
  if(offset_time_str != "")
  {
    vector<string> offsets;
    size_t comma_pos = 0;
    while ((comma_pos=offset_time_str.find_first_of(",")) != string::npos)
    {
      string offset = offset_time_str.substr(0,comma_pos);
      offsets.push_back(offset);
      offset_time_str = offset_time_str.substr(comma_pos+1);
    }
    
    // Check we have a correct number of offsets
    if(offsets.size() != path_to_input_file.size())
    {
      Cerr("**** castor-GATERootToCastor :: Unmatching number of offsets with -ot option ! "
           << offsets.size() << " have been found, while "<< path_to_input_file.size() <<"input file have been provided !" << endl);
      Exit(EXIT_FAILURE);
    }
    
    for(size_t o=0 ; o<offsets.size() ; o++)
    {
      double offset;
      if(ConvertFromString(offsets[o] , &offset) )
      {
        Cerr("**** castor-GATERootToCastor :: Error while trying to convert offset : "<< offsets[o] <<" in ms ! " << endl);
        Exit(EXIT_FAILURE);
      }
      
      // convert in ms
      offset_time_list[o] = (uint32_t)offset*1000;
    }
  }
  
  // ============================================================================================================
  // *********************************************** CONVERSION *************************************************
  // ============================================================================================================
  Cout(" --- Start conversion of datafile(s) : " << input_file <<  endl
    << "     using mac file: " << path_to_mac_file << endl
    << "     CASToR output header datafile: " << path_to_header_filename << endl
    << "     CASToR output binary datafile: " << path_to_data_filename << endl << endl);

  // ============================================================================================================
  // Variables declarartion & initialization
  // ============================================================================================================
  
  // General clock
  clock_t clock_start_main = clock();
  time_t time_start_main = time(NULL);
  
  // Scanner variables (General)
  uint32_t nCrystalsTot = 0;
  uint32_t start_time_ms = 0; // 0 as default initialization
  uint32_t duration_ms = 1000; // 1 second as default initialization
  uint32_t bins_elts = 0;
  
  // Scanner variables (PET)
  uint32_t nRsectorsAngPos = 1, nRsectorsAxial = 1;
  // index order of rsectors (transaxial/axial) for cylindricalPET depending on the GATE macro
  int      rsector_id_order = 0;
  
  // orientation of detector ordering depending on the GATE macro
  bool     invert_det_order = false;
  uint32_t nModulesTransaxial = 1, nModulesAxial = 1;
  uint32_t nSubmodulesTransaxial = 1, nSubmodulesAxial = 1;
  uint32_t nCrystalsTransaxial = 1, nCrystalsAxial = 1;
  uint32_t nBlocksLine = 1, nBlocksPerRing = 1;
  uint8_t  nLayers = 1;
  uint32_t nb_crystal_per_layer[4] = {0,0,0,0};
  
  // Scanner variables (SPECT)
  uint32_t nProjectionsByHead =1;
  uint32_t nProjectionsTot = 1;
  uint32_t nHeads =1;
  uint32_t nPixelsAxial = 1;
  uint32_t nPixelsTransaxial = 1;
  uint32_t nbSimulatedPixels = 1;
  float_t  distToDetector = 0.;
  int      headRotDirection = GEO_ROT_CW;
  float_t  head1stAngleDeg = 0;
  float_t  headAngPitchDeg = -1;
  float_t  headAngStepDeg = -1;
  float_t  crystalSizeAxl=-1., 
           crystalSizeTrs=-1.;
  FLTNB*   p_proj_spect_image=NULL;
  
  // layer repeaters
  vector<uint32_t> nLayersRptTransaxial, nLayersRptAxial; 
  
  
  // Counter for the number of events (and the number of trues, scatters and randoms for simulated data)
  uint64_t  nLORs_tot = 0;  
  uint64_t *nLORs_trues = new uint64_t[ p_ID->GetNbThreadsForProjection() ];
  uint64_t *nLORs_rdms = new uint64_t[ p_ID->GetNbThreadsForProjection() ];
  uint64_t *nLORs_scatters = new uint64_t[ p_ID->GetNbThreadsForProjection() ];
  uint64_t *nLORs_unknown = new uint64_t[ p_ID->GetNbThreadsForProjection() ];
  uint64_t *nBINs = new uint64_t[ p_ID->GetNbThreadsForProjection() ];
  
  
  // Castor variables
  uint8_t*** p_bins = NULL;//new uint8_t[ p_ID->GetNbThreadsForProjection() ];
  uint8_t*** p_scat = NULL;//new uint8_t[ p_ID->GetNbThreadsForProjection() ];
  uint8_t*** p_rdm = NULL;//new uint8_t[ p_ID->GetNbThreadsForProjection() ];
  HPFLTNB  *dt_max = new HPFLTNB[ p_ID->GetNbThreadsForProjection() ];
  HPFLTNB  *dt_min = new HPFLTNB[ p_ID->GetNbThreadsForProjection() ];
  
  if( src_img_flag )
  {
    p_src_img = new FLTNB*[ p_ID->GetNbThreadsForProjection() ];
    
    for (int th=0 ; th < p_ID->GetNbThreadsForProjection() ; th++)
    {
      p_src_img[ th ] = new FLTNB[dim_src_img[0]*dim_src_img[1]*dim_src_img[2]];
      
      for(int v=0 ; v<dim_src_img[0]*dim_src_img[1]*dim_src_img[2] ; v++)
        p_src_img[ th ][v] = 0;
    }
  }
         
  // Counter for the number of events (and the number of trues, scatters and randoms for simulated data)
  for (int th=0 ; th < p_ID->GetNbThreadsForProjection() ; th++)
  {
    nLORs_trues[ th ]    = 0;
    nLORs_rdms[ th ]     = 0;
    nLORs_scatters[ th ] = 0;
    nLORs_unknown[ th ]  = 0;
    nBINs[ th ]          = 0;

    dt_max[ th ] = 0;
    dt_min[ th ] = 1000000000000.;
  }




  #ifdef CASTOR_ROOT
  
  uint32_t *castorID1 = new uint32_t[ p_ID->GetNbThreadsForProjection() ];
  uint32_t *castorID2 = new uint32_t[ p_ID->GetNbThreadsForProjection() ];
  
  // ROOT data variables
  int32_t *crystalID1 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t *crystalID2 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  
  // cylindricalPET specific
  int32_t *rsectorID1 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t *rsectorID2 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t *moduleID1  = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t *moduleID2  = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t *submoduleID1 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t *submoduleID2 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t *layerID1 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t *layerID2 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  
  // ecat specific
  int32_t *blockID1 = new int32_t[ p_ID->GetNbThreadsForProjection() ]; 
  int32_t *blockID2 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  
  // SPECThead specific
  float_t *rotAngle = new float_t[ p_ID->GetNbThreadsForProjection() ];
  float_t *gPosX    = new float_t[ p_ID->GetNbThreadsForProjection() ];
  float_t *gPosY    = new float_t[ p_ID->GetNbThreadsForProjection() ];
  float_t *gPosZ    = new float_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t *headID   = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t *pixelID  = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  
  
  // others
  int32_t  *eventID1 = new int32_t[ p_ID->GetNbThreadsForProjection() ]; 
  int32_t  *eventID2 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t  *comptonPhantomID1 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t  *comptonPhantomID2 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t  *rayleighPhantomID1 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t  *rayleighPhantomID2 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  int32_t  *sourceID1 = new int32_t[ p_ID->GetNbThreadsForProjection() ]; 
  int32_t  *sourceID2 = new int32_t[ p_ID->GetNbThreadsForProjection() ];
  float_t  *sourcePosX1 = new float_t[ p_ID->GetNbThreadsForProjection() ];
  float_t  *sourcePosY1 = new float_t[ p_ID->GetNbThreadsForProjection() ];
  float_t  *sourcePosZ1 = new float_t[ p_ID->GetNbThreadsForProjection() ];
            
           
  for (int th=0 ; th < 1 ; th++)
  {
    castorID1[ th ] = 0;
    castorID2[ th ] = 0;
    
    // ROOT data variables
    crystalID1[ th ] = 0;
    crystalID2[ th ] = 0;
    
    // cylindricalPET specific
    rsectorID1[ th ] = 0;
    rsectorID2[ th ] = 0;
    moduleID1[ th ] = 0;
    moduleID2[ th ] = 0;
    submoduleID1[ th ] = 0;
    submoduleID2[ th ] = 0;
    layerID1[ th ] = 0;
    layerID2[ th ] = 0;
    
    // ecat specific
    blockID1[ th ] = 0;
    blockID2[ th ] = 0;
    
    // SPECThead specific
    rotAngle[ th ] = 0.;
    gPosX[ th ] = 0.;
    gPosY[ th ] = 0.;
    gPosZ[ th ] = 0.;
    headID[ th ] = 0;
    pixelID[ th ] = 0;
    
    // others
    eventID1[ th ] = 0;
    eventID2[ th ] = 0;
    comptonPhantomID1[ th ] = 0;
    comptonPhantomID2[ th ] = 0;
    rayleighPhantomID1[ th ] = 0;
    rayleighPhantomID2[ th ] = 0;
    sourceID1[ th ] = 0;
    sourceID2[ th ] = 0;
    sourcePosX1[ th ] = 0.;
    sourcePosY1[ th ] = 0.;
    sourcePosZ1[ th ] = 0.;
  }
           
  // ============================================================================================================
  // ROOT objects declaration
  // ============================================================================================================

  TApplication *Tapp = new TApplication("tapp",0,0);

  TTree** Ttree    = new TTree*[path_to_input_file.size()];
  TFile** Tfile = new TFile*[path_to_input_file.size()];
  
  uint64_t *cum_sum_events = new uint64_t[path_to_input_file.size()];
          
  if(!input_is_intf)
  {
    for (size_t iFic=0 ; iFic<path_to_input_file.size() ; iFic++)
    { 
      Tfile[iFic] = new TFile(path_to_input_file[iFic].c_str(),"READ","ROOT file with histograms");
      if(GATE_system_type == GATE_SYS_SPECT)
        Ttree[iFic] = (TTree*)Tfile[iFic]->Get("Singles");
      else
        Ttree[iFic] = (TTree*)Tfile[iFic]->Get("Coincidences");
        
      cum_sum_events[iFic] = Ttree[iFic]->GetEntries();
      if (iFic > 0 ) cum_sum_events[iFic] += cum_sum_events[iFic-1];
      
      nLORs_tot += Ttree[iFic]->GetEntries();
    }
  }
  
  for (size_t iFic=0 ; iFic<path_to_input_file.size() ; iFic++)
    delete Tfile[iFic];
  
  delete[] Tfile;
  delete[] Ttree;
                                 
  #endif
  
  // Indexes for progression output
  uint64_t printing_index = 0;
  uint64_t printing_ratio = (nLORs_tot>10000) ? 
                                        10000 : 
                                     nLORs_tot;
  
  // ============================================================================================================
  // Recover system geometric elements values
  // ============================================================================================================

  // SPECThead system
  if(GATE_system_type == GATE_SYS_SPECT)
  {
    
    if( RecoverSPECTGeometry(         input_is_intf,
                              path_to_input_file[0],
                                   path_to_mac_file,
                                     distToDetector,
                                             nHeads,
                                       nPixelsAxial,
                                  nPixelsTransaxial,
                                     crystalSizeAxl,
                                     crystalSizeTrs,
                                    nProjectionsTot,
                                 nProjectionsByHead,
                                    head1stAngleDeg,
                                    headAngPitchDeg,
                                     headAngStepDeg,
                                   headRotDirection,
                                      start_time_ms,
                                        duration_ms,
                                       nCrystalsTot,
                                 p_proj_spect_image,
                                      spect_bin_axl,
                                      spect_bin_trs,
                                  nbSimulatedPixels,
                                     histo_out_flag,
                                             p_bins,
                                             p_scat,
                                          bins_elts,
                                     scat_corr_flag,
                  p_ID->GetNbThreadsForProjection(),
                                                 vb) )
      {
        Cerr("**** castor-GATERootToCastor :: Error when during SPECT geometry data recovery !" << endl);
        Exit(EXIT_FAILURE);
      }
  
    // Specific case for interfile (not supported for now)
    if(input_is_intf)
    {
      for (size_t p=0; p<bins_elts ; p++)
        for (size_t c=0; c<nCrystalsTot ; c++)
        {
          nLORs_tot          += p_proj_spect_image[p*nCrystalsTot+c];
          nLORs_unknown[ 0 ] += p_proj_spect_image[p*nCrystalsTot+c];
        }
    }
    
  }
  
  else // PET systems
  {
    if ( RecoverPETGeometry( GATE_system_type,
                            path_to_mac_file,
                                     nLayers,
                        nb_crystal_per_layer,
                                nCrystalsTot,
                              nCrystalsAxial, 
                         nCrystalsTransaxial, 
                             nLayersRptAxial,
                        nLayersRptTransaxial,
                            nSubmodulesAxial, 
                       nSubmodulesTransaxial, 
                               nModulesAxial, 
                          nModulesTransaxial, 
                              nRsectorsAxial,
                             nRsectorsAngPos,
                                 nBlocksLine, 
                              nBlocksPerRing,
                            invert_det_order,
                            rsector_id_order,
                            pet_coinc_window,
                               start_time_ms,
                                 duration_ms,
                              histo_out_flag,
                                      p_bins,
                                      p_scat,
                                       p_rdm,
                                   bins_elts,
                              scat_corr_flag,
                               rdm_corr_flag,
           p_ID->GetNbThreadsForProjection(),
                                          vb) )
      {
        Cerr("**** castor-GATERootToCastor :: Error when reading PET mac file : "<< path_to_mac_file <<" !" << endl);
        Exit(EXIT_FAILURE);
      } 

  }  // end of PET section
  


  // ============================================================================================================
  // Histogram computation step
  // ============================================================================================================
  if(!input_is_intf) // ROOT data processing (the only exception is currently SPECT projection interfile for which data is already processed)
  {
    #ifdef CASTOR_ROOT
    
    // Histogram building : only if output is histogram AND/OR if scatter or random correction rates are estimated
    if( histo_out_flag 
    || (scat_corr_flag || rdm_corr_flag) )
    {
      // Clock start for histogram pre_computation
      clock_t clock_start_histo_pre_comp = clock();
      time_t time_start_histo_pre_comp = time(NULL);
    
      Cout(endl << "Step 1/"<< 2 <<": Computing scatter and/or random histogram..." << endl);
    
      //---------------------------------------------- /ROOT data variables
          
      // Loop on the GEvents in the current datafile
      
      // mthread loop: no need to keep track of the chronological
      //               order of the event here are we are just filling
      //               the histogram 
      
      int th=0;
      //#pragma omp parallel for private(i) schedule(static, chunk)//, default(none)
      int nb_threads = p_ID->GetNbThreadsForProjection();
      //nb_threads = 1;
      #pragma omp parallel for private(th) schedule(static, 1), num_threads(nb_threads)//, default(none)
      for (th=0; th<nb_threads ; th++)
      {
        TTree** GEvents    = new TTree*[path_to_input_file.size()];
        TFile** Tfile_root = new TFile*[path_to_input_file.size()];
        
        Troot_pet   p_Rpet;
        Troot_spect p_Rspect;
        
        // Compute the total number of LORs in the dataset
        #pragma omp critical(root_file_init)
        {
          for (size_t iFic=0 ; iFic<path_to_input_file.size() ; iFic++)
          { 
            Tfile_root[iFic] = new TFile(path_to_input_file[iFic].c_str(),"READ","ROOT file with histograms");
            
            if(GATE_system_type == GATE_SYS_SPECT)
            {
              GEvents[iFic] = (TTree*)Tfile_root[iFic]->Get("Singles");
              InitTrootSpect( &p_Rspect );
              SetROOTSPECTBranch(GEvents[iFic], &p_Rspect);
            }
            else
            {
              GEvents[iFic] = (TTree*)Tfile_root[iFic]->Get("Coincidences");
              InitTrootPet( &p_Rpet );
              p_Rpet.bool_ecat_mode = (GATE_system_type == GATE_SYS_ECAT) ? 1 : 0;
              SetROOTPETBranch(GEvents[iFic], &p_Rpet);
            }
          }
        }

        uint64_t i = 0;
        uint64_t chunk = nLORs_tot/nb_threads;
        uint64_t start = th*chunk;

        uint64_t stop = (th+1)==nb_threads ? nLORs_tot : start + chunk;
        
        for (uint64_t i=start ; i<stop ; i++)
        {
          uint16_t th_cur_fic = 0;
          
          // Get file number and entry corresponding to this event
          while( i >= cum_sum_events[ th_cur_fic ] )
            th_cur_fic++;
     
          uint64_t entry = (th_cur_fic == 0)                  ?
                           i                                  :
                           i - cum_sum_events[ th_cur_fic-1 ] ;
              
          uint8_t kind = KIND_UNKNOWN ;

          GEvents[th_cur_fic]->GetEntry(entry);
          
          if (th==0) printing_index++;
          
          // --- ID Conversions PET systems ---
          if (GATE_system_type == GATE_SYS_CYLINDRICAL
           || GATE_system_type == GATE_SYS_ECAT )
          {              
            if (ProcessPETEvent( GATE_system_type, &p_Rpet, castorID1[ th ], castorID2[ th ],
                                 nLayers, nb_crystal_per_layer,
                                 nCrystalsTot, nCrystalsAxial, nCrystalsTransaxial, nLayersRptAxial, nLayersRptTransaxial,
                                 nSubmodulesAxial,  nSubmodulesTransaxial, nModulesAxial, nModulesTransaxial, nRsectorsAxial, nRsectorsAngPos, nBlocksLine, nBlocksPerRing, invert_det_order, rsector_id_order, kind,
                                 nLORs_unknown[ th ], nLORs_trues[ th ], nLORs_scatters[ th ], nLORs_rdms[ th ], 
                                 vb ) )
            {
              Cerr("**** castor-GATERootToCastor :: Error when processing PET event !" << endl);
              Exit(EXIT_FAILURE);
            }
            
             // --- Check if some restrictions based on event type are enabled, skip next steps in this case ---
             if (kind == KIND_TRUE
             && (scat_only_flag || rdm_only_flag) )
              continue;
    
             if ((kind == KIND_SCAT || kind == KIND_MSCAT)
             && (true_only_flag || rdm_only_flag || true_rdm_flag) )
              continue;
              
             if ((kind == KIND_RDM)
             && (true_only_flag || scat_only_flag || true_scat_flag) )
              continue;
              
              
            // Fill Histograms
            uint64_t id1 = (castorID1[ th ] < castorID2[ th ]) ? castorID1[ th ] : castorID2[ th ];
            uint64_t id2 = (castorID1[ th ] < castorID2[ th ]) ? castorID2[ th ] : castorID1[ th ];
            
            if(histo_out_flag)
              p_bins[ th ][id1][id2-id1-1]++;
            
            if( scat_corr_flag
            && (kind == KIND_SCAT || kind == KIND_MSCAT) )
              p_scat[ th ][id1][id2-id1-1]++;

            if( rdm_corr_flag
            &&  kind == KIND_RDM )
              p_rdm[ th ][id1][id2-id1-1]++;
              
            
            // Source image (if no rdm event) for histogram output
            if(histo_out_flag && src_img_flag
            && p_Rpet.eventID1 == p_Rpet.eventID2
            && p_Rpet.eventID1 >0 ) // No noise event (==-2)
              {
                INTNB x,y,z;
                x = (int)(( p_Rpet.sourcePosX1 + dim_src_img[0]/2*vox_src_img[0]) / vox_src_img[0]);
                // CASToR Y axis inverted in comparison with GATE (left-handed coordinate)
                y = (int)(( p_Rpet.sourcePosY1 + dim_src_img[1]/2*vox_src_img[1]) / vox_src_img[1]);
                z = (int)(( p_Rpet.sourcePosZ1 + dim_src_img[2]/2*vox_src_img[2]) / vox_src_img[2]);
    
                if(x >= 0 && x < dim_src_img[0] &&
                   y >= 0 && y < dim_src_img[1] &&
                   z >= 0 && z < dim_src_img[2] )
                  p_src_img[ th ][z*dim_src_img[1]*dim_src_img[0] + y*dim_src_img[0] + x]++;
              }
              
          }
          
        
          // --- ID Conversions SPECT systems ---
          
          else if (GATE_system_type == GATE_SYS_SPECT)
          {
            // Compute castor ids
            if (ProcessSPECTEvent(  &p_Rspect, castorID1[ th ], castorID2[ th ],
                                    headAngPitchDeg, headAngStepDeg, nbSimulatedPixels, nPixelsTransaxial, nPixelsAxial, crystalSizeAxl, crystalSizeTrs, nProjectionsByHead, kind,
                                    nLORs_unknown[ th ], nLORs_trues[ th ], nLORs_scatters[ th ], nLORs_rdms[ th ], 
                                    vb ) )
            {
              Cerr("**** castor-GATERootToCastor :: Error when processing SPECT event !" << endl);
              Exit(EXIT_FAILURE);
            }
            
            // --- Check if some restrictions based on event type are enabled, skip next steps in this case ---

            if (kind == KIND_TRUE
            && (scat_only_flag || rdm_only_flag) )
              continue;
            
            if ((kind == KIND_SCAT || kind == KIND_MSCAT)
            && (true_only_flag || rdm_only_flag || true_rdm_flag) )
              continue;
            
            if (kind == KIND_RDM
            && (true_only_flag || scat_only_flag || true_scat_flag) )
              continue;

            // Fill Histograms
            
            if( histo_out_flag )
              p_bins[ th ][castorID1[ th ]][castorID2[ th ]]++;

            if( scat_corr_flag 
            && (kind == KIND_SCAT || kind == KIND_MSCAT) )
              p_scat[ th ][castorID1[ th ]][castorID2[ th ]]++;


            // Source image (if no rdm event) for histogram output
            if(histo_out_flag && src_img_flag
            && p_Rspect.eventID >0 ) // No noise event (==-2))
              {
                INTNB x,y,z;
                x = (int)(( p_Rspect.sourcePosX + dim_src_img[0]/2*vox_src_img[0]) / vox_src_img[0]);
                // CASToR Y axis inverted in comparison with GATE (left-handed coordinate)
                y = (int)(( p_Rspect.sourcePosY + dim_src_img[1]/2*vox_src_img[1]) / vox_src_img[1]);
                z = (int)(( p_Rspect.sourcePosZ + dim_src_img[2]/2*vox_src_img[2]) / vox_src_img[2]);
    
                if(x >= 0 && x < dim_src_img[0] &&
                   y >= 0 && y < dim_src_img[1] &&
                   z >= 0 && z < dim_src_img[2] )
                  p_src_img[ th ][z*dim_src_img[1]*dim_src_img[0] + y*dim_src_img[0] + x]++;
              }
          }
          
          // Progression
          if (th == 0 && printing_index %( nLORs_tot/printing_ratio ) == 0)
          {
            FLTNB percent = ( ((FLTNB)(printing_index+1))/((FLTNB)nLORs_tot/(FLTNB)nb_threads) ) * ((FLTNB)100);
            cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
              << percent << " %                    ";
          }
        
        } // end of loop on data events

        // Deallocation
        for (size_t iFic=0 ; iFic<path_to_input_file.size() ; iFic++)
        {
          if(Tfile_root[ iFic ]) delete Tfile_root[ iFic ];
        }
        
        if(Tfile_root) delete[] Tfile_root;
        if(GEvents)   delete[] GEvents;
        
      } // end of multithreaded loop on threads
      
      // Clock stop histogram computation
      clock_t clock_stop_histo_pre_comp = clock();
      time_t time_stop_histo_pre_comp = time(NULL);
  
      if (vb>=2) Cout ("Histogram pre_computation -> Time spent: " << time_stop_histo_pre_comp-time_start_histo_pre_comp
                         << " sec | CPU: " << (clock_stop_histo_pre_comp-clock_start_histo_pre_comp)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);
        
      printing_index = 0;
    
      // Mthread: recover histogram tables
      if(vb>=3) Cout(endl << "Reduce ... " << endl);
          
      // Reduce
      for (int th=1 ; th < p_ID->GetNbThreadsForProjection() ; th++)
      {   
        if(GATE_system_type == GATE_SYS_SPECT) 
        {
          for (size_t id1=0 ; id1<bins_elts ; id1++)
            for (size_t id2=0; id2<nCrystalsTot ; id2++)
            {
              p_bins[0][id1][id2] += p_bins[th][id1][id2];
              
              if(scat_corr_flag && !input_is_intf)
                p_scat[0][id1][id2] += p_scat[th][id1][id2];
                
            }
        }
        else // PET
        {
          for (size_t id1=0 ; id1<bins_elts ; id1++)
            for (size_t id2 = id1+1; id2<nCrystalsTot ; id2++)
            {
              p_bins[0][id1][id2-id1-1] += p_bins[th][id1][id2-id1-1];
              
              if(scat_corr_flag)
                p_scat[0][id1][id2-id1-1] += p_scat[th][id1][id2-id1-1];
                
              if(rdm_corr_flag)
                p_rdm[0][id1][id2-id1-1]  += p_rdm[th][id1][id2-id1-1];
            }
        }
      }

      if(vb>=3) Cout(endl << "Reduce ok " << endl);
      
    } // end of Histogram building step
    
    
    // ============================================================================================================
    // LIST MODE processing step ...
    // ============================================================================================================
  
    if(!histo_out_flag) // Not histogram output
    {
      // Clock start for list-mode computation
      clock_t clock_start_list_comp = clock();
      time_t time_start_list_comp = time(NULL);
      
      if(vb>=2) Cout(endl << "Step "<< 1+(scat_corr_flag || rdm_corr_flag)
                          <<"/"<< 1 +  (scat_corr_flag || rdm_corr_flag) << endl);
               
      int th=0;
      
      //#pragma omp parallel for private(i) schedule(static, chunk)//, default(none)      
      #pragma omp parallel for private(th) schedule(static, 1), num_threads(p_ID->GetNbThreadsForProjection())
      for (th=0; th<p_ID->GetNbThreadsForProjection() ; th++)
      {
        TTree** GEvents    = new TTree*[path_to_input_file.size()];
        TFile** Tfile_root = new TFile*[path_to_input_file.size()];
        
        Troot_pet   p_Rpet;
        Troot_spect p_Rspect;
        
        // Compute the total number of LORs in the dataset
        #pragma omp critical(root_file_init)
        {
          for (size_t iFic=0 ; iFic<path_to_input_file.size() ; iFic++)
          { 
            Tfile_root[iFic] = new TFile(path_to_input_file[iFic].c_str(),"READ","ROOT file with histograms");
            
            if(GATE_system_type == GATE_SYS_SPECT)
            {
              GEvents[iFic] = (TTree*)Tfile_root[iFic]->Get("Singles");
              InitTrootSpect( &p_Rspect );
              SetROOTSPECTBranch(GEvents[iFic], &p_Rspect);
            }
            else
            {
              GEvents[iFic] = (TTree*)Tfile_root[iFic]->Get("Coincidences");
              InitTrootPet( &p_Rpet );
              p_Rpet.bool_ecat_mode = (GATE_system_type == GATE_SYS_ECAT) ? 1 : 0;
              SetROOTPETBranch(GEvents[iFic], &p_Rpet);
            }
          }
        }

        uint64_t i = 0;
        uint64_t chunk = nLORs_tot/p_ID->GetNbThreadsForProjection() + 1;
        uint64_t start = th*chunk;
        uint64_t stop  = min(start + chunk, nLORs_tot);
          
        for (uint64_t i=start ; i<stop ; i++)
        {        
          uint16_t th_cur_fic = 0;
          
          // Get file number and entry corresponding to this event
          while( i >= cum_sum_events[ th_cur_fic ] )
            th_cur_fic++;
     
          uint64_t entry = (th_cur_fic == 0)                ?
                           i                                      :
                           i - cum_sum_events[ th_cur_fic-1 ] ;
              
          uint8_t kind = KIND_UNKNOWN ;
              
          GEvents[th_cur_fic]->GetEntry( entry );
          if (th==0) printing_index++;
  
          if(vb >= 4 && th==0)
            Cout("thread " << th << "File#" << th_cur_fic << ", event#" << i << endl;);
         
         
          // --- ID Conversions PET ---
          
          if (GATE_system_type == GATE_SYS_CYLINDRICAL
           || GATE_system_type == GATE_SYS_ECAT )
          {
            //Troot_pet p_Rt;
            //InitTrootPet( &p_Rt );
            
            if (ProcessPETEvent( GATE_system_type, &p_Rpet, castorID1[ th ], castorID2[ th ],
                                 nLayers, nb_crystal_per_layer, nCrystalsTot, nCrystalsAxial, nCrystalsTransaxial, nLayersRptAxial, nLayersRptTransaxial,
                                 nSubmodulesAxial,  nSubmodulesTransaxial, nModulesAxial, nModulesTransaxial, nRsectorsAxial, nRsectorsAngPos, nBlocksLine, nBlocksPerRing, invert_det_order, rsector_id_order, kind,
                                 nLORs_unknown[ th ], nLORs_trues[ th ], nLORs_scatters[ th ], nLORs_rdms[ th ], 
                                 vb ) )
            {
              Cerr("**** castor-GATERootToCastor :: Error when processing PET event !" << endl);
              Exit(EXIT_FAILURE);
            }
                        
            // --- Check if some restrictions based on event type are enabled, skip next steps in this case ---
            if (kind == KIND_TRUE
            && (scat_only_flag || rdm_only_flag) )
              continue;
    
            if ((kind == KIND_SCAT || kind == KIND_MSCAT)
            && (true_only_flag || rdm_only_flag || true_rdm_flag) )
              continue;
              
            if ((kind == KIND_RDM)
            && (true_only_flag || scat_only_flag || true_scat_flag) )
              continue;
              
            // --- Write Event ---
            iEventListPET* Event = new iEventListPET();
            
            // Write event in the datafile
            uint32_t time_in_ms = p_Rpet.time1*1000;
            int nb_lines_in_event = 1; // 1 by default for GATE root files
              
            Event->SetNbLines(nb_lines_in_event);
            Event->AllocateID();
            Event->SetID1(0, castorID1[ th ]); // 1st argument : line, 2nd argument :index
            Event->SetID2(0, castorID2[ th ]); // 1st argument : line, 2nd argument :index
            Event->SetTimeInMs(time_in_ms);

            // Estimate acf for the bin
            if(estimate_acf_flag)
            {  
              // Compute the system matrix elements for the two crystals 
              oProjectionLine* line = p_ProjectorManager->ComputeProjectionLine(Event, th);
      
              // Compute forward projection of attenuation image
              FLTNB fp = 0.;
              if (line->NotEmptyLine())
                fp = line->ForwardProject(p_ImageSpace->m4p_image[0][0][0]) * 0.1; // 0.1 -> convert in mm-1
              
              // Write atn correction factor in Event
              Event->SetAttenuationCorrectionFactor(1/std::exp(-fp));
            }
          
            // scatter correction rate
            if(scat_corr_flag)
            {
              uint64_t nb_scat = 0;
              
              uint64_t id1 = (castorID1[ th ] < castorID2[ th ]) ? castorID1[ th ] : castorID2[ th ];
              uint64_t id2 = (castorID1[ th ] < castorID2[ th ]) ? castorID2[ th ] : castorID1[ th ];
              
              nb_scat = p_scat[ 0 ][id1][id2-id1-1];

              Event->SetScatterRate(0, (FLTNB)nb_scat*1000. / duration_ms );
            }
            
            // random correction rate
            if(rdm_corr_flag)
            {
              FLTNB nb_rdm = 0.;
              
              uint64_t id1 = (castorID1[ th ] < castorID2[ th ]) ? castorID1[ th ] : castorID2[ th ];
              uint64_t id2 = (castorID1[ th ] < castorID2[ th ]) ? castorID2[ th ] : castorID1[ th ];
              
              nb_rdm = p_rdm[ 0 ][id1][id2-id1-1];
                
              Event->SetRandomRate( (FLTNB)nb_rdm*1000 / duration_ms );
            }
            
            // Compute ToF by default, only written if enabled in datafile
            double_t dt_in_ps = (p_Rpet.time1 - p_Rpet.time2)*1.0e+12 ;
            Event->SetTOFMeasurementInPs(dt_in_ps);
            dt_max[ th ] = (dt_max[ th ]<fabs(dt_in_ps) ) ? fabs(dt_in_ps) : dt_max[ th ] ; 
            dt_min[ th ] = (dt_min[ th ]>fabs(dt_in_ps) ) ? fabs(dt_in_ps) : dt_min[ th ] ;

            Event->SetKind(kind);
            Out_data_file->WriteEvent(Event, th);
            delete Event;
            
            // Source image (if no rdm event)
            if(src_img_flag  
              && p_Rpet.eventID1 >0 && p_Rpet.eventID2 >0) // No noise event (==-2)
              {
                INTNB x,y,z;
                x = (int)(( p_Rpet.sourcePosX1 + dim_src_img[0]/2*vox_src_img[0]) / vox_src_img[0]);
                // CASToR Y axis inverted in comparison with GATE (left-handed coordinate)
                y = (int)(( p_Rpet.sourcePosY1 + dim_src_img[1]/2*vox_src_img[1]) / vox_src_img[1]);
                z = (int)(( p_Rpet.sourcePosZ1 + dim_src_img[2]/2*vox_src_img[2]) / vox_src_img[2]);
    
                if(x >= 0 && x < dim_src_img[0] &&
                   y >= 0 && y < dim_src_img[1] &&
                   z >= 0 && z < dim_src_img[2] )
                  p_src_img[ th ][z*dim_src_img[1]*dim_src_img[0] + y*dim_src_img[0] + x]++;
              }
          }
          
        
          // --- ID Conversions SPECT ---
          
          else if (GATE_system_type == GATE_SYS_SPECT)
          {
            //Troot_spect p_Rt;
            //InitTrootSpect( &p_Rt );
            
            // Compute castor ids
            if (ProcessSPECTEvent(  &p_Rspect, castorID1[ th ], castorID2[ th ],
                                    headAngPitchDeg, headAngStepDeg, nbSimulatedPixels, nPixelsTransaxial, nPixelsAxial, crystalSizeAxl, crystalSizeTrs, nProjectionsByHead, kind,
                                    nLORs_unknown[ th ], nLORs_trues[ th ], nLORs_scatters[ th ], nLORs_rdms[ th ], 
                                    vb ) )
            {
              Cerr("**** castor-GATERootToCastor :: Error when processing SPECT event !" << endl);
              Exit(EXIT_FAILURE);
            }
            
            // --- Check if some restrictions based on event type are enabled, skip next steps in this case ---
            if (kind == KIND_TRUE
            && (scat_only_flag || rdm_only_flag) )
              continue;
            
            if ((kind == KIND_SCAT || kind == KIND_MSCAT)
            && (true_only_flag || rdm_only_flag || true_rdm_flag) )
              continue;
            
            if ((kind == KIND_RDM)
            && (true_only_flag || scat_only_flag || true_scat_flag) )
              continue;

          
            // --- Write Event ---
            iEventListSPECT* Event = new iEventListSPECT();
            
            // Write event in the datafile
            uint32_t time_in_ms = p_Rspect.time*1000;
            
            Event->SetNbLines(1);
            Event->AllocateID();
            Event->SetID1(0, castorID1[ th ]); // 1st argument : line, 2nd argument :index
            Event->SetID2(0, castorID2[ th ]); // 1st argument : line, 2nd argument :index
            Event->SetTimeInMs(time_in_ms);
          
            Event->SetKind(kind);
             
            // scatter correction rate
            if(scat_corr_flag)
            {
              uint64_t nb_scat = p_scat[ 0 ][castorID1[ th ]][castorID2[ th ]];                
              Event->SetScatterRate((FLTNB) nb_scat*1000 / duration_ms);
            }
             
            Out_data_file->WriteEvent(Event, th);
            delete Event;
            
            // Source image (if no rdm event)
            if(src_img_flag  
              && p_Rspect.eventID >0 ) // No noise event (==-2)
              {
                INTNB x,y,z;
                x = (int)(( p_Rspect.sourcePosX + dim_src_img[0]/2*vox_src_img[0]) / vox_src_img[0]);
                // CASToR Y axis inverted in comparison with GATE (left-handed coordinate)
                y = (int)(( p_Rspect.sourcePosY + dim_src_img[1]/2*vox_src_img[1]) / vox_src_img[1]);
                z = (int)(( p_Rspect.sourcePosZ + dim_src_img[2]/2*vox_src_img[2]) / vox_src_img[2]);
    
                if(x >= 0 && x < dim_src_img[0] &&
                   y >= 0 && y < dim_src_img[1] &&
                   z >= 0 && z < dim_src_img[2] )
                  p_src_img[ th ][z*dim_src_img[1]*dim_src_img[0] + y*dim_src_img[0] + x]++;
              }
          }
          
    

      
        
          // Progression
          if (th == 0 && printing_index %( nLORs_tot/printing_ratio ) == 0)
          {
            FLTNB percent = ( ((FLTNB)(printing_index+1))/((FLTNB)nLORs_tot/(FLTNB)p_ID->GetNbThreadsForProjection()) ) * ((FLTNB)100);
            cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
              << percent << " %                    ";
          }
          
        } // end of loop on events

      } // end of multithreaded loop
      
      // Clock stop list-mode computation
      clock_t clock_stop_list_comp = clock();
      time_t time_stop_list_comp = time(NULL);
        
      if (vb>=2) Cout ("List-mode computation -> Time spent: " << time_stop_list_comp-time_start_list_comp
                           << " sec | CPU: " << (clock_stop_list_comp-clock_start_list_comp)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);
                           
    } // end of list mode processing step
    
    #endif
  
    // Reduce
    for (int th=1 ; th < p_ID->GetNbThreadsForProjection() ; th++)
    {
      nLORs_trues[ 0 ] += nLORs_trues[ th ];
      nLORs_rdms[ 0 ] += nLORs_rdms[ th ];
      nLORs_scatters[ 0 ] += nLORs_scatters[ th ];
      nLORs_unknown[ 0 ] += nLORs_unknown[ th ];
      
      
      dt_max[ 0 ] += dt_max[ th ];
      dt_min[ 0 ] += dt_min[ th ];
      
      if(src_img_flag)
        for(int v=0 ; v<dim_src_img[0]*dim_src_img[1]*dim_src_img[2] ; v++)
          p_src_img[ 0 ][ v ] += p_src_img[ th ][ v ];
      
    }

  } // end of root data processing
  
  // Free Root objects
  #ifdef CASTOR_ROOT
  delete[] castorID1;
  delete[] castorID2;
  
  // ROOT data variables
  delete[] crystalID1;
  delete[] crystalID2;
  
  // cylindricalPET specific
  delete[] rsectorID1;
  delete[] rsectorID2;
  delete[] moduleID1;
  delete[] moduleID2;
  delete[] submoduleID1;
  delete[] submoduleID2;
  delete[] layerID1;
  delete[] layerID2;
  
  // ecat specific
  delete[] blockID1; 
  delete[] blockID2;
  
  // SPECThead specific
  delete[] rotAngle;
  delete[] gPosX;
  delete[] gPosY;
  delete[] gPosZ;
  delete[] headID;
  delete[] pixelID;
  
  
  // others
  delete[] eventID1; 
  delete[] eventID2;
  delete[] comptonPhantomID1;
  delete[] comptonPhantomID2;
  delete[] rayleighPhantomID1;
  delete[] rayleighPhantomID2;
  delete[] sourceID1; 
  delete[] sourceID2;
  delete[] sourcePosX1;
  delete[] sourcePosY1;
  delete[] sourcePosZ1;
           
  delete[] cum_sum_events;
  
  delete   Tapp;
  #endif

 
  // ============================================================================================================
  // Datafile and header writing step
  // ============================================================================================================
  
  // Clock start datafile generation
  clock_t clock_start_file_gen = clock();
  time_t time_start_file_gen = time(NULL);
  
  if(GATE_system_type == GATE_SYS_SPECT) 
  {
    if (histo_out_flag)
    {
      printing_index=0;
      
      if(vb>=2) Cout(endl << endl << "Step "<< 2
                     <<"/"<< 2
                     <<": Generate the histogram datafile" << endl);


      uint64_t nb_bins = bins_elts*nCrystalsTot;
      printing_ratio = (nb_bins>1000) ? 1000 : nb_bins/10;
      
      // Loop on the crystal IDs
      for (size_t id1=0 ; id1<bins_elts ; id1++)
      {
        size_t id2=0;
        #pragma omp parallel for private(id2) schedule(static, 1)
        for (id2=0; id2<nCrystalsTot ; id2++)
        {
          // Get the thread index
          int th = 0;
          #ifdef CASTOR_OMP
          th = omp_get_thread_num();
          #endif
          
          uint32_t nb_events_in_bin = p_bins[th][id1][id2];
          iEventHistoSPECT* Event = new iEventHistoSPECT();
          int nb_lines = 1; // 1 by default for GATE root file;
          Event->SetNbLines(nb_lines);
          Event->AllocateID();
          Event->SetID1(0, id1);  // 1st argument : line, 2nd argument :index
          Event->SetID2(0, id2);  // 1st argument : line, 2nd argument :index
          Event->SetEventValue(0, nb_events_in_bin);                  
                  
          // scatter correction rate
          if(scat_corr_flag)
            Event->SetScatterRate( (FLTNB)p_scat[th][id1][id2]*1000 / duration_ms );
          
          // Write event in DataFile
          Out_data_file->WriteEvent(Event, th);
          delete Event;
          nBINs[ th ]++;
            
          // Progression
          if (printing_index%((nb_bins)/printing_ratio) == 0)
          {
            FLTNB percent = ( ((FLTNB)(printing_index+1))/((FLTNB)nb_bins) ) * ((FLTNB)100);
            cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
              << percent << " %                    ";
          }
          
          printing_index++;
        }
      }
      // Reduce nbins
      for (int th=1 ; th < p_ID->GetNbThreadsForProjection() ; th++)
        nBINs[ 0 ] += nBINs[ th ];
      
      // Set metadata info and write header
      Out_data_file->SetStartTime((FLTNB)(start_time_ms)/1000); // start time in s
      Out_data_file->SetDuration((FLTNB)(duration_ms)/1000); // duration in s
      Out_data_file->SetNbEvents( nBINs[ 0 ] );

      Cout(endl << "The output histogram contains " << nBINs[ 0 ]  << " events." << endl);    
    }
    // Write the List-mode datafile header 
    else
    {
      Out_data_file->SetStartTime((FLTNB)(start_time_ms)/1000); // start time in s
      Out_data_file->SetDuration((FLTNB)(duration_ms)/1000); // duration in s

      uint64_t nLORs =  (true_only_flag) ? nLORs_trues[ 0 ] : 
                        (scat_only_flag) ? nLORs_scatters[ 0 ] :
                        (true_scat_flag) ? nLORs_trues[ 0 ] + nLORs_scatters[ 0 ] :
                         nLORs_tot;

      Out_data_file->SetNbEvents( nLORs );
    }
    
    FLTNB* p_angles = new FLTNB[nProjectionsTot];
    FLTNB* p_distToDetector = new FLTNB[nProjectionsTot];
    
    for(size_t h=0 ; h<nHeads ; h++)
      for(size_t p=0 ; p<nProjectionsByHead ; p++)
      {
        int idx_proj = h*nProjectionsByHead+p;
        p_angles[idx_proj] = head1stAngleDeg // initial angular position of the head
                           + p*headAngStepDeg // projection angular position
                           + h*headAngPitchDeg; // angular shift for each head
        
        // same distance for each head with GATE
        p_distToDetector[idx_proj] = distToDetector; 
      }

    ((iDataFileSPECT*)Out_data_file)->SetNbBins(nPixelsTransaxial, nPixelsAxial);
    ((iDataFileSPECT*)Out_data_file)->SetNbProjections(nProjectionsTot);
    ((iDataFileSPECT*)Out_data_file)->SetNbHeads(nHeads);
    
    if( ((iDataFileSPECT*)Out_data_file)->InitAngles(p_angles) )
    {
      Cerr("**** castor-GATERootToCastor :: Error when trying to set projection angles values !" << endl);
      Exit(EXIT_FAILURE);
    }
    
    if( ((iDataFileSPECT*)Out_data_file)->InitCorToDetectorDistance(p_distToDetector) ) 
    {
      Cerr("**** castor-GATERootToCastor :: Error when trying to set distance between center of rotation and detectors !" << endl);
      Exit(EXIT_FAILURE);
    }
    
    ((iDataFileSPECT*)Out_data_file)->SetHeadRotDirection(headRotDirection);

    delete[] p_angles;
    delete[] p_distToDetector;
  
    Out_data_file->WriteHeader();
  } // end of SPECT section
  
  else // PET 
  {
    if (histo_out_flag) // Histogram
    {
      printing_index=0;
      
      // Writing the datafile      
      if(vb>=2) Cout(endl << endl 
                     << "Step "<< 2
                     << "/"<< 2 
                     << ": Generate the histogram datafile" << endl);
                          
      uint64_t nb_bins = (nCrystalsTot*nCrystalsTot - nCrystalsTot)/2;

      printing_ratio = (nb_bins>1000) ? 1000 : nb_bins/10;
      
      // Loop on the crystal IDs
      size_t id1=0;
      //#pragma omp parallel for private(id1) schedule(static, (bins_elts/p_ID->GetNbThreadsForProjection())+1)   
      #pragma omp parallel for private(id1) schedule(static, 1)
      for (id1=0 ; id1 <bins_elts ; id1++)
        for (size_t id2 = id1+1; id2 < nCrystalsTot;id2++)
        {
          // Get the thread index
          int th = 0;
          #ifdef CASTOR_OMP
          th = omp_get_thread_num();
          #endif
          
          uint32_t nb_events_in_bin = p_bins[0][id1][id2-id1-1];
          
          iEventHistoPET* Event = new iEventHistoPET();
          int nb_lines = 1; // 1 by default for GATE root file;
          Event->SetNbLines(nb_lines);
          Event->AllocateID();
          Event->SetID1(0, id1);  // 1st argument : line, 2nd argument :index
          Event->SetID2(0, id2);  // 1st argument : line, 2nd argument :index
          Event->SetEventValue(0, nb_events_in_bin);
          Event->SetTimeInMs(start_time_ms); // Set to acquisition start time
          //cout << th << " " << id1 << " " << id2 <<  endl;
          // Estimate acf for the bin
          if(estimate_acf_flag)
          {  
            // Compute the system matrix elements for the two crystals 
            oProjectionLine* line = p_ProjectorManager->ComputeProjectionLine(Event, th);
            
            // Compute forward projection of attenuation image
            FLTNB fp = 0.;
            if (line->NotEmptyLine())
              fp = line->ForwardProject(p_ImageSpace->m4p_image[0][0][0]) * 0.1; // 0.1 -> convert in mm-1
            
            // Write atn correction factor in Event
            Event->SetAttenuationCorrectionFactor(1/std::exp(-fp));
          }
          
          // scatter correction rate
          if(scat_corr_flag)
          {
            FLTNB nb_scat = 0.;
            if (id1 < id2)
              nb_scat = p_scat[th][id1][id2-id1-1];
            else
              nb_scat = p_scat[th][id2][id1-id2-1];

            Event->SetScatterRate(0, nb_scat*1000 / duration_ms);
          }
          
          // random correction rate
          if(rdm_corr_flag)
          {
            FLTNB nb_rdm = 0.;
            if (id1 < id2)
              nb_rdm = p_rdm[th][id1][id2-id1-1];
            else
              nb_rdm = p_rdm[th][id2][id1-id2-1];
              
            Event->SetRandomRate(nb_rdm*1000 / duration_ms);
          }
          
          
          // Write event in DataFile
          Out_data_file->WriteEvent(Event, th);
          
          delete Event;
          nBINs[ th ]++;
            
          // Progression
          if (th ==0 && printing_index%((nb_bins)/printing_ratio) == 0)
          {
            FLTNB percent = ( ((FLTNB)(printing_index+1))/((FLTNB)nb_bins) ) * ((FLTNB)100);
            cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
              << percent << " %                    ";
          }
          printing_index++;
        }
  
      // Reduce nbins
      for (int th=1 ; th < p_ID->GetNbThreadsForProjection() ; th++)
        nBINs[ 0 ] += nBINs[ th ];
        
      // Set metadata info and write header
      Out_data_file->SetStartTime((FLTNB)(start_time_ms)/1000); // start time in s
      Out_data_file->SetDuration((FLTNB)(duration_ms)/1000); // duration in s
      Out_data_file->SetNbEvents( nBINs[ 0 ] );
      Out_data_file->WriteHeader(); 
      
      Cout(endl << "The output histogram contains " << nBINs[ 0 ]  << " events." << endl);    
    }
    // Write the List-mode datafile header 
    else
    {
      // If TOF enabled, compute TOF range from dt_max & dt_min recovered from the data
      // Set by user
      if (tof_range>0)
      {
        ((iDataFilePET*)Out_data_file)->SetTOFMeasurementRangeInPs(tof_range);
        if (vb>=3) Cout("Max range of TOF measurements set by user to: " << tof_range << " ps " << endl);
      }
      else if (pet_coinc_window>0)
      {
        ((iDataFilePET*)Out_data_file)->SetTOFMeasurementRangeInPs(pet_coinc_window);
        if (vb>=3) Cout("Max range of TOF measurements set to coincidence window: " << pet_coinc_window << " ps " << endl);
      }
      else
      {
        ((iDataFilePET*)Out_data_file)->SetTOFMeasurementRangeInPs(dt_max[ 0 ] - dt_min[ 0 ]);
        if (vb>=3) Cout("Max range of TOF measurements estimated from data: " << dt_max - dt_min << " ps " << endl);
      }
      
        
      Out_data_file->SetStartTime((FLTNB)(start_time_ms)/1000); // start time in s
      Out_data_file->SetDuration((FLTNB)(duration_ms)/1000); // duration in s

      Out_data_file->SetNbEvents( (true_only_flag) ? 
                                  nLORs_trues[ 0 ] : 
                                        nLORs_tot );
      uint64_t nLORs =  (true_only_flag) ? nLORs_trues[ 0 ] : 
                        (scat_only_flag) ? nLORs_scatters[ 0 ] :
                        (rdm_only_flag)  ? nLORs_rdms[ 0 ] :
                        (true_scat_flag) ? nLORs_trues[ 0 ] + nLORs_scatters[ 0 ] :
                        (true_rdm_flag)  ? nLORs_trues[ 0 ] + nLORs_rdms[ 0 ] :
                         nLORs_tot;

      Out_data_file->SetNbEvents( nLORs );
      Out_data_file->WriteHeader();
    }
  }  // end of PET section


  // Write the source image
  if(src_img_flag)
  {
    if (vb>=2) Cout("Writing source image ..." << endl);
    
    // Initialize Intf_fields object with the source image dimensions
    Intf_fields IF;
    IntfKeyInitFields(&IF);
    
    IF.mtx_size[0] = dim_src_img[0];
    IF.mtx_size[1] = dim_src_img[1];
    IF.mtx_size[2] = dim_src_img[2];
    IF.vox_size[0] = vox_src_img[0];
    IF.vox_size[1] = vox_src_img[1];
    IF.vox_size[2] = vox_src_img[2];
    IF.nb_total_imgs = dim_src_img[2];
    IF.image_start_time.push_back((FLTNB)(start_time_ms)/1000);
    IF.image_duration.push_back((FLTNB)(duration_ms)/1000);
    if (IntfWriteImgFile(path_to_src_image.append("_src"), p_src_img[ 0 ], IF, vb) )
    {
      Cerr("***** castor-GATERootToCastor :: Error writing Interfile of output image !" << endl);  
      Exit(EXIT_FAILURE);
    }
  }


  // Clock stop datafile generation
  clock_t clock_stop_file_gen = clock();
  time_t time_stop_file_gen = time(NULL);
      
  if (vb>=2) Cout ("Datafile Generation -> Time spent: " << time_stop_file_gen-time_start_file_gen
                           << " sec | CPU: " << (clock_stop_file_gen-clock_start_file_gen)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);


  Cout(endl << "The simulated dataset contained " << nLORs_tot           << " coincidences/singles, with: " << endl
            << "                                " << nLORs_trues[ 0 ]    << " trues ("    << 100.*((HPFLTNB)nLORs_trues[ 0 ])   /((HPFLTNB)nLORs_tot)    << " %), "  << endl
            << "                                " << nLORs_scatters[ 0 ] << " scatters (" << 100.*((HPFLTNB)nLORs_scatters[ 0 ])/((HPFLTNB)nLORs_tot) << " %),"  << endl
            << "                                " << nLORs_rdms[ 0 ]     << " randoms ("  << 100.*((HPFLTNB)nLORs_rdms[ 0 ])    /((HPFLTNB)nLORs_tot)   << " %),"  << endl
            << "                                " << nLORs_unknown[ 0 ]  << " unknown ("  << 100.*((HPFLTNB)nLORs_unknown[ 0 ]) /((HPFLTNB)nLORs_tot)   << " %)." << endl);
  
  
  // Merging and writing data
  if (vb>=2) Cout("Writing raw datafile ..." << endl); // todo just change name if only one datafile
  
  // Clock start datafile generation
  clock_t clock_start_data_w = clock();
  time_t time_start_data_w  = time(NULL);
  
  Out_data_file->PROJ_WriteData();  
  Out_data_file->PROJ_DeleteTmpDataFile();

  clock_t clock_stop_data_w = clock();
  time_t time_stop_data_w  = time(NULL);
  
  if (vb>=2) Cout ("Data writing -> Time spent: " << time_stop_data_w-time_start_data_w
                      << " sec | CPU: " << (clock_stop_data_w-clock_start_data_w)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);
                           
  Cout(endl << " --- Conversion completed --- " << endl);
  
  // ============================================================================================================
  // End
  // ============================================================================================================

  if (vb>=2) Cout(" Deallocating objects ..." << endl);
      
      
  // Delete objects
  if(nLORs_trues)    delete[] nLORs_trues;
  if(nLORs_rdms)     delete[] nLORs_rdms;
  if(nLORs_scatters) delete[] nLORs_scatters;
  if(nLORs_unknown)  delete[] nLORs_unknown;
  if(nBINs)          delete[] nBINs;
  
  if(dt_max) delete[] dt_max;
  if(dt_min) delete[] dt_min;
  
  if(scat_corr_flag)
  {
    for (int th=0 ; th < p_ID->GetNbThreadsForProjection() ; th++)
    {
      for(size_t b=0; b<bins_elts ; b++)
        delete[] p_scat[th][b];

      delete[] p_scat[th];
    }
    if(p_scat) delete[] p_scat;
  }

  if(rdm_corr_flag)
  {
    for (int th=0 ; th < p_ID->GetNbThreadsForProjection() ; th++)
    {
      for(size_t b=0; b<bins_elts ; b++)
        delete[] p_rdm[th][b];
  
      delete[]p_rdm[th];
    }
    if(p_rdm)  delete[] p_rdm;
  }

  delete[] offset_time_list;
  
  if(histo_out_flag)
  {
    for (int th=0 ; th < p_ID->GetNbThreadsForProjection() ; th++)
    {
      for(size_t b=0; b<bins_elts ; b++)
        if(p_bins[ th ][ b ]) delete[] p_bins[ th ][ b ];
  
      if(p_bins[ th ]) delete[] p_bins[ th ];
    }
    if(p_bins) delete[] p_bins;
    
  }

  if(Out_data_file)
  {
    if(GATE_system_type == GATE_SYS_SPECT)
      delete (iDataFileSPECT*) Out_data_file;
    else
      if(Out_data_file) delete (iDataFilePET*) Out_data_file;
  }
  
  if(src_img_flag && p_src_img)
  {
    for (int th=0 ; th < p_ID->GetNbThreadsForProjection() ; th++)
      delete[] p_src_img[ th ];
      
    delete[] p_src_img;
  }
      
  // Free objects created for analytic projection (acf estimation)
  if(estimate_acf_flag)
    p_ImageSpace->LMS_DeallocateImage();
  if(p_ImageSpace)       delete p_ImageSpace;
  if(p_ProjectorManager) delete p_ProjectorManager;
  if(p_ID)               delete p_ID;


  // Main clock computation
  clock_t clock_stop_main = clock();
  time_t time_stop_main = time(NULL);

  if (vb>=2) Cout ("Total Time spent: " << time_stop_main-time_start_main
                     << " sec | CPU: " << (clock_stop_main-clock_start_main)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);
                         
  Cout(" ---          END         --- " << endl << endl);
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InitTrootPet
  \param ap_rt : Structure to initialize
  \brief Init the variable of the root tree structure
*/
void InitTrootPet(Troot_pet* ap_rt)
{
  // cylindricalPET specific
  ap_rt->rsectorID1 = 0;
  ap_rt->rsectorID2 = 0;
  ap_rt->moduleID1 = 0;
  ap_rt->moduleID2 = 0;
  ap_rt->submoduleID1 = 0;
  ap_rt->submoduleID2 = 0;
  ap_rt->layerID1 = 0;
  ap_rt->layerID2 = 0;
  // ecat specific
  ap_rt->blockID1 = 0; 
  ap_rt->blockID2 = 0;
  // common
  ap_rt->crystalID1 = 0;
  ap_rt->crystalID2 = 0;
  ap_rt->eventID1 = 0; 
  ap_rt->eventID2 = 0;
  ap_rt->comptonPhantomID1 = 0;
  ap_rt->comptonPhantomID2 = 0;
  ap_rt->rayleighPhantomID1 = 0;
  ap_rt->rayleighPhantomID2 = 0;
  ap_rt->time1 = 0.;
  ap_rt->time2 = 0.;
  ap_rt->sourceID1 = 0; 
  ap_rt->sourceID2 = 0;
  ap_rt->sourcePosX1 = 0.;
  ap_rt->sourcePosY1 = 0.;
  ap_rt->sourcePosZ1 = 0.;
  ap_rt->bool_ecat_mode = false;
}




/*
  \fn InitTrootSpect
  \param ap_rt : Structure to initialize
  \brief Init the variable of the root tree structure
*/
void InitTrootSpect(Troot_spect* ap_rt)
{
  // SPECThead specific
  ap_rt->rotAngle = 0.;
  ap_rt->gPosX = 0.;
  ap_rt->gPosY = 0.;
  ap_rt->gPosZ = 0.;
  ap_rt->headID = 0;
  ap_rt->pixelID = 0;
  // common
  ap_rt->crystalID = 0;
  ap_rt->eventID = 0;
  ap_rt->comptonPhantomID = 0;
  ap_rt->rayleighPhantomID = 0;
  ap_rt->time = 0.;
  ap_rt->sourceID = 0;
  ap_rt->sourcePosX = 0.;
  ap_rt->sourcePosY = 0.;
  ap_rt->sourcePosZ = 0.;
}


#ifdef CASTOR_ROOT

/*
  \fn SetROOTPETBranch
  \param ap_tree : Root tree to use
  \param ap_rt : Structure to use for the branching
  \brief Set the branch of the PET root tree to the structure variables
*/
void SetROOTPETBranch(TTree* ap_tree, Troot_pet* ap_rt)
{
  ap_tree->SetBranchAddress("time1",           &ap_rt->time1);
  ap_tree->SetBranchAddress("crystalID1",      &ap_rt->crystalID1);
  ap_tree->SetBranchAddress("comptonPhantom1", &ap_rt->comptonPhantomID1);
  ap_tree->SetBranchAddress("RayleighPhantom1",&ap_rt->rayleighPhantomID1);
  ap_tree->SetBranchAddress("eventID1",        &ap_rt->eventID1);
  ap_tree->SetBranchAddress("sourceID1",       &ap_rt->sourceID1);
          
  ap_tree->SetBranchAddress("time2",           &ap_rt->time2);
  ap_tree->SetBranchAddress("crystalID2",      &ap_rt->crystalID2);    
  ap_tree->SetBranchAddress("comptonPhantom2", &ap_rt->comptonPhantomID2);
  ap_tree->SetBranchAddress("RayleighPhantom2",&ap_rt->rayleighPhantomID2);
  ap_tree->SetBranchAddress("eventID2",        &ap_rt->eventID2);
  ap_tree->SetBranchAddress("sourceID2",       &ap_rt->sourceID2);

  ap_tree->SetBranchAddress("sourcePosX1",     &ap_rt->sourcePosX1);
  ap_tree->SetBranchAddress("sourcePosY1",     &ap_rt->sourcePosY1);
  ap_tree->SetBranchAddress("sourcePosZ1",     &ap_rt->sourcePosZ1);
  
  // ECAT
  if(ap_rt->bool_ecat_mode)
  {
    ap_tree->SetBranchAddress("blockID1",     &ap_rt->blockID1); 
    ap_tree->SetBranchAddress("blockID2",     &ap_rt->blockID2); 
  }
  // CylindricalPET
  else
  {
    ap_tree->SetBranchAddress("rsectorID1",   &ap_rt->rsectorID1);
    ap_tree->SetBranchAddress("moduleID1",    &ap_rt->moduleID1);
    ap_tree->SetBranchAddress("submoduleID1", &ap_rt->submoduleID1); 
    ap_tree->SetBranchAddress("layerID1",     &ap_rt->layerID1);

    ap_tree->SetBranchAddress("rsectorID2",   &ap_rt->rsectorID2);
    ap_tree->SetBranchAddress("moduleID2",    &ap_rt->moduleID2);
    ap_tree->SetBranchAddress("submoduleID2", &ap_rt->submoduleID2); 
    ap_tree->SetBranchAddress("layerID2",     &ap_rt->layerID2);
  }
}



/*
  \fn SetROOTSPECTBranch
  \param ap_tree : Root tree to use
  \param ap_rt : Structure to use for the branching
  \brief Set the branch of the SPECT root tree to the structure variables
*/
void SetROOTSPECTBranch(TTree* ap_tree, Troot_spect* ap_rt)
{          
  ap_tree->SetBranchAddress("time",           &ap_rt->time);
  ap_tree->SetBranchAddress("headID",         &ap_rt->headID);
  ap_tree->SetBranchAddress("crystalID",      &ap_rt->crystalID);
  ap_tree->SetBranchAddress("pixelID",        &ap_rt->pixelID);
  ap_tree->SetBranchAddress("globalPosX",     &ap_rt->gPosX);
  ap_tree->SetBranchAddress("globalPosY",     &ap_rt->gPosY);
  ap_tree->SetBranchAddress("globalPosZ",     &ap_rt->gPosZ);
  ap_tree->SetBranchAddress("rotationAngle",  &ap_rt->rotAngle);
  ap_tree->SetBranchAddress("comptonPhantom", &ap_rt->comptonPhantomID);
  ap_tree->SetBranchAddress("RayleighPhantom",&ap_rt->rayleighPhantomID);
  ap_tree->SetBranchAddress("eventID",        &ap_rt->eventID);
  ap_tree->SetBranchAddress("sourcePosX",     &ap_rt->sourcePosX);
  ap_tree->SetBranchAddress("sourcePosY",     &ap_rt->sourcePosY);
  ap_tree->SetBranchAddress("sourcePosZ",     &ap_rt->sourcePosZ);
}
#endif




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn RecoverSPECTGeometry
  \param a_inputDataIntf : path to datafile type
  \param a_pathIntf : path to the interfile header
  \param a_pathMac : path to a macro file
  \param distToDetector : distance between center of rotation and detector surface
  \param nHeads : nb of cameras
  \param nPixAxl : nb of transaxial pixels
  \param nPixTrs : nb of axial pixels
  \param crystalSizeAxl : crystal axial dimension
  \param crystalSizeTrs : crystal transaxial dimension
  \param nProjectionsTot : Nb total projections
  \param nProjectionsByHead : Nb total projections by head
  \param head1stAngle : head first angle
  \param headAngPitch : angular pitch between heads
  \param headRotSpeed : angle between projection
  \param headRotDirection : rotation direction of the head (CW or CCW)
  \param start_time_ms : acquisition start time converted in ms
  \param duration_ms : acquisition duration converted in ms
  \param nCrystalsTot : total nb detectors
  \param p_proj_spect_image : array for spect projection image
  \param spect_bin_axl : gamma camera axial (virtual) pixels
  \param spect_bin_trs : gamma camera transaxial (virtual) pixels
  \param nbSimulatedPixels : total nb virtual pixels
  \param a_histoOutFlag : output mode histogram flag
  \param p_bins : prompt histogram
  \param p_scat : scatter histograms
  \param bins_elts : total nb bins elements
  \param a_scatCorrFlag : scatter correction flag
  \param a_nbThreads,
  \param vb : verbosity
  \brief Call dedicated functions to recover SPECT geometry information from the provided datafile
  \return 0 if success, positive value otherwise
*/
int RecoverSPECTGeometry(bool       a_inputDataIntf,
                        string      a_pathIntf,
                        string      a_pathMac,
                        float_t     &distToDetector,
                        uint32_t    &nHeads,
                        uint32_t    &nPixelsAxial,
                        uint32_t    &nPixelsTransaxial,
                        float_t     &crystalSizeAxl,
                        float_t     &crystalSizeTrs,
                        uint32_t    &nProjectionsTot,
                        uint32_t    &nProjectionsByHead,
                        float_t     &head1stAngleDeg,
                        float_t     &headAngPitchDeg,
                        float_t     &headAngStepDeg,
                        int         &headRotDirection,
                        uint32_t    &start_time_ms,
                        uint32_t    &duration_ms,
                        uint32_t    &nCrystalsTot,
                        FLTNB*      &p_proj_spect_image,
                        uint16_t    &spect_bin_axl,
                        uint16_t    &spect_bin_trs,
                        uint32_t    &nbSimulatedPixels,
                        bool        a_histoOutFlag,
                        uint8_t  ***&p_bins,
                        uint8_t  ***&p_scat,
                        uint32_t    &bins_elts,
                        bool        a_scatCorrFlag,
                        int         a_nbThreads,
                        int         vb)
{
  duration_ms =0.;
   
  // Data provided as an interfile projection image
  if(a_inputDataIntf)
  {
    if(ReadIntfSPECT( a_pathIntf,
                             distToDetector,
                                     nHeads,
                               nPixelsAxial,
                          nPixelsTransaxial,
                             crystalSizeAxl,
                             crystalSizeTrs,
                            nProjectionsTot,
                         nProjectionsByHead,
                            head1stAngleDeg,
                            headAngPitchDeg,
                             headAngStepDeg,
                           headRotDirection,
                              start_time_ms,
                                duration_ms,
                                         vb) )
    {
      Cerr("***** castor-GATERootToCastor::RecoverSPECTGeometry()-> Error when reading interfile : "<< a_pathIntf <<" !" << endl);
      Exit(EXIT_FAILURE);
    }
    
    nCrystalsTot = nPixelsAxial * nPixelsTransaxial;

    // Read Interfile
    Intf_fields IF_proj_spect_data;
    p_proj_spect_image = new FLTNB[nHeads*nProjectionsByHead*nCrystalsTot];

    if( IntfReadProjectionImage(a_pathIntf, p_proj_spect_image, &IF_proj_spect_data, vb, false) )
    {
      Cerr("***** castor-GATERootToCastor::RecoverSPECTGeometry()->  Error when trying to read image : "<< a_pathIntf <<" !" << endl);
      Exit(EXIT_FAILURE);
    }      
  }
  
  // Data provided as a root file
  else
  {
    if(ReadMacSPECT( a_pathMac,
                distToDetector,
                        nHeads,
                  nPixelsAxial,
             nPixelsTransaxial,
                crystalSizeAxl,
                crystalSizeTrs,
               nProjectionsTot,
            nProjectionsByHead,
               head1stAngleDeg,
               headAngPitchDeg,
                headAngStepDeg,
              headRotDirection,
                 start_time_ms,
                   duration_ms,
                            vb) )
    {
      Cerr("***** castor-GATERootToCastor::RecoverSPECTGeometry()->  Error when reading mac file : "<< a_pathMac <<" !" << endl);
      Exit(EXIT_FAILURE);
    }
    
    nbSimulatedPixels = nPixelsAxial*nPixelsTransaxial;

    if((spect_bin_trs>0 || spect_bin_axl>0)   &&   nbSimulatedPixels > 1 )
    {
      Cerr("***** castor-GATERootToCastor::RecoverSPECTGeometry()->  WARNING : Spect bins have been initialized, but the simulation already provide a specific number of pixels (="<< nPixelsAxial*nPixelsTransaxial <<") !"<< endl <<
           "                                                                   Pixel matrix used by default !" << endl);
      Exit(EXIT_FAILURE);
    }
    else // Check bins have been provided, and initialize pixel sizes with provided values
    {
      if(spect_bin_trs == 0 && spect_bin_axl == 0 && nbSimulatedPixels==1)
      {
        Cerr("***** castor-GATERootToCastor::RecoverSPECTGeometry()-> : Error : Axial and transaxial bins values expected (use option -sp_bins) !"<< endl);
        Exit(EXIT_FAILURE);
      }

      // check crystal sizes have been correctly read
      if(crystalSizeAxl<0 || crystalSizeTrs<0)
      {
        Cerr("***** castor-GATERootToCastor::RecoverSPECTGeometry()->  Crystal dimensions not correctly read in the mac files !" << endl);
        Exit(EXIT_FAILURE);
      }
      
      nPixelsTransaxial = spect_bin_trs;
      nPixelsAxial      = spect_bin_axl;
      
      if(vb>=2) Cout("castor-GATERootToCastor::RecoverSPECTGeometry()-> Transaxial/Axial nb pixels : " << nPixelsTransaxial << " , " << nPixelsAxial << endl << 
                     "castor-GATERootToCastor::RecoverSPECTGeometry()-> Transaxial/Axial pixel sizes : " << crystalSizeTrs/nPixelsTransaxial << " , " << crystalSizeAxl/nPixelsAxial << endl);
    }
  }
  
  nCrystalsTot = nPixelsAxial * nPixelsTransaxial;
  
  if(vb>=2) Cout("castor-GATERootToCastor::RecoverSPECTGeometry()-> Number of Projections: " << nProjectionsByHead << endl <<
                 "castor-GATERootToCastor::RecoverSPECTGeometry()-> Detected number of crystals in the system : " << nCrystalsTot << endl);


  // Histogram bin vector initialization using the total number of projections & pixels
  if(a_histoOutFlag)
  {
    Cout("castor-GATERootToCastor::RecoverSPECTGeometry()-> Allocating memory for histogram... " << endl);
    
    bins_elts = nHeads*nProjectionsByHead;
   
    p_bins = new uint8_t**[a_nbThreads];
   
    for(int th=0 ; th< a_nbThreads ; th++)
    {
      p_bins[th] = new uint8_t*[bins_elts];
      
      for (size_t p=0; p<bins_elts ; p++)
      {
        p_bins[th][p] = new uint8_t[nCrystalsTot];
        
        for (size_t c=0; c<nCrystalsTot ; c++)
        {
          p_bins[th][p][c] = a_inputDataIntf  && (th==0) ? 
                             (uint8_t)p_proj_spect_image[p*nCrystalsTot+c] : 
                             0 ;
          
        }
      }
    }
  }
  
  if(a_scatCorrFlag && !a_inputDataIntf)
  {
    Cout("castor-GATERootToCastor::RecoverSPECTGeometry()-> Allocating memory for scatter histogram... " << endl);
         
    bins_elts = nHeads*nProjectionsByHead;
    
    p_scat = new uint8_t**[a_nbThreads];
    
    for(int th=0 ; th< a_nbThreads ; th++)
    {
      p_scat[th] = new uint8_t*[bins_elts];
      for (size_t p=0; p<bins_elts ; p++)
      {
        p_scat[th][p] = new uint8_t[nCrystalsTot];
        
        for (size_t c=0; c<nCrystalsTot ; c++)
          p_scat[th][p][c] = 0. ;
      }
    }
    
    
    Cout("castor-GATERootToCastor::RecoverSPECTGeometry()-> Memory allocation for bins completed " << endl << endl);
  }
    
  delete[] p_proj_spect_image;
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================



/*
  \fn ProcessSPECTEvent
  \param Rt
  \param castorID1,
  \param castorID2,
  \param headAngPitchDeg,
  \param headAngStepDeg,
  \param nbSimulatedPixels,
  \param nPixelsTransaxial,
  \param nPixelsAxial,
  \param crystalSizeAxl,
  \param crystalSizeTrs,
  \param nProjectionsByHead,
  \param kind,
  \param nLORs_unknown,
  \param nLORs_trues,
  \param nLORs_scatters,
  \param nLORs_rdms,
  \param vb )
  \brief Compute CASToR IDs for a SPECT event
  \return 0 if success, positive value otherwise
*/
int ProcessSPECTEvent(  Troot_spect *Rt,
                        uint32_t   &castorID1,
                        uint32_t   &castorID2,
                        float_t     headAngPitchDeg,
                        float_t     headAngStepDeg,
                        uint32_t    nbSimulatedPixels,
                        uint32_t    nPixelsTransaxial,
                        uint32_t    nPixelsAxial,
                        float_t     crystalSizeAxl,
                        float_t     crystalSizeTrs,
                        uint32_t    nProjectionsByHead,
                        uint8_t    &kind,
                        uint64_t   &nLORs_unknown,
                        uint64_t   &nLORs_trues,
                        uint64_t   &nLORs_scatters,
                        uint64_t   &nLORs_rdms,
                        int         vb )
{
  if(vb >= 4) 
  {
    Cout("Projection ID : headID: " << Rt->headID << ", rotation angle (deg): " << Rt->rotAngle << endl;);
    Cout("Crystal ID : crystalID: " << Rt->crystalID << ", pixelID: " << Rt->pixelID << ", globalPosX " << Rt->gPosX <<  ", globalPosY " << Rt->gPosY <<  ", globalPosZ " << Rt->gPosZ << endl;);
  }
  
  // Get projection ID
  castorID1 = ConvertIDSPECTRoot1(Rt->headID, 
                                  Rt->rotAngle, 
                                  headAngStepDeg,
                                  nProjectionsByHead);
  
  // Get crystal ID
  castorID2 = ConvertIDSPECTRoot2(nbSimulatedPixels, 
                                  nPixelsTransaxial, 
                                  nPixelsAxial, 
                                  Rt->headID, 
                                  Rt->crystalID, 
                                  Rt->pixelID, 
                                  Rt->rotAngle,
                                  headAngPitchDeg,
                                  crystalSizeAxl,
                                  crystalSizeTrs,
                                  Rt->gPosX, 
                                  Rt->gPosY, 
                                  Rt->gPosZ);
 
  if(vb >= 4) 
  {
    Cout("--> castorID1: " << castorID1 << endl;);
    Cout("--> castorID2: " << castorID2 << endl;);
  }
  
  // Find out the kind of coincidence (true, scatter, random)
  kind = ComputeKindGATEEvent(Rt->eventID, Rt->eventID, Rt->comptonPhantomID, Rt->comptonPhantomID, Rt->rayleighPhantomID, Rt->rayleighPhantomID);
  
  if(vb >= 5) 
    Cout("--> Kind: " << static_cast<int16_t>(kind) << endl;);

  // Count nb LORs according to kind
  switch (kind)
  {
    case 1:
      nLORs_trues++;
      break;

    case 2:
      nLORs_scatters++;
      break;

    case 3:
      nLORs_scatters++;
      break;

    case 4:
      nLORs_rdms++;
      break;

    default:
      nLORs_unknown++;
  }   

  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ProcessPETEvent
  \param GATE_system_type,
  \param castorID1,
  \param castorID2,
  \param nLayers,
  \param nb_crystal_per_layer,
  \param nCrystalsTot,
  \param nCrystalsAxial, 
  \param nCrystalsTransaxial, 
  \param nLayersRptAxial,
  \param nLayersRptTransaxial,
  \param nSubmodulesAxial, 
  \param nSubmodulesTransaxial, 
  \param nModulesAxial, 
  \param nModulesTransaxial, 
  \param nRsectorsAxial,
  \param nRsectorsAngPos,
  \param nBlocksLine, 
  \param nBlocksPerRing,
  \param invert_det_order,
  \param rsector_id_order,
  \param kind,
  \param nLORs_unknown,
  \param nLORs_trues,
  \param nLORs_scatters,
  \param nLORs_rdms, 
  \param vb 
  \brief Compute CASToR IDs for a PET event
  \return 0 if success, positive value otherwise
*/
int ProcessPETEvent(    int         GATE_system_type,
                        Troot_pet   *Rt,
                        uint32_t   &castorID1,
                        uint32_t   &castorID2,
                        uint8_t     nLayers,
                        uint32_t    *nb_crystal_per_layer,
                        uint32_t    nCrystalsTot,
                        uint32_t    nCrystalsAxial, 
                        uint32_t    nCrystalsTransaxial, 
                 vector<uint32_t>   nLayersRptAxial,
                 vector<uint32_t>   nLayersRptTransaxial,
                        uint32_t    nSubmodulesAxial, 
                        uint32_t    nSubmodulesTransaxial, 
                        uint32_t    nModulesAxial, 
                        uint32_t    nModulesTransaxial, 
                        uint32_t    nRsectorsAxial,
                        uint32_t    nRsectorsAngPos,
                        uint32_t    nBlocksLine, 
                        uint32_t    nBlocksPerRing,
                        bool        invert_det_order,
                        int         rsector_id_order,
                        uint8_t    &kind,
                        uint64_t   &nLORs_unknown,
                        uint64_t   &nLORs_trues,
                        uint64_t   &nLORs_scatters,
                        uint64_t   &nLORs_rdms, 
                        int         vb )

{
  // ID Conversions
  if (GATE_system_type == GATE_SYS_CYLINDRICAL)
  {
    if(vb >= 4)
    {
      Cout("Crystal 1 : RsectorID: " << Rt->rsectorID1 << ", moduleID: " << Rt->moduleID1 << ", submoduleID: " << Rt->submoduleID1 << ", crystalID: " << Rt->crystalID1
       <<  ", layerID: " << Rt->layerID1 << endl;);
      Cout("Crystal 2 :RsectorID: " << Rt->rsectorID2 << ", moduleID2: " << Rt->moduleID2 << ", submoduleID: " << Rt->submoduleID2 << ", crystalID: " << Rt->crystalID2
       <<  ", layerID: " << Rt->layerID2 << endl;);
    }
    
    castorID1 = ConvertIDcylindrical(nRsectorsAngPos, nRsectorsAxial, invert_det_order, rsector_id_order,
                                     nModulesTransaxial, nModulesAxial, 
                                     nSubmodulesTransaxial, nSubmodulesAxial, 
                                     nCrystalsTransaxial, nCrystalsAxial,
                                     nLayers, nb_crystal_per_layer, 
                                     nLayersRptTransaxial, nLayersRptAxial,
                                     Rt->layerID1, Rt->crystalID1, Rt->submoduleID1, Rt->moduleID1, Rt->rsectorID1);

    castorID2 = ConvertIDcylindrical(nRsectorsAngPos, nRsectorsAxial, invert_det_order, rsector_id_order,
                                     nModulesTransaxial, nModulesAxial, 
                                     nSubmodulesTransaxial, nSubmodulesAxial, 
                                     nCrystalsTransaxial, nCrystalsAxial,
                                     nLayers, nb_crystal_per_layer, 
                                     nLayersRptTransaxial, nLayersRptAxial,
                                     Rt->layerID2, Rt->crystalID2, Rt->submoduleID2, Rt->moduleID2, Rt->rsectorID2);
                                     
    // Small check that the crystal IDs have consistent values
    if(castorID1 >= nCrystalsTot )
    {
      Cerr("***** castor-GATERootToCastor :: Issue detected during the conversion! " << endl);
      Cerr("*****                            detector ID larger than the max nb of detectors: " << nCrystalsTot << endl);
      Cerr("*****                            please email to the castor-users mailing list with this log and the mac file" << endl);
      Cerr("*****                            rsectorID1:"<< Rt->rsectorID1<< endl);
      Cerr("*****                            moduleID1:"<< Rt->moduleID1<< endl);
      Cerr("*****                            submoduleID1:"<< Rt->submoduleID1<< endl);
      Cerr("*****                            crystalID1:"<< Rt->crystalID1<< endl);
      Cerr("*****                            layerID1:"<< Rt->layerID1<< endl);
      Cerr("*****                            computed castorID:"<< castorID1<< endl);
      Exit(EXIT_FAILURE);
    }
  
    if(castorID2 >= nCrystalsTot )
    {
      Cerr("***** castor-GATERootToCastor :: Issue detected during the conversion! " << endl);
      Cerr("*****                            detector ID larger than the max nb of detectors: " << nCrystalsTot << endl);
      Cerr("*****                            please email to the castor-users mailing list with this log and the mac file" << endl);
      Cerr("*****                            rsectorID2:"<< Rt->rsectorID2<< endl);
      Cerr("*****                            moduleID2:"<< Rt->moduleID2<< endl);
      Cerr("*****                            submoduleID2:"<< Rt->submoduleID2<< endl);
      Cerr("*****                            crystalID2:"<< Rt->crystalID2<< endl);
      Cerr("*****                            layerID2:"<< Rt->layerID2<< endl);
      Cerr("*****                            computed castorID:"<< castorID2<< endl);
      Exit(EXIT_FAILURE);
    }
    
    
    if(vb >= 4)
    {
      Cout("--> castorID1: " << castorID1 << endl;);
      Cout("--> castorID2: " << castorID2 << endl;);
    }
  }
  else if (GATE_system_type == GATE_SYS_ECAT)
  {
    if(vb >= 4)
    {
      Cout("Crystal 1 : BlockID: " << Rt->blockID1 << ", crystalID: " << Rt->crystalID1 << endl;);
      Cout("Crystal 2 : BlockID: " << Rt->blockID2 << ", crystalID: " << Rt->crystalID2 << endl;);
    }
    
    castorID1 = ConvertIDecat(nBlocksPerRing, nBlocksLine, nCrystalsTransaxial, nCrystalsAxial, Rt->crystalID1, Rt->blockID1);
    castorID2 = ConvertIDecat(nBlocksPerRing, nBlocksLine, nCrystalsTransaxial, nCrystalsAxial, Rt->crystalID2, Rt->blockID2);
    
    if(vb >= 4)
    {
      Cout("--> castorID1: " << castorID1 << endl;);
      Cout("--> castorID2: " << castorID2 << endl;);
    }
    
  }
          
          
  // Find out the kind of coincidence (true, scatter, random)
  kind = ComputeKindGATEEvent(Rt->eventID1, Rt->eventID2, Rt->comptonPhantomID1, Rt->comptonPhantomID2, Rt->rayleighPhantomID1, Rt->rayleighPhantomID2);
  
  // Count nb LORs according to kind
  switch (kind)
  {
    case 1:
      nLORs_trues++;
      break;

    case 2:
      nLORs_scatters++;
      break;

    case 3:
      nLORs_scatters++;
      break;

    case 4:
      nLORs_rdms++;
      break;

    default:
      nLORs_unknown++;
  }   
      
  return 0;
}
      




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

/*
  \fn RecoverPETGeometry
  \param GATE_system_type,
  \param a_pathMac,
  \param nLayers,
  \param nb_crystal_per_layer,
  \param nCrystalsTot,
  \param nCrystalsAxial, 
  \param nCrystalsTransaxial, 
  \param nLayersRptAxial,
  \param nLayersRptTransaxial,
  \param nSubmodulesAxial, 
  \param nSubmodulesTransaxial, 
  \param nModulesAxial, 
  \param nModulesTransaxial, 
  \param nRsectorsAxial,
  \param nRsectorsAngPos,
  \param nBlocksLine, 
  \param nBlocksPerRing,
  \param invert_det_order,
  \param rsector_id_order,
  \param pet_coinc_window,
  \param start_time_ms,
  \param duration_ms,
  \param a_histoOutFlag,
  \param p_bins,
  \param p_scat,
  \param p_rdm,
  \param bins_elts,
  \param a_scatCorrFlag,
  \param a_rdmCorrFlag,
  \param a_nbThreads,
  \param vb : verbosity
  \brief Call dedicated functions to recover PET geometry information from the provided datafile
  \return 0 if success, positive value otherwise
*/
int RecoverPETGeometry( int         GATE_system_type,
                        string      a_pathMac,
                        uint8_t     &nLayers,
                        uint32_t    *nb_crystal_per_layer,
                        uint32_t    &nCrystalsTot,
                        uint32_t    &nCrystalsAxial, 
                        uint32_t    &nCrystalsTransaxial, 
                 vector<uint32_t>   &nLayersRptAxial,
                 vector<uint32_t>   &nLayersRptTransaxial,
                        uint32_t    &nSubmodulesAxial, 
                        uint32_t    &nSubmodulesTransaxial, 
                        uint32_t    &nModulesAxial, 
                        uint32_t    &nModulesTransaxial, 
                        uint32_t    &nRsectorsAxial,
                        uint32_t    &nRsectorsAngPos,
                        uint32_t    &nBlocksLine, 
                        uint32_t    &nBlocksPerRing,
                        bool        &invert_det_order,
                        int         &rsector_id_order,
                        FLTNB       &pet_coinc_window,
                        uint32_t    &start_time_ms,
                        uint32_t    &duration_ms,
                        bool        a_histoOutFlag,
                        uint8_t  ***&p_bins,
                        uint8_t  ***&p_scat,
                        uint8_t  ***&p_rdm,
                        uint32_t    &bins_elts,
                        bool        a_scatCorrFlag,
                        bool        a_rdmCorrFlag,
                        int         a_nbThreads,
                        int         vb
                        )
{
  // cylindricalPET system
  if(GATE_system_type == GATE_SYS_CYLINDRICAL)
  {
    if(ReadMacCylindrical(       a_pathMac,
                                   nLayers,
                      nb_crystal_per_layer,
                              nCrystalsTot,
                            nCrystalsAxial, 
                       nCrystalsTransaxial, 
                           nLayersRptAxial,
                      nLayersRptTransaxial,
                          nSubmodulesAxial, 
                     nSubmodulesTransaxial, 
                             nModulesAxial, 
                        nModulesTransaxial, 
                            nRsectorsAxial,
                           nRsectorsAngPos,
                          invert_det_order,
                          rsector_id_order, 
                             start_time_ms,
                               duration_ms,
                          pet_coinc_window,
                                        vb) )
    {
      Cerr("***** castor-GATERootToCastor::RecoverPETGeometry()-> Error when reading mac file : "<< a_pathMac <<" !" << endl);
      Exit(EXIT_FAILURE);
    } 
  }
  
  else // ECAT
  {
    // Reading the macro file
    if(ReadMacECAT(       a_pathMac,
                       nCrystalsTot,
                     nCrystalsAxial, 
                nCrystalsTransaxial,
                        nBlocksLine, 
                     nBlocksPerRing, 
                      start_time_ms,
                        duration_ms,
                   pet_coinc_window,
                                 vb) )
    {
      Cerr("***** castor-GATERootToCastor::RecoverPETGeometry()-> Error when reading mac file : "<< a_pathMac <<" !" << endl);
      Exit(EXIT_FAILURE);
    }
    
  }

  if(vb>=2) Cout("castor-GATERootToCastor::RecoverPETGeometry()-> Detected number of crystals in the system : " << nCrystalsTot << endl);
  
  
  // Histogram bin vector initialization using the total number of crystals
  if(a_histoOutFlag || a_scatCorrFlag)
  {
    Cout("castor-GATERootToCastor::RecoverPETGeometry()-> Allocating memory for histogram... " << endl <<
         "castor-GATERootToCastor::RecoverPETGeometry()-> Warning : this step can require huge amount of memory if the system contains a high number of crystals !" << endl);
         
    bins_elts = nCrystalsTot;
    
    p_bins = new uint8_t**[a_nbThreads];
    
    for(int th=0 ; th< a_nbThreads ; th++)
    {
      p_bins[th] = new uint8_t*[bins_elts];
      for (size_t c=0; c<bins_elts ; c++)
      {
        p_bins[th] [c] = new uint8_t[nCrystalsTot-c];
        
        for (size_t c2=0; c2<nCrystalsTot-c ; c2++)
          p_bins[th] [c][c2] = 0;
      }
    }
    
    Cout("castor-GATERootToCastor::RecoverPETGeometry()-> Memory allocation for bins completed " << endl << endl);
  }

  if(a_scatCorrFlag)
  {
    Cout("castor-GATERootToCastor::RecoverPETGeometry()-> Allocating memory for scatter histogram... " << endl <<
         "castor-GATERootToCastor::RecoverPETGeometry()-> Warning : this step can require huge amount of memory if the system contains a high number of crystals !" << endl);
         
    bins_elts = nCrystalsTot;
    
    p_scat = new uint8_t**[a_nbThreads];
        
    for(int th=0 ; th< a_nbThreads ; th++)
    {
      p_scat[th]  = new uint8_t*[bins_elts];
      for (size_t c=0; c<bins_elts ; c++)
      {
        p_scat[th] [c] = new uint8_t[nCrystalsTot-c];
        
        for (size_t c2=0; c2<nCrystalsTot-c ; c2++)
          p_scat[th] [c][c2] = 0;
      }
    }
    
    Cout("castor-GATERootToCastor::RecoverPETGeometry()-> Memory allocation for bins completed " << endl << endl);
  }
  
  if(a_rdmCorrFlag)
  {
    Cout("castor-GATERootToCastor::RecoverPETGeometry()-> Allocating memory for random histogram... " << endl <<
         "castor-GATERootToCastor::RecoverPETGeometry()-> Warning : this step can require huge amount of memory if the system contains a high number of crystals !" << endl);
         
    bins_elts = nCrystalsTot;
    
    p_rdm = new uint8_t**[a_nbThreads];
        
    for(int th=0 ; th< a_nbThreads ; th++)
    {
      p_rdm[th]  = new uint8_t*[bins_elts];
      for (size_t c=0; c<bins_elts ; c++)
      {
        p_rdm[th][c] = new uint8_t[nCrystalsTot-c];
        
        for (size_t c2=0; c2<nCrystalsTot-c ; c2++)
          p_rdm[th][c][c2] = 0;
      }
    }
    
    Cout("castor-GATERootToCastor::RecoverPETGeometry()-> Memory allocation for bins completed " << endl << endl);
  }
  
  return 0;
}



