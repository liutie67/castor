
/*!
  \file
  \ingroup main_programs
  \brief  This is the main program of the analytic projection tool for CASToR
          It reads/parses/checks the command-line options provided by the user, initialize each manager classes with the correct set of options, and launch the projection algorithm (managed by oAnalyticProjection).
*/

#include "gVariables.hh"
#include "gOptions.hh"
#include "oAnalyticProjection.hh"
#include "oImageConvolverManager.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "iDataFilePET.hh"
#include "iDataFileSPECT.hh"
#include "sOutputManager.hh"
#include "sScannerManager.hh"
#include "sRandomNumberGenerator.hh"
#include "iScannerPET.hh"


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelp()
  \brief   Display main command line options for castor-proj
*/
void ShowHelp(int a_returnCode, int a_mpiRank)
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif

  // Show help
  cout << endl;
  cout << "Usage:  castor-proj.exe  -img path_to_file.hdr  -sc name  -(f/d)out output   [settings]" << endl;
  cout << endl;
  cout << "[Main options]:" << endl;
  cout << "  -img path_to_img.hdr : Give the interfile header of the image to project (no default)." << endl;
  cout << "  -sc  scanner_alias   : Give the scanner model. It must correspond to the name of one of the file in the scanner repository (w/o file extension)" << endl;
  cout << "  -fout name           : Give the root name for all output files (no default, alternative to -dout)" << endl;
  cout << "  -dout name           : Give the name of the output directory where all output files will be written (no default, alternative to -fout)" << endl;
  cout << "  -fileType type       : Give the type of events that will be generated : ('TYPE_PET' ou '0' = PET,  'TYPE_SPECT' ou '1' = SPECT) (default : TYPE_PET)" << endl;
  cout << "  -vb                  : Give the verbosity level, from 0 (no verbose) to above 5 (at the event level) (default: 1)." << endl;
  cout << endl;
  cout << "[Specific options]:" << endl;
  cout << "  -help-in             : Print out specific help about input settings." << endl; // managed by main
  cout << "  -help-out            : Print out specific help about output settings." << endl; // managed by main
  cout << "  -help-proj           : Print out specific help about projector settings." << endl; // managed by main
  cout << "  -help-comp           : Print out specific help about computing settings." << endl; // managed by main
  cout << "  -help-misc           : Print out specific help about miscellaneous and verbose settings." << endl; // managed by main
  cout << "  -help-dynamic        : Print out specific help about dynamic options." << endl; // managed by main
  cout << "  -help-pet            : Print out specific help about PET settings for analytic projection." << endl; // managed by main
  cout << "  -help-spect          : Print out specific help about SPECT settings for analytic projection." << endl; // managed by main
  cout << endl;
  cout << "[Implemented Modules]:" << endl;
  cout << "  -help-scan           : Show the list of all scanners from the configuration directory." << endl; // managed by sScannerManager
  cout << "  -help-projm          : Show the list and description of all implemented projectors." << endl; // managed by oProjectorManager
  cout << endl;
  cout << "  --help,-h,-help      : Print out this help page." << endl; // managed by main
  cout << endl;
  #ifdef CASTOR_MPI
  cout << "  Compiled with MPI" << endl;
  #endif
  #ifdef CASTOR_OMP
  cout << "  Compiled with OpenMP" << endl;
  #endif
  #ifdef CASTOR_GPU
  cout << "  Compiled for GPU" << endl;
  #endif
  #if defined(CASTOR_OMP) || defined(CASTOR_MPI) || defined(CASTOR_GPU)
  cout << endl;
  #endif
  #ifdef BUILD_DATE
  cout << "  Build date: " << BUILD_DATE << endl;
  cout << endl;
  #endif
  #ifdef CASTOR_VERSION
  cout << "  This program is part of the CASToR release version " << CASTOR_VERSION << "." << endl;
  cout << endl;
  #endif
  Exit(a_returnCode);
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpInput()
  \brief   Display command line options related to input settings for castor-proj
*/
void ShowHelpInput()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "[Input settings]:" << endl;
  cout << endl;
  cout << "  -off x,y,z           : Give the offset of the field-of-view in each dimension, in mm (default: 0.,0.,0.)." << endl;
  cout << endl;
  cout << "  -atn path_to_img.hdr : Give the interfile header of an input attenuation image (unit has to be cm-1) (default: uniform value)." << endl;
  cout << endl;
  cout << "  -help-in             : Print out this help." << endl;
  cout << endl;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpOutput()
  \brief   Display command line options related to output settings for castor-proj
*/
void ShowHelpOutput()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "[Output settings]:" << endl;
  cout << endl;
  cout << "  -fout name           : Give the root name for all output files. All output files will be written as 'name_suffix.ext'." << endl;
  cout << "                         So the provided name should not end with '.' or '/' character. (no default, alternative to -dout)" << endl;
  cout << "  -dout name           : Give the name of the output directory where all output files will be written. All files will also" << endl;
  cout << "                         be prefixed by the name of the directory. The provided name should not end with '.' or '/' character." << endl;
  cout << "                         (no default, alternative to -fout)" << endl;
  cout << endl;
  cout << "  -olut                : If a scanner LUT (geometry information) is computed from a .geom file, it will be written on disk in the scanner repository" << endl;
  cout << endl;
  cout << "  -onbp prec            : By default, numbers are displayed using scientific format. This option allows to customize the format and precision" << endl;
  cout << "                        : The format is format,precision. f is the format (s=scientific, f=fixed), p is the precision" << endl;
  cout << "                          eg: -onbp f,5 --> fixed with 5 digits precision, -onbp -->  scientifici with max precision." << endl;
  cout << endl;
  cout << "  -help-out            : Print out this help." << endl;
  cout << endl;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpProj()
  \brief   Display command line options related to the Projector module for castor-proj
*/
void ShowHelpProj()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << "[Projector settings]:" << endl;
  cout << endl;
  cout << "  -proj  param         : Give the projector to be used for both forward and backward projections, along with a configuration file (proj:file.conf)" << endl;
  cout << "                         or the list of parameters associated to the projector (proj,param1,param2,...). If the projector only is specified, the" << endl;
  cout << "                         default configuration file is used. By default, the Siddon projector is used." << endl;
  cout << endl;
  cout << "  -proj-comp           : Give the strategy for projection line computation. Here are the three different strategies that can be used:" << endl;
  cout << "                     1 : Image-based system matrix elements storage: The voxels weights are added in a matrix representing the whole image, so" << endl;
  cout << "                         the addition of a new line to the previous ones is straightforward only by adding the weights to the corresponding voxels." << endl;
  cout << "                         As it is insanely long, it can possibly be used for example with extremely complex projectors that makes use of huge number" << endl;
  cout << "                         of ray tracings for a single event, where the list of contributions can become longer than the number of voxels in the image." << endl;
  cout << "                     2 : Fixed-size list storage of system matrix elements: The voxels are added one by one in two separated lists, one containing voxel" << endl;
  cout << "                         indices and the other voxel weights. When a voxel is added to the oProjectionLine, it is simply pilled-up to the list. The list" << endl;
  cout << "                         has a fixed size which is provided by the EstimateMaxNumberOfVoxelsPerLine() function from the vProjector class. There are no" << endl;
  cout << "                         ckecks at all for possible buffer overflows. This is the fastest strategy and default one." << endl;
  cout << "                     3 : Adaptative-size list storage of system matrix elements: This is the same as the fixed-size strategy except that the size can be" << endl;
  cout << "                         upgraded if the current number of contributing voxels exceed the list's size. The first allocated size corresponds to the diagonal" << endl;
  cout << "                         of the image. During the first iteration, this size will be upgraded until it will reach a size suitable for all events. Thus it" << endl;
  cout << "                         is a bit slower than the fixed-list strategy during the first iteration, but is optimal with respect to RAM usage." << endl;
  cout << endl;
  cout << "  -help-projm          : Print out specific help about projector settings." << endl; // managed by oProjectorManager
  cout << endl;
  cout << "  -help-proj           : Print out this help." << endl; // managed by oProjectorManager
  cout << endl;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpComputation()
  \brief   Display command line options related to the computation settings for castor-proj
*/
void ShowHelpComputation()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "[Computation settings]:" << endl;
  cout << endl;
  #ifdef CASTOR_OMP
  cout << "  -th param            : Set the number of threads for parallel computation (default: 1). If 0 is given, then the maximum number of available threads is automatically determined." << endl;
  cout << "                         Can also give two parameters separated by a comma (e.g. 16,4), to distinguish between the number of threads for projection and image operations respectively." << endl;
  cout << endl;
  #endif
  #ifdef CASTOR_GPU
  cout << "  -gpu                 : Flag to say that we want to use the GPU device (default: use the CPU only)." << endl;
  cout << endl;
  #endif
  cout << "  -conv  param         : Give an image convolver model to be used within the algorithm, along with a configuration file (conv:file.conf) or the" << endl;
  cout << "                         list of parameters associated to the convolver (conv,param1,param2,...). If the convolver only is specified, its default" << endl;
  cout << "                         configuration file is used. By default, no convolver is applied" << endl;
  cout << endl;
  cout << "  -help-conv           : Show the list and description of all implemented image convolvers." << endl; // managed by oImageConvolverManager
  cout << endl;
  cout << "  -noise poisson       : Give the noise model to be used. Only Poisson noise is implemented for the moment (default : no noise model)" << endl;
  cout << endl;
  cout << "  -rng                 : Give a seed for the random number generator (should be >=0)" << endl;
  cout << endl;
  cout << "  -no0                 : Discard zero events (not saved)." << endl;
  cout << endl;
  cout << "  -help-comp           : Print out this help." << endl;
  cout << endl;
  #ifdef CASTOR_MPI
  cout << "  Compiled with MPI" << endl;
  #endif
  #ifdef CASTOR_OMP
  cout << "  Compiled with OpenMP" << endl;
  #endif
  #ifdef CASTOR_GPU
  cout << "  Compiled for GPU" << endl;
  #endif
  #if defined(CASTOR_OMP) || defined(CASTOR_MPI) || defined(CASTOR_GPU)
  cout << endl;
  #endif
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpDynamic()
  \brief   Display command line options related to the dynamic features for castor-proj
*/
void ShowHelpDynamic()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "[Dynamic settings]:" << endl;
  cout << endl;
  cout << "  -frm list            : Give the framing for the reconstruction where 'list' is a list of all frames durations in seconds separated" << endl;
  cout << "                         by commas. To specify a dead frame, add 'df' concatenated to the frame duration. Milliseconds maximum precision." << endl;
  cout << "                         (default: read from the 'image duration (sec)' key in the interfile image header" << endl;
  cout << endl;
  cout << "  (about gated images) : An image contening several respiratory/cardiac 'gates' can be provided to generate one datafile containing histograms of these gated images." << endl;
  cout << "                         The number of gates being directly read from the 'number of respiratory gates' and 'number of cardiac gates' keys in the interfile header" << endl;
  cout << endl; 
  cout << "  -help-dynamic        : Print out this help." << endl;
  cout << endl;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpMiscellaneous()
  \brief   Display command line settings related to SPECT modality for castor-proj
*/
void ShowHelpSPECT()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << "[SPECT settings]:" << endl;
  cout << endl;
  cout << "  -SPECT_bins  par  : Give two integer parameters corresponding to the transaxial number of bins, separated by a comma (nbBinsX, nbBinsY). " << endl;
  cout << "                      Default = transaxial number of pixels in the SPECT scanner file " << endl;
  cout << endl;
  cout << "  -SPECT_ang   par  : Give the total number of projections, followed by the first angle and the angle step in float" << endl;
  cout << "                      Format : nb_projs:first_angle,step_angle" << endl;
  cout << endl;
  cout << "  -SPECT_c_ang par  : SPECT custom angles : give the total number of projections, followed by the projections angles in float" << endl;
  cout << "                      Format : nb_projs:angle1,angle2,... " << endl;
  cout << endl;
  cout << "  -SPECT_radius par : Give a distance between center of rotation to detectors which will be used for all heads" << endl;
  cout << endl;
  cout << "  -SPECT_rot par    : Give the heads rotation direction. Accept two parameters: CW (clockwise, default), CCW (counter-clockwise)" << endl;
  cout << endl;
  cout << "  -help-spect       : Print out this specific help." << endl; // managed by oProjectorManager
  cout << endl;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpMiscellaneous()
  \brief   Display command line settings related to PET modality for castor-proj
*/
void ShowHelpPET()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << "[PET settings]:" << endl;
  cout << endl;
  cout << "  -PET_maxAxialDiffmm param : Set the maximum axial difference in mm between two crystals in a LOR (default = no restriction)" << endl;
  cout << endl;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelpMiscellaneous()
  \brief   Display command line options related to miscellaneous settings for castor-proj
*/
void ShowHelpMiscellaneous()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << endl;
  cout << "[Miscellaneous settings]:" << endl;
  cout << endl;
  cout << "  -vb                  : Give the general verbosity level, from 0 (no verbose) to 5 and above (at the event level) (default: 1)." << endl;
  cout << "  -vb-algo             : Give the verbose level specific to the analytic projection algorithm (default: same as general verbose level)." << endl;
  cout << "  -vb-proj             : Give the verbose level specific to the projector (default: same as general verbose level)." << endl;
  cout << "  -vb-conv             : Give the verbose level specific to the image convolver (default: same as general verbose level)." << endl;
  cout << "  -vb-scan             : Give the verbose level specific to the scanner (default: same as general verbose level)." << endl;
  cout << "  -vb-data             : Give the verbose level specific to the data and image management (default: same as general verbose level)." << endl;
  cout << endl;
  cout << "  -conf                : Give the path to the CASToR config directory (default: located through the CASTOR_CONFIG environment variable)." << endl;
  cout << endl;
  cout << "  -help-misc           : Print out this help." << endl;
  cout << endl;
}





// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  Main program
*/
int main(int argc, char** argv)
{
  // ============================================================================================================
  // MPI management
  // ============================================================================================================
  int mpi_rank = 0;
  //int mpi_size = 1; 
  #ifdef CASTOR_MPI
  int mpi_size = 1;
  MPI_Init(&argc, &argv );
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  #endif

  // No argument, then show help
  if (argc==1) ShowHelp(0,mpi_rank);
  if (mpi_rank==0) cout << endl; // TODO : Most stuff is handled with processor 1 at first

  // ============================================================================================================
  // Parameterized variables with their default values
  // ============================================================================================================

  // Initialization of the voxel and field-of-view values in each spatial dimensions.
  // Read from the interfile header of the image
  int nb_VoxX=-1, nb_VoxY=-1, nb_VoxZ=-1;
  FLTNB fov_SizeX=-1, fov_SizeY=-1, fov_SizeZ=-1;
  FLTNB vox_SizeX=-1, vox_SizeY=-1, vox_SizeZ=-1;
  // Initialization offset values of the field-of-view in each spatial dimensions
  FLTNB offsetX = 0, offsetY = 0, offsetZ = 0;

  // Path to an image used as initialization. no default
  string path_to_initial_img = "";
  // Path to an image used for the attenuation. default : uniform image
  string path_to_attenuation_img = "";
  // Scanner name provided by the user
  string scanner_name = "";
  
  // DataFile mode
  int datafile_mode = MODE_HISTOGRAM;
  int datafile_type = TYPE_PET;

  // Frame list descriptor
  string frame_list = "";
  FLTNB acquisition_duration_sec = 0.;
  
  // Number of resp gates
  int nb_resp_gates = 1;
  // Number of card gates
  int nb_card_gates = 1;
   
  // Number of beds position in the data. default : 1
  int nb_beds = 1;

  // Maximum axial difference between the crystals in a LOR (negative value means no restriction)
  FLTNB max_axial_diff_mm = -1.;

  // Output directory name.
  string path_dout = "";
  // Or root name
  string path_fout = "";

  // Write scanner LUT generated by a geom file on disk
  bool save_LUT_flag = false;
  
  // String gathering the name of the used noise model (default : none) 
  string noise_model = "";
  
  // String gathering the name of the used projector for forward/backward projection (default : Siddon) 
  string options_projector = "incrementalSiddon";

  // Default projector computation strategy
  int projector_computation_strategy = FIXED_LIST_COMPUTATION_STRATEGY;
  
  // String vector gathering the name of the image convolvers with specific options (default: no image convolver)
  vector<string> options_image_convolver;

  // Using GPU (flag) ->NOTE : default : only CPU
  bool discard_nil_events_flag = false;
  
  // Using GPU (flag) ->NOTE : default : only CPU
  bool gpu_flag = false;
  // Size of the buffer used to read the data file (0 by default for projection as it is not used, but should be initialized to use the DataFile objects)
  int data_file_percentage_load = 0;

  // Precision for output number display
  string onb_prec = "s,0";
  
  // Verbose level
  int verbose_general = 1, 
      verbose_algo = -1, 
      verbose_proj = -1,
      verbose_conv = -1,
      verbose_scan = -1,
      verbose_data = -1;

  // Path to config directory
  string path_to_config_dir = "";
  
  // Initial seed for random number generator
  int64_t seed_RNG = -1;
  // Number of threads
  string nb_threads = "1";

  // SPECT specific parameters
  
  // SPECT number of transaxial bins
  uint16_t SPECT_nb_bins[2] = {0, 0};
  // SPECT Number of projection angles
  uint32_t SPECT_nb_projection_angles = 0;
  // SPECT projections angles
  FLTNB* SPECT_projection_angles = NULL;
  // SPECT first angle and step angle for automatic initialization
  FLTNB SPECT_first_angle= -1., SPECT_step_angle=-1.;
  // SPECT distance from center of rotation to detectors
  FLTNB SPECT_radius = -1.;
  // SPECT head rotation orientation
  string SPECT_rot = "CW";
    
  // ============================================================================================================
  // Read command-line parameters
  // ============================================================================================================

  for (int i=1; i<argc; i++)
  {
    string option = (string)argv[i];

    if (option=="-h" || option=="--help" || option=="-help") ShowHelp(0,mpi_rank);

    // Specific help for integrated scanners (managed by scanner manager)
    else if (option=="-help-scan")
    {
      if(sScannerManager::GetInstance()->ShowScannersDescription())
        Cerr("***** castor-proj() -> An error occurred when trying to output the available scanners from the scanner repository !'" << endl;);
      Exit(EXIT_SUCCESS);
    }
    

    // Specific help for projector settings (managed by projector children)
    else if (option=="-help-projm")
    {
      // Call specific showHelp function from vProjector children
      sAddonManager::GetInstance()->ShowHelpProjector();
      // Call the static showHelp function from vProjector
      vProjector::ShowCommonHelp();
      Exit(EXIT_SUCCESS);
    }
    
    // Specific help for image convolver settings (managed by image convolver children)
    else if (option=="-help-conv")
    {
      // Call specific showHelp function from vImageConvolver children
      sAddonManager::GetInstance()->ShowHelpImageConvolver();
      // Call the static showHelp function from oImageConvolverManager
      // Call the static showHelp function from oImageConvolverManager
      oImageConvolverManager::ShowCommonHelp();
      Exit(EXIT_SUCCESS);
    }
    
    // Specific help for input settings (managed by main)
    else if (option=="-help-in")
    {
      ShowHelpInput();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for output settings (managed by main)
    else if (option=="-help-out")
    {
      ShowHelpOutput();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for projector settings (managed by main)
    else if (option=="-help-proj")
    {
      ShowHelpProj();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for computation settings (managed by main)
    else if (option=="-help-comp")
    {
      ShowHelpComputation();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for dynamic settings (managed by main)
    else if (option=="-help-dynamic")
    {
      ShowHelpDynamic();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for PET projection settings (managed by main)
    else if (option=="-help-pet")
    {
      ShowHelpPET();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for SPECT projection settings (managed by main)
    else if (option=="-help-spect")
    {
      ShowHelpSPECT();
      Exit(EXIT_SUCCESS);
    }
    // Specific help for miscellaneous settings (managed by main)
    else if (option=="-help-misc")
    {
      ShowHelpMiscellaneous();
      Exit(EXIT_SUCCESS);
    }
    
    else if (option=="-off")
    {
      if (i>=argc-1)
      {
       Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
       Exit(EXIT_FAILURE);
      }
      else
      {
       FLTNB input[3];
       if (ReadStringOption(argv[i+1], input, 3, ",", option))
       {
         Cerr("***** castor-proj() -> Invalid argument " << argv[i+1] << " for option " << option << " !" << endl);
         Exit(EXIT_FAILURE);
       }
        offsetX = input[0];
        offsetY = input[1];
        offsetZ = input[2];
      }
      i++;
    }
    
    else if (option=="-frm")
    {
    if (i>=argc-1)
    {
      Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
      Exit(EXIT_FAILURE);
    }
    frame_list = (string)argv[i+1];
    i++;
    }
    
    // Scanner name
    else if (option=="-sc")
    {
    if (i>=argc-1)
    {
     Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
     Exit(EXIT_FAILURE);
    }
    scanner_name = argv[i+1];
    i++;
    }
    
    // Image for the initialisation of the algorithm : What should we do in case of multiple frames ? 
    else if (option=="-img")
    {
    if (i>=argc-1)
    {
     Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
     Exit(EXIT_FAILURE);
    }
    path_to_initial_img = argv[i+1];
    i++;
    }
    
    
    else if (option=="-fileType")
    {
      if (i>=argc-1)
      {
       Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
       Exit(EXIT_FAILURE);
      }
      
      string input_str = argv[i+1];
      
      if(!input_str.compare("TYPE_PET") )
        datafile_type = TYPE_PET;
      else if(!input_str.compare("TYPE_SPECT") )
        datafile_type = TYPE_SPECT;
      else if(!input_str.compare("TYPE_CT") )
        datafile_type = TYPE_CT;
      else // check for number
      {
        if(ConvertFromString(input_str , &datafile_type))
        {
          Cerr("***** castor-proj() -> Exception when trying to read provided entry '" << input_str << " for option: " << option << endl);
          Exit(EXIT_FAILURE);
        }
      }
  
      i++;
    }
    
    // Image for the attenuation
    else if (option=="-atn")
    {
    if (i>=argc-1)
    {
     Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
     Exit(EXIT_FAILURE);
    }
    path_to_attenuation_img = argv[i+1];
    i++;
    }
    
    // Number of beds
    else if (option=="-bed")
    {
    if (i>=argc-1)
    {
     Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
     Exit(EXIT_FAILURE);
    }
    nb_beds = atoi(argv[i+1]);
    i++;
    }
    
    // maximum axial difference between two crystals in a lor
    else if (option=="-PET_maxAxialDiffmm")
    {
    if (i>=argc-1)
    {
     Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
     Exit(EXIT_FAILURE);
    }
    max_axial_diff_mm = atoi(argv[i+1]);
    i++;
    }
    
    // Name of the output directory
    else if (option=="-dout") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_dout = argv[i+1];
      i++;
    }
    // Base name of the output files
    else if (option=="-fout") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_fout = argv[i+1];
      i++;
    }
    // Flag to say that we want to save time basis functions too
    else if (option=="-olut")
    {
      save_LUT_flag = true;
    }
    // Output number precision
    else if (option=="-onbp")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      onb_prec = argv[i+1];
      i++;
    }
    else if (option=="-noise")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      noise_model = argv[i+1];
      i++;
    }
    
    else if (option=="-proj")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      options_projector = (string)argv[i+1];
      i++;
    }

    // Projection line computation strategy
    else if (option=="-proj-comp")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      projector_computation_strategy = atoi(argv[i+1]);
      i++;
    }
    
    // Image convolver settings
    else if (option=="-conv")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      string convolver = (string)argv[i+1];
      convolver.append("::forward"); // convolution only allowed before projection
      options_image_convolver.push_back(convolver);
      i++;
    }
    
    // Image convolver settings
    else if (option=="-no0")
    {
      discard_nil_events_flag = true;
    }
    
    // --- Verbosities ---
    
    // General verbosity level
    else if (option=="-vb")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_general))
      {
        Cerr("***** castor-proj() -> Exception when trying to read provided verbosity level '" << verbose_general << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    
    // Algorithm verbosity level
    else if (option=="-vb-algo")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_algo))
      {
        Cerr("***** castor-proj() -> Exception when trying to read provided verbosity level '" << verbose_algo << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }

    // Projector verbosity level
    else if (option=="-vb-proj")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_proj))
      {
        Cerr("***** castor-proj() -> Exception when trying to read provided verbosity level '" << verbose_proj << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }

    // Image convolver verbosity level
    else if (option=="-vb-conv")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_conv))
      {
        Cerr("***** castor-proj() -> Exception when trying to read provided verbosity level '" << verbose_conv << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    
    // Scanner verbosity level
    else if (option=="-vb-scan")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_scan))
      {
        Cerr("***** castor-proj() -> Exception when trying to read provided verbosity level '" << verbose_scan << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    
    // Data and image verbosity level
    else if (option=="-vb-data")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose_data))
      {
        Cerr("***** castor-proj() -> Exception when trying to read provided verbosity level '" << verbose_data << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
        
    
    // RNG seed
    else if (option=="-rng")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if(ConvertFromString(argv[i+1], &seed_RNG))
      {
        Cerr("***** castor-proj() -> Exception when trying to read provided number '" << seed_RNG << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }

    // Path to config directory
    else if (option=="-conf")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_to_config_dir = (string)argv[i+1];
      i++;
    }
    
    #ifdef CASTOR_GPU
    else if (option=="-gpu")
    {
      gpu_flag = 1;
    }
    #endif
    #ifdef CASTOR_OMP
    else if (option=="-th")
    {
    if (i>=argc-1)
    {
     Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
     Exit(EXIT_FAILURE);
    }
    nb_threads = (string)argv[i+1];
    i++;
    }
    #endif
    
    else if (option=="-load")
    {
     if (i>=argc-1)
     {
       Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
       Exit(EXIT_FAILURE);
     }
     data_file_percentage_load = atoi(argv[i+1]);
     if (data_file_percentage_load>100 || data_file_percentage_load<0)
     {
       Cerr("***** castor-proj() -> Incorrect initialization of the size data buffer" << endl);
       Cerr("                       Number provided: " << argv[i+1] << " is not in the expected [0;100] interval" << endl);
       Exit(EXIT_FAILURE);
     }
     i++;
    }
    else if (option=="-SPECT_bins")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      
      if(ReadStringOption(argv[i+1], SPECT_nb_bins, 2, ",", option))
      {
        Cerr("***** castor-proj() -> Invalid argument " << argv[i+1] << " for option " << option << " !" << endl);
        Cerr("***** castor-proj() -> Expected two uint16 parameters separated by a comma" << endl);
        Exit(EXIT_FAILURE);
      }
      
      i++;
    }

    // Custom initialization of SPECT projection angles using a number of projections, and an equivalent number of pre-set projection angles
    else if (option=="-SPECT_c_ang")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      
      // parsing
      string input_str = argv[i+1];
      size_t colon = input_str.find_first_of(":");
      
      if (colon==string::npos)
      {
        Cerr("***** castor-proj() -> Parsing error : colon missing in option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
      {
        if(ConvertFromString(input_str.substr(0,colon), &SPECT_nb_projection_angles))
        {
          Cerr("***** castor-proj() -> An error occurred while trying to read the number of SPECT projection values " << input_str.substr(0,colon) << " for option " << option << " !" << endl);
          Exit(EXIT_FAILURE);
        }   
        
        SPECT_projection_angles = new FLTNB[SPECT_nb_projection_angles];
        
        if(ReadStringOption(input_str.substr(colon+1), SPECT_projection_angles, SPECT_nb_projection_angles, ",", option))
        {
          Cerr("***** castor-proj() -> Invalid argument " << argv[i+1] << " for option " << option << " !" << endl);
          Exit(EXIT_FAILURE);
        }
      }
      i++;
    }

    // Generic initialization of SPECT projection angles using a first angle, an angle step, and a number of projection
    else if (option=="-SPECT_ang")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }

      // parsing
      string input_str = argv[i+1];
      size_t colon = input_str.find_first_of(":");

      if (colon==string::npos)
      {
        Cerr("***** castor-proj() -> Parsing error : colon missing in option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      else
      {
        if(ConvertFromString(input_str.substr(0,colon), &SPECT_nb_projection_angles))
        {
          Cerr("***** castor-proj() -> An error occurred while trying to read the number of SPECT projection values " << input_str.substr(0,colon) << " for option " << option << " !" << endl);
          Exit(EXIT_FAILURE);
        }   
        
        FLTNB input[2];
        if(ReadStringOption(input_str.substr(colon+1), input, 2, ",", option))
        {
          Cerr("***** castor-proj() -> Invalid argument " << argv[i+1] << " for option " << option << " !" << endl);
          Exit(EXIT_FAILURE);
        }
        else
        {
          SPECT_first_angle = input[0];
          SPECT_step_angle = input[1];
        }
      }
      i++;
    }

    else if (option=="-SPECT_radius")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      
      else
      {
        if(ConvertFromString(argv[i+1], &SPECT_radius) )
        {
          Cerr("***** castor-proj() -> An error occurred while trying to read the SPECT global radius " << argv[i+1] << " for option " << option << " !" << endl);
          Exit(EXIT_FAILURE);
        }   
      }
      i++;
    }

    else if (option=="-SPECT_rot")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-proj() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      
      else
      {
        SPECT_rot = argv[i+1];
        
        if(SPECT_rot != "CW" && SPECT_rot != "CCW")
        {
          Cerr("***** castor-proj() -> '" << SPECT_rot <<"' is unknown argument for option " << option << ".");
          Cerr("                       SPECT rotation direction must be either 'CW' (clockwise) or 'CCW' (counter-clockwise) !" << endl);
          Exit(EXIT_FAILURE);
        }   
      }
      i++;
    }

    else
    {
      if (mpi_rank==0) Cerr("***** castor-proj() -> Unknown option '" << option << "' !" << endl);
      Exit(EXIT_FAILURE);
    }

  }

  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif


  // Affect specific verbose levels if not set
  if (verbose_algo==-1) verbose_algo = verbose_general;
  if (verbose_proj==-1) verbose_proj = verbose_general;
  if (verbose_conv==-1) verbose_conv = verbose_general;
  if (verbose_scan==-1) verbose_scan = verbose_general;
  if (verbose_data==-1) verbose_data = verbose_general;



  // ============================================================================================================
  // Mandatory checks: the idea is that the most likely problems are checked here, so that we do not apply a
  //                   too much secured approach in the rest of the classes in order to not overload the code
  //                   with checks everywhere, which could also possibly slows down its execution.
  // ============================================================================================================

  // Basic initialization checks (minimal initializations mandatory for the next steps)

  // scanner
  if (scanner_name.empty())
  {
    Cerr("***** castor-proj() -> Please provide a scanner name (format: Modality_Constuctor_Model) :" << endl);
    Cerr("                       -sc  scanner_alias: Give the scanner model." << endl);
    Cerr("                       It must correspond to the name of one of the file in the scanner repository (without file extension)" << endl);
    Exit(EXIT_FAILURE);
  }
  else
  {
    if (verbose_general>=2) Cout(" selected scanner name: " << scanner_name << endl);
  }
 
  // Output files
  if (path_fout.empty() && path_dout.empty())
  {
    Cerr("***** castor-proj() -> Please provide an output option for output files (-fout or -dout) !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Check that only one option has been provided
  if (!path_fout.empty() && !path_dout.empty())
  {
    Cerr("***** castor-proj() -> Please provide either output option -fout or -dout but not both !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Image source
  if (path_to_initial_img.empty())
  {
    Cerr("*****  castor-proj() -> -img path_to_img.hdr : Give an interfile header of the image to project" << endl);
    Exit(EXIT_FAILURE);
  }
  else
  {
    if (verbose_general>=2) Cout(" selected input image path: " << path_to_initial_img  << endl);
  }

  Cout(endl);




  // ============================================================================================================
  // Singletons initialization: create here all the needed singletons.
  // ============================================================================================================

  // Get user endianness (for interfile)
  GetUserEndianness();

  // OutputManager
  sOutputManager* p_outputManager; 
  p_outputManager = sOutputManager::GetInstance();
  p_outputManager->SetVerbose(verbose_general);
  
  // Set datafile name(s) in order to be able to recover them from interfile classes
  vector <string> vpath_to_initial_img;
  vpath_to_initial_img.push_back(path_to_initial_img);
  p_outputManager->SetDataFileName(vpath_to_initial_img);
  // Set output number precision
  p_outputManager->SetOutNbPrec(onb_prec);
  
  // Set path to the config directory
  if (p_outputManager->CheckConfigDir(path_to_config_dir))
  {
    Cerr("***** castor-proj() -> A problem occurred while checking for the config directory path !" << endl);
    Exit(EXIT_FAILURE);
  }
  
  // Initialize output directory and base name
  if (p_outputManager->InitOutputDirectory(path_fout, path_dout))
  {
    Cerr("***** castor-proj() -> A problem occurred while initializing output directory !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Log command line
  if (p_outputManager->LogCommandLine(argc,argv))
  {
    Cerr("***** castor-proj() -> A problem occurred while logging command line arguments !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ScannerManager
  sScannerManager* p_scannermanager = sScannerManager::GetInstance();  
  p_scannermanager->SetVerbose(verbose_scan);
  p_scannermanager->SetSaveLUTFlag(save_LUT_flag);
  
  
  
  // ============================================================================================================
  // Main part of the program: create all elements needed for the algorithm, create the algorithm, initialize it
  //                           and launch it.
  // ============================================================================================================
  

  // ============================================================================================================
  // Generate system geometry
  // ============================================================================================================
  
  scanner_name = (scanner_name.find(OS_SEP)) ? 
                  scanner_name.substr(scanner_name.find_last_of(OS_SEP)+1) :
                  scanner_name;
                  
  if(p_scannermanager->FindScannerSystem(scanner_name) )
  {
    Cerr("***** castor-proj() -> A problem occurred while searching for scanner system !" << endl);
    Exit(EXIT_FAILURE);
  }
  
  if(p_scannermanager->BuildScannerObject() ) //TODO verifier si correct d'envoyer uniquement le chemin du premier DataFile (pb si dataFile ne proviennent pas du même scan ou de la même acquisition)
  {
    Cerr("***** castor-proj() -> A problem occurred during scanner object construction !" << endl);
    Exit(EXIT_FAILURE);
  }
  
  if(p_scannermanager->InstantiateScanner() )
  {
    Cerr("***** castor-proj() -> A problem occurred while creating Scanner object !" << endl);
    Exit(EXIT_FAILURE);
  } 
 

  //Set the specific mandatory geometry information of the modality using the SetXXXSpecificParameters functions of the scannerManager
  if (p_scannermanager->GetScannerType() == SCANNER_PET)
  {
    if(p_scannermanager->PROJ_SetPETSpecificParameters(max_axial_diff_mm) )
    {
      Cerr("***** castor-proj() -> A problem occurred while setting PET geometric parameters !" << endl);
      Exit(EXIT_FAILURE);
    }
  }
  else if (p_scannermanager->GetScannerType() == SCANNER_SPECT_CONVERGENT || p_scannermanager->GetScannerType() == SCANNER_SPECT_PINHOLE )
  {
    if(p_scannermanager->PROJ_SetSPECTSpecificParameters(SPECT_nb_bins, 
                                                         SPECT_nb_projection_angles, 
                                                         SPECT_first_angle, 
                                                         SPECT_step_angle, 
                                                         SPECT_projection_angles, 
                                                         SPECT_radius,
                                                         SPECT_rot) )
    {
      Cerr("***** castor-proj() -> A problem occurred while setting SPECT geometric parameters !" << endl);
      Exit(EXIT_FAILURE);
    }
  } 


  if(p_scannermanager->BuildLUT() )
  {
    Cerr("***** castor-proj() -> A problem occurred while generating/reading the LUT !" << endl);
    Exit(EXIT_FAILURE);
  }
  

  // Check the scanner manager parameters and initialize the scanner
  if (p_scannermanager->CheckParameters())
  {
    Cerr("***** castor-proj() -> A problem occurred while checking scanner manager parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_scannermanager->Initialize())
  {
    Cerr("***** castor-proj() -> A problem occurred while initializing scanner !" << endl);
    Exit(EXIT_FAILURE);
  }


  // ============================================================================================================
  // Recover image dimensions informations from the interfile header of the image to project
  // ============================================================================================================
  
  Intf_fields IF;
  IntfKeyInitFields(&IF);
    
  if(IntfReadHeader(path_to_initial_img, &IF, verbose_data) )
  {
    Cerr("***** castor-proj() -> An error occurred while trying to read the interfile header of file reading file " << path_to_initial_img << " !" << endl);  
    Exit(EXIT_FAILURE);
  }

  if(verbose_general >= 3) IntfKeyPrintFields(IF);

  nb_VoxX = nb_VoxX<0 ? IF.mtx_size[0] : nb_VoxX ;
  nb_VoxY = nb_VoxY<0 ? IF.mtx_size[1] : nb_VoxY ;
  nb_VoxZ = nb_VoxZ<0 ? IF.mtx_size[2] : nb_VoxZ ;

  // TODO : nb_beds ?
  // If voxel sizes not provided/found, set them to the default value (1mm)
  vox_SizeX = IF.vox_size[0]<0 ? 1. : IF.vox_size[0];
  vox_SizeY = IF.vox_size[1]<0 ? 1. : IF.vox_size[1];
  vox_SizeZ = IF.vox_size[2]<0 ? 1. : IF.vox_size[2];

  // Add image offsets to user-defined offsets
  offsetX -= IF.vox_offset[0];
  offsetY -= IF.vox_offset[1];
  offsetZ -= IF.vox_offset[2];
  
   // Get the frame duration from the image in the case it is not provided by the user
  if (frame_list == "")
  {
    vector <string> image_filenames;
    if ( IntfIsMHD(path_to_initial_img, image_filenames) )
    {
      Cerr("***** castor-proj() -> Error : while trying to read the interfile header '"<< path_to_initial_img << "' !" << endl);
      Exit(EXIT_FAILURE);
    }

    // Image durations have been provided 
    if(IF.image_duration.size()>0) 
    {
      // Several files (serie of 3D images)
      // TODO : fill IF.study_duration with image durations if not init ?
      if(image_filenames.size() > 1) 
      {
        // Check nb of frames and parse framing
        if(IF.image_duration.size() != (uint32_t)IF.nb_time_frames *
                                                 IF.nb_resp_gates  *
                                                 IF.nb_card_gates )
        {
          Cerr("***** castor-proj() -> Interfile reading error of the input image :"<<endl);
          Cerr("      The number of provided image duration:  '"<< IF.image_duration.size() 
            << "'     does not match the actual number of time frames * respiratory gates * cardiac gates ) '"<< IF.nb_time_frames *
                                                                                                                 IF.nb_resp_gates  *
                                                                                                                 IF.nb_card_gates <<"' !" << endl);
          Exit(EXIT_FAILURE);
        }
            
        // If image files are splitted over a serie of 3D image files, we will have image durations for each image, including gates
        // As we need to recover the duration for each gate once, we skip gates if the data contains any
        for(int fr=0 ; fr<IF.nb_time_frames ; fr++)
        {
          // Manage start time of the frame, add dead frame if delay
          /*
          if(fr==0 && IF.image_start_time.at(0) > 0 )
          {
            stringstream ss;
            ss << IF.image_start_time.at(fr);
            //frame_list.append(ss.str()).append("df,");
            frame_list.append(ss.str());
          }
          // If the next frame starts after the previous frame, add frame duration
          else if(IF.image_start_time.at(fr) > ( IF.image_start_time.at(fr-1) + IF.image_duration.at(fr-1) )  )
          {
            stringstream ss;
            ss << IF.image_start_time.at(fr);
            frame_list.append(ss.str()).append("df,");
          }
          */
          // First frame : just add frame start
          if( fr==0  )
          {
            stringstream ss;
            ss << IF.image_start_time.at(fr);
            frame_list.append(ss.str());
          }
          // If current frame starts after the previous frame
          // -> add ':' and frame duration of the previous frame
          // -> add start time of the actual frame
          else if(IF.image_start_time.at(fr) > ( IF.image_start_time.at(fr-1) + IF.image_duration.at(fr-1) )  )
          {
            stringstream ss;
            ss << IF.image_duration.at(fr-1);
            frame_list.append(":").append(ss.str());
            ss.str("");
            ss << IF.image_start_time.at(fr);
            frame_list.append(",").append(ss.str());
          }
          // Just add frame start otherwise
          else
          {
            stringstream ss;
            ss << IF.image_start_time.at(fr);
            frame_list.append(",").append(ss.str());
          }
          

          
          // TODO TM: Check if the following is still required
          
          // Recover image duration and dead frames (if any) in the frame list
          /*
          stringstream ss_id, ss_df;
          // Go to the next frame (skip gate if any)
          ss_id << IF.image_duration.at(fr * IF.nb_resp_gates* IF.nb_card_gates);
          string frame_duration_str = ss_id.str();
          
          // Check if there is any dead frame (pause between frame groups)
          if(IF.frame_group_pause.size()>0)
          {
            FLTNB dead_frame_sec = IF.frame_group_pause.at(fr * IF.nb_resp_gates* IF.nb_card_gates);
            if(dead_frame_sec > 0.)
            {
              ss_df << dead_frame_sec;
              frame_duration_str.append(",").append(ss_df.str()).append("df");
            }
          }
          */
          //frame_list.append(frame_duration_str);
                    
          /*
          // add comma if not last frame
          if(fr+1 < IF.nb_time_frames)
            frame_list.append(",");
          */
            
          // If last frame : add frame duration
          if(fr+1 >= IF.nb_time_frames)
          {
            stringstream ss;
            ss << IF.image_duration.at(fr);
            frame_list.append(":").append(ss.str());
          }

        }
      }
      else
      {
        if(IF.image_duration.size() != IF.nb_time_frames)
          {
            Cerr("***** castor-proj() -> Interfile reading error of the input image :"<<endl);
            Cerr("      The number of provided image duration:  '"<< IF.image_duration.size() 
              << "'     does not match the actual number of time frames * respiratory gates * cardiac gates ) '"<< IF.nb_time_frames <<"' !" << endl);
            Exit(EXIT_FAILURE);
          }  
          
        for(int fr=0 ; fr<IF.nb_time_frames ; fr++)
        {
          /*
          if(fr==0 && IF.image_start_time.at(0) > 0 )
          {
            stringstream ss;
            ss << IF.image_start_time.at(fr);
            frame_list.append(ss.str()).append("df,");
          }
          else if(IF.image_start_time.at(fr) > ( IF.image_start_time.at(fr-1) + IF.image_duration.at(fr-1) )  )
          {
            stringstream ss;
            ss << IF.image_start_time.at(fr);
            frame_list.append(ss.str()).append("df,");
          }
          
          stringstream ss;
          ss << IF.image_duration.at(fr);
          frame_list.append(ss.str());
          // add comma if not last frame
          if(fr+1 < IF.nb_time_frames)
            frame_list.append(",");
          */
          
          
          // First frame : just add frame start
          if( fr==0  )
          {
            stringstream ss;
            ss << IF.image_start_time.at(fr);
            frame_list.append(ss.str());
          }
          // If current frame starts after the previous frame
          // -> add ':' and frame duration of the previous frame
          // -> add start time of the actual frame
          else if(IF.image_start_time.at(fr) > ( IF.image_start_time.at(fr-1) + IF.image_duration.at(fr-1) )  )
          {
            stringstream ss;
            ss << IF.image_duration.at(fr-1);
            frame_list.append(":").append(ss.str());
            ss.str("");
            ss << IF.image_start_time.at(fr);
            frame_list.append(",").append(ss.str());
          }
          // Just add frame start otherwise
          else
          {
            stringstream ss;
            ss << IF.image_start_time.at(fr);
            frame_list.append(",").append(ss.str());
          }
          
          // If last frame : add frame duration
          if(fr+1 >= IF.nb_time_frames)
          {
            stringstream ss;
            ss << IF.image_duration.at(fr);
            frame_list.append(":").append(ss.str());
          }
        } // end of loop on frames
      }  // end of condition on dynamic image type (one image or several)
    } // end of condition on image duration tag (present or not)
    // Image durations have not been provided,
    // use study duration which is initialized with 1s (default)
    else
    {
      stringstream ss;
      ss << IF.study_duration;
      frame_list.append("0:").append(ss.str());
    }
  }


  // Check nb gating
  (IF.nb_resp_gates >1) ? nb_resp_gates = IF.nb_resp_gates : 1 ;
  (IF.nb_card_gates >1) ? nb_card_gates = IF.nb_card_gates : 1 ;

  //TODO : offsets ?

  if(nb_VoxX<0 || nb_VoxY<0  || nb_VoxZ<0)
  {
    ReadDataASCIIFile(p_scannermanager->GetPathToScannerFile(), "voxels number transaxial", &nb_VoxX, 1, KEYWORD_MANDATORY);
    ReadDataASCIIFile(p_scannermanager->GetPathToScannerFile(), "voxels number transaxial", &nb_VoxY, 1, KEYWORD_MANDATORY);
    ReadDataASCIIFile(p_scannermanager->GetPathToScannerFile(), "voxels number axial", &nb_VoxZ, 1, KEYWORD_MANDATORY);
  }

  if ( (fov_SizeX<0 || fov_SizeY<0 || fov_SizeZ<0) && (vox_SizeX<0 || vox_SizeY<0 || vox_SizeZ<0) )
  {
    ReadDataASCIIFile(p_scannermanager->GetPathToScannerFile(), "field of view transaxial", &fov_SizeX, 1, KEYWORD_MANDATORY);
    ReadDataASCIIFile(p_scannermanager->GetPathToScannerFile(), "field of view transaxial", &fov_SizeY, 1, KEYWORD_MANDATORY);
    ReadDataASCIIFile(p_scannermanager->GetPathToScannerFile(), "field of view axial", &fov_SizeZ, 1, KEYWORD_MANDATORY);
  }
/*
  else if ( (fov_SizeX<0 || fov_SizeY<0 || fov_SizeZ<0) && (vox_SizeX>=0 || vox_SizeY>=0 || vox_SizeZ>=0))
  {
    fov_SizeX = vox_SizeX*(FLTNB)nb_VoxX;
    fov_SizeY = vox_SizeY*(FLTNB)nb_VoxY;
    fov_SizeZ = vox_SizeZ*(FLTNB)nb_VoxZ;
    cout << vox_SizeX << " " << fov_SizeX << endl;
  }
  else
  {
    vox_SizeX = fov_SizeX/(FLTNB)nb_VoxX;
    vox_SizeY = fov_SizeY/(FLTNB)nb_VoxY;
    vox_SizeZ = fov_SizeZ/(FLTNB)nb_VoxZ;
  }
*/
  if (verbose_general>=2)
  {
    Cout(" Image dimensions for analytical projections: " << endl);
    Cout(" voxels number (X,Y,Z): " << nb_VoxX << "," << nb_VoxY << "," << nb_VoxZ << endl);
    Cout(" voxels size (X,Y,Z): "     << vox_SizeX << "," << vox_SizeY << "," << vox_SizeZ << endl);
    Cout(" FOV (X,Y,Z): "     << fov_SizeX << "," << fov_SizeY << "," << fov_SizeZ << endl);
  }

  // ============================================================================================================
  // oImageDimensionsAndQuantification settings
  // ============================================================================================================
  
  // Create oImageDimensionsAndQuantification
  oImageDimensionsAndQuantification* p_ImageDimensionsAndQuantification = new oImageDimensionsAndQuantification(); 
  p_ImageDimensionsAndQuantification->SetNbVoxX(nb_VoxX);
  p_ImageDimensionsAndQuantification->SetNbVoxY(nb_VoxY);
  p_ImageDimensionsAndQuantification->SetNbVoxZ(nb_VoxZ);
  p_ImageDimensionsAndQuantification->SetNbThreads(nb_threads);
  p_ImageDimensionsAndQuantification->SetNbBeds(nb_beds);
  p_ImageDimensionsAndQuantification->SetVoxSizeX(vox_SizeX);
  p_ImageDimensionsAndQuantification->SetVoxSizeY(vox_SizeY);
  p_ImageDimensionsAndQuantification->SetVoxSizeZ(vox_SizeZ);
  p_ImageDimensionsAndQuantification->SetFOVSizeX(fov_SizeX);
  p_ImageDimensionsAndQuantification->SetFOVSizeY(fov_SizeY);
  p_ImageDimensionsAndQuantification->SetFOVSizeZ(fov_SizeZ);
  p_ImageDimensionsAndQuantification->SetFOVOutMasking(0.,  0);
  p_ImageDimensionsAndQuantification->SetOffsetX(offsetX);
  p_ImageDimensionsAndQuantification->SetOffsetY(offsetY);
  p_ImageDimensionsAndQuantification->SetOffsetZ(offsetZ);
  p_ImageDimensionsAndQuantification->SetVerbose(verbose_data);
  p_ImageDimensionsAndQuantification->SetNbRespGates(nb_resp_gates);
  p_ImageDimensionsAndQuantification->SetNbCardGates(nb_card_gates);
  p_ImageDimensionsAndQuantification->SetFrames(frame_list);  
  if (p_ImageDimensionsAndQuantification->CheckParameters())
  {
    Cerr("***** castor-proj() -> A problem occurred while checking image dimensions parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ImageDimensionsAndQuantification->Initialize())
  {
    Cerr("***** castor-proj() -> A problem occurred while initializing image dimensions !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Initialization of DynamicDataManager class, related 4D data splitting management 
  if (p_ImageDimensionsAndQuantification->InitDynamicData( "", 0, 0, 0, nb_resp_gates, nb_card_gates ) )
  {
    Cerr("***** castor-proj() -> A problem occurred while initializing Dynamic data manager's class !" << endl);
    Exit(EXIT_FAILURE);
  }
  
  
  // ----------------------------------------------------------------------------------------
  //  Random Number Generator initialization: (we first required to know the number of threads to use from p_ImageDimensionsAndQuantification)
  // ----------------------------------------------------------------------------------------
  
  sRandomNumberGenerator* p_RNG = sRandomNumberGenerator::GetInstance(); 
  p_RNG->SetVerbose(verbose_general);
  // Use a user-provided seed to initialize the RNG if one has been provided. Use random number otherwise.
  (seed_RNG>=0 )?
    p_RNG->Initialize(seed_RNG, p_ImageDimensionsAndQuantification->GetNbThreadsForProjection(),0):
    p_RNG->Initialize(p_ImageDimensionsAndQuantification->GetNbThreadsForProjection(),0);


  // ============================================================================================================
  // DataFile initialization
  // ============================================================================================================
  
  vDataFile** p_DataFile = new vDataFile*[nb_beds];

  if (p_scannermanager->GetScannerType() < 0)
  {
    Cerr("***** castor-proj() -> Unknown scanner type !" << endl);
    Exit(EXIT_FAILURE);
  }
  else if (p_scannermanager->GetScannerType() == SCANNER_PET)
  {
    if (datafile_type!=TYPE_PET)
    {
      Cerr("***** castor-proj() -> Only PET type events are allowed while using a PET scanner !" << endl);
      Cerr("***** castor-proj() -> (use '-fileType' to define modality (0=PET, 1=SPECT, 2=CT) )" << endl);
      Exit(EXIT_FAILURE);
    }
    for (int i=0 ; i<nb_beds ; i++)
    {
      p_DataFile[i] = new iDataFilePET();
      if(path_to_attenuation_img != "") // Enable correction factor for attenuation
        (dynamic_cast<iDataFilePET*>(p_DataFile[i]))->SetAtnCorrectionFlagOn();
    }
  }
  else if (p_scannermanager->GetScannerType() == SCANNER_SPECT_CONVERGENT)
  {
    if (datafile_type!=TYPE_SPECT)
    {
      Cerr("***** castor-proj() -> Only SPECT type events are allowed while using a SPECT scanner !" << endl);
      Cerr("***** castor-proj() -> (use '-fileType' to define modality (0=PET, 1=SPECT, 2=CT) )" << endl);
      Exit(EXIT_FAILURE);
    }
    for (int i=0 ; i<nb_beds ; i++)
    {
      p_DataFile[i] = new iDataFileSPECT(); 
    }
  }
  else
  {
    Cerr("***** castor-proj() -> Unsupported scanner type !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Compute total acquisition duration in seconds
  for (int fr=0 ; fr<p_ImageDimensionsAndQuantification->GetNbTimeFrames() ; fr++)
    acquisition_duration_sec += p_ImageDimensionsAndQuantification->GetFrameDurationInSec(0, fr);
    
  // Load raw data in memory and do other stuff if needed.
  // TODO TM : Even if the loop is there, no support for multi-bed projection as for now (case we have to project an image as in a multi-bad acquisition)
  // It would first require some informations about the FOV of each bed (get this from interfile header of image or for command line options ?)
  // In the loop, we need to provide the bed in argument of PROJ_InitFile() in order to correctly initialize the beds
  for (int bed=0 ; bed<nb_beds ; bed++)
  {
    p_DataFile[bed]->SetBedIndex(bed);
    p_DataFile[bed]->SetVerbose(verbose_data);
    p_DataFile[bed]->SetDataMode(datafile_mode);
    p_DataFile[bed]->SetDataType(datafile_type);
    p_DataFile[bed]->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
    p_DataFile[bed]->PROJ_GetScannerSpecificParameters();
    p_DataFile[bed]->PROJ_InitFile();
    p_DataFile[bed]->ComputeSizeEvent();
    p_DataFile[bed]->SetDuration(acquisition_duration_sec);
    // The Check parameters function is dedicated to reconstruction, as it checks things such as consistenties between size of the datafile and user-provided information
    if(p_DataFile[bed]->CheckParameters())
    {
      Cerr("***** castor-proj() -> A problem occurred while checking datafile parameters ! Abort." << endl);
      Exit(EXIT_FAILURE);
    }
    p_DataFile[bed]->PrepareDataFile();
    
    /*
    if(p_DataFile[bed]->SetScannerSpecificParameters())
    {
      Cerr("***** A problem occurred while setting scanner parameters from the datafile !" << endl);
      Exit(EXIT_FAILURE);
    }
    */
    
    //p_DataFile[bed]->InitDynamicData(p_ImageDimensionsAndQuantification->GetNbTimeFrames(), nb_resp_gates, nb_card_gates, "", 0, 0, bed); 

    /*
    if (p_DataFile[bed]->InitDynamicData(nb_resp_gates, nb_card_gates, "", 0, 0, 0, p_DataFile[bed]) )
    {
      Cerr("***** A problem occurred while initializing Dynamic data !" << endl);
      Exit(EXIT_FAILURE);
    }
    */
  }


  // ============================================================================================================
  // Projector manager initialization
  // ============================================================================================================

  // Create object
  oProjectorManager* p_ProjectorManager = new oProjectorManager();

  // Set all parameters
  p_ProjectorManager->SetScanner(p_scannermanager->GetScannerObject());
  p_ProjectorManager->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_ProjectorManager->SetDataFile(p_DataFile[0]);
  p_ProjectorManager->SetComputationStrategy(projector_computation_strategy);
  p_ProjectorManager->SetOptionsForward(options_projector);
  p_ProjectorManager->SetOptionsBackward(options_projector);
  p_ProjectorManager->SetOptionsCommon("");
  p_ProjectorManager->SetVerbose(verbose_proj);

  // Check parameters
  if (p_ProjectorManager->CheckParameters())
  {
    Cerr("***** castor-proj() -> A problem occurred while checking projector manager's parameters !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Initialize projector manager
  if (p_ProjectorManager->Initialize())
  {
    Cerr("***** castor-proj() -> A problem occurred while initializing projector manager !" << endl);
    Exit(EXIT_FAILURE);
  }




  // ============================================================================================================
  // Image Convolver manager initialization
  // ============================================================================================================

  // Create object
  oImageConvolverManager* p_ImageConvolverManager = new oImageConvolverManager();
  // Set all parameters
  p_ImageConvolverManager->SetVerbose(verbose_conv);
  p_ImageConvolverManager->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_ImageConvolverManager->SetOptions(options_image_convolver);
  
  // Check parameters
  if (p_ImageConvolverManager->CheckParameters())
  {
    Cerr("***** castor-proj() -> A problem occurred while checking optimizer manager's parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize optimizer manager
  if (p_ImageConvolverManager->Initialize())
  {
    Cerr("***** castor-proj() -> A problem occurred while initializing optimizer manager !" << endl);
    Exit(EXIT_FAILURE);
  }



  // ============================================================================================================
  // Image Space initialization
  // ============================================================================================================

  oImageSpace* p_ImageSpace = new oImageSpace();
  p_ImageSpace->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_ImageSpace->SetVerbose(verbose_data);
  
  
  
  // ============================================================================================================
  // Init alogorithm and proceed to analytical projection
  // ============================================================================================================

  oAnalyticProjection* p_Projection = new oAnalyticProjection(); 
  p_Projection->InitOptimizer(p_ImageDimensionsAndQuantification);
  if (p_Projection->InitNoiseModel(noise_model) )
  {
    Cerr("***** castor-proj() -> A problem occurred while initializing noise model !" << endl);
    Exit(EXIT_FAILURE);
  }

  p_Projection->SetImageDimensionsAndQuantification(p_ImageDimensionsAndQuantification);
  p_Projection->SetProjectorManager(p_ProjectorManager);
  p_Projection->SetImageConvolverManager(p_ImageConvolverManager);
  p_Projection->SetImageSpace(p_ImageSpace);
  p_Projection->SetDataFile(p_DataFile);
  p_Projection->SetGPUflag(gpu_flag);
  p_Projection->SetVerbose(verbose_algo);
  p_Projection->SetNbBeds(nb_beds);
  p_Projection->SetPathInitImage(path_to_initial_img);
  p_Projection->SetPathAtnImage(path_to_attenuation_img);
  p_Projection->SetScanner(p_scannermanager->GetScannerObject());
  p_Projection->SetNoZeroEvent(discard_nil_events_flag);

  // Launch analytical projection
  if (p_Projection->Launch())
  {
    Cerr("***** castor-proj() ->Error occurred during Projections" << endl;);
  }


  // ============================================================================================================
  // End
  // ============================================================================================================

  if(SPECT_projection_angles) delete SPECT_projection_angles;

  if (verbose_general>=5) Cout("----- Deleting CASToR objects ... -----" << endl);
  delete p_Projection;
  delete p_ImageSpace;
  delete p_ImageConvolverManager;
  delete p_ProjectorManager;

  
  // Delete objects
  for (int i=0 ; i<nb_beds ; i++) delete p_DataFile[i];
  delete[] p_DataFile;
  delete p_ImageDimensionsAndQuantification;


  if (verbose_general>=5) Cout("----- CASToR objects successfully deleted -----" << endl);
    
  // Ending
  if (verbose_general>=1) Cout(endl);
  #ifdef CASTOR_MPI
  MPI_Finalize();
  #endif
  return EXIT_SUCCESS;
}
