
/*!
  \file
  \ingroup utils
  \brief  This program is used to merge a group of images corresponding to multiple bed positions reconstructed separately.
*/

#include "gVariables.hh"
#include "gOptions.hh"
#include "sOutputManager.hh"
#include "sScannerManager.hh"
#include "oInterfileIO.hh"


// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================
//                                                        H E L P     F U N C T I O N S
// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================


/*!
  \fn      ShowHelp()
  \brief   Display main command line options for castor-mergeBedPositions
*/
void ShowHelp()
{
  // Show help
  cout << endl;
  cout << "Usage:  castor-mergeBedPositions  -img image1.hdr  -img image2.hdr  -(f/d)out output  [settings]" << endl;
  cout << endl;
  cout << "This program can be used to merge several images corresponding to multiple bed positions of the same acquisition. The slices of any image" << endl;
  cout << "are weighted during merging. By default, these weights are uniform (1 for all slices of all images). The -sens option can be used to provide" << endl;
  cout << "a sensitivity image for each image to be merged, where the voxels' value of the sensitivity images are used as weights." << endl;
  cout << "The bed relative positions are recovered from the image headers as 'bed relative position (mm)'." << endl;
  cout << "If not found, the default bed displacement between two successive bed positions from the scanner is used." << endl;
  cout << "Important note: this program assumes the provided images exactly match the axial FOV of the scanner; if bigger or smaller results will be wrong." << endl;
  cout << "Important note: this program only works with 3D images (without time dimension)." << endl;
  cout << endl;
  cout << "[Mandatory parameters]:" << endl;
  cout << "  -img imageX.hdr      : Give an image corresponding to a bed position (must give at least two images; the provided order defines the order" << endl;
  cout << "                         in which bed positions are merged)" << endl;
  cout << "  -fout name           : Give the root name for all output files (no default, alternative to -dout)" << endl;
  cout << "  -dout name           : Give the name of the output directory where all output files will be written (no default, alternative to -fout)" << endl;
  cout << endl;
  cout << "[Options]:" << endl;
  cout << "  -sens sensX.hdr      : Give the sensitivity image to weight the imageX.hdr while merging the bed positions (must provide as many sensitivity" << endl;
  cout << "                         images as images to be merged; the association of sensitivity images to the images is done with respect to the order" << endl;
  cout << "                         both images are provided)" << endl;
  cout << "  -oweight             : Also save the global weights" << endl;
  cout << "  -inv                 : Inverse the order of the provided images (and sensitivity images)" << endl;
  cout << "  -vb value            : Give the verbosity level, from 0 (no verbose) to 2 (default: 1)" << endl;
  cout << "  -conf value          : Give the path to the CASToR config directory (default: located through the CASTOR_CONFIG environment variable)." << endl;
  cout << "  -tolerance value     : Give the tolerance of bed positions with respect to the output image slices, which will be interpreted as a" << endl;
  cout << "                         percentage of the slice thickness [default: 1.]." << endl;
  cout << "  -flipZ               : Flip the output image, once merged." << endl;
  cout << "  --help,-h,-help      : Print out this help page." << endl; // managed by main
  cout << endl;
  #ifdef CASTOR_OMP
  cout << "  Compiled with OpenMP (will use as many cores as available)" << endl;
  #endif
  #ifdef BUILD_DATE
  cout << "  Build date: " << BUILD_DATE << endl;
  cout << endl;
  #endif
  #ifdef CASTOR_VERSION
  cout << "  This program is part of the CASToR release version " << CASTOR_VERSION << "." << endl;
  cout << endl;
  #endif
}

// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================
//                                                        M A I N     P R O G R A M
// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================

int main(int argc, char** argv)
{
  // ============================================================================================================
  // MPI stuff (we make all instances but the first one returning 0 directly)
  // ============================================================================================================
  int mpi_rank = 0;
  #ifdef CASTOR_MPI
  int mpi_size = 1;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_rank!=0) return 0;
  #endif

  // No argument, then show help
  if (argc==1)
  {
    ShowHelp();
    Exit(EXIT_SUCCESS);
  }

  // ============================================================================================================
  // Parameterized variables with their default values
  // ============================================================================================================

  // --------------------------------------------------------------------------------
  // Input settings
  // --------------------------------------------------------------------------------

  // Number of bed position 
  uint32_t nb_beds = 0;
  // Vector containing string pointing to the images to be merged
  vector<string> path_to_image_filename;
  // Vector containing string pointing to the sensitivity images used as weights
  vector<string> path_to_sens_filename;
  // Invert bed datafile names order
  bool invert_order_flag = false;

  // --------------------------------------------------------------------------------
  // Output settings
  // --------------------------------------------------------------------------------

  // Output directory name.
  string path_dout = "";
  // Or root name
  string path_fout = "";
  // Save global weight
  bool save_weights = false;
  // Flip output image axially
  bool flip_output_axial = false;

  // --------------------------------------------------------------------------------
  // Miscellaneous settings
  // --------------------------------------------------------------------------------

  // Tolerance of bed positions with respect to image slices, defined as a percentage of the slice thickness (default is 1%)
  FLTNB tolerance_as_percent_slice_thickness = 1.;
  // Verbose level
  int verbose = 1;
  // Path to config directory
  string path_to_config_dir = "";

  // ============================================================================================================
  // Read command-line parameters
  // ============================================================================================================

  // Must manually increment the option index when an argument is needed after an option
  for (int i=1; i<argc; i++)
  {
    // Get the option as a string
    string option = (string)argv[i];

    // --------------------------------------------------------------------------------
    // Miscellaneous settings
    // --------------------------------------------------------------------------------

    // Show help
    if (option=="-h" || option=="--help" || option=="-help")
    {
      ShowHelp();
      Exit(EXIT_SUCCESS);
    }

    // Tolerance of the bed positions with respect to the output image slices, interpreted as a percentage of the slice thickness
    else if (option=="-tolerance")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-mergeBedPositions() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &tolerance_as_percent_slice_thickness))
      {
        Cerr("***** castor-mergeBedPositions() -> Exception when trying to read provided verbosity level '" << verbose << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // General verbosity level
    else if (option=="-vb")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-mergeBedPositions() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose))
      {
        Cerr("***** castor-mergeBedPositions() -> Exception when trying to read provided verbosity level '" << verbose << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // Path to config directory
    else if (option=="-conf")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-mergeBedPositions() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_to_config_dir = (string)argv[i+1];
      i++;
    }

    // --------------------------------------------------------------------------------
    // Input settings
    // --------------------------------------------------------------------------------

    // Images
    else if (option=="-img") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-mergeBedPositions() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      string file_name = (string)argv[i+1];
      path_to_image_filename.push_back(file_name);
      nb_beds++;
      i++;
    }
    // Sensitivity images
    else if (option=="-sens")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-mergeBedPositions() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      string file_name = (string)argv[i+1];
      path_to_sens_filename.push_back(file_name);
      i++;
    }
    // Invert bed positions order
    else if (option=="-inv")
    {
      invert_order_flag = true;
    }

    // --------------------------------------------------------------------------------
    // Output settings
    // --------------------------------------------------------------------------------

    // Name of the output directory
    else if (option=="-dout") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-mergeBedPositions() -> Argument missing for option: " << option << endl);
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
        Cerr("***** castor-mergeBedPositions() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_fout = argv[i+1];
      i++;
    }
    // Save global weights
    else if (option=="-oweight")
    {
      save_weights = true;
    }
    // Flip output image axially
    else if (option=="-flipZ")
    {
      flip_output_axial = true;
    }

    // --------------------------------------------------------------------------------
    // Unknown option!
    // --------------------------------------------------------------------------------

    else
    {
      Cerr("***** castor-mergeBedPositions() -> Unknown option '" << option << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  // ============================================================================================================
  // Some checks
  // ============================================================================================================

  // Data files
  if (nb_beds < 2)
  {
    Cerr("***** castor-mergeBedPositions() -> Please provide at least two image filename !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Output files
  if (path_fout.empty() && path_dout.empty())
  {
    Cerr("***** castor-mergeBedPositions() -> Please provide an output option for output files (-fout or -dout) !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Check that only one option has been provided
  if (!path_fout.empty() && !path_dout.empty())
  {
    Cerr("***** castor-mergeBedPositions() -> Please provide either output option -fout or -dout but not both !" << endl);
    Exit(EXIT_FAILURE);
  }
  // If sensitivity is provided, then check that there are as many files as input images
  if (path_to_sens_filename.size()>0 && path_to_sens_filename.size()!=path_to_image_filename.size())
  {
    Cerr("***** castor-mergeBedPositions() -> If using sensitivity images as weights, please provide as many sensitivity images as input images !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Set maximum number of threads if compiled with OpenMP
  // ============================================================================================================

  #ifdef CASTOR_OMP
  omp_set_num_threads(omp_get_max_threads());
  #endif

  // ============================================================================================================
  // Create sOutputManager
  // ============================================================================================================
  sOutputManager* p_OutputManager = sOutputManager::GetInstance();  
  // Set verbose level
  p_OutputManager->SetVerbose(verbose);
  // Set MPI rank
  p_OutputManager->SetMPIRank(mpi_rank);
  // Set path to the config directory
  if (p_OutputManager->CheckConfigDir(path_to_config_dir))
  {
    Cerr("***** castor-mergeBedPositions() -> A problem occurred while checking for the config directory path !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize output directory and base name
  if (p_OutputManager->InitOutputDirectory(path_fout, path_dout))
  {
    Cerr("***** castor-mergeBedPositions() -> A problem occurred while initializing output directory !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Log command line
  if (p_OutputManager->LogCommandLine(argc,argv))
  {
    Cerr("***** castor-mergeBedPositions() -> A problem occurred while logging command line arguments !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Create sScannerManager
  // ============================================================================================================
  sScannerManager* p_ScannerManager = sScannerManager::GetInstance();  
  p_ScannerManager->SetVerbose(verbose);

  // ============================================================================================================
  // Read all input images
  // ============================================================================================================

  // Verbose
  if (verbose>=VERBOSE_LIGHT) Cout("castor-mergeBedPositions() -> Load input images" << endl);

  // Get user endianness (interfile I/O)
  GetUserEndianness();

  // Create interfile image fields and image pointers
  Intf_fields** p_image_fields = (Intf_fields**)malloc(nb_beds*sizeof(Intf_fields*));
  FLTNB** p_images = (FLTNB**)malloc(nb_beds*sizeof(FLTNB*));

  // Read interfile images
  for (uint32_t bed=0; bed<nb_beds; bed++)
  {
    // Create the image fields
    p_image_fields[bed] = new Intf_fields;
    // Read the image
    if (invert_order_flag) p_images[bed] = IntfLoadImageFromScratch(path_to_image_filename[nb_beds-1-bed],p_image_fields[bed],verbose);
    else p_images[bed] = IntfLoadImageFromScratch(path_to_image_filename[bed],p_image_fields[bed],verbose);
    // Check reading
    if (p_images[bed]==NULL)
    {
      Cerr("***** castor-mergeBedPositions() -> A problem occurred while reading image for bed " << bed+1 << " !" << endl);
      return 1;
    }
  }

  // Check dimensions consistency between all image fields, as well as originating system and bed displacement
  for (uint32_t bed=1; bed<nb_beds; bed++)
  {
    if (IntfCheckDimensionsConsistency(*p_image_fields[0],*p_image_fields[bed]))
    {
      Cerr("***** castor-mergeBedPositions() -> Inconsistency between image dimensions !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  // Verbose
  if (verbose>=VERBOSE_NORMAL)
  {
    Cout("  --> Image dimensions: [ " << p_image_fields[0]->mtx_size[0] << " : "
                                      << p_image_fields[0]->mtx_size[1] << " : "
                                      << p_image_fields[0]->mtx_size[2] << " ]" << endl);
    Cout("  --> Voxel sizes: [ " << p_image_fields[0]->vox_size[0] << " : "
                                 << p_image_fields[0]->vox_size[1] << " : "
                                 << p_image_fields[0]->vox_size[2] << " ]" << endl);
  }

  // Get dimensions locally
  INTNB dim_x = p_image_fields[0]->mtx_size[0];
  INTNB dim_y = p_image_fields[0]->mtx_size[1];
  INTNB dim_z = p_image_fields[0]->mtx_size[2];
  INTNB dim_xy = dim_x * dim_y;
  INTNB dim_xyz = dim_xy * dim_z;

  // ============================================================================================================
  // Check if at least one or all image headers have the bed relative position field
  // ============================================================================================================

  // Analyze if one is not zero and if all are not zero
  bool one_position_not_zero = false;
  bool all_positions_not_zero = true;
  for (uint32_t bed=0; bed<nb_beds; bed++)
  {
    // Here, we compare the bed relative position to FLT_MIN as it is the default value of the interfile field
    one_position_not_zero = (one_position_not_zero || p_image_fields[0]->bed_relative_position!=FLT_MIN);
    all_positions_not_zero = (all_positions_not_zero && p_image_fields[0]->bed_relative_position!=FLT_MIN);
  }

  // If one is not zero but not all, then throw a warning
  if (one_position_not_zero && !all_positions_not_zero)
  {
    Cerr("!!!!! castor-mergeBedPositions() -> The bed relative position is provided in some but not all interfile headers !" << endl);
  }

  // ============================================================================================================
  // If the provided relative positions are used, then recenter the mass to 0.
  // ============================================================================================================

  if (all_positions_not_zero)
  {
    // Verbose
    if (verbose>=1) Cout("castor-mergeBedPositions() -> Recenter relative bed positions" << endl);
    // We have to compute the center of mass of the provided bed positions, and move it to the 0-center of the Z axis
    FLTNB center = 0.;
    // Compute the center of mass
    for (uint32_t bed=0; bed<nb_beds; bed++) center += p_image_fields[bed]->bed_relative_position;
    center /= ((FLTNB)nb_beds);
    // Compute the new bed positions from their relative values to recenter the mass at 0
    for (uint32_t bed=0; bed<nb_beds; bed++) p_image_fields[bed]->bed_relative_position = p_image_fields[bed]->bed_relative_position - center;
  }

  // ============================================================================================================
  // If no relative positions provided, then compute them from the default scanner displacement
  // ============================================================================================================

  if (!all_positions_not_zero)
  {
    // Verbose
    if (verbose>=1) Cout("castor-mergeBedPositions() -> Recover default bed displacement from scanner configuration file" << endl);
    // Find scanner configuration file
    if (p_ScannerManager->FindScannerSystem(p_image_fields[0]->originating_system) )
    {
      Cerr("***** castor-mergeBedPositions() -> A problem occurred while searching for scanner system '" << p_image_fields[0]->originating_system << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
    // Get the bed displacement from the scanner configuration file
    FLTNB default_bed_displacement = -1.;
    if (ReadDataASCIIFile(p_ScannerManager->GetPathToScannerFile(), "multiple bed displacement", &default_bed_displacement, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** castor-mergeBedPositions() -> Multiple bed displacement field not found in the scanner configuration file '" <<
           p_ScannerManager->GetPathToScannerFile() << " !" << endl);
      Exit(EXIT_FAILURE);
    }
    // Check that the default bed displacement is more than 0. in the scanner
    if (default_bed_displacement<=0.)
    {
      Cerr("***** castor-mergeBedPositions() -> Default bed displacement from the scanner configuration file must be strictly positive !" << endl);
      Exit(EXIT_FAILURE);
    }
    // Verbose
    if (verbose>=1) Cout("  --> Bed displacement of " << default_bed_displacement << " mm" << endl);
    // Loop on the bed positions
    for (uint32_t bed=0; bed<nb_beds; bed++)
    {
      // Compute bed offset (remember here that the 0. on the Z axis is at the center of the whole image)
      FLTNB bed_offset = 0.;
      // For odd numbers of bed positions
      if (nb_beds%2==1) bed_offset = ((FLTNB)( ((int32_t)(bed))-((int32_t)(nb_beds))/2 )) * default_bed_displacement;
      // For even numbers of bed positions
      else bed_offset = (((FLTNB)( ((int32_t)(bed))-((int32_t)(nb_beds))/2 )) + 0.5) * default_bed_displacement;
      // Record the value
      p_image_fields[bed]->bed_relative_position = bed_offset;
    }
  }

  // ============================================================================================================
  // Find the bed with the minimum position and compute the displacement of each bed with respect to that
  // ============================================================================================================

  // Verbose
  if (verbose>=1)
  {
    Cout("castor-mergeBedPositions() -> Use following relative bed positions:" << endl);
    for (uint32_t bed=0; bed<nb_beds; bed++) Cout("  --> Bed " << bed << " | Relative position (mm): " << p_image_fields[bed]->bed_relative_position << endl);
  }

  // Find the minimum position that will be used as the reference
  FLTNB reference_minimum_position = p_image_fields[0]->bed_relative_position;
  for (uint32_t bed=1; bed<nb_beds; bed++)
    if (p_image_fields[bed]->bed_relative_position<reference_minimum_position)
      reference_minimum_position = p_image_fields[bed]->bed_relative_position;

  // Compute bed displacement of each bed position from the reference
  FLTNB* p_displacement_from_reference = (FLTNB*)malloc(nb_beds*sizeof(FLTNB));
  for (uint32_t bed=0; bed<nb_beds; bed++)
    p_displacement_from_reference[bed] = p_image_fields[bed]->bed_relative_position - reference_minimum_position;

  // ============================================================================================================
  // Check that the displacement of each bed is compatible with the slice thickness
  // ============================================================================================================

  // Define the tolerance associated to the rounding from the slice thickness and the tolerance provided as a percent of the slice thickness
  FLTNB tolerance = p_image_fields[0]->vox_size[2] * tolerance_as_percent_slice_thickness / 100.;
  // The number of slices defining the displacement of each bed with respect to the reference
  INTNB* p_displacement_slices = (INTNB*)calloc(nb_beds,sizeof(INTNB));

  // Loop on beds
  for (uint32_t bed=0; bed<nb_beds; bed++)
  {
    // Compute floating point number of slices
    FLTNB nb_slices_fltnb = p_displacement_from_reference[bed] / p_image_fields[0]->vox_size[2];
    // Compute integer number of slices
    FLTNB nb_slices_intnb = ((FLTNB)((INTNB)nb_slices_fltnb));
    // Case where the floating number minus the integer number is less than 0.5; this means that
    // the closest integer to the floating number of slices is below it.
    if (nb_slices_fltnb-nb_slices_intnb<0.5)
    {
      // If the difference is below the tolerance, we are ok taking the 'integerized' number of slices as the displacement
      if (nb_slices_fltnb-nb_slices_intnb<tolerance) p_displacement_slices[bed] = ((INTNB)nb_slices_fltnb);
      // Otherwise, the displacement is not compatible with the slice thickness
      else
      {
        Cerr("***** castor-mergeBedPositions() -> Bed position of bed " << bed << " is not compatible with the slice thickness of " << p_image_fields[0]->vox_size[2] << " mm !" << endl);
        Exit(EXIT_FAILURE);
      }
    }
    // Case where the floating number minus the integer number is more than 0.5; this means that
    // the closest integer to the floating number of slices is above it.
    else
    {
      // If minus the difference + 1 is below the tolerance, we are ok taking the 'integerized' number of slices + 1 as the displacement
      if (nb_slices_intnb+1.-nb_slices_fltnb<tolerance) p_displacement_slices[bed] = ((INTNB)nb_slices_fltnb) + 1;
      // Otherwise, the displacement is not compatible with the slice thickness
      else
      {
        Cerr("***** castor-mergeBedPositions() -> Bed position of bed " << bed << " is not compatible with the slice thickness of " << p_image_fields[0]->vox_size[2] << " mm !" << endl);
        Exit(EXIT_FAILURE);
      }
    }
  }

  // ============================================================================================================
  // Compute the dimension of the output image based on the reference bed position and the maximum displacement
  // ============================================================================================================

  // Search for the maximum displacement
  INTNB max_displacement_slices = p_displacement_slices[0];
  for (uint32_t bed=1; bed<nb_beds; bed++) if (p_displacement_slices[bed]>max_displacement_slices) max_displacement_slices = p_displacement_slices[bed];

  // The output image dimension will be the size of the input image plus the maximum displacement in number of slices
  INTNB whole_dim_z = max_displacement_slices + dim_z;
  INTNB whole_dim_xyz = dim_xy * whole_dim_z;

  // Verbose
  if (verbose>=1) Cout("castor-mergeBedPositiosn() -> Output merged image has " << whole_dim_z << " slices" << endl);

  // ============================================================================================================
  // If sensitivity, then load it as weights, otherwise set everything to 1
  // ============================================================================================================

  // Create the weights matrix
  FLTNB** p_weights = (FLTNB**)malloc(nb_beds*sizeof(FLTNB*));
  for (uint32_t bed=0; bed<nb_beds; bed++) p_weights[bed] = NULL;

  // If we have sensitivity
  if (path_to_sens_filename.size()>0)
  {
    // Verbose
    if (verbose>=1) Cout("castor-mergeBedPositions() -> Load sensitivity images" << endl);
    // Create interfile image fields and image pointers for sensitivity
    Intf_fields** p_sens_fields = (Intf_fields**)malloc(nb_beds*sizeof(Intf_fields*));
    // Read interfile images
    for (uint32_t bed=0; bed<nb_beds; bed++)
    {
      // Create the image fields
      p_sens_fields[bed] = new Intf_fields;
      // Read the image
      if (invert_order_flag) p_weights[bed] = IntfLoadImageFromScratch(path_to_sens_filename[nb_beds-1-bed],p_sens_fields[bed],verbose);
      else p_weights[bed] = IntfLoadImageFromScratch(path_to_sens_filename[bed],p_sens_fields[bed],verbose);
      // Check reading
      if (p_weights[bed]==NULL)
      {
        Cerr("***** castor-mergeBedPositions() -> A problem occurred while reading sensitivity for bed " << bed+1 << " !" << endl);
        return 1;
      }
    }
    // Check dimensions consistency between all image fields
    for (uint32_t bed=1; bed<nb_beds; bed++)
    {
      if (IntfCheckDimensionsConsistency(*p_sens_fields[0],*p_sens_fields[bed]))
      {
        Cerr("***** castor-mergeBedPositions() -> Inconsistency between sensitivity dimensions !" << endl);
        Exit(EXIT_FAILURE);
      }
    }
    // Check dimensions consistency between input images and sensitivity images
    for (uint32_t bed=0; bed<nb_beds; bed++)
    {
      if (IntfCheckDimensionsConsistency(*p_image_fields[bed],*p_sens_fields[bed]))
      {
        Cerr("***** castor-mergeBedPositions() -> Inconsistency between input image and sensitivity dimensions for bed " << bed+1 << " !" << endl);
        Exit(EXIT_FAILURE);
      }
    }
  }
  // If we do not have sensitivity
  else
  {
    // Verbose
    if (verbose>=1) Cout("castor-mergeBedPositions() -> Initialize weights to 1." << endl);
    // We set everything to 1.
    for (uint32_t bed=0; bed<nb_beds; bed++)
    {
      p_weights[bed] = (FLTNB*)malloc(dim_xyz*sizeof(FLTNB));
      for (INTNB v=0; v<dim_xyz; v++) p_weights[bed][v] = 1.;
    }
  }

  // ============================================================================================================
  // Compute the whole body image
  // ============================================================================================================

  // Verbose
  if (verbose>=1) Cout("castor-mergeBedPositiosn() -> Compute whole image" << endl);
  // Allocate whole image
  FLTNB* p_whole_image = (FLTNB*)malloc(whole_dim_xyz*sizeof(FLTNB));
  for (INTNB v=0; v<whole_dim_xyz; v++) p_whole_image[v] = 0.;
  // Allocate whole weights
  FLTNB* p_whole_weight = (FLTNB*)malloc(whole_dim_xyz*sizeof(FLTNB));
  for (INTNB v=0; v<whole_dim_xyz; v++) p_whole_weight[v] = 0.;

  // Loop on all beds
  for (uint32_t bed=0; bed<nb_beds; bed++)
  {
    // Loop on all slices
    INTNB z;
    #pragma omp parallel for private(z) schedule(guided)
    for (z=0; z<dim_z; z++)
    {
      // Compute whole z index
      INTNB whole_z = z + p_displacement_slices[bed];
      // For efficiency
      INTNB base_z = z * dim_xy;
      INTNB whole_base_z = whole_z * dim_xy;
      // Loop on all voxels inside the slice
      for (INTNB xy=0; xy<dim_xy; xy++)
      {
        // Add the contribution of this bed to the whole image
        p_whole_image[whole_base_z+xy] += p_images[bed][base_z+xy] * p_weights[bed][base_z+xy];
        // Add the weights to the global weights
        p_whole_weight[whole_base_z+xy] += p_weights[bed][base_z+xy];
      }
    }
  }
  // Normalize the whole image with respect to the weights
  for (INTNB v=0; v<whole_dim_xyz; v++)
  {
    if (p_whole_weight[v]>0.) p_whole_image[v] /= p_whole_weight[v];
    else p_whole_image[v] = 0.;
  }

  // ============================================================================================================
  // Flip output image axially
  // ============================================================================================================

  if (flip_output_axial)
  {
    // Verbose
    if (verbose>=1) Cout("castor-mergeBedPositions() -> Flip output images axially" << endl);
    // Half loop over Z
    for (INTNB z_1=0; z_1<whole_dim_z/2; z_1++)
    {
      // Compute opposite Z
      INTNB z_2 = whole_dim_z - 1 - z_1;
      // For efficiency
      INTNB base_z1 = z_1 * dim_xy;
      INTNB base_z2 = z_2 * dim_xy;
      // Loop over Y and X
      for (INTNB xy=0; xy<dim_xy; xy++)
      {
        // Compute both indices
        INTNB indice_1 = base_z1 + xy;
        INTNB indice_2 = base_z2 + xy;
        // Switch voxels for the image
        FLTNB buffer = p_whole_image[indice_1];
        p_whole_image[indice_1] = p_whole_image[indice_2];
        p_whole_image[indice_2] = buffer;
        // Switch voxels for the sensitivity
        buffer = p_whole_weight[indice_1];
        p_whole_weight[indice_1] = p_whole_weight[indice_2];
        p_whole_weight[indice_2] = buffer;
      }
    }
  }

  // ============================================================================================================
  // Save output image
  // ============================================================================================================

  // Verbose
  if (verbose>=1) Cout("castor-mergeBedPositions() -> Write output image(s)" << endl);

  // Build interfile fields (we use first one of the input images and modify the relevant fields)
  // Set header metadata using Image Dimensions object
  p_image_fields[0]->mtx_size[2] = whole_dim_z;
  p_image_fields[0]->nb_total_imgs = whole_dim_z;
  p_image_fields[0]->nb_bytes_pixel = sizeof(FLTNB);
  /* TODO: UPDATE TIME START AND DURATION BUT WILL HAVE TO WRITE THEM IN THE HEADERS OF IMAGES WRITTEN DURING ITERATIONS
  ap_IF->study_duration = ap_ID->GetFinalTimeStopInSec(0) -
                          ap_ID->GetFrameTimeStartInSec(0,0);
  for(int fr=0 ; fr<ap_ID->GetNbTimeFrames() ; fr++)
  {
    ap_IF->image_duration.push_back(ap_ID->GetFrameDurationInSec(0, fr));
    ap_IF->frame_group_pause.push_back((fr == 0) ? 0 
                                                 : ap_ID->GetFrameTimeStartInSec(0,fr) - ap_ID->GetFrameTimeStopInSec(0,fr-1));
  }
  */
  // File name
  string whole_image_filename = p_OutputManager->GetPathName() + p_OutputManager->GetBaseName() + "_whole";
  // Write the image
  if (IntfWriteImageFromIntfFields(whole_image_filename, p_whole_image, *p_image_fields[0], verbose))
  {
    Cerr("***** castor-mergeBedPositions() -> A problem occurred while writing output whole image !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Save global weights if asked for
  if (save_weights)
  {
    // File name
    string whole_weight_filename = p_OutputManager->GetPathName() + p_OutputManager->GetBaseName() + "_weights";
    // Write the image
    if (IntfWriteImageFromIntfFields(whole_weight_filename, p_whole_weight, *p_image_fields[0], verbose))
    {
      Cerr("***** castor-mergeBedPositions() -> A problem occurred while writing global weights image !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  // ============================================================================================================
  // Free memory properly
  // ============================================================================================================

  // Input images
  if (p_images)
  {
    for (uint32_t bed=0; bed<nb_beds; bed++) if (p_images[bed]) free(p_images[bed]);
    free(p_images);
  }
  // Weights
  if (p_weights)
  {
    for (uint32_t bed=0; bed<nb_beds; bed++) if (p_weights[bed]) free(p_weights[bed]);
    free(p_weights);
  }
  // Whole image
  if (p_whole_image) free(p_whole_image);
  // Whole weights
  if (p_whole_weight) free(p_whole_weight);

  // Ending
  #ifdef CASTOR_MPI
  MPI_Finalize();
  #endif
  return EXIT_SUCCESS;
}
