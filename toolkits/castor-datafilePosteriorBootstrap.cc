
/*!
  \file
  \ingroup utils
  \brief  This program is used to create a resampled histogram datafile for the posterior bootstrap reconstruction
          see "Reconstruction, analysis and interpretation of posterior probability distributions of PET images, using the posterior bootstrap" by Marina Filipovic et al., PMB 2021.
*/

#include "gVariables.hh"
#include "gOptions.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "vDataFile.hh"
#include "iDataFilePET.hh"
#include "iDataFileSPECT.hh"
#include "iDataFileCT.hh"
#include "iEventHistoPET.hh"
#include "iEventHistoSPECT.hh"
#include "iEventHistoCT.hh"
#include "iEventListPET.hh"
#include "iEventListSPECT.hh"
#include "iEventListCT.hh"
#include "sOutputManager.hh"
#include "sRandomNumberGenerator.hh"

// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================
//                                                        H E L P     F U N C T I O N S
// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================


/*!
  \fn      ShowHelp()
  \brief   Display main command line options for castor-datafilePosteriorBootstrap
*/
void ShowHelp()
{
  // Show help
  cout << endl;
  cout << "Usage:  castor-datafilePosteriorBootstrap  -df datafile.cdh  -(f/d)out output" << endl;
  cout << endl;
  cout << "This program can be used to create a resampled histogram datafile for the posterior bootstrap reconstruction." << endl;
  cout << "See the related paper 'Reconstruction, analysis and interpretation of posterior probability distributions of PET images, using the posterior bootstrap' by Marina Filipovic et al., PMB 2021." << endl;
  cout << endl;
  cout << "[Mandatory parameters]:" << endl;
  cout << "  -df datafile.cdh     : Give a CASToR histogram datafile." << endl;
  cout << "  -fout name           : Give the root name for all output files (no default, alternative to -dout)" << endl;
  cout << "  -dout name           : Give the name of the output directory where all output files will be written (no default, alternative to -fout)" << endl;
  cout << endl;
  cout << "[Options]:" << endl;
  cout << "  -rng                 : Give a seed for the random number generator (should be >=0)" << endl;
  cout << "  --norm-w             : Set an additional normalization step to maintain the same number of counts: implies more RAM usage and forces the random weights to form a probability distribution" << endl;
  cout << "                         -> by default set to false" << endl;
  cout << "  -vb value            : Give the verbosity level, from 0 (no verbose) to 2 (default: 1)" << endl;
  cout << "  --help,-h,-help      : Print out this help page." << endl; // managed by main
  cout << endl;
  #ifdef CASTOR_OMP
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
  // MPI stuff (we make all instances except the first one returning 0 directly)
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

  // Input datafile
  string datafile = "";

  // --------------------------------------------------------------------------------
  // Output settings
  // --------------------------------------------------------------------------------

  // Output directory name.
  string path_dout = "";
  // Or root name
  string path_fout = "";

  // --------------------------------------------------------------------------------
  // Miscellaneous settings
  // --------------------------------------------------------------------------------

  // Verbose level
  int verbose = 1;
  // Initial seed for random number generator
  int64_t random_generator_seed = -1;
  // normalize the histogram after resampling so that the number of counts remains identical to the number of counts in the original histogram
  // forces the random weights to form a probability distribution (sum to 1)
  bool norm_w = false;

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
    // RNG seed
    else if (option=="-rng")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-datafileBootstrap() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &random_generator_seed))
      {
        Cerr("***** castor-datafileBootstrap() -> Exception when trying to read provided number '" << random_generator_seed << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // normalize so that the number of counts remains the same
    else if (option=="--norm-w")
    {
      norm_w = true;
      i++;
    }
    // General verbosity level
    else if (option=="-vb")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-datafileBootstrap() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose))
      {
        Cerr("***** castor-datafileBootstrap() -> Exception when trying to read provided verbosity level '" << verbose << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }

    // --------------------------------------------------------------------------------
    // Input settings
    // --------------------------------------------------------------------------------

    // Images
    else if (option=="-df") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-datafileBootstrap() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      datafile = (string)argv[i+1];
      i++;
    }

    // --------------------------------------------------------------------------------
    // Output settings
    // --------------------------------------------------------------------------------

    // Name of the output directory
    else if (option=="-dout") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-datafileBootstrap() -> Argument missing for option: " << option << endl);
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
        Cerr("***** castor-datafileBootstrap() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_fout = argv[i+1];
      i++;
    }
    // --------------------------------------------------------------------------------
    // Unknown option!
    // --------------------------------------------------------------------------------

    else
    {
      Cerr("***** castor-datafileBootstrap() -> Unknown option '" << option << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  // ============================================================================================================
  // Some checks
  // ============================================================================================================

  // Data file
  if (datafile=="")
  {
    Cerr("***** castor-datafileBootstrap() -> Please provide an input datafile !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Output files
  if (path_fout.empty() && path_dout.empty())
  {
    Cerr("***** castor-datafileBootstrap() -> Please provide an output option for output files (-fout or -dout) !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Check that only one option has been provided
  if (!path_fout.empty() && !path_dout.empty())
  {
    Cerr("***** castor-datafileBootstrap() -> Please provide either output option -fout or -dout but not both !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Create sOutputManager
  // ============================================================================================================
  sOutputManager* p_OutputManager = sOutputManager::GetInstance();  
  // Set verbose level
  p_OutputManager->SetVerbose(verbose);
  // Set MPI rank
  p_OutputManager->SetMPIRank(mpi_rank);
  // Initialize output directory and base name
  if (p_OutputManager->InitOutputDirectory(path_fout, path_dout))
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred while initializing output directory !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Log command line
  if (p_OutputManager->LogCommandLine(argc,argv))
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred while logging command line arguments !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Random number generator
  // ============================================================================================================

  sRandomNumberGenerator* p_RandomNumberGenerator = sRandomNumberGenerator::GetInstance(); 
  p_RandomNumberGenerator->SetVerbose(verbose);
  // Use a user-provided seed to initialize the RNG if one has been provided. Use random number otherwise.
  int extra_gen = 0;
  int nb_threads = 1;
  if (random_generator_seed>=0) p_RandomNumberGenerator->Initialize(random_generator_seed, nb_threads, extra_gen);
  else p_RandomNumberGenerator->Initialize(nb_threads, extra_gen);

  // ============================================================================================================
  // Create sScannerManager (in order to get the datafile type)
  // ============================================================================================================
  sScannerManager* p_ScannerManager = sScannerManager::GetInstance();  
  p_ScannerManager->SetVerbose(verbose);
  // Get system name from the dataFile
  string scanner_name = "";
  if (ReadDataASCIIFile(datafile, "Scanner name", &scanner_name, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred while trying to find the system name in the datafile header !" << endl);
    Exit(EXIT_FAILURE);
  } 
  if (p_ScannerManager->FindScannerSystem(scanner_name) )
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred while searching for scanner system !" << endl);
    Exit(EXIT_FAILURE);
  } 
  if (p_ScannerManager->BuildScannerObject() )
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred during scanner object construction ! !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->InstantiateScanner() )
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred while creating Scanner object !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->GetGeometricInfoFromDataFile(datafile))
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred while retrieving scanner fields from the datafile header !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->BuildLUT() )
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred while generating/reading the LUT !" << endl);
    Exit(EXIT_FAILURE);
  } 

  // ============================================================================================================
  // Create the input datafile
  // ============================================================================================================

  // Create a default image dimensions and quantification object
  oImageDimensionsAndQuantification* p_ID = new oImageDimensionsAndQuantification();
  p_ID->SetDefault();
  p_ID->SetVerbose(verbose);
  // Create the datafile based on the data type
  vDataFile* p_DataFile = NULL;
  if (p_ScannerManager->GetScannerType() == SCANNER_PET) p_DataFile = new iDataFilePET();
  else if (p_ScannerManager->GetScannerType() == SCANNER_SPECT_CONVERGENT) p_DataFile = new iDataFileSPECT(); 
  else if (p_ScannerManager->GetScannerType() == SCANNER_CT) p_DataFile = new iDataFileCT(); 
  else
  {
    Cerr("***** castor-datafileBootstrap() -> Unknown scanner type (" << p_ScannerManager->GetScannerType() << ") for datafile construction ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  p_DataFile->SetImageDimensionsAndQuantification(p_ID);
  p_DataFile->SetHeaderDataFileName(datafile);
  p_DataFile->SetVerbose(verbose);
  p_DataFile->SetBedIndex(0);
  bool do_not_affect_quantification = false;
  if (p_DataFile->ReadInfoInHeader(do_not_affect_quantification))
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred during datafile header reading ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_DataFile->CheckParameters())
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred while checking datafile parameters ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_DataFile->ComputeSizeEvent())
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred in datafile initialization ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_DataFile->InitializeMappedFile())
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred in datafile initialization ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_DataFile->PrepareDataFile())
  {
    Cerr("***** castor-datafileBootstrap() -> A problem occurred in datafile preparation ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  // Check if datafile is a histogram, otherwise, throw an error
  if (p_DataFile->GetDataMode()!=MODE_HISTOGRAM)
  {
    Cerr("***** castor-datafileBootstrap() -> The input datafile is not a histogram, this program is only suitable to histogram files !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Create the output datafile
  // ============================================================================================================

  // Create output datafile object
  vDataFile* p_OutputDataFile = NULL;
  if (p_ScannerManager->GetScannerType() == SCANNER_PET) p_OutputDataFile = new iDataFilePET();
  else if (p_ScannerManager->GetScannerType() == SCANNER_SPECT_CONVERGENT) p_OutputDataFile = new iDataFileSPECT(); 
  else if (p_ScannerManager->GetScannerType() == SCANNER_CT) p_OutputDataFile = new iDataFileCT(); 
  // Build output data file from the input one
  if (p_OutputDataFile->SetParametersFrom(p_DataFile))
  {
    Cerr("***** castor-datafileBootstrap() -> An error occurred while setting parameters of output file from input file !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Open output file
  if (p_OutputDataFile->OpenFileForWriting())
  {
    Cerr("***** castor-datafileBootstrap() -> An error occurred while opening file for writing !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Compute posterior bootstrap values for each event bin
  // ============================================================================================================

  // Get the number of events
  int64_t nb_events = p_DataFile->GetSize();
  // Verbose
  if (verbose>=1) Cout("castor-datafileBootstrap() -> Draw posterior bootstrap values for all the " << nb_events << " events" << endl);

  sRandomNumberGenerator::Engine& rgen = sRandomNumberGenerator::GetInstance()->GetGenerator();

  // number of counts for the original and the resampled histograms
  HPFLTNB nb_counts_new = 0.;
  HPFLTNB nb_counts_old = 0.;
  // renormalization factor to ensure that the number of counts remains the same
  HPFLTNB norm_factor = 1.;
  // variable needed for the renormalization step, used only if renormalization required
  HPFLTNB** vals = NULL;
  // allocate memory if required
  if (norm_w) vals = new HPFLTNB*[nb_events];
  int64_t printing_index = 0;
  int64_t i = 0;
  for (i=0; i<nb_events; i++)
  {
    // Verbose
    if (verbose>=1)
    {
      if (printing_index%5000==0)
      {
        int percent = ((int)( ((FLTNB)(i)) * 100. / ((FLTNB)(nb_events/nb_threads)) ));
        cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
             << percent << " %                    " << flush;
      }
      printing_index++;
    }

    // get current event
    vEvent* event = p_DataFile->GetEvent(i);
    if (event==NULL)
    {
      Cerr("***** castor-datafileBootstrap() -> An error occurred while getting the event from index " << i << " !" << endl);
      if (p_OutputDataFile->CloseFile())
        Cerr("***** castor-datafileBootstrap() -> An error occurred while closing file during writing !" << endl);
      Exit(EXIT_FAILURE);
    }

    if (norm_w) vals[i] = new HPFLTNB[event->GetNbValueBins()];
    // loop over bins
    for (int b=0; b<event->GetNbValueBins(); b++)
    {
      // build a Gamma distribution of shape parameter equal to the observed bin value and of scale parameter equal to 1
      // only if there are some detected counts, otherwise keep the 0.
      FLTNB val = event->GetEventValue(b);
      if (norm_w) vals[i][b] = val;
      if (val>0.)
      {
        nb_counts_old += val;
        // build the Gamma distribution for this bin
        gamma_distribution<FLTNB> gamma_distrib(val, 1.);
        // draw the new posterior bootstrap value
        val = gamma_distrib(rgen);
        nb_counts_new += val;
        // save the randomized value or write it into the event
        if (norm_w) vals[i][b] = val;
        else event->SetEventValue(b, (FLTNBDATA)val);
      }
    }
    // Write the event
    if (!norm_w) p_OutputDataFile->WriteEvent(event);

  }

  // If required, normalize the randomized values and actually write the output file
  if (norm_w)
  {
    printing_index = 0;
    norm_factor = nb_counts_old/nb_counts_new;
    nb_counts_new = 0.;
    for (i=0; i<nb_events; i++)
    {
      // Verbose
      if (verbose>=1)
      {
        if (printing_index%5000==0)
        {
          int percent = ((int)( ((FLTNB)(i)) * 100. / ((FLTNB)(nb_events/nb_threads)) ));
          cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
               << percent << " %                    " << flush;
        }
        printing_index++;
      }

      // get current event
      vEvent* event = p_DataFile->GetEvent(i);
      if (event==NULL)
      {
        Cerr("***** castor-datafileBootstrap() -> An error occurred while getting the event from index " << i << " !" << endl);
        if (p_OutputDataFile->CloseFile())
          Cerr("***** castor-datafileBootstrap() -> An error occurred while closing file during writing !" << endl);
        Exit(EXIT_FAILURE);
      }

      // loop over bins
      for (int b=0; b<event->GetNbValueBins(); b++)
      {
        if (vals[i][b]>0.)
        {
          // apply the renormalization factor
          vals[i][b] *= norm_factor;
          nb_counts_new += vals[i][b];
          // modify the event bin accordingly
          event->SetEventValue(b, (FLTNBDATA)(vals[i][b]));
        }
      }
      // free memory
      delete vals[i];

      // Write the event
      p_OutputDataFile->WriteEvent(event);
    }
    // free memory
    delete vals;
  }


  // Verbose
  if (verbose>=2) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                       << "  --> 100 %                        " << endl;
                       cout<< "number of counts (original histogram) "<<nb_counts_old<<endl;
                       cout<< "number of counts (resampled histogram) "<<nb_counts_new<<endl;
                       if (norm_w) cout<< "renormalization factor "<<norm_factor<<endl;
  // Close file
  if (p_OutputDataFile->CloseFile())
  {
    Cerr("***** castor-datafileBootstrap() -> An error occurred while closing file after writing !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Set number of events
  p_OutputDataFile->SetNbEvents(nb_events);
  // Write header
  if (p_OutputDataFile->WriteHeader())
  {
    Cerr("***** castor-datafileBootstrap() -> An error occurred while writing output header file !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Exit
  // ============================================================================================================

  // Ending
  #ifdef CASTOR_MPI
  MPI_Finalize();
  #endif
  return EXIT_SUCCESS;
}
