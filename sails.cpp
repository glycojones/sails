
// SAILS
// Software for the Automated Identification of Linked Sugars
// 2013-2019 by Mihaela Atanasova, Kevin Cowtan & Jon Agirre
// Department of Chemistry, University of York, UK
// For funding details see README.md
// mailto: jon.agirre@york.ac.uk
//        -

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>
#include <stdlib.h>

#include "sails-lib.h"

int main( int argc, char** argv )
{

    clipper::String program_version = "MKI";
    CCP4Program prog( "Sails", program_version.c_str(), "$Date:07/06/2019" );

    // defaults
    clipper::String title;
    clipper::String ipmtz =     "NONE";
    clipper::String ipcol_fo =  "NONE";
    clipper::String ipcol_hl =  "NONE";
    clipper::String ipcol_pw =  "NONE";
    clipper::String ipcol_fc =  "NONE";
    clipper::String ipcol_fr =  "NONE";
    clipper::String ippdb =     "NONE";
    clipper::String oppdb =     "xyzout.pdb";
    clipper::String opmap =     "NONE";
    clipper::String ipcodes =   "NONE";
    std::vector<clipper::String> input_codes;
    double res_in = 2.0;         // Resolution limit
    int nhit = 10;
    int verbose = 0;
    double step = 10; // default value, works fine for pyranoses
    bool useMap = false;
    bool usePDB = false; // for using an input model
    clipper::CCP4MAPfile mapFile;
    clipper::Xmap<float> xwrk;
    clipper::String building_options = ""; // can be 'all', or { 'nglycans' 'oglycans' 'ligands' }
    sails::data::build_options options = { false, false, false };

    std::cout << std::endl << "Copyright 2013-2019 Mihaela Atanasova, Kevin Cowtan and Jon Agirre." << std::endl;

    // command input
    CCP4CommandInput args( argc, argv, true );
    int arg = 0;

    while ( ++arg < args.size() )
    {
        if ( args[arg] == "-title" )
        {
            if ( ++arg < args.size() )
                title = args[arg];
        }
        else if ( args[arg] == "-mtzin" )
        {
            if ( ++arg < args.size() )
                ipmtz = args[arg];
        }
        else if ( args[arg] == "-mapin" )
        {
            if ( ++arg < args.size() )
            {
                useMap = true;
                mapFile.open_read(args[arg]);
                mapFile.import_xmap(xwrk);
            }
        }
        else if ( args[arg] == "-pdbin" )
        {
            if ( ++arg < args.size() )
            {
                ippdb = args[arg];
                usePDB = true;
            }
        }
        else if ( args[arg] == "-pdbout" )
        {
            if ( ++arg < args.size() )
                oppdb = args[arg];
        }
        else if ( args[arg] == "-mapout" )
        {
            if ( ++arg < args.size() )
                opmap  = args[arg];
        }
        else if ( args[arg] == "-colin-fo" )
        {
            if ( ++arg < args.size() )
                ipcol_fo = args[arg];
        }
        else if ( args[arg] == "-colin-hl" )
        {
            if ( ++arg < args.size() )
                ipcol_hl = args[arg];
        }
        else if ( args[arg] == "-colin-phifom" )
        {
            if ( ++arg < args.size() )
                ipcol_pw = args[arg];
        }
        else if ( args[arg] == "-colin-fc" )
        {
            if ( ++arg < args.size() )
                ipcol_fc = args[arg];
        }
        else if ( args[arg] == "-colin-free" )
        {
            if ( ++arg < args.size() )
                ipcol_fr = args[arg];
        }
        else if ( args[arg] == "-hits" )
        {
            if ( ++arg < args.size() )
                nhit = clipper::String(args[arg]).i();
        }
        else if ( args[arg] == "-step" )
        {
            if ( ++arg < args.size() )
                step = clipper::String(args[arg]).i();
        }
        else if ( args[arg] == "-resolution" )
        {
            if ( ++arg < args.size() )
                res_in = clipper::String(args[arg]).f();
                if ( res_in > 6.0 || res_in < 0.5 ) {
                  std::cout << "Error: Sails is unlikely to work at the specified resolution. Exiting... " << std::endl;
                  prog.set_termination_message( "Failed" );
                  return(1);
                }
        }
        else if ( args[arg] == "-searchfor" ) {
          if ( ++arg < args.size() ) {
            ipcodes = args[arg];
            sails::get_input_codes (ipcodes, input_codes);
          }
        }
        else if ( args[arg] == "-build" )
        {
            if ( ++arg < args.size() )
                building_options = args[arg];
        }
        else if ( args[arg] == "-verbose" )
        {
            if ( ++arg < args.size() )
                verbose = clipper::String(args[arg]).i();
        }
        else
        {
            std::cout << "\nUnrecognised:\t" << args[arg] << std::endl;
            args.clear();
        }
    }

    if ( args.size() <= 1 || ( !useMap && ( ipmtz == "NONE" || ipcol_fo == "NONE" || ipcol_hl == "NONE")) ) // print help and exit
    {
        std::cout << "Usage: sails\n\n"
                  << "\t-mtzin <filename>\t\tCOMPULSORY\n"
                  << "\t-colin-fo <colpath>\t\tspecify F,SIGF\n"
                  << "\t-colin-hl <colpath>\t\tphase probabilities (ABCD)\n"
                  << "\t-colin-phifom <colpath>\t\tphase probabilities (Phi/Fom)\n"
                  << "\t-colin-fc <colpath>\t\tinput map coefficients\n"
                  << "\t-colin-free <colpath>\t\tfree R column\n"
                  << "\t-pdbin <filename>\t\tinitial model to be extended\n"
                  << "\t-pdbout <filename>\t\toutput atomic model\n"
                  << "\t-build \t\t\t\t{all,nglycans,oglycans,ligands}\n\t\t\t\t\tDefault is nglycans\n"
                  << "\t-resolution <resolution/A>\t(def: 2.0 angstroems)\n"
                  << "\t-step <degrees>\t\t\tangular sampling (def: 15 deg)\n"
                  << "\t-mapin <filename>\t\tuse precomputed electron density\n"
                  << "\t\t\t\t\tignores -mtzin and -colin" << std::endl;

        prog.set_termination_message( "Failed" );
        return(1);
    }

    sails::process_building_options ( building_options, options );

    std::cout << "Detecting:\n"
              << "N-linked glycans: ";
    if ( options.nglycans )
        std::cout << "yes\n";
    else
        std::cout << "no\n";

    std::cout << "O-linked glycans: ";
    if ( options.oglycans )
        std::cout << "yes\n";
    else
        std::cout << "no\n";

    std::cout << "Ligands: ";
    if ( options.ligands )
        std::cout << "yes\n" << std::endl;
    else
        std::cout << "no\n" << std::endl;


    using clipper::data32::Compute_fphi_from_fsigf_phifom;
    clipper::Resolution resol;
    clipper::CCP4MTZfile mtzfile;

    if (!useMap)
    {
        // other initialisations
        mtzfile.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
        mtzfile.open_read( ipmtz );
        double res = clipper::Util::max( mtzfile.resolution().limit(), res_in );
        std::cout << std::endl << "Using reflections up to " << res << " Angstroms" << std::endl;
        mtzfile.close_read();
        resol = clipper::Resolution( res );

        // Get work reflection data
        clipper::HKL_info hkls;
        mtzfile.open_read( ipmtz );
        hkls.init( mtzfile.spacegroup(), mtzfile.cell(), resol, true );

        clipper::HKL_data<clipper::data32::F_sigF>  wrk_f ( hkls );
        clipper::HKL_data<clipper::data32::ABCD>    wrk_hl( hkls );
        clipper::HKL_data<clipper::data32::Phi_fom> wrk_pw( hkls );
        clipper::HKL_data<clipper::data32::F_phi>   fphi( hkls );
        clipper::HKL_data<clipper::data32::Flag>    flag( hkls );
        if ( ipcol_fo != "NONE" ) mtzfile.import_hkl_data( wrk_f ,ipcol_fo );
        if ( ipcol_hl != "NONE" ) mtzfile.import_hkl_data( wrk_hl,ipcol_hl );
        if ( ipcol_pw != "NONE" ) mtzfile.import_hkl_data( wrk_pw,ipcol_pw );
        if ( ipcol_fc != "NONE" ) mtzfile.import_hkl_data( fphi,  ipcol_fc );
        if ( ipcol_fr != "NONE" ) mtzfile.import_hkl_data( flag,  ipcol_fr );
        mtzfile.close_read();

        // apply free flag
        clipper::HKL_data<clipper::data32::F_sigF> wrk_f1 = wrk_f;

        //wrk_f1.mask( flag != 0 );
        for ( clipper::HKL_data_base::HKL_reference_index ih = hkls.first(); !ih.last(); ih.next() ) if ( flag[ih].flag() == 0 ) wrk_f1[ih] = clipper::data32::F_sigF();
        // and fill in hl
        clipper::Spacegroup cspg = hkls.spacegroup();
        clipper::Cell cxtl = hkls.cell();
        clipper::Grid_sampling grid = clipper::Grid_sampling( cspg, cxtl, hkls.resolution() );
        xwrk = clipper::Xmap<float>( cspg, cxtl, grid );

        // work map
        if ( ipcol_hl == "NONE" )
            wrk_hl.compute( wrk_pw, clipper::data32::Compute_abcd_from_phifom() );

        if ( ipcol_pw == "NONE" )
            wrk_pw.compute( wrk_hl, clipper::data32::Compute_phifom_from_abcd() );

        if ( ipcol_fc == "NONE" )
            fphi.compute( wrk_f1, wrk_pw, Compute_fphi_from_fsigf_phifom() );

            xwrk.fft_from( fphi );
    }

    clipper::MiniMol tmpmol;
    clipper::MMDBManager mmdb;
    clipper::MMDBfile mfile;

    const int mmdbflags = mmdb::MMDBF_IgnoreBlankLines | mmdb::MMDBF_IgnoreDuplSeqNum | mmdb::MMDBF_IgnoreNonCoorPDBErrors | mmdb::MMDBF_IgnoreRemarks | mmdb::MMDBF_EnforceUniqueChainID;

    if ( usePDB )
    {
        // load partial input model
        mfile.SetFlag( mmdbflags );
        mfile.read_file( ippdb );
        mfile.import_minimol( tmpmol );

        std::cout << "Using a partially-built structure as input: " << ippdb << std::endl;
        std::cout << tmpmol.spacegroup().symbol_hm() << " " << tmpmol.cell().format() << " " << tmpmol.atom_list().size() << std::endl;
    }

		clipper::MiniMol mol_new = sails::build_sugars (xwrk, options, step, nhit);

    clipper::MMDBfile pdbfile;
    pdbfile.export_minimol( mol_new );
    pdbfile.write_file( oppdb );

    prog.set_termination_message( "Normal termination" );
}
