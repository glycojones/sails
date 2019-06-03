
// Clipper app to extract fingerprints from pyranose rings
// 2013 Jon Agirre & Kevin Cowtan @ University of York
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk
//

#include <string.h>
#include <sails-lib.h>
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-contrib.h>
#include <clipper/minimol/minimol_utils.h>
#include <privateer/clipper-glyco.h>
#include <algorithm>
#include <stdlib.h>


int main( int argc, char** argv )
{
    CCP4Program prog( "Backend - minimax", "0.1", "$Date: 2013/05/21" );
    prog.set_termination_message( "Failed" );

    // defaults
    clipper::String title;
    clipper::String ipmtz = "NONE";
    clipper::String ippdb = "NONE";
    clipper::String ipcol_fc = "/*/*/FC_BEST";
    clipper::String atoms = " C1 , C3 , C5 "; // default atoms for aldopyranoses
    clipper::String atomt = " C2 , O5 , C4 ";
    clipper::String atomr = ""; // we're going to load up this vector with all the atoms from the input carbohydrate model
    clipper::String sugar = "";

    double res_in = 2.0;
    double cutoff = 0.0;
    bool verbose = false;
    clipper::String validation_string = "";
    int uiso = 9999;
    int dump = 0;

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
        else if ( args[arg] == "-pdbin" )
        {
            if ( ++arg < args.size() )
                ippdb = args[arg];
        }
        else if ( args[arg] == "-colin-fc" )
        {
            if ( ++arg < args.size() )
                ipcol_fc = args[arg];
        }
        else if ( args[arg] == "-validation" )
        {
            if ( ++arg < args.size() )
                validation_string = args[arg];
        }
        else if ( args[arg] == "-resolution" )
        {
            if ( ++arg < args.size() )
                res_in = clipper::String(args[arg]).f();
        }
        else if ( args[arg] == "-verbose" )
        {
            verbose = true;
        }
        else if ( args[arg] == "-cutoff" )
        {
            if ( ++arg < args.size() )
                cutoff = clipper::String(args[arg]).i();
        }
        else if ( args[arg] == "-sugar" )
        {
            if ( ++arg < args.size() )
            {
                sugar = args[arg];
                if (!clipper::MSugar::search_database ( sugar.trim().c_str() ) )
                {
                    std::cout << std::endl << "ERROR: " << sugar.trim() << " is not a supported sugar" << std::endl << std::endl;
                    exit(1);
                }
            }
        }
        else if ( args[arg] == "-uiso" )
        {
            if ( ++arg < args.size() )
                uiso = clipper::String(args[arg]).i();
        }
        else if ( args[arg] == "-dump" )
        {
            dump = 1;
        }
        else
        {
            std::cout << "\nUnrecognised:\t" << args[arg] << std::endl;
            args.clear();
        }
    }

    // end command input
    //

    if ( args.size() <= 1 )
    {
        std::cout << "\nUsage: minimax\n\t-pdbin <filename> COMPULSORY\n\t-mtzin <filename> Fc map will be calculated if no map coeffs are supplied.\n\t";
        std::cout << "-validation all,none,anomer,geometry,handedness,conformation\n\t-colin-fc <colpath>\n\t";
        std::cout << "-resolution <resolution/A>\n\t-sugar <3-Letter code>\n\t-cutoff <percentage>\n\t-uiso <int> optional smoothing\n" << std::endl << std::endl;
        exit(1);
    }

    sails::data::validation_flags flags;
    flags.validate_anomer = false;
    flags.validate_handedness = false;
    flags.validate_conformation = false;
    flags.validate_geometry = false;

    sails::process_validation_options ( validation_string, flags );

    flags.validate_geometry ? std::cout << "\ngeometry \ton\n" : std::cout << "geometry \toff\n";
    flags.validate_conformation ? std::cout << "conformation \ton\n" : std::cout << "conformation \toff\n";
    flags.validate_anomer ? std::cout << "anomer \t\ton\n" : std::cout << "anomer \t\toff\n";
    flags.validate_handedness ? std::cout << "handedness \ton\n" << std::endl : std::cout << "handedness \toff\n" << std::endl;

    // other initialisations
    using clipper::data32::Compute_fphi_from_fsigf_phifom;
    clipper::CCP4MTZfile mtzfile;

    std::vector<clipper::String> atomv; // create a vector with each of the comma-separated atoms coming from a String
    std::vector<clipper::String> atomvt;
    std::vector<clipper::String> atomvr;

    // Get reference model
    clipper::MiniMol mol;

    clipper::MMDBManager mmdb;
    const int mmdbflags = mmdb::MMDBF_IgnoreBlankLines | mmdb::MMDBF_IgnoreDuplSeqNum | mmdb::MMDBF_IgnoreNonCoorPDBErrors | mmdb::MMDBF_IgnoreRemarks | mmdb::MMDBF_EnforceUniqueChainID;
    clipper::MMDBfile mfile;
    mfile.SetFlag( mmdbflags );
    mfile.read_file( ippdb );
    mfile.import_minimol( mol );

    std::cout << mol.spacegroup().symbol_hm() << " " << mol.cell().format() << " " << mol.atom_list().size() << std::endl;

    const clipper::MAtomNonBond nb(mol, 5.0);

    clipper::HKL_data<clipper::data32::F_phi> fphi;

    if ( ipmtz != "NONE" ) // experimental data
    {
        // Get resolution for calculation
        mtzfile.open_read( ipmtz );
        mtzfile.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
        double res = clipper::Util::max( mtzfile.resolution().limit(), res_in );
        mtzfile.close_read();
        clipper::Resolution resol = clipper::Resolution( res );

        // Get reference reflection data
        mtzfile.open_read( ipmtz );
        clipper::HKL_sampling hkls( mtzfile.cell(), resol );
        fphi.init( mtzfile.spacegroup(), mtzfile.cell(), hkls );
        mtzfile.import_hkl_data( fphi, ipcol_fc );
        mtzfile.close_read();
    }
    else     // use calculated data: preferred way as of 2015
    {
        clipper::Resolution resol = clipper::Resolution( res_in );
        clipper::HKL_sampling hkls( mol.cell(), resol );
        fphi.init( mol.spacegroup(), mol.cell(), hkls );
        clipper::SFcalc_aniso_fft<float> sfc( 2.5, 1.5, 0.25 );
        sfc( fphi, mol.atom_list() );
        std::cout << fphi.num_obs() << std::endl;  // calculates structure factors for the whole molecule
        std::cout << fphi.hkl_info().num_reflections() << std::endl;
    }

    // scale the data when asked for it. We use a temperature factor U. B = 8 Pi ^^2 * U

    if (uiso != 9999)
    {
        clipper::data32::Compute_scale_u_iso_fphi compute_uiso( 1.0, uiso );
        fphi.compute( fphi, compute_uiso );
    }

    // make map
    clipper::Grid_sampling grid( fphi.hkl_info().spacegroup(), fphi.hkl_info().cell(), fphi.hkl_info().resolution() );
    clipper::Xmap<float> xref( fphi.hkl_info().spacegroup(), fphi.hkl_info().cell(), grid );
    xref.fft_from( fphi ); // calculates an Fc map containing all the loaded fragments

    if (dump == 1)
    {
        clipper::CCP4MAPfile mapOut;
        mapOut.open_write("mapDump.map");      // dump internal map
        mapOut.export_xmap( xref );
        mapOut.close_write();
    }

    // stats
    clipper::Map_stats ms( xref );
    std::cout << "Mean: " << ms.mean() << " ; Standard deviation: " << ms.std_dev() << std::endl;

    // make model
    clipper::Cell_descr cd( 100.0, 100.0, 100.0 );
    clipper::Cell cell( cd );
    clipper::MiniMol mol_tmp[59];

    for (int i =0; i < 59 ; i++)
        mol_tmp[i] = clipper::MiniMol( clipper::Spacegroup::p1(), cell );
        // a variable to store the transformed coordinates for each conformer
        // was 41 but now incorporates furanoses and should be 59

    // make rtops
    std::vector<clipper::RTop_orth> rtops[59];

    const clipper::MAtomNonBond& manb = clipper::MAtomNonBond( mol, 5.0 );

    for ( int chn = 0; chn < mol.size(); chn++ )  // for each chain
    {
        for ( int res = 0; res < mol[chn].size(); res++ ) // and for each residue
        {
            clipper::MMonomer mm = mol[chn][res];
            int an[3];

            if (mm.type() == sugar)
            {
                clipper::MSugar ms(mol, mm, manb);
                clipper::MAtom m_anomeric = ms.anomeric_substituent ( );
//                m_anomeric.set_id("AN");
                mm.insert ( m_anomeric );
                std::cout << " m_anomeric is " << m_anomeric.id() << " " << m_anomeric.id() << std::endl;
                std::vector<clipper::MAtom> ring_atoms = ms.ring_members();
                atomvr.clear();

                bool worthy = true;

                if ( flags.validate_anomer && ( ! ms.ok_with_anomer() ) )
                    worthy = false;

                if ( flags.validate_handedness && ( ! ms.ok_with_chirality() ) )
                    worthy = false;

                if ( flags.validate_conformation && ( ! ms.ok_with_conformation() ) )
                    worthy = false;

                if ( flags.validate_geometry && ((! ms.ok_with_bonds_rmsd()) || (! ms.ok_with_angles_rmsd() ) ) )
                    worthy = false;

                if ( worthy )
                {

                    if ( ring_atoms.size() == 5 ) // furanose sugar
                    {
                        atomv.clear();
                        atomv.push_back( ring_atoms[1].name().trim() ); // first carbon
                        atomv.push_back( ring_atoms[3].name().trim() ); // third carbon
                        atomv.push_back( ring_atoms[4].name().trim() ); // fourth carbon, need three atoms

                        atomvt.clear();
                        atomvt.push_back( ring_atoms[0].name().trim() ); // ring oxygen
                        atomvt.push_back( ring_atoms[2].name().trim() ); // second carbon
                        atomvt.push_back( ring_atoms[4].name().trim() ); // fourth carbon, need three atoms as well
                    }
                    else // we've got a pyranose in our hands
                    {
                        atomv.clear();
                        atomv.push_back( ring_atoms[1].name().trim() ); // first carbon
                        atomv.push_back( ring_atoms[3].name().trim() ); // third carbon
                        atomv.push_back( ring_atoms[5].name().trim() ); // fourth carbon, need three atoms

                        atomvt.clear();
                        atomvt.push_back( ring_atoms[0].name().trim() ); // ring oxygen
                        atomvt.push_back( ring_atoms[2].name().trim() ); // second carbon
                        atomvt.push_back( ring_atoms[4].name().trim() ); // fourth carbon, need three atoms as well
                    }

                    for ( int i = 0; i < atomv.size(); i++ )
                        an[i] = mm.lookup( atomv[i], clipper::MM::ANY );

                    if ( an[0] >= 0 && an[1] >= 0 && an[2] >= 0 )
                    {
                        atomr = mm[0].name().trim();
                        atomvr.push_back ( mm[0].name().trim() );

                        for ( int index = 1; index < mm.size() ; index++ )
                        {
                            atomr += ", " + mm[index].name().trim();
                            atomvr.push_back ( mm[index].name().trim() );
                        }

                        // std::cout << "an[0]=" << an[0] << " an[1]=" << an[1] << " an[2]=" << an[2] << std::endl;
                        // std::cout << "atomr=" << atomr << std::endl;

                        // for ( int j = 0 ; j < atomvt.size() ; j++ )
                        //    std::cout << "atomvt[" << j << "]=" << atomvt[j] << " ";

                        // std::cout << std::endl;

                        // for ( int j = 0 ; j < atomv.size() ; j++ )
                        //    std::cout << "atomv[" << j << "]=" << atomv[j] << " ";

                        // make rtop
                        clipper::Coord_orth c1 = mm[an[0]].coord_orth()-mm[an[1]].coord_orth();
                        clipper::Coord_orth c2 = mm[an[2]].coord_orth()-mm[an[1]].coord_orth();
                        clipper::Vec3<> v1( (c1.unit()+c2.unit()).unit() );
                        clipper::Vec3<> v2( clipper::Vec3<>::cross(c1,c2).unit() );
                        clipper::Vec3<> v3( clipper::Vec3<>::cross(v1,v2).unit() );
                        clipper::RTop_orth rtop( clipper::Mat33<>( v1[0], v2[0], v3[0],
                                                                   v1[1], v2[1], v3[1],
                                                                   v1[2], v2[2], v3[2] ), mm[an[1]].coord_orth() );

                        rtops[ms.conformation_code()].push_back( rtop );

                        // and model
                        clipper::MPolymer mp;
                        mm = mol[chn][res];
                        mm.insert ( m_anomeric );
                        mm.transform( rtop.inverse() );
                        mp.insert( mm );
                        mol_tmp[ms.conformation_code()].insert( mp );
                    }
                }
            }
        }
    }

    std::vector<int> resvt( atomvt.size(), 0 );
    std::vector<int> resvr( atomvr.size(), 0 );

    for ( int i = 0; i <  atomvt.size(); i++ )
    {
        if ( atomvt[i].length() > 4 )
        {
            if ( atomvt[i][4] == '+' ) resvt[i]++;
            if ( atomvt[i][4] == '-' ) resvt[i]--;
            atomvt[i] = atomvt[i].substr(0,4); // crop the name to 0..3
        }
        //std::cout << atomvt[i] << "[" << resvt[i] << "]";
    }

    std::cout << std::endl;

    for ( int i = 0; i <  atomvr.size(); i++ )
    {
        if ( atomvr[i].length() > 4 )
        {
            if ( atomvr[i][4] == '+' ) resvr[i]++;
            if ( atomvr[i][4] == '-' ) resvr[i]--;
            atomvr[i] = atomvr[i].substr(0,4); // crop the name to 0..3
        }
        //std::cout << atomvr[i] << "[" << resvr[i] << "]";
    }

    int nfound = 0;

    for (int i = 0; i < 59 ; i++)
        if (rtops[i].size() != 0)
        {
            nfound += rtops[i].size();
            std::cout << "Found " << rtops[i].size() << " in conformation " << clipper::data::conformational_landscape[i] << std::endl;
        }

    if (nfound == 0)
    {
        system("touch not_in_scope");
        std::cout << std::endl << "Not in scope - exiting..." << std::endl;
        return 0;
    }


    // now check for outliers
    std::vector<bool> flag[59];
    bool any_results[59] = { false };

    for (int conf = 0; conf < 59 ; conf++) // filter out incomplete models
    {
        flag[conf] = std::vector<bool>( mol_tmp[conf].size(), true );

        for ( int c = 0; c < mol_tmp[conf].size(); c++ )
        {
            for ( int index = 0; index < atomvt.size(); index++ )
            {
                int a = mol_tmp[conf][c][0].lookup( atomvt[index], clipper::MM::ANY );

                if ( a < 0 ) // a=-1 in case that minimol didn't find it
                {
                    flag[conf][c] = false;
                }
                else if ( mol_tmp[conf][c][0][a].occupancy() < 1.0 ) flag[conf][c] = false;
            }
        }

        double d2cut = 999.0;  // filter out most different models

        while ( atomvt.size() > 0 && int(d2cut) > cutoff && mol_tmp[conf].size() > 0) // was inflooping because ~0 (double) > 0 (int)
        {
            std::vector<double> d2( mol_tmp[conf].size(), 0.0 );

            for ( int index = 0; index < atomvt.size(); index++ )
            {
                int r = resvt[index];
                clipper::Coord_orth ca(0.0,0.0,0.0);
                double na = 0.0;

                for ( int c = 0; c < mol_tmp[conf].size(); c++ )
                {

                    if ( flag[conf][c] )
                    {
                        int a = mol_tmp[conf][c][r].lookup( atomvt[index], clipper::MM::ANY );
                        ca += mol_tmp[conf][c][r][a].coord_orth();
                        na += 1.0;
                    }
                }
                ca = (1.0/na) * ca; ///////////////////////////////////// calc centre

                for ( int c = 0; c < mol_tmp[conf].size(); c++ )
                {
                    if ( flag[conf][c] )
                    {
                        int a = mol_tmp[conf][c][r].lookup( atomvt[index], clipper::MM::ANY );
                        d2[c] += ( mol_tmp[conf][c][r][a].coord_orth() - ca ).lengthsq();
                    }
                }
            }

            for ( int c = 0; c < mol_tmp[conf].size(); c++ ) d2[c] /= atomvt.size();

            std::vector<double> s2 = d2;
            std::sort( s2.begin(), s2.end() );
            d2cut = s2[ 99*s2.size()/100 ]; // was s2[ (100-int(cutoff))*s2.size()/100 ];

            for ( int c = 0; c < mol_tmp[conf].size(); c++ )
                if ( d2[c] > d2cut )
                    flag[conf][c] = false;

            int nf = 0;

            for ( int c = 0; c < mol_tmp[conf].size(); c++ )
                if ( flag[conf][c] )
                    nf++;

            std::cout << "Accepted " << nf << " with cutoff distance " << d2cut << std::endl;

            if (nf != 0) any_results[conf] = true; // prevent dumping the core

        }
    }


    // make maps
    const double radius = 6.0;
    const double sampling = 0.25;
    double extent = radius + 2.0 * sampling;  // model box size
    double ng = rint( extent / sampling ); // number of grid points, round to integer

    clipper::Grid nxgrid( 2*int(ng)+1, 2*int(ng)+1, 2*int(ng)+1 );
    clipper::RTop<> rtop( clipper::Mat33<>( 1.0/sampling, 0.0, 0.0,
                                            0.0, 1.0/sampling, 0.0,
                                            0.0, 0.0, 1.0/sampling ),
                                            clipper::Vec3<>( ng, ng, ng ) );

///////////////////////////////////////////////////////


    clipper::NXmap<float> rhomin[59];
    clipper::NXmap<float> rhomax[59];
    clipper::NXmap<float> rhoavg[59];

    for (int i = 0; i < 59 ; i++)   // THIS IS QUITE A BIG 'FOR'
    {
        if ( mol_tmp[i].size() > 0 && any_results[i] )
        {
            rhomin[i] = clipper::NXmap<float>( nxgrid, rtop );
            rhomax[i] = clipper::NXmap<float>( nxgrid, rtop );
            rhoavg[i] = clipper::NXmap<float>( nxgrid, rtop );

            rhomin[i] = 1.0e12;
            rhomax[i] = 1.0e12;
            rhoavg[i] = 0.0;

            std::cout << std::endl << "Setting up maps for conformer " << clipper::data::conformational_landscape[i] << std::endl;

            double nop = 0.0;  // number of points

            typedef clipper::NXmap<float>::Map_reference_index MRI;

            for ( int r = 0; r < rtops[i].size(); r++ ) // for each fragment
            {
                if ( flag[i][r] )
                {
                    for ( MRI ix = rhomin[i].first(); !ix.last(); ix.next() )  // run across the map
                    {
                        float rho = xref.interp<clipper::Interp_cubic>( (rtops[i][r]*ix.coord_orth()).coord_frac(xref.cell()) ); // get map value for each fragment from the big map
                        rhomin[i][ix] = clipper::Util::min( rhomin[i][ix], rho );
                        rhomax[i][ix] = clipper::Util::min( rhomax[i][ix], -rho );
                        rhoavg[i][ix] += rho;
                    }
                    nop += 1.0;
                }
            }


            for ( MRI ix = rhomin[i].first(); !ix.last(); ix.next() )
                rhoavg[i][ix] *= 1.0/nop;

            // get list of peaks
            clipper::Coord_grid cg;
            clipper::Coord_grid cgu( 1, 0, 0 ), cgv( 0, 1, 0 ), cgw( 0, 0, 1 );
            float rho0;
            std::vector<std::pair<float,int> > peakmax, peakmin; // beware of these ones

            for ( cg.w() = 1; cg.w() < 2*ng; cg.w()++ )
                for ( cg.v() = 1; cg.v() < 2*ng; cg.v()++ )
                    for ( cg.u() = 1; cg.u() < 2*ng; cg.u()++ )
                    {
                        rho0 = rhomin[i].get_data( cg );

                        if ( rho0 > rhomin[i].get_data( cg - cgu ) &&
                             rho0 > rhomin[i].get_data( cg + cgu ) &&
                             rho0 > rhomin[i].get_data( cg - cgv ) &&
                             rho0 > rhomin[i].get_data( cg + cgv ) &&
                             rho0 > rhomin[i].get_data( cg - cgw ) &&
                             rho0 > rhomin[i].get_data( cg + cgw ) )
                                peakmin.push_back( std::pair<float,int>( rho0, nxgrid.index(cg) ) );

                        rho0 = rhomax[i].get_data( cg );

                        if ( rho0 > rhomax[i].get_data( cg - cgu ) &&
                             rho0 > rhomax[i].get_data( cg + cgu ) &&
                             rho0 > rhomax[i].get_data( cg - cgv ) &&
                             rho0 > rhomax[i].get_data( cg + cgv ) &&
                             rho0 > rhomax[i].get_data( cg - cgw ) &&
                             rho0 > rhomax[i].get_data( cg + cgw ) )
                                peakmax.push_back( std::pair<float,int>( rho0, nxgrid.index(cg) ) );
                    }

                    std::cout << peakmin.size() << " " << peakmax.size() << std::endl;

                    std::sort( peakmin.begin(), peakmin.end() );
                    std::sort( peakmax.begin(), peakmax.end() );
                    std::reverse( peakmin.begin(), peakmin.end() );
                    std::reverse( peakmax.begin(), peakmax.end() );

                    for ( int index = 0; index < 10; index++ )  // output the first 10 min peaks
                        std::cout << peakmin[index].first << " " << rhomin[i].coord_orth(nxgrid.deindex(peakmin[index].second).coord_map()).format() << std::endl;

                    for ( int index = 0; index < 10; index++ )  // output the first 10 max peaks
                        std::cout << peakmax[index].first << " " << rhomax[i].coord_orth(nxgrid.deindex(peakmax[index].second).coord_map()).format() << std::endl;

                    // export the maps
                    clipper::CCP4MAPfile mapfile;
                    std::string namePDB = ippdb.substr(0, ippdb.find(".pdb"));
                    mapfile.open_write( clipper::String(namePDB + "-" + sugar + "-" +  clipper::data::conformational_landscape[i] + "_min.map") );
                    mapfile.set_cell( cell );
                    mapfile.export_nxmap( rhomin[i] );
                    mapfile.close_write();
                    mapfile.open_write( clipper::String(namePDB + "-" + sugar + "-" +  clipper::data::conformational_landscape[i] + "_max.map") );
                    mapfile.set_cell( cell );
                    mapfile.export_nxmap( rhomax[i] );
                    mapfile.close_write();
                    mapfile.open_write( clipper::String(namePDB + "-" + sugar + "-" +  clipper::data::conformational_landscape[i] + "_avg.map") );
                    mapfile.set_cell( cell );
                    mapfile.export_nxmap( rhoavg[i] );
                    mapfile.close_write();

                    // model
                    clipper::MiniMol mol_new( mol_tmp[i].spacegroup(), mol_tmp[i].cell() );

                    clipper::MPolymer mp1; mp1.set_id("A"); mol_new.insert(mp1); // chain A: accepted
                    clipper::MPolymer mp2; mp2.set_id("B"); mol_new.insert(mp2); // chain B: turned down

                    for ( int c = 0; c < mol_tmp[i].size(); c++ )
                    {
                        int nc = flag[i][c] ? 0 : 1;
                        int nr = 0;
                        if ( mol_new[nc].size() > 0 )
                            nr = mol_new[nc][mol_new[nc].size()-1].seqnum() + 2;

                        for ( int r = 0; r < 1; r++ )
                        {
                            clipper::MMonomer mm = mol_tmp[i][c][r];
                            mm.set_seqnum( nr++ );
                            mol_new[nc].insert( mm );
                        }
                    }

                    clipper::MMDBfile pdbfile;
                    pdbfile.export_minimol( mol_new );
                    pdbfile.write_file(clipper::String(namePDB + "-" + sugar + "-" + clipper::data::conformational_landscape[i] + ".pdb") );

                    // now get most representative fragment
                    std::vector<clipper::Coord_orth> cavg( atomvr.size(), clipper::Coord_orth( 0.0, 0.0, 0.0 ) );
                    double navg = 1.0;

                    for ( int c = 0; c < mol_tmp[i].size(); c++ )
                        if ( flag[i][c] )
                        {
                            for ( int irep = 0; irep < atomvr.size(); irep++ ) // for each representative atom
                            {

                                int r = resvr[irep];
                                int a;

                                for ( a = 0; a < mol_tmp[i][c][r].size(); a++ )
                                    if ( mol_tmp[i][c][r][a].id().trim() == atomvr[irep].trim() )
                                        break;

                                if ( a == mol_tmp[i][c][r].size() )
                                {
                                    std::cout << "resvr[i]: " << r << " a: " << a << " ChainID: " << mol_tmp[i][c][r].id() << std::endl;
                                    clipper::Message::message(  clipper::Message_fatal( "Missing " + atomvr[irep]) );
                                }

                                cavg[irep] += mol_tmp[i][c][r][a].coord_orth();
                                navg += 1.0;
                            }
                        }

                    for ( int index = 0; index < atomvr.size(); index++ )
                        cavg[index] = (1.0/navg) * cavg[index];

                    double cmin = 0;
                    double dmin = 1.0e20;

                    for ( int c = 0; c < mol_tmp[i].size(); c++ )
                        if ( flag[i][c] )
                        {
                            double d = 0.0;
                            for ( int index = 0; index < atomvr.size(); index++ )
                            {
                                int r = resvr[index];
                                int a;

                                for ( a = 0; a < mol_tmp[i][c][r].size(); a++ )
                                    if ( mol_tmp[i][c][r][a].id().trim() == atomvr[index].trim() )
                                        break;

                                if ( a == mol_tmp[i][c][r].size() )
                                {
                                    std::cout << "resvr[i]: " << r << " a: " << a << " ChainID: " << mol_tmp[i][c][r].id() << std::endl;
                                }

                                d += ( cavg[index] - mol_tmp[i][c][r][a].coord_orth() ).lengthsq();
                            }

                            if ( d < dmin )
                            {
                                dmin = d;
                                cmin = c;
                            }
                        }

                    clipper::MiniMol mol_r( mol_tmp[i].spacegroup(), mol_tmp[i].cell() );
                    clipper::MPolymer mp;
                    clipper::MMonomer mm;

                    mm.set_seqnum( 0 );
                    mm.set_type( "SUG" ); // was U
                    mp.insert( mm );

                    for ( int index = 0; index < atomvr.size(); index++ )
                    {
                        int r = resvr[index];
                        for ( int a = 0; a < mol_tmp[i][cmin][r].size(); a++ )
                            if ( mol_tmp[i][cmin][r][a].id() == atomvr[index] )
                                mp[r].insert( mol_tmp[i][cmin][r][a] );
                    }

                    if (atomvr.size() > 0)
                    {
                        mp.set_id( "H" );
                        mol_r.insert( mp );
                        clipper::MMDBfile pdbfile1;
                        pdbfile1.export_minimol( mol_r );
                        pdbfile1.write_file( clipper::String(namePDB + "-" + sugar + "-" +  clipper::data::conformational_landscape[i] + "_repr.pdb") );
                    }
                }
            } ///////// END OF THE BIG 'FOR'

    prog.set_termination_message( "Normal termination" );
    system("touch fingerprinted");
}
