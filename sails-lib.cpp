/* \file sails-lib.cpp
    A set of utilities that helps the Privateer do his job */
// version  0.1.0
// 2013 Jon Agirre & Kevin Cowtan, The University of York
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk
//
// Compile with 'make DEBUG=-DDUMP' to enable debugging messages
//
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA

#include "sails-lib.h"

void sails::insert_coot_prologue_scheme ( std::fstream& output )
{
    output  << "; This script has been created by Sails (Atanasova, Cowtan and Agirre, 2013-19)\n"
            << "(set-graphics-window-size 1873 968)\n"
            << "(set-graphics-window-position 0 0)\n"
            << "(set-go-to-atom-window-position 0 19)\n"
            << "(set-display-control-dialog-position 366 20)\n"
            << "(vt-surface 2)\n"
            << "(set-clipping-front  0.00)\n"
            << "(set-clipping-back  0.00)\n"
            << "(set-map-radius 10.00)\n"
            << "(set-iso-level-increment  0.0500)\n"
            << "(set-diff-map-iso-level-increment  0.0050)\n"
            << "(set-colour-map-rotation-on-read-pdb 21.00)\n"
            << "(set-colour-map-rotation-on-read-pdb-flag 1)\n"
            << "(set-colour-map-rotation-on-read-pdb-c-only-flag 1)\n"
            << "(set-swap-difference-map-colours 0)\n"
            << "(set-background-colour  0.00  0.00  0.00)\n"
            << "(set-symmetry-size 13.00)\n"
            << "(set-symmetry-colour-merge  0.50)\n"
            << "(set-symmetry-colour  0.10  0.20  0.80)\n"
            << "(set-symmetry-atom-labels-expanded 0)\n"
            << "(set-active-map-drag-flag 1)\n"
            << "(set-show-aniso 0)\n"
            << "(set-aniso-probability 50.00)\n"
            << "(set-smooth-scroll-steps 40)\n"
            << "(set-smooth-scroll-limit 10.00)\n"
            << "(set-font-size 2)\n"
            << "(set-rotation-centre-size  0.10)\n"
            << "(set-default-bond-thickness 5)\n"
            << "(scale-zoom  0.20)\n"
            << "(set-nomenclature-errors-on-read \"auto-correct\")"
            << "(set-run-state-file-status 0)\n";
}

void sails::insert_coot_files_loadup_scheme ( std::fstream& output, const clipper::String& pdb, const clipper::String& mapbest, const clipper::String& mapdiff, const clipper::String& mapomit, bool mode )
{
    if (!mode) output << "(handle-read-draw-molecule \"" << pdb << "\")\n";

    if ( mapbest == "" ) // no map output
    {
        output << "(set-last-map-colour 1.00  0.13  0.89)\n"
           << "(interesting-things-gui \"Report from Sails\"\n\t(list\n\t\t";
    }
    else
    {
        if (!mode)
            output << "(handle-read-ccp4-map \"" << mapbest << "\" 0)\n" << "(handle-read-ccp4-map \"" << mapomit << "\" 1)\n";

        output << "(set-last-map-colour 1.00  0.13  0.89)\n"
           << "(interesting-things-gui \"Report from Sails\"\n\t(list\n\t\t";
    }
}

void sails::insert_coot_files_loadup_python ( std::fstream& output, const clipper::String& pdb, const clipper::String& mapbest, const clipper::String& mapdiff, const clipper::String& mapomit, bool mode )
{
    if (!mode) output  << "handle_read_draw_molecule (\"" << pdb << "\")\n";

        if ( mapbest == "" ) // no map output
    {
        output << "set_last_map_colour  (1.00,  0.13,  0.89)\n"
               << "interesting_things_gui (\"Report from Sails\",[\n";
    }
    else
    {
        if (!mode)
                output << "handle_read_ccp4_map (\"" << mapbest << "\", 0)\n" << "handle_read_ccp4_map (\"" << mapomit << "\", 1)\n";

            output << "set_last_map_colour  (1.00,  0.13,  0.89)\n"
           << "interesting_things_gui (\"Report from Sails\",[\n";
    }
}

void sails::insert_coot_epilogue_scheme ( std::fstream& output )
{
    output  << "\n\n))\n(set-scroll-wheel-map 3)\n"
            << "(set-matrix 60.00)\n"
                        << "(set-refine-with-torsion-restraints 1)\n"
            << "(set-show-symmetry-master 0)\n";
}

void sails::insert_coot_prologue_python ( std::fstream& output )
{

    output  << "# This script has been created by Sails\n"
            << "set_graphics_window_size (1873, 968)\n"
        << "set_graphics_window_position (0, 0)\n"
        << "set_go_to_atom_window_position (0, 19)\n"
        << "vt_surface (2)\n"
        << "set_clipping_front  (0.00)\n"
        << "set_clipping_back  (0.00)\n"
        << "set_map_radius (10.00)\n"
        << "set_iso_level_increment  (0.0500)\n"
        << "set_diff_map_iso_level_increment  (0.0050)\n"
        << "set_colour_map_rotation_on_read_pdb (21.00)\n"
        << "set_colour_map_rotation_on_read_pdb_flag (1)\n"
        << "set_colour_map_rotation_on_read_pdb_c_only_flag (1)\n"
        << "set_swap_difference_map_colours (0)\n"
        << "set_background_colour  (0.00,  0.00,  0.00)\n"
        << "set_symmetry_size (13.00)\n"
        << "set_symmetry_colour_merge  (0.50)\n"
        << "set_symmetry_colour  (0.10,  0.20,  0.80)\n"
        << "set_symmetry_atom_labels_expanded (0)\n"
        << "set_active_map_drag_flag (1)\n"
        << "set_show_aniso (0)\n"
        << "set_aniso_probability (50.00)\n"
        << "set_smooth_scroll_steps (40)\n"
        << "set_smooth_scroll_limit (10.00)\n"
        << "set_font_size (2)\n"
        << "set_rotation_centre_size (0.10)\n"
        << "set_default_bond_thickness (4)\n"
        << "scale_zoom (0.20)\n"
        << "set_nomenclature_errors_on_read (\"auto-correct\")\n"
        << "set_run_state_file_status (0)\n"
            << "toggle_idle_spin_function\n";
}

void sails::insert_coot_epilogue_python ( std::fstream& output )
{
    output  << "\n\n])\nset_scroll_wheel_map (3)\n"
            << "set_matrix (60.00)\n"
                        << "set_refine_with_torsion_restraints (1)\n"
            << "set_show_symmetry_master (0)\n";
}

void sails::insert_coot_go_to_sugar_scheme ( std::fstream& output, const clipper::Coord_orth& sugar_centre, const clipper::String& diagnostic )
{
    output  << "\t(list\t\"" << diagnostic << "\"\t" << sugar_centre.x() << "\t" << sugar_centre.y() << "\t" << sugar_centre.z() << ")\n";
}

void sails::insert_coot_go_to_sugar_python ( std::fstream& output, const clipper::Coord_orth& sugar_centre, const clipper::String& diagnostic )
{
    output  << "\t[\"" << diagnostic << "\",\t" << sugar_centre.x() << ",\t" << sugar_centre.y() << ",\t" << sugar_centre.z() << "],\n";
}

void sails::insert_coot_statusbar_text_scheme ( std::fstream& output, clipper::String& text)
{
    output  << "(add-status-bar-text \"" << text << "\")" ;
}

void sails::insert_coot_statusbar_text_python ( std::fstream& output, clipper::String& text )
{
    output  << "add_status_bar_text (\"" << text << "\")" ;
}

clipper::ftype sails::real_space_correlation ( const clipper::Xmap<float>& map1, const clipper::Xmap<float>& map2 )
{
    return 0.0;
}

clipper::MMonomer sails::get_ideal_monomer ( const sails::data::fingerprint& fp )
{
    clipper::MMonomer tmp_mon;
    tmp_mon.set_type ( fp.name_short );
    tmp_mon.set_id ( 0 );

    for ( int index = 0; index < fp.num_control_points ; index++ )
    {
        clipper::MAtom tmp_atm;
        clipper::Coord_orth coords(fp.atoms[index].x, fp.atoms[index].y, fp.atoms[index].z );
        tmp_atm.set_coord_orth ( coords );
        tmp_atm.set_name ( fp.atoms[index].atom_name );
        tmp_atm.set_element ( fp.atoms[index].atom_name[0] );
        tmp_atm.set_id ( 1 );
        tmp_atm.set_occupancy ( 1.0 );
        tmp_atm.set_u_iso ( 0.25 ); // aiming for an iso-B of ~20.0
        tmp_mon.insert ( tmp_atm );
    }
    return tmp_mon;
}

clipper::MMonomer sails::get_peak_monomer ( const sails::data::fingerprint& fp )
{
    clipper::MMonomer tmp_mon;
    tmp_mon.set_type ( fp.name_short );
    tmp_mon.set_id ( 0 );

    for ( int index = 0; index < fp.num_control_points ; index++ )
    {
        clipper::MAtom tmp_atm;
        clipper::Coord_orth coords(fp.peaks[index].x, fp.peaks[index].y, fp.peaks[index].z );
        tmp_atm.set_coord_orth ( coords );
        tmp_atm.set_name ( fp.peaks[index].atom_name );
        tmp_atm.set_element ( fp.peaks[index].atom_name[0] );
        tmp_atm.set_id ( 1 );
        tmp_atm.set_occupancy ( 1.0 );
        tmp_atm.set_u_iso ( 0.25 ); // aiming for an iso-B of ~20.0
        tmp_mon.insert ( tmp_atm );
    }
    return tmp_mon;
}

clipper::MMonomer sails::get_void_monomer ( const sails::data::fingerprint& fp )
{
    clipper::MMonomer tmp_mon;
    tmp_mon.set_type ( fp.name_short );
    tmp_mon.set_id ( 0 );

    for ( int index = 0; index < fp.num_control_points ; index++ )
    {
        clipper::MAtom tmp_atm;
        clipper::Coord_orth coords(fp.voids[index].x, fp.voids[index].y, fp.voids[index].z );
        tmp_atm.set_coord_orth ( coords );
        tmp_atm.set_name ( fp.voids[index].atom_name );
        tmp_atm.set_element ( fp.voids[index].atom_name[0] );
        tmp_atm.set_id ( 1 );
        tmp_atm.set_occupancy ( 1.0 );
        tmp_atm.set_u_iso ( 0.25 ); // aiming for an iso-B of ~20.0
        tmp_mon.insert ( tmp_atm );
    }
    return tmp_mon;
}

void sails::process_building_options ( clipper::String building_options, sails::data::build_options &flags )
{
    std::vector < clipper::String > buffer = building_options.split( "," );

    for ( int i = 0; i < buffer.size(); i++ )
    {
        if ( buffer[i].trim() == "all" )
        {
            flags.nglycans = true;
            flags.oglycans = true;
            flags.ligands = true;
            return;
        }
        else if ( buffer[i].trim() == "nglycans" )
        {
            flags.nglycans = true;
        }
        else if ( buffer[i].trim() == "oglycans" )
        {
            flags.oglycans = true;
        }
        else if ( buffer[i].trim() == "ligands" )
        {
            flags.ligands = true;
        }
    }
    return;
}

void sails::process_validation_options ( clipper::String validation_string, sails::data::validation_flags &flags )
{
    std::vector < clipper::String > buffer = validation_string.split( "," );

    for ( int i = 0; i < buffer.size(); i++ )
    {
        if ( buffer[i].trim() == "all" )
        {
            flags.validate_anomer = true;
            flags.validate_handedness = true;
            flags.validate_conformation = true;
            flags.validate_geometry = true;
            return;
        }
        else if ( buffer[i].trim() == "none" )
        {
            flags.validate_anomer = false;
            flags.validate_handedness = false;
            flags.validate_conformation = false;
            flags.validate_geometry = false;
            return;
        }
        else if ( buffer[i].trim() == "anomer" )
        {
            flags.validate_anomer = true;
        }
        else if ( buffer[i].trim() == "geometry" )
        {
            flags.validate_geometry = true;
        }
        else if ( buffer[i].trim() == "handedness" )
        {
            flags.validate_handedness = true;
        }
        else if ( buffer[i].trim() == "conformation" )
        {
            flags.validate_conformation = true;
        }
    }
    return;
}


clipper::MiniMol sails::build_sugars ( clipper::Xmap<float>& xwrk, sails::data::build_options& options, double step, int nhit )
{
    double rad, sigcut;
    clipper::MPolymer mprep;
    typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Pair_coord;
    std::vector<Pair_coord> all_co;


    // get cutoff (for optimisation)
    clipper::Map_stats stats( xwrk );
    sigcut = stats.mean() + 0.5*stats.std_dev();

    clipper::MAtom peak_atom, void_atom;
    clipper::MMonomer peaks, voids, ideal;

    ideal = sails::get_ideal_monomer ( sails::data::fingerprint_list[0] );
    peaks = sails::get_peak_monomer  ( sails::data::fingerprint_list[0] );
    voids = sails::get_void_monomer  ( sails::data::fingerprint_list[0] );

    mprep.insert ( ideal, -1 );

    // set up targets
    for ( int r = 0; r < peaks.size(); r++ )
        all_co.push_back( Pair_coord( peaks[r].coord_orth(), voids[r].coord_orth() ) );

    // get map radius
    clipper::Atom_list atoms = voids.atom_list();
    double r2 = 0.0;

    for ( int a = 0; a < atoms.size(); a++ )
    {
        double d2 = atoms[a].coord_orth().lengthsq();
        if ( d2 > r2 )
            r2 = d2;
    }

    rad = sqrt( r2 ) + 1.0;

    // make a list of rotations
    std::vector<clipper::RTop_orth> rots;

    // make a list of rotation ops to try
    float glim = 360.0;  // gamma
    float blim = 180.0;  // beta
    float alim = 360.0;  // alpha

    // do a uniformly sampled search of orientation space
    float anglim = clipper::Util::min( alim, glim );

    for ( float bdeg=step/2; bdeg < 180.0; bdeg += step )
    {
        float beta = clipper::Util::d2rad(bdeg);
        float spl = anglim/clipper::Util::intf(cos(0.5*beta)*anglim/step+1);
        float smi = anglim/clipper::Util::intf(sin(0.5*beta)*anglim/step+1);
        for ( float thpl=spl/2; thpl < 720.0; thpl += spl )
            for ( float thmi=smi/2; thmi < 360.0; thmi += smi )
            {
                float adeg = clipper::Util::mod(0.5*(thpl+thmi),360.0);
                float gdeg = clipper::Util::mod(0.5*(thpl-thmi),360.0);

                if ( adeg <= alim && bdeg <= blim && gdeg <= glim )
                {
                    float alpha = clipper::Util::d2rad(adeg);
                    float gamma = clipper::Util::d2rad(gdeg);
                    clipper::Euler_ccp4 euler( alpha, beta, gamma );
                    rots.push_back(clipper::RTop_orth(clipper::Rotation(euler).matrix()));
                }
            }
    }

    // feature search
    SSfind ssfind;
    ssfind.prep_xmap( xwrk, rad );
    ssfind.prep_search( xwrk );
    std::vector<SearchResult> results =

    ssfind.search( all_co, rots, sigcut, 0.0 );

    std::sort( results.begin(), results.end() );
    std::reverse( results.begin(), results.end() );
    std::cout << results.size() << std::endl;

    for ( int i = 0; i < clipper::Util::min( int(results.size()), nhit ); i++ )
        std::cout << i << " " << results[i].score << " " << results[i].rot << " " << results[i].trn << std::endl;

    // output
    clipper::MiniMol mol_new( xwrk.spacegroup(), xwrk.cell() );
    const clipper::String chainid1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const clipper::String chainid2 = "abcdefghijklmnopqrstuvwxyz";

    clipper::Grid_sampling grid = xwrk.grid_sampling();

    for ( int i = 0; i < clipper::Util::min( int(results.size()), nhit ); i++ )
    {
        clipper::String id; int roff;

        if ( i < 26 )
        {
            id = chainid1.substr( i, 1 );
            roff = 0;
        }
        else
        {
            id = chainid2.substr( (i-26)/100, 1 );
            roff = 10*(i%100);
        }

        clipper::MPolymer mprot = mprep;
        int ir = results[i].rot;
        int it = results[i].trn;
        clipper::RTop_orth rtop( rots[ir].rot(),
                                 xwrk.coord_orth( grid.deindex(it).coord_map() ) );
        mprot.transform( rtop );
        mprot.set_id( id );

        for ( int r = 0; r < mprot.size(); r++ )
            mprot[r].set_seqnum( roff+mprot[r].seqnum() );

        mol_new.insert( mprot );

    }

    return mol_new;
}

const sails::data::fingerprint sails::data::fingerprint_list[] =
{
     
    {
        "ARA", "ligand", "A", "L", "alpha-L-Arabinose", 10,
        { { "AO5", 2.826, 0.443, 0.013 }, { "AC5", 2.182, 0.000,-1.193 }, { "AC4", 0.772, 0.418,-1.248 },
          { "AO4", 0.631, 1.605,-1.628 }, { "AC1", 2.206, 0.000, 1.183 }, { "AO1", 2.875, 0.588, 2.298 },
          { "AC2", 0.731, 0.512, 1.249 }, { "AO2", 0.131, 0.300, 2.188 }, { "AC3", 0.000, 0.000, 0.000 }, { "AO3",-1.125, 0.610, 0.003 } } ,
        { { "V01", 2.758, 2.629, 3.794 }, { "V02", 2.490, 5.743,-5.634 }, { "V03", 0.870, 1.227, 5.015 },
          { "V04",-5.376, 2.851, 3.476 }, { "V05",-3.243, 3.161, 4.321 }, { "V06",-1.485,-2.860, 2.153 },
          { "V07",-3.001,-0.086,-6.081 }, { "V08", 2.412,-2.704,-0.215 }, { "V09", 3.461, 0.901,-3.155 }, { "V10",-0.204, 4.505,-2.072 } } ,
        { { "C1",  2.178, 0.000, 1.206 }, { "C2",  0.733, 0.508, 1.236 }, { "C3",  0.000, 0.000, 0.000 },
          { "C4",  0.728, 0.462,-1.271 }, { "C5",  2.172, 0.000,-1.203 }, { "O1",  2.874, 0.596, 2.298 },
          { "O2",  0.088, 0.115, 2.439 }, { "O3", -1.357, 0.457,-0.031 }, { "O4",  0.635, 1.868,-1.478 }, { "O5", 2.817, 0.460, 0.003 } }
    },  
     {
        "NAG", "nglycan", "B", "D", "N-Acetyl-D-Glucosamine", 15,
        { { "PK1", 0, 0, 0 }, { "PK1", 0.738, 0.455, -1.276 }, { "PK1", 2.222, 0, -1.237 },
          { "PK1", 3.031, 0.678, -2.389 }, { "PK1", -0.699, 0.854, 3.364 }, { "PK1", -1.386, 0.144, 4.516 },
          { "PK1", -0.002, -0.007, 2.517 }, { "PK1", -1.338, 0.547, -0.016 }, { "PK1", 0.061, -0.095, -2.424 },
          { "PK1", 3.139, -0.092, -3.582 }, { "PK1", -0.76, 2.089, 3.2 }, { "PK1", 3.064, 0.507, 2.434 } } ,
        { { "VD1", 5.25, 1.25, 3.5 }, { "VD1", 5, -0.5, -1.25 }, { "VD1", -2.75, -8.88178e-16, 2 },
          { "VD1", 3.5, -5, 1 }, { "VD1", 1, 1.5, -4.25 }, { "VD1", -3.25, -3.5, -3 },
          { "VD1", -1, 5, 2.5 }, { "VD1", 8.88178e-16, -2, 5.75 }, { "VD1", -4.5, 3.25, -2.75 },
          { "VD1", 6.5, -2.5, 3 }, { "VD1", 2.5, 5.75, -2 }, { "VD1", -2.25, -5.75, 1.5 } } ,
        { { "ATM", 2.278, 0, 1.268 }, { "ATM", 0.751, 0.417, 1.309 }, { "ATM", 0, 0, 0 },
          { "ATM", 0.738, 0.455, -1.276 }, { "ATM", 2.222, 0, -1.237 }, { "ATM", 3.031, 0.678, -2.389 },
          { "ATM", -0.699, 0.854, 3.364 }, { "ATM", -1.386, 0.144, 4.516 }, { "ATM", -0.002, -0.007, 2.517 },
          { "ATM", -1.338, 0.547, -0.016 }, { "ATM", 0.061, -0.095, -2.424 }, { "ATM", 2.887, 0.396, 0.013 },
          { "ATM", 3.139, -0.092, -3.582 }, { "ATM", -0.76, 2.089, 3.2 }, { "ATM", 3.064, 0.507, 2.434 } }


    },

    {
        "BMA", "nglycan", "B", "D", "beta-D-mannopyranose", 3,
        { { "P", 0.0, 0.0, 0.0 }, { "P", 0.0, 0.0, 0.0 }, { "P", 0.0, 0.0, 0.0 } } ,
        { { "V", 0.0, 0.0, 0.0 }, { "V", 0.0, 0.0, 0.0 }, { "V", 0.0, 0.0, 0.0 } } ,
        { { "C", 0.0, 0.0, 0.0 }, { "C", 0.0, 0.0, 0.0 }, { "C", 0.0, 0.0, 0.0 } }
    } /* , 

         
   {
        "NAG", "nglycan" , "B", "D", "N-Acetyl-D-Glucosamine", 6,
		{ { "PK1", -2.405, 16.081, 13.767 }, { "PK1", -1.235, 16.595, 12.999 }, { "PK1", 0.011, 16.688, 13.863 },
          { "PK1", -0.059, 16.986, 15.056 } } ,
        { { "VD1", -0.75, 1.5, 1.25 }, { "VD1", 1.25, 1.5, -3 }, { "VD1", 1.25, 1.25, 3 }, 
          { "VD1", 4.25, -0.25, -0.75 } } ,
        { { "ATM", -2.405, 16.081, 13.767 }, { "ATM", -1.235, 16.595, 12.999 }, { "ATM", 0.011, 16.688, 13.863 }, 
          { "ATM", -0.059, 16.986, 15.056 }, { "ATM", -1.54, 17.971, 12.396 }, { "ATM", -0.385, 18.548, 11.802 } }

 
    {
		"NAG", "nglycan", "B", "D", "N-Acetyl-D-Glucosamine", 15,
		{ { "PK1", 2.229, 0, 1.191 }, { "PK1", 0.732, 0.427, 1.267 }, { "PK1", 0, 0, 0 },
          { "PK1", 0.754, 0.418, -1.258 }, { "PK1", 2.222, -0, -1.187 }, { "PK1", 3.038, 0.523, -2.343 },
          { "PK1", -0.356, 0.655, 3.478 }, { "PK1", 0.04, -0.086, 2.429 }, { "PK1", -1.304, 0.567, 0.022 },
          { "PK1", 0.156, -0.215, -2.388 }, { "PK1", 2.815, 0.523, 0.011 }, { "PK1", -0.12, 1.85, 3.556 },
          { "PK1", 3.024, 0.447, 2.302 } } ,
        { { "VD1", 8.88178e-16, 3, -1.5 }, { "VD1", 2.25, -2.25, 0 }, { "VD1", -2.25, -2.5, -4 },
          { "VD1", 3, 1.75, 5 }, { "VD1", 4, -2.25, 1.25 }, { "VD1", -4, 1, 0 },
          { "VD1", 1.5, -2.5, 5 }, { "VD1", 0.25, 2.75, -4 }, { "VD1", 1.75, -2.5, -2.25 },
          { "VD1", -4.5, -2.5, 2.5 }, { "VD1", 1.25, 3, 1.25 }, { "VD1", -4.25, 4.25, -0.75 }, 
          { "VD1", 4.5, 2.75, 2.5 } } ,
        { { "ATM", 2.229, 0, 1.191 }, { "ATM", 0.732, 0.427, 1.267 }, { "ATM", 0, 0, 0 }, 
          { "ATM", 0.754, 0.418, -1.258 }, { "ATM", 2.222, -0, -1.187 }, { "ATM", 3.038, 0.523, -2.343 }, 
          { "ATM", -0.356, 0.655, 3.478 }, { "ATM", -1.107, -0.066, 4.553 },{ "ATM", 0.04, -0.086, 2.429 }, 
          { "ATM", -1.304, 0.567, 0.022 }, { "ATM", 0.156, -0.215, -2.388 }, { "ATM", 2.815, 0.523, 0.011 }, 
          { "ATM", 3.09, 1.943, -2.323 }, { "ATM", -0.12, 1.85, 3.556 }, { "ATM", 3.024, 0.447, 2.302 } }  
   }, 

    {
        
		"MAN", "nglycan" , "B", "D", "alpha-D-mannopyranose", 12,
		{ { "PK1", 2.196, 0, 1.204 }, { "PK1", 0.737, 0.447, 1.258 }, { "PK1", 0, 0, 0 }, 
          { "PK1", 0.736, 0.434, -1.262 }, { "PK1", 2.201, -0, -1.206 }, { "PK1", 0.619, 1.859, 1.401 }, 
          { "PK1", -1.313, 0.549, 0.002 }, { "PK1", 0.131, -0.151, -2.412 }, { "PK1", 2.808, 0.499, 0.009 }, 
          { "PK1", 2.149, -1.464, 1.199 } },
        { { "VD1", 4.75, -0.5, 1.75 }, { "VD1", -1.25, -2, -1 }, { "VD1", 8.88178e-16, 2.5, -3 }, 
          { "VD1", 2, -0.25, 3.75 }, { "VD1", 5.5, -1.5, -0.5 }, { "VD1", -1, -1.25, 2.25 },
          { "VD1", -0.5, -1.75, -4.5 }, { "VD1", 1.25, 3.25, -1 }, { "VD1", -2.25, 3, 1.75 }, 
          { "VD1", 2, -2.5, -2.25 } }, 
        { { "ATM", 2.196, 0, 1.204 }, { "ATM", 0.737, 0.447, 1.258 }, { "ATM", 0, 0, 0 }, 
          { "ATM", 0.736, 0.434, -1.262 }, { "ATM", 2.201, -0, -1.206 }, { "ATM", 3.019, 0.53, -2.371 }, 
          { "ATM", 0.619, 1.859, 1.401 }, { "ATM", -1.313, 0.549, 0.002 }, { "ATM", 0.131, -0.151, -2.412 }, 
          { "ATM", 2.808, 0.499, 0.009 }, { "ATM", 2.602, 1.847, -2.729 }, { "ATM", 2.149, -1.464, 1.199 } } 

        
        "MAN", "nglycan" , "B", "D", "alpha-D-mannopyranose", 12,
		{ { "PK1", 2.3, 0, 1.171 }, { "PK1", 0.851, 0.508, 1.118 }, { "PK1", 0, 0, 0 }, 
          { "PK1", 0.79, 0.177, -1.262 }, { "PK1", 2.31, 0, -1.176 }, { "PK1", 0.878, 1.969, 1.252 }, 
          { "PK1", -1.19, 0.73, 0.174 }, { "PK1", 0.059, -0.054, -2.442 }, { "PK1", 2.903, 0.293, 0.022 }, 
          { "PK1", 2.259, -1.38, 1.058 } } , 
        { { "VD1", -5, 2.5, 3 }, { "VD1", -5, -4.75, -2.75 }, { "VD1", -0.25, 4.25, 1.75 }, 
          { "VD1", -4.5, -3, 2.25 }, { "VD1", -2.25, 3.75, 3 }, { "VD1", -2, -5.75, 0.25 }, 
          { "VD1", -4, 1, -2.25 }, { "VD1", -2.75, -0.25, 4.25 }, { "VD1", -1, -3, 6 }, 
          { "VD1", 3, -4, -1.25 } } ,
        { { "ATM", 2.3, 0, 1.171 }, { "ATM", 0.851, 0.508, 1.118 }, { "ATM", 0, 0, 0 }, 
          { "ATM", 0.79, 0.177, -1.262 }, { "ATM", 2.31, 0, -1.176 }, { "ATM", 3.259, 1.188, -1.408 }, 
          { "ATM", 0.878, 1.969, 1.252 }, { "ATM", -1.19, 0.73, 0.174 }, { "ATM", 0.059, -0.054, -2.442 }, 
          { "ATM", 2.903, 0.293, 0.022 }, { "ATM", 4.276, 0.924, -2.342 }, { "ATM", 2.259, -1.38, 1.058 } } */

  /*   {
        
		"MAN", "nglycan" , "B", "D", "Alpha-D-Mannopyranose", 12,         ////4C1 //////////////////////////////////////////////////
		{ { "PK1", 2.198, 0, 1.218 }, { "PK1", 0.747, 0.45, 1.245 }, { "PK1", 0, 0, 0 },
          { "PK1", 0.712, 0.474, -1.254 }, { "PK1", 2.17, -0, -1.202 }, { "PK1", 2.989, 0.567, -2.34 }, 
          { "PK1", 0.709, 1.872, 1.383 }, { "PK1", -1.363, 0.473, 0.019 }, { "PK1", 0.145, -0.059, -2.437 }, 
          { "PK1", 2.809, 0.433, 0.008 }, { "PK1", 2.222, -1.414, 1.317 } } ,
        { { "VD1", 2.5, -1.75, -4.75 }, { "VD1", 0.5, -0.5, 6.5 }, { "VD1", 2.5, 1.75, 3.25 },
          { "VD1", -1, 2.25, -2 }, { "VD1", -1.75, -5.25, 1.75 }, { "VD1", 1.75, -3.75, -4.75 },
          { "VD1", -3, 5, -0.5 }, { "VD1", 5, 0.75, 1.25 }, { "VD1", 2.25, 1.5, 5.75 }, 
          { "VD1", -2.75, 2.25, 2 }, { "VD1", 5.25, -0.75, 3.25 } } ,
        { { "ATM", 2.198, 0, 1.218 }, { "ATM", 0.747, 0.45, 1.245 }, { "ATM", 0, 0, 0 }, 
          { "ATM", 0.712, 0.474, -1.254 }, { "ATM", 2.17, -0, -1.202 }, { "ATM", 2.989, 0.567, -2.34 },
          { "ATM", 0.709, 1.872, 1.383 }, { "ATM", -1.363, 0.473, 0.019 }, { "ATM", 0.145, -0.059, -2.437 }, 
          { "ATM", 2.809, 0.433, 0.008 }, { "ATM", 2.981, 2.005, -2.359 }, { "ATM", 2.222, -1.414, 1.317 } }


        
    }

    {
		"BMA", "nglycan" , "B", "D", "beta-D-mannopyranose", 12,
		{ { "PK1", 2.201, 0, 1.208 }, { "PK1", 0.734, 0.434, 1.265 }, { "PK1", 0, 0, 0 }, 
          { "PK1", 0.741, 0.484, -1.238 }, { "PK1", 2.19, 0, -1.202 }, { "PK1", 3.021, 0.511, -2.372 }, 
          { "PK1", -1.334, 0.51, 0.008 }, { "PK1", 2.823, 0.465, 0.004 }, { "PK1", 3.075, 1.936, -2.385 }, 
          { "PK1", 2.856, 0.535, 2.391 } },
        { { "VD1", 1, 3, -1 }, { "VD1", 0.75, 8.88178e-16, 4 }, { "VD1", 5, -1.75, 1.75 },
          { "VD1", 3, -0.5, -4.75 }, { "VD1", 5, 1.75, -0.25 }, { "VD1", 0.5, -2.75, -0.75 }, 
          { "VD1", -3.25, 1.5, -2.25 }, { "VD1", 4, 1.5, 4.5 }, { "VD1", -2.25, -2.5, 5 }, 
          { "VD1", 0.75, 4.5, 2 } },
        { { "ATM", 2.201, 0, 1.208 }, { "ATM", 0.734, 0.434, 1.265 }, { "ATM", 0, 0, 0 }, 
          { "ATM", 0.741, 0.484, -1.238 }, { "ATM", 2.19, 0, -1.202 }, { "ATM", 3.021, 0.511, -2.372 }, 
          { "ATM", 0.648, 1.845, 1.439 }, { "ATM", -1.334, 0.51, 0.008 }, { "ATM", 0.102, -0.003, -2.417 }, 
          { "ATM", 2.823, 0.465, 0.004 }, { "ATM", 3.075, 1.936, -2.385 }, { "ATM", 2.856, 0.535, 2.391 } } 
    }  */
}; 

int sails::data::count_number_of_fingerprints () {
    return sizeof(sails::data::fingerprint_list) / sizeof(sails::data::fingerprint_list[0]);
} 
    