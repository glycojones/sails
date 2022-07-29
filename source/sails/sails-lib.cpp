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
  // This needs rewriting to get something from mon_lib using Privateer
  // and superpose it onto the peaks
    clipper::MMonomer tmp_mon;
    tmp_mon.set_type ( fp.name_short );
    tmp_mon.set_id ( 0 );

    privateer::restraints::CarbohydrateDictionary dict;
    dict.read_from_monlib ( fp.name_short );

    for ( int index = 0; index < fp.num_control_points ; index++ )
    {
        clipper::MAtom tmp_atm;      
        clipper::Coord_orth coords(fp.atoms[index].x, fp.atoms[index].y, fp.atoms[index].z );
        tmp_atm.set_coord_orth ( coords );
        tmp_atm.set_name ( fp.atoms[index].atom_name );
        tmp_atm.set_element ( fp.atoms[index].atom_name[0] );
        tmp_atm.set_id ( fp.peaks[index].atom_name );
        //tmp_atm.set_id ( 1 );
        tmp_atm.set_occupancy ( 1.0 );
        tmp_atm.set_u_iso ( 0.25 ); // aiming for an iso-B of ~20.0
        tmp_mon.insert ( tmp_atm );
    }
    return tmp_mon;
}

clipper::MMonomer sails::get_monomer_dictionary( clipper::MMonomer fp )
{ 
  from_monlib = true;
//   std::string path_to_lib = privateer::restraints::check_monlib_access ( );
  std::string path_to_lib = "/Users/mihaelaatanasova/sails/dependencies/privateer/dependencies/lib/data/monomers/";

  std::stringstream str;
  std::locale loc;
  
  std::string ccd_id;
  ccd_id = fp.type() ;

  std::vector<size_t> atom_indices;
  clipper::MMonomer tmp_mon;

  tmp_mon.set_type ( ccd_id );

  char initial = std::tolower(ccd_id[0],loc);

    if (!path_to_lib.empty()) {
        str << path_to_lib << initial << "/" << ccd_id << ".cif";
        // std::cout << str.str() << std::endl;
        cif_document = gemmi::cif::read_file( str.str() );
        for (gemmi::cif::Block& block : cif_document.blocks){
            chemical_component = gemmi::make_chemcomp_from_block(block);
            gemmi::cif::Table chem_comp_atom = block.find("_chem_comp_atom.",
                                                            {"comp_id", "atom_id", "x", "y", "z"});
            for (auto atom_row : chem_comp_atom) {
                clipper::Coord_orth atom_coords ( std::stof(atom_row[2]), 
                                            std::stof(atom_row[3]), 
                                            std::stof(atom_row[4]) ); 
                clipper::MAtom tmp_atm; 
                tmp_atm.set_coord_orth ( atom_coords );
                tmp_atm.set_name ( atom_row[1] );  //atom_id C1
                tmp_atm.set_element ( atom_row[1][0] ); //atom type = first letter of atom id C
                tmp_atm.set_id ( atom_row[1] ); // atom_id C1
                tmp_atm.set_occupancy ( 1.0 );
                tmp_atm.set_u_iso ( 0.25 ); // aiming for an iso-B of ~20.0
                tmp_mon.insert(tmp_atm); 
            } 
        }
    }
    return tmp_mon;
}

clipper::MMonomer sails::get_peak_monomer ( const sails::data::fingerprint& fp )
{
    clipper::MMonomer tmp_mon;
    tmp_mon.set_type ( fp.name_short );
    tmp_mon.set_id ( 0 );

    for ( int index = 0; index < fp.num_peaks ; index++ )
    //for ( int index = 0; index < fp.num_control_points ; index++ )
    {
        clipper::MAtom tmp_atm;
        clipper::Coord_orth coords(fp.peaks[index].x, fp.peaks[index].y, fp.peaks[index].z );
        tmp_atm.set_coord_orth ( coords );
        tmp_atm.set_name ( fp.peaks[index].atom_name );
        tmp_atm.set_element ( fp.peaks[index].atom_name[0] );
        tmp_atm.set_id ( fp.peaks[index].atom_name );
        //tmp_atm.set_id ( 1 );
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

    for ( int index = 0; index < fp.num_voids ; index++ )
    //for ( int index = 0; index < fp.num_control_points ; index++ )
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

void sails::get_input_codes ( clipper::String input_string, std::vector<clipper::String> &codes) {
  codes = input_string.split ( "," );
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


// clipper::MiniMol sails::build_sugars ( clipper::Xmap<float>& xwrk, sails::data::build_options& options, double step, int nhit )
// {
//     double rad, sigcut;
//     typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Pair_coord;

//     // get cutoff (for optimisation)
//     clipper::Map_stats stats( xwrk );
//     sigcut = stats.mean() + 0.5*stats.std_dev();

//     clipper::MAtom peak_atom, void_atom;
//     clipper::MMonomer peaks, voids, ideal;

//     // to do: get the relevant set of fingerprints
//     // maybe use a std::vector ?
//     ////TODO: add a "build all" option
//     std::string context="ligand";
//     if ( options.nglycans )
//       context = "nglycan";
//     else if ( options.oglycans )
//       context = "oglycan";
 
//     // make a list of rotations
//     std::vector<clipper::RTop_orth> rots;

//     // make a list of rotation ops to try
//     float glim = 360.0;  // gamma
//     float blim = 180.0;  // beta
//     float alim = 360.0;  // alpha

//     // do a uniformly sampled search of orientation space
//     float anglim = clipper::Util::min( alim, glim );

//     for ( float bdeg=step/2; bdeg < 180.0; bdeg += step )
//     {
//         float beta = clipper::Util::d2rad(bdeg);
//         float spl = anglim/clipper::Util::intf(cos(0.5*beta)*anglim/step+1);
//         float smi = anglim/clipper::Util::intf(sin(0.5*beta)*anglim/step+1);
//         for ( float thpl=spl/2; thpl < 720.0; thpl += spl )
//             for ( float thmi=smi/2; thmi < 360.0; thmi += smi )
//             {
//                 float adeg = clipper::Util::mod(0.5*(thpl+thmi),360.0);
//                 float gdeg = clipper::Util::mod(0.5*(thpl-thmi),360.0);

//                 if ( adeg <= alim && bdeg <= blim && gdeg <= glim )
//                 {
//                     float alpha = clipper::Util::d2rad(adeg);
//                     float gamma = clipper::Util::d2rad(gdeg);
//                     clipper::Euler_ccp4 euler( alpha, beta, gamma );
//                     rots.push_back(clipper::RTop_orth(clipper::Rotation(euler).matrix()));
//                 }
//             }
//     }

//     const clipper::String chainid1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
//     const clipper::String chainid2 = "abcdefghijklmnopqrstuvwxyz";
    
//     clipper::MiniMol mol_new( xwrk.spacegroup(), xwrk.cell() );
//     clipper::Grid_sampling grid = xwrk.grid_sampling();

//     std::vector < sails::data::fingerprint > fp_data = sails::data::get_fingerprints_by_context(context);
//     std::cout << "Number of " << context << " fingerprints: " << fp_data.size() << std::endl;

//     for (int i=0; i < fp_data.size(); i++){
//         clipper::MPolymer mprep;
//         std::vector<Pair_coord> all_co;

//         std::cout << "Checking for " << fp_data[i].name_short << std::endl;
//         ideal = sails::get_ideal_monomer ( fp_data.at(i) );
//         peaks = sails::get_peak_monomer  ( fp_data.at(i) );
//         voids = sails::get_void_monomer  ( fp_data.at(i) );
//         mprep.insert ( ideal, -1 );

//         for (int j = 0; j < ideal.size(); j++) {
//             clipper::Coord_orth coord = ideal[j].coord_orth();
//             std::cout << "atom:" << ideal[j].id() << " " << coord.x() << " " << coord.y() << " " << coord.z() << std::endl;
//         }

//         for (int k = 0; k < peaks.size(); k++) {
//             clipper::Coord_orth coord = peaks[k].coord_orth();
//             std::cout << "peak:" << peaks[k].id() << " " << coord.x() << " " << coord.y() << " " << coord.z() << std::endl;
//         }

//         for (int l = 0; l < voids.size(); l++) {
//             clipper::Coord_orth coord = voids[l].coord_orth();
//             std::cout << "void:" << voids[l].id() << " " << coord.x() << " " << coord.y() << " " << coord.z() << std::endl;
//         }

//         // set up targets
//         for ( int r = 0; r < peaks.size(); r++ ) 
//             all_co.push_back( Pair_coord( peaks[r].coord_orth(), voids[r].coord_orth() ) );

//         // get map radius
//         clipper::Atom_list atoms = voids.atom_list();
//         double r2 = 0.0;

//         for ( int a = 0; a < atoms.size(); a++ )
//         {
//             double d2 = atoms[a].coord_orth().lengthsq();
//             if ( d2 > r2 )
//                 r2 = d2;
//         }

//         rad = sqrt( r2 ) + 1.0;

//         SSfind ssfind;
//         ssfind.prep_xmap( xwrk, rad );
//         ssfind.prep_search( xwrk );
//         std::vector<SearchResult> results = ssfind.search( all_co, rots, sigcut, 0.0 );

//         std::sort( results.begin(), results.end() );
//         std::reverse( results.begin(), results.end() );
//         std::cout << results.size() << std::endl;

//         //for ( int j = 0; j < int(results.size()); j++ ) /////////////////////////////
//         //for ( int j = 0; j < clipper::Util::min( int(results.size()), nhit ); j++ )
//             //std::cout << j << " " << results[j].score << " " << results[j].rot << " " << results[j].trn << std::endl; //////////////////
        
//         // initialise vector of placed coorinates
//         std::vector<clipper::Coord_orth> placed_coordinates; 

//         for ( int k = 0; k <  int(results.size()); k++ )
//         //for ( int k = 0; k < clipper::Util::min( int(results.size()), nhit ); k++ )
//         {
//         clipper::String id; int roff;

//         if ( k < 26 )
//         {
//             id = chainid1.substr( k, 1 );
//             roff = 0;
//         }
//         else
//         {
//             id = chainid2.substr( (k-26)/100, 1 );
//             roff = 10*(i%100);
//         }

//         clipper::MPolymer mprot = mprep;
//         int ir = results[k].rot;
//         int it = results[k].trn;

//         // Get coord_orth
//         clipper::Coord_orth coord_result = xwrk.coord_orth(grid.deindex(it).coord_map());

//         // std::cout << "coord_result:" << coord_result.x() << " " << coord_result.y() << " " << coord_result.z() << std::endl;

//         double cutoff=0;
//         std::vector<double> distance_vector;

//         if (k == 0){
//             placed_coordinates.push_back(coord_result);
//             std::cout << "coord_result1:" << coord_result.x() << " " << coord_result.y() << " " << coord_result.z() << std::endl;//////////
//         }
//         // std::cout << "coord_result:" << coord_result.x() << " " << coord_result.y() << " " << coord_result.z() << std::endl;
//         else
//         {
//             // If it's close to any in placed coord vector then continue the results loop        
//             for (int s = 0; s < int(placed_coordinates.size()); s++)
//             {
//             double distance = (coord_result - placed_coordinates[s]).lengthsq();
//             distance_vector.push_back(distance);
//             }
//             if (std::any_of( distance_vector.begin(), distance_vector.end(), [] (int t){return t<0.1;}) )
//             continue;
//             // break;
//             else
//             // If not add it to the placed coord vector and do all stuff below
//             {
//                 placed_coordinates.push_back(coord_result);
//                 std::cout << "coord_result2:" << coord_result.x() << " " << coord_result.y() << " " << coord_result.z() << std::endl;//////////
//             }
//         }
    
//         clipper::RTop_orth rtop( rots[ir].rot(),
//                 placed_coordinates[k]);  

//         std::cout << "coord_result3:" << placed_coordinates[k].x() << " " << placed_coordinates[k].y() << " " << placed_coordinates[k].z() << std::endl;//////////
    
              
//         mprot.transform( rtop );
//         mprot.set_id( id );
        
//         for ( int r = 0; r < mprot.size(); r++ )
//         {
//             mprot[r].set_seqnum( roff+mprot[r].seqnum() );
//         }

//         for (int m = 0; m < mprot.size(); m++) 
//         {   
//             clipper::MMonomer mon = mprot.operator[](m);

//             for ( int i = 0; i < mon.size(); i++ )
//             {
//                 std::cout << "mprot atom:" << mon[i].coord_orth().x() << " " << mon[i].coord_orth().y() << " " << mon[i].coord_orth().z() << std::endl;
//             }

//         }

//         mol_new.insert( mprot );
//         }


//     }
//     return mol_new;
// }

// need to change this to return results, then make another function to turn that into minimol
std::vector<SearchResult> sails::sugar_results ( clipper::Xmap<float>& xwrk, sails::data::build_options& options, double step, int nhit )
{
    double rad, sigcut;
    typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Pair_coord;

    // get cutoff (for optimisation)
    clipper::Map_stats stats( xwrk );
    sigcut = stats.mean() + 0.5*stats.std_dev();

    clipper::MAtom peak_atom, void_atom;
    clipper::MMonomer peaks, voids, ideal;

    // to do: get the relevant set of fingerprints
    // maybe use a std::vector ?
    ////TODO: add a "build all" option
    std::string context="ligand";
    if ( options.nglycans )
      context = "nglycan";
    else if ( options.oglycans )
      context = "oglycan";

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

    std::vector < sails::data::fingerprint > fp_data = sails::data::get_fingerprints_by_context(context);
    std::cout << "Number of " << context << " fingerprints: " << fp_data.size() << std::endl;

    std::vector<SearchResult>  results;

    for (int i=0; i < fp_data.size(); i++)
    {
        clipper::MPolymer mprep;
        std::vector<Pair_coord> all_co;

        std::cout << "Checking for " << fp_data[i].name_short << std::endl;
        ideal = sails::get_ideal_monomer ( fp_data.at(i) );
        peaks = sails::get_peak_monomer  ( fp_data.at(i) );
        voids = sails::get_void_monomer  ( fp_data.at(i) );
        mprep.insert ( ideal, -1 );

        // for (int j = 0; j < ideal.size(); j++) {
        //     clipper::Coord_orth coord = ideal[j].coord_orth();
        //     std::cout << "atom:" << ideal[j].id() << " " << coord.x() << " " << coord.y() << " " << coord.z() << std::endl;
        // }

        // for (int k = 0; k < peaks.size(); k++) {
        //     clipper::Coord_orth coord = peaks[k].coord_orth();
        //     std::cout << "peak:" << peaks[k].id() << " " << coord.x() << " " << coord.y() << " " << coord.z() << std::endl;
        // }

        // for (int l = 0; l < voids.size(); l++) {
        //     clipper::Coord_orth coord = voids[l].coord_orth();
        //     std::cout << "void:" << voids[l].id() << " " << coord.x() << " " << coord.y() << " " << coord.z() << std::endl;
        // }

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
    
        SSfind ssfind;
        ssfind.prep_xmap( xwrk, rad );
        ssfind.prep_search( xwrk );
    
        std::string name_sugar = fp_data[i].name_short;
        float bestcut = fp_data[i].bestcut;

        std::vector<SearchResult>  temp_result;   
        temp_result = ssfind.search( all_co, rots, sigcut, 0.0, name_sugar, bestcut );

        for ( int k = 0; k < int(temp_result.size()); k++ ) 
        {
            int ir = temp_result[k].rot;
            temp_result[k].sugar_rot = rots[ir].rot(); 
        }

        std::sort( temp_result.begin(), temp_result.end());
        std::reverse( temp_result.begin(), temp_result.end() );

        results.insert(results.end(), temp_result.begin(), temp_result.end());
    }

    // std::sort( results.begin(), results.end());
    // std::reverse( results.begin(), results.end() );
    std::cout << "Number of results: " << results.size() << std::endl;

    for ( int j = 0; j < int(results.size()); j++ ) 
    {
        //for ( int j = 0; j < clipper::Util::min( int(results.size()), nhit ); j++ )
        std::cout << j << " " << results[j].score << " " << results[j].rot << " " << results[j].trn << " " << results[j].name_short << std::endl; 
    }
    return results;
        
}
        
clipper::MiniMol sails::build_sugars ( std::vector<SearchResult> results, clipper::Xmap<float>& xwrk, sails::data::build_options& options, double step, int nhit ) {
    // initialise vector of placed coorinates
    std::vector<clipper::Coord_orth> placed_coordinates; 
    clipper::MiniMol mol_new( xwrk.spacegroup(), xwrk.cell() );
    clipper::Grid_sampling grid = xwrk.grid_sampling();

    const clipper::String chainid1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const clipper::String chainid2 = "abcdefghijklmnopqrstuvwxyz";

    std::string context="ligand";
    if ( options.nglycans )
      context = "nglycan";
    else if ( options.oglycans )
      context = "oglycan";
    
    std::vector < sails::data::fingerprint > fp_data = sails::data::get_fingerprints_by_context(context);
    
    clipper::MPolymer mprot; 

    for ( int k = 0; k < int(results.size()); k++ ) 
    {
        clipper::String id; int roff;

        std::cout << "k = " << k << std::endl;
        
        if ( k < 26 )
        {
            id = chainid1.substr( k, 1 );
            roff = 1;
        }
        else
        {
            id = chainid2.substr( (k-26)/100, 1 );
            roff = 1;//10*(i%100);
        } 

        // clipper::MPolymer mprot; 
        clipper::MMonomer placed_sugar;
    
        int it = results[k].trn;
        std::string name_sugar = results[k].name_short;

        // Get coord_orth

        clipper::Coord_orth coord_result = xwrk.coord_orth(grid.deindex(it).coord_map());
        std::vector<double> distance_vector;
        
        
        if (k == 0)
        {
            placed_coordinates.push_back(coord_result);
            clipper::RTop_orth rtop( results[k].sugar_rot, coord_result); 
            for (int p = 0; p < int(fp_data.size()); p++)
            {
                if (fp_data[p].name_short == name_sugar)
                {
                    placed_sugar = get_ideal_monomer(fp_data.at(p));
                    mprot.insert ( placed_sugar, -1 );
                    mprot.operator[](mprot.size() - 1).transform( rtop );   
                    std::cout << name_sugar << " placed" << std::endl;
                    break;
                }
            }
        }
        else 
        {
            // If it's close to any in placed coord vector then continue the results loop        
            for (int s = 0; s < int(placed_coordinates.size()); s++)
            {
                double distance = (coord_result - placed_coordinates[s]).lengthsq();
                distance_vector.push_back(distance);        
            }

            // if (std::any_of( distance_vector.begin(), distance_vector.end(), [] (int t){return t<0.1;}) ) 
            if (std::any_of( distance_vector.begin(), distance_vector.end(), [] (int t){return t<0.1;}) ) 
            {
                std::cout << "Too close" << std::endl;
                // break;
            }
            else
            // If not add it to the placed coord vector and do all stuff below
            {
                std::cout << "Not too close" << std::endl;
                placed_coordinates.push_back(coord_result);
                clipper::RTop_orth rtop( results[k].sugar_rot, coord_result);  
                for (int q = 0; q < int(fp_data.size()); q++)
                {
                    if (fp_data[q].name_short == name_sugar)
                    {
                        placed_sugar = get_ideal_monomer(fp_data.at(q));
                        mprot.insert ( placed_sugar, -1 );
                        mprot.operator[](mprot.size() - 1).transform( rtop ); 
                        std::cout << name_sugar << " placed" << std::endl;
                        break;
                    }
                }
            } 
        }
   
        // for ( int r = 0; r < mprot.size(); r++ )
        // {
        //     mprot[r].set_seqnum( roff+mprot[r].seqnum() );
        // }
        
    }

    // getting the dictinary for each sugar
    clipper::MPolymer mol_dict;
    mol_dict.set_id( mprot.id() );

    for (int m = 0; m < mprot.size(); m++) 
    {   
        std::cout << "m = " << m << std::endl;

        std::vector<clipper::Coord_orth> mol_coord;
        std::vector<clipper::Coord_orth> mol_dict_coord;

        clipper::MMonomer mon = mprot.operator[](m);
        clipper::MMonomer mon_dict = sails::get_monomer_dictionary(mon);

        clipper::MiniMol mol( xwrk.spacegroup(), xwrk.cell() );
        const clipper::MAtomNonBond& manb = clipper::MAtomNonBond( mol, 5.0 );

        clipper::MSugar ms(mol, mon, manb);
        clipper::MSugar ms_dict(mol, mon_dict, manb);

        int an[3];
        int an_dict[3];

        std::vector<clipper::MAtom> ring_atoms = ms.ring_members();
        std::vector<clipper::MAtom> ring_atoms_dict = ms_dict.ring_members();

        std::vector<clipper::String> atomv; // create a vector of atoms 
        std::vector<clipper::String> atomv_dict;

        atomv.clear();
        atomv.push_back( ring_atoms[1].name().trim() ); // first carbon
        atomv.push_back( ring_atoms[3].name().trim() ); // third carbon
        atomv.push_back( ring_atoms[5].name().trim() ); // fifth carbon

        atomv_dict.clear();
        atomv_dict.push_back( ring_atoms_dict[1].name().trim() ); // first carbon
        atomv_dict.push_back( ring_atoms_dict[3].name().trim() ); // third carbon
        atomv_dict.push_back( ring_atoms_dict[5].name().trim() ); // fifth carbon

        for ( int i = 0; i < atomv.size(); i++ )
        {
            an[i] = mon.lookup( atomv[i], clipper::MM::ANY );
            mol_coord.push_back(mon[an[i]].coord_orth());
            std::cout << "mol_coord:" << mon[an[i]].coord_orth().x() << " " << mon[an[i]].coord_orth().y() << " " << mon[an[i]].coord_orth().z() << std::endl;
            std::cout << "mol_coord.size:" << mol_coord.size() << std::endl;
        }

        for ( int j = 0; j < atomv_dict.size(); j++ )
        {
            an_dict[j] = mon_dict.lookup( atomv_dict[j], clipper::MM::ANY ); 
            mol_dict_coord.push_back(mon_dict[an_dict[j]].coord_orth());
            std::cout << "mol_dict_coord:" << mon_dict[an_dict[j]].coord_orth().x() << " " << mon_dict[an_dict[j]].coord_orth().y() << " " << mon_dict[an_dict[j]].coord_orth().z() << std::endl;
            std::cout << "mol_dict_coord.size:" << mol_dict_coord.size() << std::endl;
        }
            
        clipper::RTop_orth rtop_dict(mol_dict_coord, mol_coord);       
        mon_dict.transform(rtop_dict);
        mol_dict.insert ( mon_dict, -1 );  
        std::cout << "mol_dict.size: " << mol_dict.size() << std::endl;

        for ( int f = 0; f < mol_dict.size(); f++ )
        {
            mol_dict[f].set_seqnum( f );
            // mol_dict[f].set_seqnum( mprot[f].seqnum() );
        }  
  
    }
        
    mol_new.insert( mol_dict );   
    // mol_new.insert( mprot );

    return mol_new;
}

void sails::initialise_fingerprints () {
  sails::data::fingerprint_list = {    
    {
      "ARA", "ligand", 10, 10, 10, 0,
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

    // {
  	//   "BMA", "nglycan" , 12, {
    //     { "PK1", 2.201, 0, 1.208 }, { "PK1", 0.734, 0.434, 1.265 }, { "PK1", 0, 0, 0 },
    //     { "PK1", 0.741, 0.484, -1.238 }, { "PK1", 2.19, 0, -1.202 }, { "PK1", 3.021, 0.511, -2.372 },
    //     { "PK1", -1.334, 0.51, 0.008 }, { "PK1", 2.823, 0.465, 0.004 }, { "PK1", 3.075, 1.936, -2.385 },
    //     { "PK1", 2.856, 0.535, 2.391 }
    //   }, {
    //     { "VD1", 1, 3, -1 }, { "VD1", 0.75, 8.88178e-16, 4 }, { "VD1", 5, -1.75, 1.75 },
    //     { "VD1", 3, -0.5, -4.75 }, { "VD1", 5, 1.75, -0.25 }, { "VD1", 0.5, -2.75, -0.75 },
    //     { "VD1", -3.25, 1.5, -2.25 }, { "VD1", 4, 1.5, 4.5 }, { "VD1", -2.25, -2.5, 5 },
    //     { "VD1", 0.75, 4.5, 2 }
    //   }, {
    //     { "ATM", 2.201, 0, 1.208 }, { "ATM", 0.734, 0.434, 1.265 }, { "ATM", 0, 0, 0 },
    //     { "ATM", 0.741, 0.484, -1.238 }, { "ATM", 2.19, 0, -1.202 }, { "ATM", 3.021, 0.511, -2.372 },
    //     { "ATM", 0.648, 1.845, 1.439 }, { "ATM", -1.334, 0.51, 0.008 }, { "ATM", 0.102, -0.003, -2.417 },
    //     { "ATM", 2.823, 0.465, 0.004 }, { "ATM", 3.075, 1.936, -2.385 }, { "ATM", 2.856, 0.535, 2.391 }
    //   },
    // },


{
		"NAG", "nglycan", 15, 15,15,0.2,//0.2, //3bwh maskradius 	2.3   mapradius 	9.5
		{
{ "C8", -0.856, -0.832, 4.587 }, 
{ "C7", -0.283, 0.076, 3.541 }, 
{ "O7", -0.159, 1.284, 3.743 }, 
{ "N2", 0.078, -0.479, 2.369 }, 
{ "C2", 0.721, 0.291, 1.313 }, 
{ "C3", 0, 0, 0 }, 
{ "C4", 0.764, 0.473, -1.23 }, 
{ "C5", 2.215, 0, -1.185 }, 
{ "O5", 2.812, 0.538, 0.013 }, 
{ "C1", 2.223, 0, 1.189 }, 
{ "C6", 3.053, 0.454, -2.36 }, 
{ "O6", 3.012, 1.885, -2.474 }, 
{ "O4", 0.147, -0.095, -2.389 }, 
{ "O3", -1.285, 0.628, 0.047 }, 
{ "O4", 2.849, 0.423, 2.453 }}, 
{{ "VD1", 4.5, -5, 3.5 }, 
{ "VD1", 4.75, -2.75, 4.5 }, 
{ "VD1", 1.5, -1.75, -3.25 }, 
{ "VD1", -2.25, 4.25, -4.5 }, 
{ "VD1", 1.25, 3.75, -1.25 }, 
{ "VD1", 2, -1.75, 3 }, 
{ "VD1", -1.75, -1.5, -0.75 }, 
{ "VD1", -2.5, 2, -1.5 }, 
{ "VD1", 6.5, -1.75, 1.25 }, 
{ "VD1", 4.75, 1.75, 0.5 }, 
{ "VD1", 2.5, -5, -0.5 }, 
{ "VD1", 4.5, -3.75, -4.75 }, 
{ "VD1", 3.5, -6.5, 2.25 }, 
{ "VD1", 0.5, 2.25, -2.75 }, 
{ "VD1", 0.75, 0.25, -6.5 }}, 
{{ "C8", -0.856, -0.832, 4.587 }, 
{ "C7", -0.283, 0.076, 3.541 }, 
{ "O7", -0.159, 1.284, 3.743 }, 
{ "N2", 0.078, -0.479, 2.369 }, 
{ "C2", 0.721, 0.291, 1.313 }, 
{ "C3", 0, 0, 0 }, 
{ "C4", 0.764, 0.473, -1.23 }, 
{ "C5", 2.215, 0, -1.185 }, 
{ "O5", 2.812, 0.538, 0.013 }, 
{ "C1", 2.223, 0, 1.189 }, 
{ "C6", 3.053, 0.454, -2.36 }, 
{ "O6", 3.012, 1.885, -2.474 }, 
{ "O4", 0.147, -0.095, -2.389 }, 
{ "O3", -1.285, 0.628, 0.047 }, 
{ "O4", 2.849, 0.423, 2.453 }}}, 


	{
		"FUL", "nglycan", 11,11,11,0.7,//0.3,  //7c38_refmac-FUL-1c4.pdb maskradius 	2.5 mapradius 	11
		{
{ "C6", 3.051, -0.42, -2.307 }, 
{ "C5", 2.148, -0, -1.183 }, 
{ "C4", 0.728, -0.498, -1.265 }, 
{ "C3", 0.748, -0.441, 1.248 }, //
{ "C2", 2.176, -0, 1.198 }, //
{ "C1", 2.794, -0.501, 0.009 }, //
{ "O5", 2.871, -0.457, 2.294 }, //
{ "O1", 0.113, 0.162, 2.393 }, //
{ "O2", -1.362, -0.444, 0.033 }, 
{ "O3", 0.81, -1.933, -1.344 }}, 
// { "O1", 2.871, -0.457, 2.294 }, 
{{ "VD1", 0.5, 2.5, -6.25 }, 
{ "VD1", -6.25, -2.5, 3.5 }, 
{ "VD1", -3.25, -3.75, -0.75 }, 
{ "VD1", 2, -5, 1 }, 
{ "VD1", -2, -6.5, 1.5 }, 
{ "VD1", 1, 3.75, -1.5 }, 
{ "VD1", -1, -4, 1.75 }, 
{ "VD1", -5.5, 0.25, -2.25 }, 
{ "VD1", 5.75, -4, 0.5 }, 
{ "VD1", -4.75, -4.25, 3.25 }, 
{ "VD1", -1.25, 5.75, -5.75 }}, 
{{ "C6", 3.051, -0.42, -2.307 }, 
{ "C5", 2.148, -0, -1.183 }, 
{ "C4", 0.728, -0.498, -1.265 }, 
{ "O4", 0, 0, 0 }, 
{ "C3", 0.748, -0.441, 1.248 }, 
{ "C2", 2.176, -0, 1.198 }, 
{ "C1", 2.794, -0.501, 0.009 }, 
{ "O5", 2.871, -0.457, 2.294 }, 
{ "O1", 0.113, 0.162, 2.393 }, 
{ "O2", -1.362, -0.444, 0.033 }, 
{ "O3", 0.81, -1.933, -1.344 }}}, 
// { "O1", 2.871, -0.457, 2.294 }, 


{
		"NDG", "nglycan", 15,11,11,0.7,//0.4, //1q4g maskradius 	2.5 mapradius 	9.5 SHIT
		{
{ "C7", 0.611, 0.787, 3.635 }, 
{ "N2", 0.048, 0.456, 2.46 }, 
{ "C2", 0.738, 0.594, 1.192 }, 
{ "C3", 0, 0, 0 }, 
{ "C4", 0.721, 0.422, -1.28 }, 
{ "C5", 2.192, 0, -1.217 }, 
{ "O5", 2.798, 0.456, 0.013 }, 
{ "C1", 2.148, 0, 1.193 }, 
{ "C6", 3.027, 0.578, -2.344 }, 
{ "O4", 0.133, -0.188, -2.426 }, 
{ "O3", -1.365, 0.391, -0.001 }}, 
{{ "VD1", 0.25, -2.75, 5 }, 
{ "VD1", 3, -8.88178e-16, 5 }, 
{ "VD1", 0.25, -2.75, 0.25 }, 
{ "VD1", 5, -0.25, 1 }, 
{ "VD1", -2, 2, 5.5 }, 
{ "VD1", 3.25, -5.25, 0.75 }, 
{ "VD1", -5.25, -8.88178e-16, 1.25 }, 
{ "VD1", 4.5, 1.5, -4.5 }, 
{ "VD1", -1.25, 4.75, -1.75 }, 
{ "VD1", -2.75, -2.75, 3 }, 
{ "VD1", 3.25, -3.5, 5 }}, 
{{ "C8", 0.248, -0.053, 4.822 }, 
{ "C7", 0.611, 0.787, 3.635 }, 
{ "O7", 1.4, 1.719, 3.748 }, 
{ "N2", 0.048, 0.456, 2.46 }, 
{ "C2", 0.738, 0.594, 1.192 }, 
{ "C3", 0, 0, 0 }, 
{ "C4", 0.721, 0.422, -1.28 }, 
{ "C5", 2.192, 0, -1.217 }, 
{ "O5", 2.798, 0.456, 0.013 }, 
{ "C1", 2.148, 0, 1.193 }, 
{ "C6", 3.027, 0.578, -2.344 }, 
{ "O6", 3.019, -0.268, -3.49 }, 
{ "O4", 0.133, -0.188, -2.426 }, 
{ "O3", -1.365, 0.391, -0.001 }, 
{ "O4", 2.233, -1.424, 1.113 }}}, 



{
		"FUC", "nglycan", 11,10,10,0.7,//0.45, //7c38 maskradius 	3 mapradius 	9.5
		{
{ "C6", 2.996, -0.475, -2.351 }, 
{ "C5", 2.15, 0, -1.197 }, 
{ "C4", 0, 0, 0 }, ///
{ "C3", 2.186, 0, 1.217 },  /////
{ "C2", 2.815, -0.418, 0.007 }, ///
{ "C1", 2.228, 1.388, 1.338 }, ////
{ "O5", 0.111, 0.052, 2.425 },  ////
{ "O1", -1.37, -0.422, 0.045 }, ///
{ "O2", 0.824, -1.974, -1.179 }}, //
// { "O1", 2.228, 1.388, 1.338 }, 
{{ "VD1", 2.5, -4.5, 1.25 }, 
{ "VD1", -0.5, -4.25, 0.75 }, 
{ "VD1", 6, -4, 0 }, 
{ "VD1", 0.75, -5.25, 4 }, 
{ "VD1", 3.75, -5.5, -1 }, 
{ "VD1", -4, 1, -5.75 }, 
{ "VD1", -2, -5, -1.75 }, 
{ "VD1", -4.75, -1, 2.75 }, 
{ "VD1", -5, 0.75, 0.5 }, 
{ "VD1", 2.25, -3.25, 3.5 }}, 
{{ "C6", 2.996, -0.475, -2.351 }, 
{ "C5", 2.15, 0, -1.197 }, 
{ "O3", 0.727, -0.539, -1.216 }, 
{ "O4", 0, 0, 0 }, 
{ "C3", 0.752, -0.418, 1.24 }, 
{ "C3", 2.186, 0, 1.217 }, 
{ "C2", 2.815, -0.418, 0.007 }, 
{ "C1", 2.228, 1.388, 1.338 }, 
{ "O5", 0.111, 0.052, 2.425 }, 
{ "O1", -1.37, -0.422, 0.045 }, 
{ "O2", 0.824, -1.974, -1.179 }}}, 
// { "O1", 2.228, 1.388, 1.338 }, 



{
		"GLC", "nglycan", 12, 11, 11, 0.7,//0.45,//maskradius 	2.3 mapradius 	6 3weo score 0.5 OK
		{
{ "C6", 3.038, 0.481, -2.34 }, 
{ "C5", 2.156, 0, -1.189 }, 
{ "C4", 0.704, 0.32, -1.314 }, 
{ "C3", 0, 0, 0 }, 
{ "C2", 0.791, 0.482, 1.225 }, 
{ "C1", 2.166, 0, 1.194 }, 
{ "O5", 2.617, 0.679, -0.001 }, 
{ "O2", 0.248, 0.084, 2.475 }, 
{ "O3", -1.241, 0.628, -0.086 }, 
{ "O4", 0.117, -0.401, -2.426 }, 
{ "O4", 2.207, -1.411, 1.53 }}, 
{{ "VD1", 0.5, -2.25, -1 }, 
{ "VD1", 0.75, 3.25, -2.5 }, 
{ "VD1", -3, 1, -2.25 }, 
{ "VD1", 0.25, 2.5, 0 }, 
{ "VD1", 1, -3.5, 1.5 }, 
{ "VD1", -0.75, 2.25, 2 }, 
{ "VD1", -1.25, -2, 0.25 }, 
{ "VD1", -1, -1, 4.25 }, 
{ "VD1", 2, -1.75, -2.75 }, 
{ "VD1", -3.25, -0.75, 2.25 }, 
{ "VD1", 2.25, 2.25, 1.75 }}, 
{{ "C6", 3.038, 0.481, -2.34 }, 
{ "C5", 2.156, 0, -1.189 }, 
{ "C4", 0.704, 0.32, -1.314 }, 
{ "C3", 0, 0, 0 }, 
{ "C2", 0.791, 0.482, 1.225 }, 
{ "C1", 2.166, 0, 1.194 }, 
{ "O5", 2.617, 0.679, -0.001 }, 
{ "O2", 0.248, 0.084, 2.475 }, 
{ "O3", -1.241, 0.628, -0.086 }, 
{ "O4", 0.117, -0.401, -2.426 }, 
{ "O6", 2.585, 1.786, -2.86 }, 
{ "O4", 2.207, -1.411, 1.53 }}},


{
		"MAN", "nglycan", 12,11,11, 0.5,//score 0.5 maskradius 2.4, mapradius 9 5o2x
		{
{ "C6", 2.991, 0.556, -2.349 }, 
{ "C5", 2.193, -0, -1.205 }, 
{ "C4", 0.728, 0.459, -1.257 }, 
{ "C3", 0, 0, 0 }, 
{ "C2", 0.756, 0.452, 1.233 }, 
{ "C1", 2.208, 0, 1.214 }, 
{ "O5", 2.818, 0.439, 0.001 }, 
{ "O2", 0.71, 1.874, 1.363 }, 
{ "O3", -1.361, 0.476, 0.008 }, 
{ "O4", 0.166, -0.076, -2.441 }, 
{ "OG", 2.237, -1.415, 1.326 }}, 
{{ "VD1", 2.25, 1.75, 3.25 }, 
{ "VD1", -1, 2, -2 }, 
{ "VD1", 2.5, -2, -4.75 }, 
{ "VD1", -3, 4.75, -0.25 }, 
{ "VD1", 4.5, 0.75, 1.75 }, 
{ "VD1", -0.75, -1.5, 1.75 }, 
{ "VD1", 2.25, 1.5, 5.75 }, 
{ "VD1", 2.75, -2, -2.5 }, 
{ "VD1", -2, -5.25, 2 }, 
{ "VD1", 0.25, -0.5, 6.25 }, 
{ "VD1", -2.75, -0.5, -3 }}, 
{{ "C6", 2.991, 0.556, -2.349 }, 
{ "C5", 2.193, -0, -1.205 }, 
{ "C4", 0.728, 0.459, -1.257 }, 
{ "C3", 0, 0, 0 }, 
{ "C2", 0.756, 0.452, 1.233 }, 
{ "C1", 2.208, 0, 1.214 }, 
{ "O5", 2.818, 0.439, 0.001 }, 
{ "O2", 0.71, 1.874, 1.363 }, 
{ "O3", -1.361, 0.476, 0.008 }, 
{ "O4", 0.166, -0.076, -2.441 }, 
{ "O6", 2.993, 1.987, -2.371 }, 
{ "OG", 2.237, -1.415, 1.326 }}},




{
		"BMA", "nglycan", 12, 11, 11, 0.5,// maskradius 2.4 mapradius 9.5 1pmh 
		{
{ "C6", 3.069, 0.471, -2.304 }, 
{ "C5", 2.228, 0, -1.17 }, 
{ "C4", 0.746, 0.454, -1.201 }, 
{ "C3", 0, 0, 0 }, 
{ "C2", 0.717, 0.391, 1.263 }, 
{ "C1", 2.161, 0, 1.135 }, 
{ "O5", 2.807, 0.488, 0.062 }, 
{ "O2", 0.609, 1.807, 1.481 }, 
{ "O3", -1.367, 0.452, 0.022 }, 
{ "O4", 0.09, -0.152, -2.316 }, 
{ "O4", 2.776, 0.365, 2.355 }}, 
{{ "VD1", 0.25, -2.5, 1 }, 
{ "VD1", -1.75, -1.75, -1 }, 
{ "VD1", 1, -2.75, 4.25 }, 
{ "VD1", 3.75, -1.75, 0.5 }, 
{ "VD1", 0.5, 2.5, -2.75 }, 
{ "VD1", 0.75, 0.25, 3.75 }, 
{ "VD1", -1, 4.5, 4 }, 
{ "VD1", 3.5, 2.5, 3.5 }, 
{ "VD1", -4.75, 3, -2.75 }, 
{ "VD1", 1.75, 0.5, -4.5 }, 
{ "VD1", 1, -5.25, 1 }}, 
{{ "C6", 3.069, 0.471, -2.304 }, 
{ "C5", 2.228, 0, -1.17 }, 
{ "C4", 0.746, 0.454, -1.201 }, 
{ "C3", 0, 0, 0 }, 
{ "C2", 0.717, 0.391, 1.263 }, 
{ "C1", 2.161, 0, 1.135 }, 
{ "O5", 2.807, 0.488, 0.062 }, 
{ "O2", 0.609, 1.807, 1.481 }, 
{ "O3", -1.367, 0.452, 0.022 }, 
{ "O4", 0.09, -0.152, -2.316 }, 
{ "O6", 2.99, 1.907, -2.417 }, 
{ "O4", 2.776, 0.365, 2.355 }}}, 

	{
		"GAL", "nglycan", 12, 11,11,0.7,//0.5, //5elb maskradius 	2.5  mapradius 	9.5
		{
{ "C6", 2.961, 0.557, -2.352 }, 
{ "C5", 2.174, 0, -1.198 }, 
{ "C4", 0.701, 0.42, -1.251 }, 
{ "C3", 0, 0, 0 }, 
{ "C2", 0.775, 0.438, 1.209 }, 
{ "C1", 2.225, 0, 1.227 }, 
{ "O5", 2.848, 0.476, 0.038 }, 
{ "O2", 0.059, -0.03, 2.34 }, 
{ "O3", -1.279, 0.633, 0.057 }, 
{ "O4", 0.564, 1.853, -1.329 }, 
{ "O4", 2.858, 0.773, 2.223 }}, 
{{ "VD1", 8.88178e-16, 5.75, 0.5 }, 
{ "VD1", 2.5, 3.75, 4.75 }, 
{ "VD1", -0.25, 4.5, -4.75 }, 
{ "VD1", 3, 3.75, -2.75 }, 
{ "VD1", -2, 2.75, 1.25 }, 
{ "VD1", -4.5, 5.25, 0.25 }, 
{ "VD1", -3.25, 1.25, -4.25 }, 
{ "VD1", 0.25, 4, 4.5 }, 
{ "VD1", 5.25, 0.25, 1 }, 
{ "VD1", -0.75, 4, -2 }, 
{ "VD1", 1, -2.25, 1.25 }}, 
{{ "C6", 2.961, 0.557, -2.352 }, 
{ "C5", 2.174, 0, -1.198 }, 
{ "C4", 0.701, 0.42, -1.251 }, 
{ "C3", 0, 0, 0 }, 
{ "C2", 0.775, 0.438, 1.209 }, 
{ "C1", 2.225, 0, 1.227 }, 
{ "O5", 2.848, 0.476, 0.038 }, 
{ "O2", 0.059, -0.03, 2.34 }, 
{ "O3", -1.279, 0.633, 0.057 }, 
{ "O4", 0.564, 1.853, -1.329 }, 
{ "O6", 4.247, -0.072, -2.338 }, 
{ "O4", 2.858, 0.773, 2.223 }}}, 


{
		"BGC", "nglycan", 12,12,12,0.7,//0.5, //3pfz maskradius 	2.5 mapradius 	9.5
		{
{ "C6", 3.037, 0.523, -2.328 }, 
{ "C5", 2.197, 0, -1.193 }, 
{ "C4", 0.744, 0.43, -1.249 }, 
{ "C3", 0, 0, 0 }, 
{ "C2", 0.71, 0.402, 1.242 }, 
{ "C1", 2.156, 0, 1.17 }, 
{ "O5", 2.77, 0.532, 0.025 }, 
{ "O1", 2.773, 0.411, 2.335 }, 
{ "O2", 0.028, -0.222, 2.337 }, 
{ "O3", -1.333, 0.513, 0.037 }, 
{ "O4", 0.107, -0.188, -2.339 }}, 
// { "O1", 2.773, 0.411, 2.335 }, 
{{ "VD1", 0.25, 6, -0.5 }, 
{ "VD1", 8.88178e-16, 3.75, -1.75 }, 
{ "VD1", -3, 2.5, -3.25 }, 
{ "VD1", 1.25, 2, -5.75 }, 
{ "VD1", -3.25, 3, 0.5 }, 
{ "VD1", -2.75, 1.75, 4.5 }, 
{ "VD1", 1.25, -1.5, 4.75 }, 
{ "VD1", -1.75, -1.75, -1 }, 
{ "VD1", 6.5, 2, 1 }, 
{ "VD1", 1.25, -0.5, -6.5 }, 
{ "VD1", 3.5, 3, 3 }, 
{ "VD1", 4.25, -3.75, -4.25 }}, 
{{ "C6", 3.037, 0.523, -2.328 }, 
{ "C5", 2.197, 0, -1.193 }, 
{ "C4", 0.744, 0.43, -1.249 }, 
{ "C3", 0, 0, 0 }, 
{ "C2", 0.71, 0.402, 1.242 }, 
{ "C1", 2.156, 0, 1.17 }, 
{ "O5", 2.77, 0.532, 0.025 }, 
{ "O1", 2.773, 0.411, 2.335 }, 
{ "O2", 0.028, -0.222, 2.337 }, 
{ "O3", -1.333, 0.513, 0.037 }, 
{ "O4", 0.107, -0.188, -2.339 }, 
{ "O6", 2.858, 1.962, -2.512 }}}, 
// { "O1", 2.773, 0.411, 2.335 }, 





{ 
		"XYP", "nglycan", 10,9,9,1,//0.65, //5lal++ maskradius 	2 mapradius 	8
		{
{ "C5", 2.177, -0, -1.197 }, 
{ "C4", 0.734, 0.449, -1.254 }, 
{ "C3", 0, 0, 0 }, 
{ "C2", 0.741, 0.432, 1.261 }, 
{ "C1", 2.185, -0, 1.201 }, 
{ "O5", 2.815, 0.493, -0.011 }, 
{ "O2", 0.112, -0.127, 2.411 }, 
{ "O4", 0.104, -0.158, -2.386 }, 
{ "O2", 2.921, 0.446, 2.168 }}, 
{{ "VD1", 0.75, 2.5, 0.75 }, 
{ "VD1", -1.75, -2.5, 2.5 }, 
{ "VD1", -4, -0.25, 1.25 }, 
{ "VD1", -1.25, 0.75, 5 }, 
{ "VD1", 3, -2.25, -0.5 }, 
{ "VD1", -3.25, 2.75, -1.5 }, 
{ "VD1", -0.5, -2, -0.5 }, 
{ "VD1", 0.75, -0.75, -4.25 }, 
{ "VD1", 0.5, -4.5, 2.25 }}, 
{{ "C5", 2.177, -0, -1.197 }, 
{ "C4", 0.734, 0.449, -1.254 }, 
{ "C3", 0, 0, 0 }, 
{ "C2", 0.741, 0.432, 1.261 }, 
{ "C1", 2.185, -0, 1.201 }, 
{ "O5", 2.815, 0.493, -0.011 }, 
{ "O2", 0.112, -0.127, 2.411 }, 
{ "O3", -1.32, 0.522, -0.024 }, 
{ "O4", 0.104, -0.158, -2.386 }, 
{ "O2", 2.921, 0.446, 2.168 }}}, 

    };
}


int sails::data::number_of_fingerprints () {
    return fingerprint_list.size();
}
