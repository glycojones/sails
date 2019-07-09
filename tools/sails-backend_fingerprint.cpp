
// Clipper app that builds a fingerprint from min/max maps and a target model
// 2013 Jon Agirre & Kevin Cowtan @ University of York
// Contains portions of the GNU Scientific Library Exponential fit by non-linear least squares
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk
//

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit_nlin.h>
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-contrib.h>

struct data
{
  size_t n;
  float * distance;
  float * density_diff;
};

void print_state (size_t iter, gsl_multifit_fdfsolver * s);      // partial definitions, full code is after main()
int expb_f (const gsl_vector * x, void *data, gsl_vector * f);
int expb_df (const gsl_vector * x, void *data, gsl_matrix * J);
int expb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J);
double expofit (size_t n_obs, size_t n_params, data& input_data );
float bin_distance (float distance); // returns a the closest bin for the input distance

int main( int argc, char** argv )
{
    CCP4Program prog( "Fingerprinting engine", "0.1", "$Date: 2013/08/29" );
    prog.set_termination_message( "Failed" );


    // command input
    CCP4CommandInput args( argc, argv, true );

    float ipradius = 2.5;
    float ipmapradius = 9.5;
    int arg, dump;
    arg = dump = 0;
    clipper::String title, ippdb, ipmaxmap, ipminmap, oppdb,sugar;
    ippdb = ipmaxmap = ipminmap = oppdb = sugar = "NONE";
    title = "FINGERPRINT";
    typedef clipper::NXmap<float>::Map_reference_index MRI;
    std::ofstream dataFile;

    while ( ++arg < args.size() )
    {
        if ( args[arg] == "-title" )
        {
            if ( ++arg < args.size() )
                title = args[arg];
        }
        else if ( args[arg] == "-mapin1" )
        {
            if ( ++arg < args.size() )
                ipmaxmap = args[arg];
        }
        else if ( args[arg] == "-mapin2" )
        {
            if ( ++arg < args.size() )
                ipminmap = args[arg];
        }
        else if ( args[arg] == "-pdbin" )
        {
        if ( ++arg < args.size() )
            ippdb = args[arg];
        }
        else if ( args[arg] == "-maskradius" )
        {
            if ( ++arg < args.size() )
                ipradius = clipper::String(args[arg]).f();
        }
        else if ( args[arg] == "-mapradius" )
        {
            if ( ++arg < args.size() )
                ipmapradius = clipper::String(args[arg]).f();
        }
        else if ( args[arg] == "-pdbout" )
        {
            if ( ++arg < args.size() )
                oppdb = args[arg];
        }
        else if ( args[arg] == "-dump" )
        {
            dump = 1;
        }
        else if ( args[arg] == "-sugar" )
        {
            if ( ++arg < args.size() ) sugar = args[arg];
        }
        else
        {
            std::cout << "\nUnrecognised:\t" << args[arg] << std::endl;
            args.clear();
        }
    }

    if ( args.size() <= 1 )
    {
        std::cout << "\nUsage: fingerprint\n\t-mapin1 <max map>\n\t-mapin2 <min map>\t\t"
                  << "COMPULSORY\n\t-pdbin <filename>\t\tCOMPULSORY\n\t-sugar <3-letter code>\t\t"
                  << "COMPULSORY\n\t-pdbout <filename>\t\tCOMPULSORY\n\t-maskradius <A>\t\t\t"
                  << "default value: 2.5\n\t-mapradius <A>\t\t\tdefault value: 9.5\n\t-dump\n\n";
        exit(1);
    }

    std::string namePDB = ippdb.substr(0, ippdb.find(".pdb"));

    if (dump == 1)
    {
        std::cout << std::endl << "Dumping internal maps and masks for debug..." << std::endl << std::endl;
    }

    // Get reference model
    clipper::MiniMol mol;
    clipper::MMDBManager mmdb;
    clipper::MMDBfile mfile;
    const int mmdbflags = mmdb::MMDBF_IgnoreBlankLines | mmdb::MMDBF_IgnoreDuplSeqNum | mmdb::MMDBF_IgnoreNonCoorPDBErrors | mmdb::MMDBF_IgnoreRemarks;// | MMDBF_EnforceUniqueChainID;
    mfile.SetFlag( mmdbflags );
    mfile.read_file( ippdb );
    mfile.import_minimol( mol );

    std::cout << mol.spacegroup().symbol_hm() << " " << mol.cell().format() << " " << mol.atom_list().size() << std::endl;


    if ( sugar == "NONE" ) {
      if ( mol.model().size() > 0 && mol.model()[0].size() > 0 ) {
        sugar = mol[0][0].type();
        std::cout << "As no three-letter code was provided (-sugar), " << sugar
                  << " has been inferred from the content of the PDB file" << std::endl;
      }
      else {
        // no code and nothing available from the PDB file? we shall not proceed
        prog.set_termination_message( "Failed" );
        return(1);
      }
    }

    clipper::NXmap<float> max_map;
    clipper::CCP4MAPfile mapfile;
    mapfile.open_read(ipmaxmap);
    mapfile.import_nxmap(max_map);
    mapfile.close_read();

    clipper::NXmap<float> min_map;
    clipper::CCP4MAPfile mapfile2;
    mapfile2.open_read(ipminmap);
    mapfile2.import_nxmap(min_map);
    mapfile2.close_read();

    const double radius = 6.0;
    const double sampling = 0.25; // was 0.25;
    double extent = radius + 2.0 * sampling;  // model box size
    double ng = rint( extent / sampling ); // number of grid points, round to integer

    clipper::Grid nxgrid( 2*int(ng)+1, 2*int(ng)+1, 2*int(ng)+1 );
    clipper::RTop<> rtop( clipper::Mat33<>( 1.0/sampling, 0.0, 0.0,
                                            0.0, 1.0/sampling, 0.0,
                                            0.0, 0.0, 1.0/sampling ),
                          clipper::Vec3<>( ng, ng, ng ) );

    clipper::NXmap<float> mask( nxgrid, rtop );

    clipper::EDcalc_mask<float> masker( ipradius );
    masker(mask, mol[0].atom_list());

    if (dump == 1)
    {
        mapfile.open_write( namePDB + ".mask" );
        mapfile.export_nxmap( mask );
        mapfile.close_write();
    }

    int n_points = 0;
    float mean = 0.0;
    float std_dev = 0.0;

    for ( MRI ix = max_map.first(); !ix.last(); ix.next() )  // masked region mean calculation
    {
        if (mask[ix] == 1.0)
        {
            mean += (max_map[ix]);
            n_points++;
        }
    }

    mean /= n_points;

    for ( MRI ix = max_map.first(); !ix.last(); ix.next() )  // std_dev calculation for the masked region
    {
        if (mask[ix] == 1.0)
        {
            std_dev += pow(max_map[ix] - mean, 2);
        }
    }

    std_dev = sqrt(std_dev/n_points);

    std::cout << "Mean (masked region): " << mean << " ; Standard deviation (masked region): " << std_dev << std::endl;

    for ( MRI ix = max_map.first(); !ix.last(); ix.next() )
    {
        if (mask[ix] == 1.0)
            max_map[ix] = ((max_map[ix] - mean) / std_dev); // normalise values inside the mask
    }

    /////////////////////// check ///////////////////////

    n_points = 0;
    mean = 0.0;
    std_dev = 0.0;

    for ( MRI ix = max_map.first(); !ix.last(); ix.next() )  // masked region mean calculation
    {
        if (mask[ix] == 1.0)
        {
            mean += (max_map[ix]);
            n_points++;
        }
    }

    mean /= n_points;

    for ( MRI ix = max_map.first(); !ix.last(); ix.next() )  // std_dev calculation for the masked region
    {
        if (mask[ix] == 1.0)
        {
            std_dev += pow(max_map[ix] - mean, 2);
        }
    }

    std_dev = sqrt(std_dev/n_points);

    std::cout << "\nMean (masked region, normalised): " << mean << " ; Standard deviation (masked region, normalised): " << std_dev << std::endl;

    ////////////////// end of check /////////////////

    if (dump == 1)
    {
        mapfile.open_write( namePDB + "-norm-masked.map" );
        mapfile.export_nxmap( max_map );
        mapfile.close_write();
    }

    int n_distances;

    for ( n_distances = n_points; n_distances > 0 ; n_points += --n_distances ); // calculate the number of distances

    std::cout << "Number of points: " << n_points << std::endl;

    float *distances, *covariances;
    int *number_of_occurrences;
    struct data exp_fit_data;

    if ( dump == 1)
    {
        distances   = (float *)malloc(n_points * sizeof(float));
        covariances = (float *)malloc(n_points * sizeof(float));

        if (( distances == NULL ) || ( covariances == NULL ))
        {
            std::cout << std::endl << "Unable to allocate enough memory. Exiting..." << std::endl;
            prog.set_termination_message( "Failed" );
        }

        exp_fit_data.n =  n_points;
        exp_fit_data.distance = distances;
        exp_fit_data.density_diff = covariances;
    }
    else
    {
        /////////// create data arrays here
        distances = (float *)malloc(43*sizeof(float));
        covariances = (float *)malloc(43*sizeof(float));
        number_of_occurrences = (int *)malloc(43*sizeof(int));

        for (int i = 0; i < 43 ; i++)  // initialise our data arrays
        {
            distances[i] = 0.0;
            covariances[i] = 0.0;
            number_of_occurrences[i]=0;
        }

        exp_fit_data.n = 43;
        exp_fit_data.distance = distances;
        exp_fit_data.density_diff = covariances;
    }

    int index_exp_fit_data = 0;

    if (dump == 1)
    {
        dataFile.open ("covariance_data");
    }

    int max_index = 0;

    for ( MRI ix = max_map.first(); !ix.last(); ix.next() )  // run across the map
    {
        //ix.next();ix.next();ix.next();ix.next();
        if ( mask[ix] == 1.0 ) // only execute the inner loop if the first point is inside the mask
        {
            for ( MRI iy = ix; !iy.last(); iy.next() )
            {
                if ( mask[iy] == 1.0 )
                {
                    if (dump == 1)
                    {
                        distances[index_exp_fit_data] = bin_distance(clipper::Coord_orth::length(ix.coord_orth(), iy.coord_orth()));
                        covariances[index_exp_fit_data] = 1-( pow(max_map[ix] - max_map[iy], 2) /2.0);
                        dataFile << distances[index_exp_fit_data] << "\t" << covariances[index_exp_fit_data] << std::endl;
                        index_exp_fit_data++;
                    }
                    else
                    {
                        int bin_index = ceil(bin_distance(clipper::Coord_orth::length(ix.coord_orth(), iy.coord_orth()))*3.0);

                        distances[bin_index] = bin_distance(clipper::Coord_orth::length(ix.coord_orth(), iy.coord_orth()));
                        covariances[bin_index] += 1-( pow(max_map[ix] - max_map[iy], 2) /2.0);
                        number_of_occurrences[bin_index]++;

                        if (bin_index > max_index)
                            max_index = bin_index;

                    }
                }
            }
        }
    }

    if (dump == 1)
    {
        dataFile.close();
        exp_fit_data.n = index_exp_fit_data;
    }
    else
    {
        // set the exact number of data points and extract the means
        exp_fit_data.n = max_index;

        for (int i = 0 ; i <= max_index ; i++)
        {
            covariances[i] /= number_of_occurrences[i];
        }
    }

    std::cout << "Number of data (exp fit): " << exp_fit_data.n << "\n";

    double r0 = expofit(exp_fit_data.n, 1, exp_fit_data);

    std::ofstream means;
    means.open("means");

    if (dump == 1)
    {
        for ( int i = 1 ; i < 50 ; i += 2)
        {
            float means_bin = 0.0;
            int n_data = 0;

            for (int j = 0 ; j < index_exp_fit_data ; j++)
            {
                if ( (abs((distances[j] * 6) - i) < 0.1 ) )
                {
                    means_bin += covariances[j];
                    n_data++;
                }
            }

            if (n_data > 0)
                means << i/6.0 << "\t" << means_bin/n_data << std::endl;
        }

        means.close();
    }

    // do the mapcrawling thing here

    std::ofstream coord_out;
    coord_out.open(oppdb.c_str());


    coord_out << "REMARK   0 THIS FILE SERVES AS A FINGERPRINT FOR DETECTING '" << sugar << "' MOLECULES" << std::endl;
    coord_out << "REMARK   0 IN ELECTRON DENSITY MAPS. IT HAS BEEN AUTOMATICALLY GENERATED BASED ON" << std::endl;
    coord_out << "REMARK   0 EXPERIMENTAL DATA AND IS DISTRIBUTED AS PART OF THE PRIVATEER SOFTWARE" << std::endl;
    coord_out << "REMARK   0 " << std::endl;
    coord_out << "REMARK   0 (C) 2013 JON AGIRRE & KEVIN COWTAN (YSBL), UNIVERSITY OF YORK, ENGLAND" << std::endl;
    coord_out << "REMARK   0 " << std::endl;
    coord_out << "REMARK   0 IF YOU WERE TO RE-DISTRIBUTE A MODIFIED VERSION OF THIS FILE, PLEASE" << std::endl;
    coord_out << "REMARK   0 PRESERVE THIS HEADER AND APPEND A SUMMARY OF THE ALTERATIONS USING" << std::endl;
    coord_out << "REMARK   0 THE 'REMARK 0' SECTION" << std::endl;
    coord_out << "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1       " << std::endl;

    int total_atom_count = 0;

    std::ios_base::fmtflags original_flags = std::cout.flags();

    // output the peaks (atoms)

    int npeaks = 0;
    float mean_minmap, stdev_minmap;
    mean_minmap = stdev_minmap = 0.0;

    std::cout << "\t{" << std::endl << "\t\t\"" << sugar << "\", \"nglycan\", " << mol[0][0].size() << ",\n\t\t{" << std::endl;

    for ( int natm = 1 ; natm <= mol[0][0].size() ; natm++) // mean calculation at the atom positions
    {
        mean_minmap += min_map.interp<clipper::Interp_cubic>(min_map.coord_map(mol[0][0][natm-1].coord_orth()));
    }

    mean_minmap /= mol[0][0].size();

    for ( int natm = 1 ; natm <= mol[0][0].size() ; natm++) // standard deviation calculation at the atom positions
    {
        stdev_minmap += pow(min_map.interp<clipper::Interp_cubic>(min_map.coord_map(mol[0][0][natm-1].coord_orth())) - mean_minmap, 2);
    }

    stdev_minmap = sqrt(stdev_minmap / mol[0][0].size());

    for ( int natm = 1 ; natm <= mol[0][0].size() ; natm++)
    {
        if ( min_map.interp<clipper::Interp_cubic>(min_map.coord_map(mol[0][0][natm-1].coord_orth())) > mean_minmap - (stdev_minmap))
        {
            coord_out.flags(original_flags);

            total_atom_count++;

            coord_out << "ATOM";
            coord_out.width(7);
            coord_out.setf( std::ios::right);
            coord_out << total_atom_count;
            coord_out.width(4);
            coord_out.setf( std::ios::right);
            coord_out << mol[0][0][natm-1].id().trim();

            coord_out.width(5);
            coord_out.setf( std::ios::right);
            coord_out << mol[0][0][natm-1].id().trim() << " A ";

            coord_out.width(3);
            coord_out << natm;

            coord_out.setf( std::ios::right);
            coord_out.setf( std::ios::fixed, std::ios::floatfield );
            coord_out.width(12);
            coord_out.precision(3);
            coord_out << mol[0][0][natm-1].coord_orth().x();
            coord_out.setf( std::ios::right);
            coord_out.setf( std::ios::fixed, std::ios::floatfield );
            coord_out.width(8);
            coord_out.precision(3);
            coord_out << mol[0][0][natm-1].coord_orth().y();
            coord_out.setf( std::ios::right);
            coord_out.setf( std::ios::fixed, std::ios::floatfield );
            coord_out.width(8);
            coord_out.precision(3);
            coord_out << mol[0][0][natm-1].coord_orth().z() << "  1.00 30.00       PEAKS" << std::endl;

            // std::cout << "Peak set at " << mol[0][0][natm-1].coord_orth().format() << " with density "
            //          << min_map.interp<clipper::Interp_cubic>(min_map.coord_map(mol[0][0][natm-1].coord_orth()))
            //          << ", assigning number " << natm << std::endl;

            std::cout << "{ \"PK1\", " << mol[0][0][natm-1].coord_orth().x() << " "
                      << mol[0][0][natm-1].coord_orth().y() << " "
                      << mol[0][0][natm-1].coord_orth().z() << " }" << std::endl;

            npeaks++;
        }
        else {}
            // std::cout << "Peak rejected at " << mol[0][0][natm-1].coord_orth().format()
            //          << " with density " << min_map.interp<clipper::Interp_cubic>(min_map.coord_map(mol[0][0][natm-1].coord_orth()))
            //          << ". Mean: " << mean_minmap << " Stdev: " << stdev_minmap << std::endl;

    }

    coord_out << "TER" << std::endl;

    // find and output the voids

    for ( MRI ix = max_map.first(); !ix.last(); ix.next() )  // run across the map, invert values,
    {                                                        // for it was constructed as a difference map.
        if (mask[ix] != 1)                                   // also upscale values far from the center
        {                                                    // to prevent picking artificial peaks
            max_map[ix] *= -1;
        }

        max_map[ix] += pow(clipper::Coord_orth::length(ix.coord_orth(), clipper::Coord_orth( 0.0, 0.0, 0.0 )) / ipmapradius, 6) ;
    }

    float minimum_density = 9999999.9;
    clipper::Coord_orth minimum_spot;
    minimum_spot = clipper::Coord_orth( 0.0, 0.0, 0.0 );

    for ( int natm = 1 ; natm <= npeaks ; natm++)
    {
        minimum_density = 9999999.9;
        minimum_spot = clipper::Coord_orth( 0.0, 0.0, 0.0 );

        for ( MRI ix = max_map.first(); !ix.last(); ix.next() )  // run across the map
        {
            if (mask[ix] != 1)
            {
                if (max_map[ix] < minimum_density )
                {
                    minimum_density = max_map[ix];
                    minimum_spot = ix.coord_orth();
                }
            }
        }

        // std::cout << std::endl << minimum_density << "\t" << minimum_density + exp(-0/r0) << std::endl;

        for ( MRI ix = max_map.first(); !ix.last(); ix.next() )  // run across the map
        {
            if (mask[ix] != 1) // downweight every value by 1 - f(r)
            {
                float dist = clipper::Coord_orth::length(ix.coord_orth(),minimum_spot);
                //std::cout << ix.coord_orth().format() << "\t" << max_map[ix] << "\t" << exp(-dist/r0) << "\t";

                max_map[ix] = max_map[ix] + exp(-dist/r0);
                //std::cout << max_map[ix] << std::endl;
            }
        }

        total_atom_count++;

        //std::cout << "Void found at " << minimum_spot.format() << " with density "
        //          << minimum_density << ", assigning number " << natm << std::endl;

        coord_out.flags(original_flags);

        coord_out << "ATOM";

        coord_out.width(7);
        coord_out.setf( std::ios::right);
        coord_out << total_atom_count;

        coord_out.width(3);
        coord_out.setf( std::ios::right);
        coord_out << natm << "X" ;

        coord_out.width(4);
        coord_out.setf( std::ios::right);
        coord_out << natm << "X" << " B ";

        coord_out.width(3);
        coord_out.setf( std::ios::right);
        coord_out << natm;

        coord_out.setf( std::ios::right);
        coord_out.setf( std::ios::fixed, std::ios::floatfield );
        coord_out.width(12);
        coord_out.precision(3);
        coord_out << minimum_spot.x();
        coord_out.setf( std::ios::right);
        coord_out.setf( std::ios::fixed, std::ios::floatfield );
        coord_out.width(8);
        coord_out.precision(3);
        coord_out << minimum_spot.y();
        coord_out.setf( std::ios::right);
        coord_out.setf( std::ios::fixed, std::ios::floatfield );
        coord_out.width(8);
        coord_out.precision(3);
        coord_out << minimum_spot.z() << "  1.00 30.00       VOIDS" << std::endl;


        std::cout << "{ \"VD1\", " << minimum_spot.x() << " "
                  << minimum_spot.y() << " "
                  << minimum_spot.z() << " }" << std::endl;
    }

    coord_out << "TER" << std::endl;

    // output the model
    // To do: get monomer from library, superpose
    // Or maybe skip this part altogether and superpose at building stage

    int resnumber = total_atom_count+1;

    for ( int natm = 1 ; natm <= mol[0][0].size() ; natm++)
    {
        coord_out.flags(original_flags);
        total_atom_count++;

        coord_out << "ATOM";
        coord_out.width(7);
        coord_out.setf( std::ios::right);
        coord_out << total_atom_count;
        coord_out.width(4);
        coord_out.setf( std::ios::right);
        coord_out << mol[0][0][natm-1].id().trim();

        coord_out.width(5);
        coord_out.setf( std::ios::right);
        coord_out << sugar.c_str() << " C ";

        coord_out.width(3);
        coord_out << resnumber;

        coord_out.setf( std::ios::right);
        coord_out.setf( std::ios::fixed, std::ios::floatfield );
        coord_out.width(12);
        coord_out.precision(3);
        coord_out << mol[0][0][natm-1].coord_orth().x();
        coord_out.setf( std::ios::right);
        coord_out.setf( std::ios::fixed, std::ios::floatfield );
        coord_out.width(8);
        coord_out.precision(3);
        coord_out << mol[0][0][natm-1].coord_orth().y();
        coord_out.setf( std::ios::right);
        coord_out.setf( std::ios::fixed, std::ios::floatfield );
        coord_out.width(8);
        coord_out.precision(3);
        coord_out << mol[0][0][natm-1].coord_orth().z() << "  1.00 30.00       MODEL" << std::endl;

        std::cout << "{ \"ATM\", " << mol[0][0][natm-1].coord_orth().x() << " "
                  << mol[0][0][natm-1].coord_orth().y() << " "
                  << mol[0][0][natm-1].coord_orth().z() << " }" << std::endl;

    }

    coord_out << "END" << std::endl;

    coord_out.close();
    prog.set_termination_message( "Normal termination" );
}


int expb_f (const gsl_vector * param_vector, void *data, gsl_vector * f) // setup the function we want to fit
{ // vector x contains the parameters for the exponential

    size_t n = ((struct data *)data)->n;
    float *distance = ((struct data *)data)->distance;
    float *density_diff = ((struct data *) data)->density_diff;

    double r0 = gsl_vector_get (param_vector, 0);

    #pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        /* Model Yi = exp(-r / r0) */
        double r = (double)distance[i];
        double Yi = exp (-r / r0);
        gsl_vector_set (f, i, (Yi - density_diff[i]));
    }

    return GSL_SUCCESS;
}

int expb_df (const gsl_vector * param_vector, void *data, gsl_matrix * jacobian) // setup the jacobian
{
    size_t n = ((struct data *)data)->n;
    float *density_diff = ((struct data *) data)->density_diff;
    float *distance = ((struct data *)data)->distance;

    float r0 = gsl_vector_get (param_vector, 0);

    #pragma omp parallel for
    for (size_t i = 0; i < n; i++)
    {
        /* Jacobian "matrix" df / dr0, */
        /* where fi = (Yi - yi),       */
        /*       Yi = exp(-r / r0)     */
        /* and r0 is the parameter     */

        double r = (double)distance[i];
        double s = density_diff[i];
        double e = exp(-r / r0);
        gsl_matrix_set (jacobian, i, 0, e * (r/pow(r0,2))); // just one parameter, so just one column in the jacobian
    }

    return GSL_SUCCESS;
}

int expb_fdf (const gsl_vector * param_vector, void *data, gsl_vector * function, gsl_matrix * jacobian) // f * df
{
    expb_f (param_vector, data, function);
    expb_df (param_vector, data, jacobian);

    return GSL_SUCCESS;
}

double expofit (size_t n_obs, size_t n_params, data& input_data )
{
    const gsl_multifit_fdfsolver_type *solver_type;
    gsl_multifit_fdfsolver *solver;
    int status;
    unsigned int i, iter = 0;

    double x_init[1] = { 1.0 };
    gsl_vector_view x = gsl_vector_view_array (x_init, n_params); // x is the initial guess
    const gsl_rng_type * type;

    gsl_multifit_function_fdf f;

    f.f = &expb_f;
    f.df = &expb_df;
    f.fdf = &expb_fdf;
    f.n = n_obs;
    f.p = n_params;
    f.params = &input_data;

    solver_type = gsl_multifit_fdfsolver_lmsder;
    solver = gsl_multifit_fdfsolver_alloc (solver_type, n_obs, n_params);
    gsl_multifit_fdfsolver_set (solver, &f, &x.vector);

    print_state (iter, solver);

    do
    {
        iter++;
        status = gsl_multifit_fdfsolver_iterate (solver);

        printf ("status = %s\n", gsl_strerror (status));

        print_state (iter, solver);

        if (status) break;
            status = gsl_multifit_test_delta (solver->dx, solver->x, 1e-5, 1e-5);

    } while (status == GSL_CONTINUE && iter < 500);

    double chi = gsl_blas_dnrm2(solver->f);
    double dof = n_obs - n_params;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof));

    printf ("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
    printf ("r0      = %.5f", gsl_vector_get(solver->x, 0));
    printf ("\nstatus = %s\n\n", gsl_strerror (status));

    double result = gsl_vector_get(solver->x, 0);

    gsl_multifit_fdfsolver_free (solver);

    return result;
}

void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
    int iteration_number = (int)iter;
    printf ("iteration: %3u r0 = % 15.8f "
            "|f(x)| = %g\n",
            iteration_number,
            gsl_vector_get (s->x, 0),
            gsl_blas_dnrm2 (s->f)
           );
}

float bin_distance(float distance)
{
    int numerator;

    if (distance < 0.3)
        return 0.0;

    for ( numerator = 1; distance >= (numerator+2)/6.0; numerator += 2 );
    return numerator/6.0;
}
