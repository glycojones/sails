// Clipper app to compare two sets of min, max & avg maps
// 2013 Jon Agirre & Kevin Cowtan @ University of York
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk
//

#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <clipper/clipper.h>
#include <clipper/clipper-cif.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/contrib/sfcalc_obs.h>


inline bool fileMissing (const std::string&);

using namespace std;


int main(int argc, char* argv[])
{

    if (argc != 3)
    {
        std::cout << "\nUsage: map_combine <reference-prefix> <input-prefix> \n\n";
        return 0;
    }

    clipper::CCP4MAPfile referenceMap_min, referenceMap_max, referenceMap_avg, inputMap_min, inputMap_max, inputMap_avg;

    if (fileMissing( clipper::String(argv[1] + clipper::String("_min.map") ))) // particular case for the first run for each sugar+conformation
    {                                                                         // we will compare the maps with themselves and thus generate a copy
        referenceMap_min.open_read( clipper::String(argv[2] + clipper::String("_min.map") ));
        referenceMap_max.open_read( clipper::String(argv[2] + clipper::String("_max.map") ));
        referenceMap_avg.open_read( clipper::String(argv[2] + clipper::String("_avg.map") ));
    }
    else
    {
        referenceMap_min.open_read( clipper::String(argv[1] + clipper::String("_min.map") ));
        referenceMap_max.open_read( clipper::String(argv[1] + clipper::String("_max.map") ));
        referenceMap_avg.open_read( clipper::String(argv[1] + clipper::String("_avg.map") ));
    }

    // read in the reference maps
    clipper::NXmap<float> referenceNXMap_min;
    referenceMap_min.import_nxmap(referenceNXMap_min);
    referenceMap_min.close_read();

    clipper::NXmap<float> referenceNXMap_max;
    referenceMap_max.import_nxmap( referenceNXMap_max);
    referenceMap_max.close_read();

    clipper::NXmap<float> referenceNXMap_avg;
    referenceMap_avg.import_nxmap(referenceNXMap_avg);
    referenceMap_avg.close_read();

    // we read & close first the reference map since it could be the same as the input one

    inputMap_min.open_read( clipper::String(argv[2] + clipper::String("_min.map") ));
    inputMap_max.open_read( clipper::String(argv[2] + clipper::String("_max.map") ));
    inputMap_avg.open_read( clipper::String(argv[2] + clipper::String("_avg.map") ));

    // read in the input maps
    clipper::NXmap<float> inputNXMap_min;
    inputMap_min.import_nxmap(inputNXMap_min);
    inputMap_min.close_read();

    clipper::NXmap<float> inputNXMap_max;
    inputMap_max.import_nxmap(inputNXMap_max);
    inputMap_max.close_read();

    clipper::NXmap<float> inputNXMap_avg;
    inputMap_avg.import_nxmap(inputNXMap_avg);
    inputMap_avg.close_read();


    typedef clipper::NXmap<float>::Map_reference_index MRI;

    for ( MRI ix = inputNXMap_min.first(); !ix.last(); ix.next() )  // run across the map
    {
        referenceNXMap_min[ix] = clipper::Util::min( referenceNXMap_min[ix], inputNXMap_min[ix] );
        referenceNXMap_max[ix] = clipper::Util::min( referenceNXMap_max[ix], -inputNXMap_max[ix] );
        referenceNXMap_avg[ix] = (referenceNXMap_avg[ix] + inputNXMap_avg[ix]) / 2;
    }

    // get list of peaks
    const double radius = 6.0;
    const double sampling = 0.25;
    double extent = radius + 2.0 * sampling;  // model box size
    double ng = rint( extent / sampling ); // number of grid points, round to integer
    clipper::Grid nxgrid( 2*int(ng)+1, 2*int(ng)+1, 2*int(ng)+1 );
    clipper::Coord_grid cg;
    clipper::Coord_grid cgu( 1, 0, 0 ), cgv( 0, 1, 0 ), cgw( 0, 0, 1 );


    clipper::NXmap<float> rhomin, rhomax;

    float rho0;
    std::vector<std::pair<float,int> > peakmax, peakmin; // beware of these ones

    for ( cg.w() = 1; cg.w() < 2*ng; cg.w()++ )
        for ( cg.v() = 1; cg.v() < 2*ng; cg.v()++ )
            for ( cg.u() = 1; cg.u() < 2*ng; cg.u()++ )
            {
                rho0 = referenceNXMap_min.get_data( cg );

                if ( rho0 > referenceNXMap_min.get_data( cg - cgu ) &&
                     rho0 > referenceNXMap_min.get_data( cg + cgu ) &&
                     rho0 > referenceNXMap_min.get_data( cg - cgv ) &&
                     rho0 > referenceNXMap_min.get_data( cg + cgv ) &&
                     rho0 > referenceNXMap_min.get_data( cg - cgw ) &&
                     rho0 > referenceNXMap_min.get_data( cg + cgw ) ) peakmin.push_back( std::pair<float,int>( rho0, nxgrid.index(cg) ) );

                     rho0 = referenceNXMap_min.get_data( cg );

                if ( rho0 > referenceNXMap_max.get_data( cg - cgu ) &&
                     rho0 > referenceNXMap_max.get_data( cg + cgu ) &&
                     rho0 > referenceNXMap_max.get_data( cg - cgv ) &&
                     rho0 > referenceNXMap_max.get_data( cg + cgv ) &&
                     rho0 > referenceNXMap_max.get_data( cg - cgw ) &&
                     rho0 > referenceNXMap_max.get_data( cg + cgw ) ) peakmax.push_back( std::pair<float,int>( rho0, nxgrid.index(cg) ) );
            }

  std::cout << peakmin.size() << " " << peakmax.size() << std::endl;

  std::sort( peakmin.begin(), peakmin.end() );
  std::sort( peakmax.begin(), peakmax.end() );
  std::reverse( peakmin.begin(), peakmin.end() );
  std::reverse( peakmax.begin(), peakmax.end() );

  for ( int index = 0; index < 10; index++ )  // output the first 10 min peaks
    std::cout << peakmin[index].first << " " << rhomin.coord_orth(nxgrid.deindex(peakmin[index].second).coord_map()).format() << std::endl;

  for ( int index = 0; index < 10; index++ )  // output the first 10 max peaks
    std::cout << peakmax[index].first << " " << rhomax.coord_orth(nxgrid.deindex(peakmax[index].second).coord_map()).format() << std::endl;

    referenceMap_min.open_write( clipper::String(argv[1] + clipper::String("_min.map") ));
    referenceMap_min.export_nxmap(referenceNXMap_min);
    referenceMap_min.close_write();

    referenceMap_max.open_write( clipper::String(argv[1] + clipper::String("_max.map") ));
    referenceMap_max.export_nxmap(referenceNXMap_max);
    referenceMap_max.close_write();

    referenceMap_avg.open_write( clipper::String(argv[1] + clipper::String("_avg.map") ));
    referenceMap_avg.export_nxmap(referenceNXMap_avg);
    referenceMap_avg.close_write();

    return 0;

}

inline bool fileMissing (const std::string& filename)
{
    ifstream f (filename.c_str());

    if (f.good())
    {
        f.close();
        return false;
    }
    else
    {
        f.close();
        return true;
    }
}

// // Clipper app to compare two sets of min, max & avg maps
// // 2013 Jon Agirre & Kevin Cowtan @ University of York
// // mailto: jon.agirre@york.ac.uk
// // mailto: kevin.cowtan@york.ac.uk
// //

// #include <fstream>
// #include <algorithm>
// #include <stdio.h>
// #include <string.h>
// #include <stdlib.h>
// #include <clipper/clipper.h>
// #include <clipper/clipper-cif.h>
// #include <clipper/clipper-mmdb.h>
// #include <clipper/clipper-ccp4.h>
// #include <clipper/clipper-contrib.h>
// #include <clipper/clipper-minimol.h>
// #include <clipper/contrib/sfcalc_obs.h>


// inline bool fileMissing (const std::string&);

// using namespace std;


// int main(int argc, char* argv[])
// {

//     if (argc != 3)
//     {
//         std::cout << "\nUsage: map_combine <reference-prefix> <input-prefix> \n\n";
//         return 0;
//     }

//     clipper::CCP4MAPfile referenceMap_min, referenceMap_max, referenceMap_avg, inputMap_min, inputMap_max, inputMap_avg;

//     if (fileMissing( clipper::String(argv[1] + clipper::String("_min.map") ))) // particular case for the first run for each sugar+conformation
//     {                                                                         // we will compare the maps with themselves and thus generate a copy
//         referenceMap_min.open_read( clipper::String(argv[2] + clipper::String("_min.map") ));
//         referenceMap_max.open_read( clipper::String(argv[2] + clipper::String("_max.map") ));
//         referenceMap_avg.open_read( clipper::String(argv[2] + clipper::String("_avg.map") ));
//     }
//     else
//     {
//         referenceMap_min.open_read( clipper::String(argv[1] + clipper::String("_min.map") ));
//         referenceMap_max.open_read( clipper::String(argv[1] + clipper::String("_max.map") ));
//         referenceMap_avg.open_read( clipper::String(argv[1] + clipper::String("_avg.map") ));
//     }

//     // read in the reference maps
//     clipper::NXmap<float> referenceNXMap_min;
//     referenceMap_min.import_nxmap(referenceNXMap_min);
//     referenceMap_min.close_read();

//     clipper::NXmap<float> referenceNXMap_max;
//     referenceMap_max.import_nxmap( referenceNXMap_max);
//     referenceMap_max.close_read();

//     clipper::NXmap<float> referenceNXMap_avg;
//     referenceMap_avg.import_nxmap(referenceNXMap_avg);
//     referenceMap_avg.close_read();

//     // we read & close first the reference map since it could be the same as the input one

//     inputMap_min.open_read( clipper::String(argv[2] + clipper::String("_min.map") ));
//     inputMap_max.open_read( clipper::String(argv[2] + clipper::String("_max.map") ));
//     inputMap_avg.open_read( clipper::String(argv[2] + clipper::String("_avg.map") ));

//     // read in the input maps
//     clipper::NXmap<float> inputNXMap_min;
//     inputMap_min.import_nxmap(inputNXMap_min);
//     inputMap_min.close_read();

//     clipper::NXmap<float> inputNXMap_max;
//     inputMap_max.import_nxmap(inputNXMap_max);
//     inputMap_max.close_read();

//     clipper::NXmap<float> inputNXMap_avg;
//     inputMap_avg.import_nxmap(inputNXMap_avg);
//     inputMap_avg.close_read();


//     typedef clipper::NXmap<float>::Map_reference_index MRI;

//     for ( MRI ix = inputNXMap_min.first(); !ix.last(); ix.next() )  // run across the map
//     {
//         referenceNXMap_min[ix] = clipper::Util::min( referenceNXMap_min[ix], inputNXMap_min[ix] );
//         referenceNXMap_max[ix] = clipper::Util::min( referenceNXMap_max[ix], -inputNXMap_max[ix] );
//         referenceNXMap_avg[ix] = (referenceNXMap_avg[ix] + inputNXMap_avg[ix]) / 2;
//     }

//     // get list of peaks
//     const double radius = 6.0;
//     const double sampling = 0.25;
//     double extent = radius + 2.0 * sampling;  // model box size
//     double ng = rint( extent / sampling ); // number of grid points, round to integer
//     clipper::Grid nxgrid( 2*int(ng)+1, 2*int(ng)+1, 2*int(ng)+1 );
//     clipper::Coord_grid cg;
//     clipper::Coord_grid cgu( 1, 0, 0 ), cgv( 0, 1, 0 ), cgw( 0, 0, 1 );


//     clipper::NXmap<float> rhomin, rhomax;

//     float rho0;
//     std::vector<std::pair<float,int> > peakmax, peakmin; // beware of these ones

//     for ( cg.w() = 1; cg.w() < 2*ng; cg.w()++ )
//         for ( cg.v() = 1; cg.v() < 2*ng; cg.v()++ )
//             for ( cg.u() = 1; cg.u() < 2*ng; cg.u()++ )
//             {
//                 rho0 = referenceNXMap_min.get_data( cg );

//                 if ( rho0 > referenceNXMap_min.get_data( cg - cgu ) &&
//                      rho0 > referenceNXMap_min.get_data( cg + cgu ) &&
//                      rho0 > referenceNXMap_min.get_data( cg - cgv ) &&
//                      rho0 > referenceNXMap_min.get_data( cg + cgv ) &&
//                      rho0 > referenceNXMap_min.get_data( cg - cgw ) &&
//                      rho0 > referenceNXMap_min.get_data( cg + cgw ) ) peakmin.push_back( std::pair<float,int>( rho0, nxgrid.index(cg) ) );

//                      rho0 = referenceNXMap_min.get_data( cg );

//                 if ( rho0 > referenceNXMap_max.get_data( cg - cgu ) &&
//                      rho0 > referenceNXMap_max.get_data( cg + cgu ) &&
//                      rho0 > referenceNXMap_max.get_data( cg - cgv ) &&
//                      rho0 > referenceNXMap_max.get_data( cg + cgv ) &&
//                      rho0 > referenceNXMap_max.get_data( cg - cgw ) &&
//                      rho0 > referenceNXMap_max.get_data( cg + cgw ) ) peakmax.push_back( std::pair<float,int>( rho0, nxgrid.index(cg) ) );
//             }

//   std::cout << peakmin.size() << " " << peakmax.size() << std::endl;

//   std::sort( peakmin.begin(), peakmin.end() );
//   std::sort( peakmax.begin(), peakmax.end() );
//   std::reverse( peakmin.begin(), peakmin.end() );
//   std::reverse( peakmax.begin(), peakmax.end() );

//   for ( int index = 0; index < 10; index++ )  // output the first 10 min peaks
//     std::cout << peakmin[index].first << " " << rhomin.coord_orth(nxgrid.deindex(peakmin[index].second).coord_map()).format() << std::endl;

//   for ( int index = 0; index < 10; index++ )  // output the first 10 max peaks
//     std::cout << peakmax[index].first << " " << rhomax.coord_orth(nxgrid.deindex(peakmax[index].second).coord_map()).format() << std::endl;

//   for ( MRI ix = inputNXMap_min.first(); !ix.last(); ix.next() )  // run across the map
//    for (int i = 0; i<20; i++)  //what is i supposed to be ?
//    {
//         referenceNXMap_min[ix] = peakmin[i].first;
//    }

//   for ( MRI ix = inputNXMap_max.first(); !ix.last(); ix.next() )  // run across the map
//     for (int i = 0; i<20; i++)
//     {
//         referenceNXMap_max[ix] = peakmax[i].first;
//     }
   
//     referenceMap_min.open_write( clipper::String(argv[1] + clipper::String("_min.map") ));
//     referenceMap_min.export_nxmap(referenceNXMap_min);
//     referenceMap_min.close_write();

//     referenceMap_max.open_write( clipper::String(argv[1] + clipper::String("_max.map") ));
//     referenceMap_max.export_nxmap(referenceNXMap_max);
//     referenceMap_max.close_write();

//     referenceMap_avg.open_write( clipper::String(argv[1] + clipper::String("_avg.map") ));
//     referenceMap_avg.export_nxmap(referenceNXMap_avg);
//     referenceMap_avg.close_write();

//   // for ( int index = 0; index < 10; index++ )  // output the first 10 max peaks
//     // std::cout << peakmax[index].first << " " << rhomax.coord_orth(nxgrid.deindex(peakmax[index].second).coord_map()).format() << std::endl;
    
//     // referenceMap_min.open_write( clipper::String(argv[1] + clipper::String("_min.map") ));
//     // referenceMap_min.export_nxmap(referenceNXMap_min);
//     // referenceMap_min.close_write();

//     // referenceMap_max.open_write( clipper::String(argv[1] + clipper::String("_max.map") ));
//     // referenceMap_max.export_nxmap(referenceNXMap_max);
//     // referenceMap_max.close_write();

//     // referenceMap_avg.open_write( clipper::String(argv[1] + clipper::String("_avg.map") ));
//     // referenceMap_avg.export_nxmap(referenceNXMap_avg);
//     // referenceMap_avg.close_write();

//     return 0;

// }

// inline bool fileMissing (const std::string& filename)
// {
//     ifstream f (filename.c_str());

//     if (f.good())
//     {
//         f.close();
//         return false;
//     }
//     else
//     {
//         f.close();
//         return true;
//     }
// }
