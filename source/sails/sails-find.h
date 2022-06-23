/*! \file ss-find.h buccaneer library */
/* (C) 2011-2019 Mihaela Atanasova, Kevin Cowtan & Jon Agirre */


#ifndef SS_FIND_H
#define SS_FIND_H

#include <clipper/clipper-contrib.h>
#include <privateer/cpp/privateer-lib.h>


/*! Result class */
class SearchResult {
  public:
    float score; int rot; int trn; std::string name_short; clipper::Mat33<> sugar_rot;
    bool operator <( SearchResult other ) const {
      return score < other.score;
    }
};

class SSfind {
  public:
    typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Pair_coord;

    void prep_xmap( const clipper::Xmap<float>& xmap,
                    const double radius );
    void prep_search( const clipper::Xmap<float>& xmap );
    void prep_search( const clipper::Xmap<float>& xmap,
                      const double rhocut,
                      const double radcut,
                      const clipper::Coord_orth centre );
    std::vector<SearchResult> search( const std::vector<Pair_coord>& target_cs,
                                      const std::vector<clipper::RTop_orth>& ops,
                                      const double rhocut,
                                      const double frccut = 0.0,
                                      std::string name_short = "") const; 

  private:
    std::vector<float> mapbox;
    std::vector<int> srctrn;
    clipper::Grid grid;
    clipper::Grid_range mxgr;
    clipper::Mat33<> grrot;
};

#endif
