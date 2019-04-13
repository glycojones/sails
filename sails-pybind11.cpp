
// Pybind11 Python bindings for the YSBL program Sails
// (Software for the Automated Identification of Linked Sugars)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2019 Mihaela Atanasova, Kevin Cowtan and Jon Agirre
// York Structural Biology Laboratory
// The University of York
// mailto: jon.agirre@york.ac.uk
//

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "sails-lib.h"
using namespace pybind11::literals;


// pybind11 module definition
//
PYBIND11_MODULE(privateer, m)
{
  m.doc() = "Sails's Python interface.\nVersion history:\n- 2019-present MKI (pybind11)"; // docstring

}
