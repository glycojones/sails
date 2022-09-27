# Sails

Sails - Software for the Automated Identification of Linked Sugars

- C++11 and Python 3 (pybind11 interface) codebase.
- CMAKE build system
- Py.test framework for pre-commit and pre-release tests and profiling
- Mainly based on CCP4 libraries: Clipper, Gemmi, Privateer

Lead authors: Mihaela Atanasova, Kevin Cowtan and Jon Agirre (YSBL, University of York)

## **Installation instructions:**

**Operating systems supported** - **MacOS**(tested on "Monterey"(aka 12.1.X).

**Requirements:** 

**bzr** 

**cmake** (minimum version required 3.12)

**Homebrew**

**gfortran**


1.) git clone https://github.com/glycojones/sails.git sails

2.) cd sails

3.) git checkout master

4.) cd dependencies/privateer

4.) git clone https://github.com/glycojones/privateer.git 

5.) cd privateer

6.) git checkout master 

7.) git submodule update --init --recursive

8.) source ccp4.envsetup-sh

9.) pip3 install -r requirements.txt

10.) python3 setup.py install

11.) cd ../..

13.) python3 setup.py install


#### Funding 
* 2013-2015 BBSRC grant awarded to Keith Wilson (York)
* 2015-2017 BBSRC grant awarded to CCP4/Kevin Cowtan (York)
* 2017-     Royal Society University Research Fellowship awarded to Jon Agirre (York)
* 2018-2021 EPSRC DTP studentship allocated to Jon Agirre & awarded to Mihaela Atanasova (York)
* 2019-2023 Royal Society Research Grant awarded to Jon Agirre (York)
