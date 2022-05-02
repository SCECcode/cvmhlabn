# CVMH15-1-Los-Angeles-Basin (cvmhlabn)

<a href="https://github.com/sceccode/cvmhlabn.git"><img src="https://github.com/sceccode/cvmhlabn/wiki/images/cvmhlabn_logo.png"></a>

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![GitHub repo size](https://img.shields.io/github/repo-size/sceccode/cvmhlabn)
[![cvmhlabn-ci Actions Status](https://github.com/SCECcode/cvmhlabn/workflows/cvmhlabn-ci/badge.svg)](https://github.com/SCECcode/cvmhlabn/actions)
[![cvmhlabn-ucvm-ci Actions Status](https://github.com/SCECcode/cvmhlabn/workflows/cvmhlabn-ucvm-ci/badge.svg)](https://github.com/SCECcode/cvmhlabn/actions)


## Description

CVMH Los Angeles Basin

## Table of Contents
1. [Software Documentation](https://github.com/SCECcode/cvmhlabn/wiki)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Contributing](#contributing)
5. [Credits](#credit)
6. [License](#license)

## Installation
This package is intended to be installed as part of the UCVM framework,
version 21.12.0 or higher. 

This package can also be build as a standalone program

<pre>
git submodule init
git submodule update

aclocal
automake --add-missing
autoconf
./configure --prefix=/path/to/install
cd data; ./make_data_files.py 
make
make install
</pre>

## Usage

### UCVM

As part of [UCVM](https://github.com/SCECcode/ucvm) installation, use 'cvmhlabn' as the model.

### vx_lite_cvmhlabn

A command line program accepts Geographic Coordinates or UTM Zone 11 to extract velocity values
from CVMHLABN.

## Support
Support for CVMHLABN is provided by the Southern California Earthquake Center
(SCEC) Research Computing Group.  Users can report issues and feature requests 
using CVMHLABN's github-based issue tracking link below. Developers will also 
respond to emails sent to the SCEC software contact listed below.
1. [CVMHLABN Github Issue Tracker](https://github.com/SCECcode/cvmhlabn/issues)
2. Email Contact: software@scec.usc.edu

## Contributing
We welcome contributions to the CVMHLABN, please contact us at software@scec.usc.edu.

## Credits
* Andreas Plesch <andreas_plesch@harvard.edu>
* John Shaw <shaw@eps.harvard.edu>

## License
This software is distributed under the BSD 3-Clause open-source license.
Please see the [LICENSE.txt](LICENSE.txt) file for more information.

