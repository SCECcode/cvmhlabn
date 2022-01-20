# CVMH15-1-Los-Angeles-Basin (cvmhlabn)
Description of CVMHLABN

## Table of Content
1. [Software Documentation](https://github.com/SCECcode/cvmhlabn/wiki)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Contributing](#contributing)
5. [Credits](#credit)
6. [License](#license)

### Installation
This package is intended to be installed as part of the UCVM framework,
version 21.10.0 or higher. 

This package can also be tested in standalone mode by,

```
aclocal
automake --add-missing
autoconf
./configure --prefix=$UCVM_INSTALL_PATH/model/cvmhlabn 

cd data; ./make_data_files.py 
make
make install
make run_test
```

### Usage


### Support
Support for CVMHLABN is provided by the Southern California Earthquake Center
(SCEC) Research Computing Group.  Users can report issues and feature requests 
using CVMHLABN's github-based issue tracking link below. Developers will also 
respond to emails sent to the SCEC software contact listed below.
1. [CVMHLABN Github Issue Tracker](https://github.com/SCECcode/cvmhlabn/issues)
2. Email Contact: software@scec.usc.edu

### Contributing
We welcome contributions to the CVMHLABN, please contact us at software@scec.usc.edu.

### Credits
* Andreas Plesch <andreas_plesch@harvard.edu>
* John Shaw <shaw@eps.harvard.edu>

