---
title: 'EDL Shock & Detonation Toolbox - Cantera 2.1 Version'
...

Shock & Detonation Toolbox - Cantera 2.1
========================================

------------------------------------------------------------------------


------------------------------------------------------------------------

### About the SD Toolbox

The Shock & Detonation Toolbox is a collection of numerical routines
that enables the solution of standard problems for gas-phase explosions
using realistic thermochemistry and detailed chemical kinetics. The SD
Toolbox employs Dave Goodwin's [Cantera](http://www.cantera.org/)
software for the chemistry functionality and uses either MATLAB or
Python (and related libraries) for scripting. The Cantera package
provides conversion utilities from legacy formats in order to make use
of existing databases of [chemical
kinetics](http://shepherd.caltech.edu/EDL/public/mechanisms.html) and
[thermochemistry](http://shepherd.caltech.edu/EDL/public/thermo.html).
The focus of the SD Toolbox is not computational fluid dynamics, which
is outside the scope of our project at the present time.\
\
 The SD Toolbox includes numerical routines for the computation of:

-   CJ detonation speed and post-detonation state
-   Postshock gas state for frozen composition
-   Postshock gas state for equilibrium composition
-   Frozen and equilibrium Hugoniot curves
-   Constant-volume explosion structure
-   ZND detonation structure (need additional routines)
-   Effective activation energies and chemical time scales from detailed
    reaction mechanisms
-   Extrapolation of low temperature thermodynamic polynomial fits to
    higher temperatures.

The SD Toolbox has been implemented in both MATLAB and Python. There are
different releases of the toolbox as described below for each
environment. Because of differences in Cantera within MATLAB and Python,
not all of the releases offer the same functions. The SD Toolbox
programs are available free of charge, but are subject to the same
[license
agreement](http://shepherd.caltech.edu/EDL/public/cantera/html/SD_Toolbox/License.txt)
as Cantera.\
\
 For the theoretical background, details on the algorithms, and
implementations for the SD Toolbox, please read

-   [Numerical Solution Methods for Shock and Detonation Jump
    Conditions](http://shepherd.caltech.edu/EDL/public/cantera/doc/tex/ShockDetonation/ShockDetonation.pdf)
    This report contains both the theoretical background and information
    on the algorithms used by the SD\_Toolbox. This is a hypertext
    document with links to the scripts that are used to generate the
    results for example problems.\
    \
-   [Numerical Solution Methods for Control Volume Explosions and ZND
    Detonation
    Structure](http://shepherd.caltech.edu/EDL/public/cantera/doc/tex/CVZND/CVZND.pdf)
    (Draft Version.)

\

------------------------------------------------------------------------

### Installation of Cantera

The SD Toolbox requires Cantera along with either MATLAB or Python. The
recommended version of Cantera is 2.1 or higher. Instructions for
installing and using Cantera under Windows in either the MATLAB or
Python environments can be found
[here](http://cantera.org/docs/sphinx/html/install.html#windows).
Instructions for building and installing Cantera from source code under
Windows, Mac OS X, and Linux are available
[here](http://cantera.org/docs/sphinx/html/compiling/index.html).\
\

------------------------------------------------------------------------

### Installation of SD Toolbox under Python

 The Python interface to Cantera was changed considerably in the
transition from Cantera 2.0 to Cantera 2.1. Many of the changes are
documented
[here](http://cantera.org/docs/sphinx/html/cython/migrating.html);
these changes include different naming conventions for constants and
thermodynamic properties as well as new syntax for reporting and setting
the thermodynamic state. The Python interface to the SD Toolbox has been
updated accordingly. However, because of the extensive changes to the
Cantera-Python interface, the new SD Toolbox is NOT backward compatible
with Cantera 2.0 and earlier.\
\
 The Python version of the SD Toolbox and a collection of demo scripts
can downloaded from the following links:

-   Shock and Detonation Toolbox
    ([zip](http://shepherd.caltech.edu/EDL/public/cantera/2.1/python/SDToolbox/SDToolbox_Python2.7_Cantera2.1.zip))
    ([tgz](http://shepherd.caltech.edu/EDL/public/cantera/2.1/python/SDToolbox/SDToolbox_Python2.7_Cantera2.1.tgz))
-   Python demo scripts
    ([zip](http://shepherd.caltech.edu/EDL/public/cantera/2.1/python/SDToolbox/demos_Python2.7_Cantera2.1.zip))
    ([tgz](http://shepherd.caltech.edu/EDL/public/cantera/2.1/python/SDToolbox/demos_Python2.7_Cantera2.1.tgz))

After successfully installing both Cantera and Python, the SD Toolbox is
installed using the following procedure:

> ```bash 
> pip install SDToolbox
> ```
1.  ~~Download and install [Scipy](http://www.scipy.org/install.html).~~
2.  ~~Create a directory named SDToolbox on the path
    *PYTHONPATH/Lib/site-packages/SDToolbox*, where *PYTHONPATH* is the
    location of your Python installation. For Windows, this is usually
    *C:/Python27* and for Linux it might be */usr/local/lib/python2.7*.~~
3.  ~~Download the SDToolbox from the link above and unzip the files into
    the SDToolbox directory just created.~~
4.  ~~For Windows only, go up one directory to
    *PYTHONPATH/Lib/site-packages* and create a text file named
    "SDToolbox.pth" containing the single line "SDToolbox" (no quotes).~~
5.  ~~For Windows only, go up one more directory to *PYTHONPATH/Lib* and
    run the python script 'site.py'~~

\
 The SD Toolbox is now installed. To test the installation, open a
terminal window and run python. In the python environment, execute
`from SDToolbox import *` to load the toolbox and then use the call\
\

```python 
[cj_speed,_] = CJspeed(101325,300,'H2:2 O2:1','gri30.cti',0)
```

\
\
 This will compute the CJ speed for a stoichiometric hydrogen-oxygen
detonation at ambient temperature and pressure. If this call executes
with no errors, the installation was successful.\
\

------------------------------------------------------------------------

### Installation of SD Toolbox under MATLAB

The MATLAB version of the SD Toolbox and a collection of demo scripts
can be downloaded here:

-   Shock and Detonation Toolbox
    ([zip](http://shepherd.caltech.edu/EDL/public/cantera/2.1/matlab/SDToolbox/SDToolbox_MATLAB_Cantera2.1.zip))
    ([tgz](http://shepherd.caltech.edu/EDL/public/cantera/2.1/matlab/SDToolbox/SDToolbox_MATLAB_Cantera2.1.tgz))
-   MATLAB demo scripts
    ([zip](http://shepherd.caltech.edu/EDL/public/cantera/2.1/matlab/SDToolbox/demos_MATLAB_Cantera2.1.zip))
    ([tgz](http://shepherd.caltech.edu/EDL/public/cantera/2.1/matlab/SDToolbox/demos_MATLAB_Cantera2.1.tgz))

After successfully installing Cantera and MATLAB, the SD Toolbox can be
installed using the following procedure:

1.  Navigate to your MATLAB installation and create the folder
    *R20xxy/toolbox/SDToolbox*, where *R20xxy* is your MATLAB release,
    e.g., *R2014a*.
2.  Unpack the SD Toolbox *zip* or *tar* file downloaded above into the
    folder just created.
3.  Start MATLAB. Go to *File--&gt;Set Path* and click *Add with
    Subfolders*. Select the SDToolbox folder created in step (1).
4.  Click *Save* and *Close*

The SD Toolbox is now installed. You can test it by calling one of the
functions, for example,\
\

```python 
U = CJspeed(101325, 300, 'H2:2 O2:1', 'gri30.cti', 1)
```

\
 which will compute the Chapman-Jouguet speed for a stoichiometric
hydrogen-oxygen detonation.\
\

------------------------------------------------------------------------

### ZND Calculations

 ZND calculations have been implemented in the MATLAB version of the
toolbox above, but are not available in the python version. However, on
Linux platforms the ZND calculations are implemented both in C++ and in
python scripts that call on a C++ executable. The source code, makefile,
and python scripts are available in
[zip](http://shepherd.caltech.edu/EDL/public/cantera/2.1/znd/znd_Python2.7_Cantera2.1.zip)
and
[tgz](http://shepherd.caltech.edu/EDL/public/cantera/2.1/znd/znd_Python2.7_Cantera2.1.tgz)
archives. These scripts have been built and tested on Ubuntu 12.04 with
GCC-4.6.3 and Python 2.7. The ZND code can be built using the following
steps.

1.  Obtain working installations of Cantera, Python (version 2.7
    recommended), GCC, and pkg-config. The last three of these can be
    installed using the package manager on most linux systems.
2.  Download and unpack the source code using the command
    ` tar -zxvf znd_Cantera2.1_Python2.7.tgz` or
    ` unzip znd_Cantera2.1_Python2.7.zip`.
3.  The unpacked directory tree contains three folders, *bin*, *build*,
    and *src*. Navigate to the *build* directory and execute
    ` make -f znd.make ` to build the executable. The program pkg-config
    should automatically link to the correct Cantera libraries.

The executable is placed in the *bin* directory. To run the ZND
calculation, execute `./znd`. You will then be asked to supply the name
of an input file which specifies the flow conditions. Several examples
of input files are available in the *bin* directory. The *bin* directory
also contains several python scripts which demonstrate how to call the
znd program from python. The scripts *znd\_fseries.py*,
*znd\_phiseries.py*, and *znd\_Pseries.py* perform ZND calculations for
sequences of overdrive values, equivalence ratios, and pressures.\
\

------------------------------------------------------------------------

### Mechanism Files for Cantera

Cantera requires a mechanism (.cti)
file with thermodynamic and reaction rate data for the species of
interest. Some files are included when Cantera is installed, but we
provide additional files below that have been modified by our group for
application to shock and detonation problems. The species thermodynamic
data in these files have coefficients that are valid to much higher
temperatures (5000-6000 K in most cases) than the original data sets
supplied with Cantera. However, it is to be noted that these high
temperature data sets were computed by extrapolation from low
temperature data and have not been extensively validated. Individual
mechanism files are available from the following links:

-   [GRI Mech
    3.0](http://shepherd.caltech.edu/EDL/public/cantera/mechs/cti/web/gri30_highT.cti)
    (for CH4-C2H4-air mixtures)
-   [Hydrogen Oxygen
    Mechanism](http://shepherd.caltech.edu/EDL/public/cantera/mechs/cti/web/h2o2_highT.cti)
-   [Hydrogen Air
    Mechanism](http://shepherd.caltech.edu/EDL/public/cantera/mechs/cti/web/h2air_highT.cti)
-   [Hydrogen Bromine
    Mechanism](http://shepherd.caltech.edu/EDL/public/cantera/mechs/cti/web/h2br2_highT.cti)
-   [Octane
    Mechanism](http://shepherd.caltech.edu/EDL/public/cantera/mechs/cti/web/octane_highT.cti)
    (from LLNL)
-   [NASA
    Data](http://shepherd.caltech.edu/EDL/public/cantera/mechs/cti/web/nasa.cti)
-   [Hai Wang
    Mechanism](http://shepherd.caltech.edu/EDL/public/cantera/mechs/cti/web/Wang_highT.cti)
-   [Hydrogen-Nitrous Oxide
    Mechanism](http://shepherd.caltech.edu/EDL/public/cantera/mechs/cti/web/h2-n2o_highT.cti)

The entire collection of mechanism files is also available as a
([zip](http://shepherd.caltech.edu/EDL/public/cantera/mechs/cti/high_T_cti_files.zip))
([tar](http://shepherd.caltech.edu/EDL/public/cantera/mechs/cti/high_T_cti_files.tar))
archive.\
\
 Please note that these files are only included for use with the demo
programs - you will need to either add these files to the
\\ProgramFiles\\Cantera\\data\\ directory (PC installation) or place
them in the working directory for the demo programs. These mechanisms
are not intended to be representative of the state of the practice in
chemical kinetics or thermodynamics. Many of these are derived from
older compilations and we make no guarrantees about the accuracy,
particularly in regard to the rate constants. The thermodynamic data is
reasonably reliable when used within the intended temperature range
although care should be taken at the upper end. The NASA 9-term
polynomial fits should be used if very high (&gt;6000 K) temperature
thermodynamic properties are required.\
\
 There are many newer sets of chemical reaction data and themodynamics
which are available for download on the www or as part of supplemental
data from archival publications. Users should seek out these newer data
sets and use them for any quantitative work that depend on the details
of the reaction mechanism. Many sets of chemistry data are stored in
CHEMKIN format instead of the cti style used by Cantera. These legacy
CHEMKIN data sets can be converted to cti format using the program
` ck2cti` that is included with Cantera. Instructions for making this
conversion are available from the [Cantera web
page](http://www.cantera.org/docs/sphinx/html/matlab/input-tutorial.html).
Note that the python converter `ck2cti.py ` is the now the preferred
tool for this rather than the previously used compiled program.\
\

------------------------------------------------------------------------

### Legacy SD Toolbox

Legacy versions of the SD Toolbox which are compatible with Cantera
1.7-2.0 remain available
[here](http://shepherd.caltech.edu/EDL/public/cantera/html/SD_Toolbox/Cantera_legacy.html),
but these versions will no longer be supported. Note that the MATLAB
interface to Cantera and the SD Toolbox was unaffected by the transition
to Cantera 2.1. The legacy webpages are still useful as they have an
extensive set of documentation with screen shots of the input and output
to demo programs. There are also discussions of how to fit and
extrapolate thermodynamic data; some users may find this useful.

### Notes and Acknowledgments

The Toolbox was updated in 2014 (Sept 9, 2014 ) for compatability with
Python 2.7 and Cantera 2.1. Significant changes to the Cantera-Python
interface demanded that the Python interface to the SDToolbox be updated
accordingly. The changes are not backward compatible. Byran Schmidt and
Neal Bitter carried out the conversion and testing of the new routines.\
\
 There have been substantial developments in Cantera since the last
major update of the SD\_Toolbox. The orginal author of Cantera, Dave
Goodwin, lost his battle with cancer and Parkinson's in 2012. Although
Dave was not able to actively contribute to development of Cantera for a
number of years due to his illness, Dave's vision of an open source tool
for energy research continues to be realized through the efforts of the
many volunteers that have contributed to the Cantera project. We
appreciate all the work by others that has gone into maintaining and
extending Cantera so that we can continue to rely on this as the
software engine underneath the SD\_Toolbox.

------------------------------------------------------------------------

\
\
