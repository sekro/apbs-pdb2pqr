Please see http://www.poissonboltzmann.org/apbs/release-history for the
complete change log.

# These are notes for APBS version 1.4.1-binary-release.

* This release coincides with the new web site, hosted using github
  pages, but is still found here:  http://www.poissonboltzmann.org/
* Git source repository has been moved from SourceForge to GitHub:
  (https://github.com/Electrostatics/apbs-pdb2pqr).  The binary releases are
  still to be found at SourceForge.

New Features:
* New Geometric Flow
    - Preliminary implementation with a limited subset of tweakable options.
    - Electrostatic energies returned by the current implementation tend to be
      a bit on the high side. We suspect that there is a bug in the code, and
      plan of resolving the issue in a future release.  The full array of
      algorithm options will be surfaced in a future release as well.
* Finite element method support has been re-enabled for Linux and OS X platforms
  via the Finite Element ToolKit (http://www.fetk.org).
* Added a progress indicator to mergedx2 before starting a really long process.
* Application Bundle distribution for Apple OS X, in addition to the tarball.

APBS defect fixes:
* A multitude of bugs related to crossing the 32-bit to 64-bit boundary have
  been dealt with.  The files involved touched nearly every program in the APBS
  suite, and thus these changes are far reaching.
* Many of the bugs fixed were in the multigrid code.  Previously more than 2^32
  grid points would crash APBS.
* The membrane potential boundary condition is once again functional, thanks to
  Frank Marcoline.
* 'Forked' the FeTK into a local branch so that it would not spit out silly
  errors when reading PQR files.  For more details, see: https://github.com/Electrostatics/apbs-pdb2pqr/blob/1.4.1-binary-release/apbs/fetk-1.5_APBS/README_APBS.md
* Fixed the code that reads PQR files so that it no longer requires a hack to
  MALOC.
* Fixed mergedx2's woefully broken handling of command line arguments.

Build System defects addressed:
* Minor library ordering changes to fix build failure on Ubuntu using a newish
  "--as-needed" linker flag.
* Removed hard-coded dependency on the GNU Fortran library.
* Included psize.py in Windows installer for PyMOL love.

Test Suite changes:
* Fixed a test that was failing because it got munged in the move to CMake.
  The test, examples/FKBP/1d7h-dmso-mol.in, now passes as expected.
* Fixed an error in apbs/examples/born/ion.pqr that was causing test failures.
* Fixed a bug in the test validator, apbs_check_results.py, that caused it to
  print negative percent errors.
* Added FEM tests back into the test suite.
* All README.html files, located in examples/<test> have been updated with results
  from this release.