README
  1 TITLE
    2D Frequency-domain Full Waveform Inversion Environment in Matlab
  2 DESCRIPTION
    Various SPOT operators
    opDFTR       - FFT for real signals
    opExtension  - pads input with zeros or constant values
    opLInterp1D  - 1D cubic Lagrange interpolation
    opLInterp2D  - 2D linear Lagrange interpolation
    opSmooth     - Smoothing by convolution with triangular kernel
    opSpline1D   - 1D cubic spline evaluation
    opSpline2D   - 2D cubic spline evaluation
    test         - unit tests
  3 COPYRIGHT
    You may use this code only under the conditions and terms of the
    license contained in the file LICENSE or COPYING.txt provided with
    this source code. If you do not agree to these terms you may not
    use this software.
  4 PREREQUISITES
    All prerequisites, except for MATLAB, are provided with the
    software release and should be installed before using any of
    SINBAD's software.
  5 INSTALLATION
    Follow the instructions in the INSTALLATION file (located in the
    root directory of this software release) to install all 3-rd party
    software (except for MATLAB) and SINBAD's software.
  6 DOCUMENTATION
    Documentation for each of the functions can be accessed by typing `help <function>' in Matlab. 
    Examples are provided in applications/FWI 
  7 RUNNING
    The functions can be called directly from Matlab. 
  8 NOTES
  9 SUPPORT
    You may contact developers of SINBAD software by means of:
    9.1 Mailing list
      Subscribe to SINBAD software mailing list at
      http://slim.eos.ubc.ca/mailman/listinfo/slimsoft and e-mail your
      question to the mailing list.
    9.2 Direct mail
      Contact SLIM administrator at nadmin@slim.eos.ubc.ca with any
      questions related to the SINBAD software release.