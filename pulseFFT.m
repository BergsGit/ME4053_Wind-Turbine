function [ fourierTermMag, fourierTermPhase  ] = ...
   pulseFFT( freq, pulseWidth, pulseHeight, nHarmonics )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: pulseFFT
%
%  PURPOSE
%    Determine the frequency components of a harmonic pulse
%    having constant magnitude
%
%  INPUT
%    freq: Frequency of the input pulses (Hz)
%    nHarmonics: # of harmonics to include in the Fourier series
%                NOTE 1: The fundamental frequency is counted as
%                        the "first" harmonic
%                NOTE 2: nHarmonics must be < 500 due to the
%                        sampling rate assumed in this function
%    pulseHeight: Magnitude of the pulse 
%                 NOTE: Units must be consistent with forceMag
%    pulseWidth: Width of the pulse (sec)
%                NOTE: pulseWidth must be less than 1/freq!
%
%  OUTPUT
%    fourierTermMag: Magnitudes of the harmonic terms 
%                    (same units as pulseHeight)
%            NOTE 1: Units must be consistent with pulseHeight
%            NOTE 2: forceMag(1) = DC component
%                    forceMag(2) = magnitude of fundamental freq
%                    forceMag(3) = magnitude of 1st harmonic
%                                  (2X fundamental freq)
%                    forceMag(i) = magnitude of (i-2)nd harmonic
%                                  ((i-2)X fundamental freq)
%   fourierTermPhase: Phase angles of the harmonic terms (radians)
%                     (see Note 2 above for meaning of each term)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: F. Kelso & T. Chase
%  DATE: 10/25/24
%
%  DESCRIPTION OF LOCAL VARIABLES
%    ampPhase: Complex number representation of the raw
%              amplitude and phase of one harmonic component
%              of the input pulse
%    denom: Scaling factor to convert the raw digital representation
%           of the discrete Fourier transform of the input signal
%           to a usable magnitude and phase angle
%    index: Loop counter index
%    nPts: Number of points used to digitally sample the pulse
%    nPulse: Number of points of the digital signal 
%            that the pulse is "on" (between 0 and nPts)
%    pulseFrac: Fraction of the signal period that the pulse is "on"
%    samplePeriod: Period of the digital sample
%    signal: Raw discrete Fourier transform of the input pulse
%    sigPeriod: Period of the input signal (s)
%    sqrPulse: Array representing the pulse as a function of time
%              over the period of the input signal
%
%  FUNCTIONS CALLED
%    fft: Standard MATLAB function for performing a
%         discrete Fourier transform of an input vector
%    floor: Standard MATLAB function for rounding down
%           to the nearest integer toward infinity of a real number
%    round: standard MATLAB function for rounding a real number
%           to the closest integer
%    zeros: standard MATLAB function for creating an array of zeros
%
%  START OF EXECUTABLE CODE
%
%  Check for a legal number of harmonics
%
   if nHarmonics >= 500
%
%  Too many harmonics requested.
%  All return arguments must be set to something to avoid a MATLAB error.
%
      fourierTermMag = NaN;
      fourierTermPhase = NaN;
      return
   end
%
%  Determine period of cyclical pulse
%
   sigPeriod = 1.0 / freq;
%
%  Check for a legal pulse width
%
   if pulseWidth >= sigPeriod
%
%  Input is not a pulse!
%  All return arguments must be set to something to avoid a MATLAB error.
%
      fourierTermMag = NaN;
      fourierTermPhase = NaN;
      return
   end
%
%  Assume that signal is sampled using 1000 data points
%
   nPts = 1000;  
%
%  Determine period of sample
%
   samplePeriod = sigPeriod / nPts;
%
%  Determine fraction of the period where pulse is "on"
%
   pulseFrac = pulseWidth / sigPeriod; 
%
%  Define an array defining the pulse as a function of time...
%
%  ...Start with "off"
%
   sqrPulse = zeros( 1, nPts );
%
%  ...Now add the "on"
%     (+1 accounts for the pulse being "on" at zero time)
%
   nPulse = round( pulseFrac * nPts ) + 1; 
   sqrPulse( 1:nPulse ) = pulseHeight;
%
%  Perform a discrete Fourier transform on the input pulse

   signal = fft( sqrPulse ); %(complex) components of the input pulse
%
%  Convert the raw complex components of the input pulse
%  to magnitudes and angles...
%  (A reference needs to be added here to explain the theory behind
%   what is being done)
%
%  ...Reserve space for the output arrays
%
   fourierTermMag = zeros( 1, nHarmonics );
   fourierTermPhase = zeros( 1, nHarmonics );
%
%  ...Determine the scale factor
%
   denom = floor( nPts );
%
%  ...Set the magnitude and phase angle (0) of the DC component
%     of the signal
%
   fourierTermMag( 1 ) = signal( 1, 1 ) / denom;
   fourierTermPhase( 1 ) = 0.0;
%
%  ...Calculate all harmonics requested by the user
%
   for index = 2:nHarmonics
%
%  ...Capture the complex # representation of one harmonic
%
      ampPhase = 2.0 * signal( index ) / denom;
%
%  ...Convert to magnitude and phase angle
%
      fourierTermMag( index ) = abs( ampPhase );
      fourierTermPhase( index ) = angle( ampPhase );
   end
end
