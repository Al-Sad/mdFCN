% Copyright (c) 2024 Mohammad Al-Sa'd
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% Email: mohammad.al-sad@helsinki.fi, alsad.mohamed@gmail.com
%
% The following reference should be cited whenever this script is used:
% Al‐Sa'd, M., Vanhatalo, S. and Tokariev, A., 2024. Multiplex dynamic
% networks in the newborn brain disclose latent links with neurobehavioral
% phenotypes. Human Brain Mapping, 45(2), https://doi.org/10.1002/hbm.26610
%
% Last Modification: 12-February-2024
%
% Description:
% It generates Bandpass filters via cascaded zero-phase low-pass and high-pass 
% infinite impulse response filters.

function flt = iir_filters(Fs)
Apass = [1 1 1 1 1];
Astop_highpass = [10 10 10 10 10];
Astop_lowpass  = [10 10 10 10 10];
Fstop_highpass = [0.3 1.3 3.8 7.7 12.7];
Fpass_highpass = [0.4 1.5 4 8 13];
Fstop_lowpass  = [1.7 4.3 8.2 13.3 22.4];
Fpass_lowpass  = [1.5 4 8 13 22];
flt = cell(2,length(Fstop_lowpass));
for i = 1:length(Fstop_highpass)
    flt{1,i} = iir_highpass(Fstop_highpass(i), Fpass_highpass(i),...
        Astop_highpass(i), Apass(i), Fs);
end
for i = 1:length(Fstop_lowpass)
    flt{2,i} = iir_lowpass(Fstop_lowpass(i), Fpass_lowpass(i),...
        Astop_lowpass(i), Apass(i), Fs);
end
function h = iir_lowpass(Fstop, Fpass, Astop, Apass, Fs)
h = designfilt('lowpassiir','PassbandFrequency',Fpass, ...
         'StopbandFrequency',Fstop,'PassbandRipple',Apass, ...
         'StopbandAttenuation',Astop,'SampleRate',Fs,...
         'MatchExactly','stopband','DesignMethod','butter');
end
function h = iir_highpass(Fstop, Fpass, Astop, Apass, Fs)
h = designfilt('highpassiir','PassbandFrequency',Fpass, ...
         'StopbandFrequency',Fstop,'PassbandRipple',Apass, ...
         'StopbandAttenuation',Astop,'SampleRate',Fs,...
         'MatchExactly','stopband','DesignMethod','butter');
end
end