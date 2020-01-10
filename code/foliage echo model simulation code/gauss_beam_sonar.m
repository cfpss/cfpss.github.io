function [ Zgain ] = gauss_beam_sonar(amp,freqorbmwdth,freqorbmwdthVal,beamcenter,varargin)
%GAUSS_BEAM_SONAR: Gaussian approximation of actual bat beampatterns mainlobe
%Inputs: amp: Peak Amplitude of the beampattern
%        freqorbmwdth: 1) 'incident_freq' for entering frequency of incidence
%        b/w 60-80 kHz in kHz to estimate spread(sigma) of gauss beampattern
%        based on horseshoe data, 2) 'beamwidth' for directly entering
%        beamwidth in degrees and estimate spread(sigma) of gauss
%        beampattern for the same
%        freqorbmwdthVal: frequency (kHz) b/w 60-80 kHz for 'incident_freq' in freqorbmwdth
%        and beamwidth (degrees) for 'beamwidth' in freqorbmwdth.
%        center_x & center_y : coordinates of gaussian beam center
%        (default- 0,0), beamcenter
%    
%        varagrin : add amplitude of gaussian white noise if needed
%Outputs: Zgain : gain matrix

% Make a 2D Gaussian Kernel in cartesian coordinates
A = amp;
x0 = beamcenter(1); y0 = beamcenter(2);

% Predict gaussian spread (standard deviation) based on horseshoe bat data
sigma_x = predict_gauss_sigma(freqorbmwdth,freqorbmwdthVal); % get sigma for beamwidth corresponding to frequency of incidence
sigma_y = predict_gauss_sigma(freqorbmwdth,freqorbmwdthVal); % get sigma for beamwidth corresponding to frequency of incidence

[X, Y] = meshgrid(linspace(-pi,pi,361),linspace(-pi/2,pi/2,181));

theta = 0;
a = ((cos(theta))^2)/(2*(sigma_x^2)) + ((sin(theta))^2)/(2*(sigma_y)^2);
b = (-sin(2*theta))/(4*(sigma_x)^2) + (sin(2*theta))/(4*(sigma_y)^2);
c = ((sin(theta))^2)/(2*(sigma_x)^2) + ((cos(theta))^2)/(2*(sigma_y)^2);
Z = A*exp( -((a*(X-x0).^2) + (2*b*(X-x0).*(Y-y0)) + (c*(Y-y0).^2)));
disp(size(varargin))
    if(~isempty(cell2mat(varargin)))
        %%% Add gaussian noise (white noise) to beampattern
        Znoise = (varargin{1})*randn(size(Z));
        Zgain = Z+Znoise;
    else
        Zgain = Z;
    end    
end


