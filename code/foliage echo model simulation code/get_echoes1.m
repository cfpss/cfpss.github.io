% This program is to produce echoes from leaves
% Inputs: 
%         frequency of the incident wave, frq, in kHz.
%
%         % LEAF INPUTS:
%         a, radius of the leaf
%         theta0, incident angle
%         Sonargain, the sound gain of each leaf from sonar
%         R_sonar_leaf, the distance between each leaf and sonar
%
%         % SONAR INPUTS:
%         frq, frequencies
%         no_bins, the number of bins in our interesting frequency range
%         [60, 80]kHz
% Outputs:
%         total_amp & total_phase, another interpretation of frequency domain
% Created by Chen Ming, Date: 06/01/2015

function [total_amp, total_phase] = get_echoes1(frq,a,theta0, no_bins,Sonargain,R_sonar_leaf)

%%% Calculate the total gain

% leaf beampatterns

[leaf_amp, leaf_phase] = leaf_beampattern( frq, a, theta0, no_bins);

%%% Calculate total amplitude gain
speed = 340;
lambda = speed./(frq*1E3);
k = 2*pi./lambda;
k = repmat(k',1,length(a));
dist_s_l = repmat(transpose(R_sonar_leaf(:)),no_bins,1);
sonar_gain = repmat(transpose(Sonargain(:)),no_bins,1);

total_amp = leaf_amp./(k.*dist_s_l).* sonar_gain.*0.2./dist_s_l;

%%% Calculate the total phase
total_phase = -2.*dist_s_l.*k - leaf_phase;

end