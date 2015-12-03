%% 6.438 FALL 2015
%% project_part_II.m
%   main file for part 2 of project

clc; clear; close all;

% / / / / / / / / / / / / / / / / / / / / / / / / /
% / / do not modify below for final submission / / 
load('test_file'); % call test_file created by create_test_file.m
addpath(genpath([pwd '/supplementary_functions_part_II'])); %
% / / do not modify above for final submission / / 
% / / / / / / / / / / / / / / / / / / / / / / / / /


% Write your code here for choosing
%   - r_dope: doping rate (scalar)
%     There will be penalization in choosing high doping rates
%     in the final calculation of compression ratio.
%   - temp (optional): row vector storing nonnegative integers
%     e.g. temp = [1,10,100] that you may use for your decoder. 
%     There will be penalization in using this extra information 
%     in the final calculation of compression ratio. 
%     Leave it empty if your decoder doesn't use this

% r_dope = -0.25 * r + 0.205; % some doping rate defined here 
                % (can depend on some variables in test_file)
                
% r_dope = 1 / (2 * (r + 0.8)) - 0.307;
% r_dope = -0.22 * r + 0.17 + 1 / (100 * (r - 0.1));  
r_dope = -0.1 * r + 0.06 + 1 / (50 * (r - 0.18));  

if r < 0.3 % this is a pretty silly case, but whatever
    r_dope = 0.6;
end
r_dope = max(r_dope, 0.01);
% r_dope = min(r_dope, 0.2);

% r_dope = 0.5;

disp(['r_dope before information adjustment: ', num2str(r_dope)])


temp = build_temp(s_to_shotgun(s, num_shotgun, l_shotgun));
temp(:, 1) = temp(:, 1) / 1;
temp = round(temp);

% increase r_dope if temp doesn't have much information
r_dope = r_dope + 1.5*r_dope * (1 - mean(temp(:, 1) / 100 .* (l_shotgun - temp(:, 3)) / l_shotgun));
r_dope = min(r_dope, 0.6);

temp = [temp(:); num_shotgun];

disp(['r_dope: ', num2str(r_dope)])

% / / / / / / / / / / / / / / / / / / / / / / / / /
% / / do not modify below for final submission / / 
[phi_source,phi_code] = dope_source(s,r_dope); % dope potentials
tic
[s_hat] = project_part_II_decoder(...
    x,H,phi_source,phi_code,psi_source,temp ,s); % call decoder
toc
if(isempty(temp)); sc = 0; else sc = numel(de2bi(temp)); end
[rH,cH] = size(H);
r_normalized = (rH+cH*r_dope+sc)/cH; % compute normalized rate
fprintf(['* normalized rate = ' num2str(r_normalized) '\n']); % print rate
if(sum(sum(s-s_hat))==0) % lossless?
   fprintf('lossless compression successful\n'); 
   flag_success = 1;
else
   fprintf('compression is not successful... \n'); 
   flag_success = 0;
end
