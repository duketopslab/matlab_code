%% Phantom study script for ps performed 10FEB2015, BRP NIH
%last modified: 10FEB2015
clear all;
clc;
sample_path = ['C:\Dropbox\Research\Data\2015_02_10_phantom_study_brp\'];
lambda = [470 480 490 500 510 560 580 600]'; %in nm
cd(sample_path);
%% Std deviation of measurements
sample_std = [];
%Start, before cleaning
clc;
clear sample
sample_nums = [1:2];

ITs = [50, 50, 50, 50, 50, 70, 70, 70, 85, 100, 120, 140];

num_phantoms = numel(ITs);
num_samps = num_phantoms;

mpix = 'c';
pixel_number = 4; %e is 10, c is 4, s is 

for i = 1:num_phantoms
    sample(i).base_name = ['phantom_study'];
    sample(i).num = sample_nums;
    sample(i).IT = ITs(i);
end

%% Set Spectralon

clc;
spec.path = sample_path;
spec.base_name = ['phantom_study_puck_am'];
spec.start_index = 1;
spec.reflectance = 99; %in percent
spec.IT = 50;

%% Data Correction Script
clc;
repeats = 10;
offset = 2;
cal_factor = 2.5; %this is basically a fudge factor because the reflectance of some phantoms is higher than the puck. 

puck = load([spec.path spec.base_name '_' num2str(spec.start_index) '_' num2str(spec.IT) 'us_int_channel_1.txt']);
puck_bgd = load([spec.path spec.base_name '_' num2str(spec.start_index) '_' num2str(spec.IT) 'us_int_bgd_channel_1.txt']);

puck = compute_array_averages_gen2( puck(:,2:17) , offset, repeats, '');
puck = puck(:,2:17).*cal_factor;
puck_bgd = mean(puck_bgd(:,2:17));
puck_bgd = ones(8,1)*puck_bgd;
puck_bgd_corrected = puck-puck_bgd;


for n = 1:num_samps
    for i = 1:numel(sample(n).num);
        name = [sample(i).base_name '_' num2str(n) '_' mpix '_' num2str(sample(n).num(i)) '_' num2str(sample(n).IT) 'us'];

        intensity = load([name '_int_channel_1.txt']);
        intensity = compute_array_averages_gen2( intensity(:,2:17) , offset, repeats, '');
        intensity = intensity(:,2:17);
        bgd = load([name '_int_bgd_channel_1.txt']);
        bgd = mean(bgd(:,2:17));
        bgd = ones(8,1)*bgd;
        sample(n).repeat(i).name = name;
        sample(n).repeat(i).int = intensity;
        sample(n).repeat(i).bgd = bgd;
        bgd_corrected = (intensity-bgd);
        sample(n).repeat(i).bgd_corrected = bgd_corrected;
        sample(n).repeat(i).reflectance = (bgd_corrected./puck_bgd_corrected)*(spec.IT./sample(n).IT);%.*(100/spec.reflectance));
    end
    
end

%% Look at reflectance values
%verify they appear in the right range

for n = 1:num_samps
    i=1;
    figure(1);
    plot(lambda, sample(n).repeat(i).reflectance(:,pixel_number));
    ylabel('R')
    xlabel('\lambda (nm)')
    hold on;
end




 
 %% Testing plots
 
 
i = 1; 
n = 1;

imagesc( reshape(mean(sample(n).repeat(i).reflectance(:,:), 1),4,4)') 
 
%% Optical Property Generation-------------------------------------------------------------------------
%% Constituent Volumes
clc;
lambda_vector = 450:1:600; 

%Enter actual mass values as measured on scale;

%stock absorber solution
stock_absorber_masses = [60.7];
stock_water_masses = [3011.4];
stock_water_volumes = stock_water_masses./1000;
stock_concentrations = stock_absorber_masses./stock_water_volumes;

%Cary measurement volumes (in the cuvette)
medium_masses = [3371.5];%typically water
stock_masses = [239.4];


%phantoms
absorber_values = [0.1 39.8 80.4 122.6 163.7 204.6 250.4 289.8 332.2 374.7 412.2 458.5];

water_values = [1792.9 1755.1 1723.7 1695.6 1670.0 1631.5 1605.2 1573.0 1549.8 1512.2 1479.2 1452.1];

beads_values = [209.2 195.2 187.4 175.7 164.4 153.1 143.2 131.9 122.4 111.3 101.1 88.3]; 

total_volumes = (water_values + absorber_values + beads_values)./1000;
%% Beads Properties
clc;
% mie_directory =  ['C:\Dropbox\Code\MATLAB\mc_model_src\mc_model_src-2.3\mc_model_src\mcinversion\'];
inversion_model_path = ['C:\Dropbox\Code\MATLAB\mc_model_src\mc_model_src-2.4\mc_model_src\mcinversion'];
cd(inversion_model_path);
load mie_gp;
% load([mie_directory 'mie_gp.mat'])

%Percent solid of beads
bead_properties.percent_solids = 2.62;
%Actual bead size - diameter in um
bead_properties.particle_diameter = .99;

%scattering properties
[ mie_properties ] = scattering_cross_section( lambda_vector, bead_properties.particle_diameter);

bead_properties.sigma_s = mie_properties.sigma_s;
bead_properties.sigma_s_prime = mie_properties.sigma_s_prime;
bead_properties.g = mie_properties.g; 

%[mus_lambda, mus_p_lambda ] = measurement_based_scattering_coefficient( beads_values, total_volumes,...
%    bead_properties,  lambda_vector);

mie_object = mie_poly_spheres_variable_refr_index;
[mus_lambda, mus_p_lambda, g ]= measurement_based_scattering_coefficient_mie_object( beads_values, total_volumes, bead_properties,  lambda_vector, mie_object);
bead_properties.g = g; 

%% Absorption properties
clc;
absorption_directory =  ['C:\Dropbox\Research\absorption_data\'];
titles = {'hemoglobin'};
load([absorption_directory 'abs_data'])

clc;
num_repeats = 1;
absorber_levels = 1;

% abs_var = load([absorption_directory 'Hb_10FEB2015_ps.txt']);
% 
%  
% [ extinction ] = compute_extinction( titles , abs_var, num_repeats, absorber_levels,...
%     stock_masses, medium_masses, stock_concentrations );
% 
concentrations = ((absorber_values./1000)*stock_concentrations)./total_volumes;
phantom_type = 'liquid';

extinction.lambda = HbO2(:,1);
extinction.molar_extinction_value = HbO2(:,2)./64500;
[mua_lambda] = measurement_based_absorption_coefficient( extinction, concentrations, lambda_vector, phantom_type);


%% Make phantom spec_data objects
clc;
pixels = 16;
clear phantom
% spectra = meas{:,i};
for n = 1:num_phantoms
    for i = 1:pixels
        ref_index = n+1;
        phantom(n).pixel(i) = spec_data; 
        phantom(n).pixel(i).reflectance.measured{1} = [lambda sample(n).repeat(1).reflectance(:,i)];
        phantom(n).pixel(i).optical_properties.expected_g = [lambda interp1(lambda_vector, g, lambda)];
        phantom(n).pixel(i).optical_properties.expected_mus = [lambda interp1(lambda_vector, mus_lambda(:,ref_index), lambda)];
        phantom(n).pixel(i).optical_properties.expected_mus_prime = [lambda interp1(lambda_vector, mus_p_lambda(:,ref_index), lambda)];
        phantom(n).pixel(i).optical_properties.expected_mua = [lambda interp1(lambda_vector, mua_lambda(:,ref_index), lambda)];
        
        figure(1);
        subplot(4,4,i);
        hold on
        %plot(phantom(n).pixel(i).optical_properties.expected_mus_prime(:,1), phantom(n).pixel(i).optical_properties.expected_mus_prime(:,2:end))
        plot(phantom(n).pixel(i).reflectance.measured{1}(:,1), phantom(n).pixel(i).reflectance.measured{1}(:,2))
        title(['Pixel: ' num2str(i) ]);
        
    end
end


%% Load inversion files
clc;
% dir = ['C:\Users\' user '\Dropbox\Code\MATLAB\Inversions'];
% refDir = ['C:\Users\' user '\Dropbox\Research\BRP\BRP Data\31JUL2012 Test_Data\'];
% cd([dir])
% inversionfp = '\MC inversion files\';
% addpath([dir inversionfp]);
absorption_directory =  ['C:\Dropbox\Research\absorption_data\'];

inversion_object_directory = ['W:\brp_scaling_obj\05may2015_workspace\'];
load mie_gp;
tic;
load([inversion_object_directory 'mcs'])
toc;


%% Setup phantom study data
clear mci

mcs_ref = 1;
mci = mc_inversion;
mie_object = mie_poly_spheres_variable_refr_index;
% one collection, one illumination ring
mci.use_rings = 1;%set to correct pixel number
% fill in reference
mci.reference.mc_scaling = mcs;
mci.target.mc_scaling = mcs;
% fill in target
mci.target.mie = mie_object;

mci.reflectance_settings.num_fits=100;

mci.reflectance_settings.preset_absorbers = [];
% mci.reflectance_settings.fitting_wavelengths = [470 480 490 500 510 560 580]';
mci.reflectance_settings.fitting_wavelengths = [470 480 490 500 510 560 580 600]';

%Add custom absorber
stockmua(:,1) = extinction.lambda;
stockmua(:,2) = extinction.molar_extinction_value.*log(10);
stockmua_interp(:,1) = lambda;
stockmua_interp(:,2) = interp1(stockmua(:,1), stockmua(:,2), lambda', 'spline', NaN);
mci.reflectance_settings.custom_absorber{1} = [stockmua_interp(:,1) stockmua_interp(:,2)];
figure;
plot(stockmua_interp(:,1), stockmua_interp(:,2))
mci.reflectance_settings.custom_absorber_label{1} = 'Hemoglobin 10FEB2015';

%% Fitting loop

ref = 1:num_phantoms;
for i = pixel_number %1:pixels    
    for j =  1:length(ref)
        for k = 1:length(ref)
            
            mci.reference.mc_scaling{i} = mcs;
            mci.target.mc_scaling{i} = mcs;
            
            samp = phantom(k).pixel(i);
            mci.reference.phantom = phantom(ref(j)).pixel(i);
            phantom(k).invResults.reference(ref(j)).inversion(i) = refl_inversion(mci,samp, [], [], 1);
            display([' pixel: ' num2str(i) ' target: ' num2str(k) ' ref: ' num2str(ref(j))])
        end
    end   
end


%%
error_threshold = 30;
subplot_length = sqrt(pixels);
for i = pixel_number
    for j=1:num_phantoms
        for k = 1:num_phantoms
        avg_mua_error(j,ref(k),i) = phantom(j).invResults.reference(ref(k)).inversion(i).reflectance.model_info.average_percent_error_mua;
        avg_mus_prime_error(j,ref(k),i) = phantom(j).invResults.reference(ref(k)).inversion(i).reflectance.model_info.average_percent_error_mus_prime;
        mua_error(j,ref(k),i)=  phantom(j).invResults.reference(ref(k)).inversion(i).reflectance.model_info.percent_rmse_mua;
        mus_prime_error(j,ref(k),i)= phantom(j).invResults.reference(ref(k)).inversion(i).reflectance.model_info.percent_rmse_mus_prime;

        end
    end
    figure(12);
   % subplot(subplot_length,subplot_length,i)    
    imagesc(avg_mua_error(:,:,i), [0 error_threshold])
    xlabel('Reference')
    ylabel('Target')
    title('\mu_a %Error')

    
    figure(13);
   % subplot(subplot_length,subplot_length,i)
    imagesc(avg_mus_prime_error(:,:,i), [0 error_threshold])
    xlabel('Reference')
    ylabel('Target')
    title('\mu_s'' %Error')
    
    

end


%%
clc;

r=1;
figure(16);
for i = pixel_number
    for j = 1:num_phantoms
        for k = 1:num_phantoms
            subplot(num_phantoms,num_phantoms,r)
            plot(phantom(j).invResults.reference(ref(k)).inversion(i).reflectance.measured{1}(:,1), phantom(j).invResults.reference(ref(k)).inversion(i).reflectance.measured{1}(:,2))
            hold on
            plot(phantom(j).invResults.reference(ref(k)).inversion(i).reflectance.fit{1}(:,1), phantom(j).invResults.reference(ref(k)).inversion(i).reflectance.fit{1}(:,2),'r')
           title(['p' num2str(j) 'r' num2str(k)])
            r=r+1;
        end
    end
end

%% Absorption Fits
clc;
r=1;
figure(17);
for i = pixel_number
    for j = 1:num_phantoms
        for k = 1:num_phantoms
            subplot(num_phantoms,num_phantoms,r)
            plot(phantom(j).invResults.reference(ref(k)).inversion(i).optical_properties.expected_mua(:,1), phantom(j).invResults.reference(ref(k)).inversion(i).optical_properties.expected_mua(:,2))
            hold on
            plot(phantom(j).invResults.reference(ref(k)).inversion(i).optical_properties.extracted_mua(:,1), phantom(j).invResults.reference(ref(k)).inversion(i).optical_properties.extracted_mua(:,2),'r')
           title(['p' num2str(j) 'r' num2str(k)])
            r=r+1;
        end
    end
end

%% Scattering Fits
clc;
r=1;
figure(18);
for i = pixel_number
    for j = 1:num_phantoms
        for k = 1:num_phantoms
            subplot(num_phantoms,num_phantoms,r)
            plot(phantom(j).invResults.reference(ref(k)).inversion(i).optical_properties.expected_mus_prime(:,1), phantom(j).invResults.reference(ref(k)).inversion(i).optical_properties.expected_mus_prime(:,2))
            hold on
            plot(phantom(j).invResults.reference(ref(k)).inversion(i).optical_properties.extracted_mus_prime(:,1), phantom(j).invResults.reference(ref(k)).inversion(i).optical_properties.extracted_mus_prime(:,2),'r')
           title(['p' num2str(j) 'r' num2str(k)])
            r=r+1;
        end
    end
end    
    
    
    
    
    
    
    