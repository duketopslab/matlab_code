function [mua] = measurement_based_absorption_coefficient( extinction, concentrations, lambda_vector, phantom_type)

user = getenv('UserName');
dye_path = ['C:\Dropbox\Research\absorption_data\'];
load([dye_path 'mua_gelatin_base']);
num_absorbers = numel(concentrations);


for i = 1:num_absorbers
    mua(i,:) = log(10).*concentrations(i).*interp1(extinction.lambda, extinction.molar_extinction_value, lambda_vector);
    switch phantom_type
        case 'solid'
            mua(i,:)  = mua(i,:)  + interp1(mua_gelatin_base(:,1), mua_gelatin_base(:,2).*65, lambda_vector);
        otherwise
    end
end

mua = [lambda_vector; mua]';

figure;
plot(lambda_vector, mua(:,2:end))
xlabel('\lambda (nm)');
ylabel('\mu_a (cm^-^1)')
title('Expected \mu_a')

