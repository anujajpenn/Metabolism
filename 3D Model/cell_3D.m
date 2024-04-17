
AR = 1:1:15;
Density = round(logspace(-0.5,1,30),3);

%% LINK COMSOL
import com.comsol.model.util.*;
ModelUtil.showProgress(true);

Energy = zeros(length(Density),length(AR));
Rho = zeros(length(Density),length(AR));
Stressvol = zeros(length(Density),length(AR));
Sigma11 = zeros(length(Density),length(AR));
Sigma22 = zeros(length(Density),length(AR));
Sigma33 = zeros(length(Density),length(AR));
Umech = zeros(length(Density),length(AR));
Cell_Umech = zeros(length(Density),length(AR));
Mat_Umech = zeros(length(Density),length(AR));
Ubinding = zeros(length(Density),length(AR));
UATP = zeros(length(Density), length(AR));
Umotor = zeros(length(Density),length(AR));
Cell_Strain = zeros(length(Density),length(AR));
Mat_Strain = zeros(length(Density),length(AR));
SPol = zeros(length(Density),length(AR));
CPol = zeros(length(Density),length(AR));

% Run the model for every specified combination of cell aspect ratio and
% matrix density
for i = 1:length(AR)
    for j = 1:length(Density)
        if ~exist('model', 'var')
            disp('Loading model...');
            model = mphopen('3Dcell.mph');
        end
        % Set Geometry and model parameters
        model.param.set('AR', AR(i));
        model.param.set('Rho0', 350);
        model.param.set('CellE', 700);
        model.param.set('CellNu', 0.4);
        model.param.set('Ef0', 10);
        model.param.set('LC', 1.001);
        model.param.set('T', 0.5);
        model.param.set('IntE', 30); 
        model.param.set('IntNu', 0.4);
        model.param.set('NucE', 1000);
        model.param.set('NucNu', 0.4);
        model.param.set('NesE', 850);
        model.param.set('NesNu', 0.4);
        model.param.set('n', 10);
        model.param.set('m', 15);
        model.param.set('DensityECM',Density(j));
        % Compute
        disp(['Running simulation for aspect ratio = ',num2str(AR(i)),', Density = ',num2str(Density(j)),' ...']);
        model.study('std1').run;
        % Calculate energies, stresses, and strains
        disp('Calculating Energy...');
        Energy(j, i) = mphint2(model, 'Utotal', 2, 'selection', 'all', 'intvolume', 'on')';
        Rho(j, i) = mphint2(model, 'RhoKK', 2, 'selection', 3, 'intvolume', 'on', 'Outersolnum', 'all')';
        Stressvol(j, i) = mphint2(model, 'SKK', 2, 'selection', 3, 'intvolume', 'on', 'Outersolnum', 'all')';
        Sigma11(j, i) = mphint2(model, 'solid.sl11', 2, 'selection', 3, 'intvolume', 'on', 'Outersolnum', 'all')';
        Sigma22(j, i) = mphint2(model, 'solid.sl22', 2, 'selection', 3, 'intvolume', 'on', 'Outersolnum', 'all')';
        Sigma33(j, i) = mphint2(model, 'solid.sl33', 2, 'selection', 3, 'intvolume', 'on', 'Outersolnum', 'all')';
        Mat_Strain(j, i) = mphint2(model, 'EKK2', 2, 'selection', 1, 'intvolume', 'on', 'Outersolnum', 'all')';
        Cell_Strain(j, i) = mphint2(model, 'EKK', 2, 'selection', 3, 'intvolume', 'on', 'Outersolnum', 'all')';
        Umech(j, i) = mphint2(model, 'Umech', 2, 'selection', 'all', 'intvolume', 'on', 'Outersolnum', 'all')';
        Cell_Umech(j, i) = mphint2(model, 'Umech', 2, 'selection', 3, 'intvolume', 'on', 'Outersolnum', 'all')';
        Mat_Umech(j, i) = mphint2(model, 'Umech', 2, 'selection', 1, 'intvolume', 'on', 'Outersolnum', 'all')';
        Ubinding(j, i) = mphint2(model, 'Ubinding', 2, 'selection', 3, 'intvolume', 'on', 'Outersolnum', 'all')';
        UATP(j, i) = mphint2(model, 'UATP', 2, 'selection', 3, 'intvolume', 'on', 'Outersolnum', 'all')';
        Umotor(j, i) = mphint2(model, 'Umotor', 2, 'selection', 3, 'intvolume', 'on', 'Outersolnum', 'all')';
        SPol(j, i) = mphint2(model, 'Spol', 2, 'selection', 3, 'intvolume', 'on', 'Outersolnum', 'all')';
        CPol(j, i) = mphint2(model, 'Cpol', 2, 'selection', 3, 'intvolume', 'on', 'Outersolnum', 'all')';
    end
end