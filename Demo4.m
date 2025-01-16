% Demo #4: Demonstrating pore size control in Berea sandstone reconstruction
% Raw data from Figshare: Micro-CT image of Berea sandstone and extracted networks
% (Dong and Blunt, 2009, DOI: http://dx.doi.org/10.1103/PhysRevE.80.036307)

clear; close all;
addpath('funcs');

%% Data Loading and Parameter Setup
load('B.mat');
tissue_volume = logical(A > 0);

% Define base parameters
params = struct();
params.SampleSize = [48, 48, 48];
params.VolumeSize = [115, 115, 200];
params.RotationAngle = 8;
params.MemoryAllowed = 0.5;
params.DownsampleFactor = 8;
params.RandomSeed = 42;
params.TargetPorosity = 0.5;
params.TargetSurfaceArea = 0.5;

% Define pore size controller values (gamma)
gamma_values = [0.4, 0.5, 0.6];
reconstructions = cell(1, length(gamma_values));

%% Generate reconstructions with different pore size controllers
fprintf('Generating reconstructions with different pore size controllers...\n');

for i = 1:length(gamma_values)
    fprintf('\nGenerating reconstruction %d/%d (γ = %.2f)\n', ...
        i, length(gamma_values), gamma_values(i));
    
    % Set pore size controller
    params.TargetPoreSize = gamma_values(i);
    
    % Convert struct to name-value pairs
    paramCell = [fieldnames(params)'; struct2cell(params)'];
    paramCell = paramCell(:)';
    
    % Generate reconstruction
    reconstructions{i} = TDSYN('InputVolume', tissue_volume, paramCell{:});
    
    % Calculate pore size value
    pore_space = ~reconstructions{i};
    dist_transform = bwdist(~pore_space);
    max_pore_size = max(dist_transform(:));
    fprintf('Max pore size: %.3f\n', max_pore_size);
end

%% Visualization
figure('Position', [100 100 1200 400]);
sgtitle('Berea Sandstone Reconstructions with Different Pore Size Controllers (γ)', ...
    'FontSize', 14);

% Reconstructions
for i = 1:length(gamma_values)
    subplot(1, 3, i);
    surfacex(reconstructions{i});
    pore_space = ~reconstructions{i};
    dist_transform = bwdist(~pore_space);
    max_pore_size = max(dist_transform(:));
    title(sprintf('γ = %.1f\nMax Pore Size = %.3f', gamma_values(i), max_pore_size));
    xlabel('X'); ylabel('Y'); zlabel('Z');
end

% Adjust figure properties
colormap('parula');
set(gcf, 'Color', 'white');
print('img/PoreSize.png','-dpng','-r200');

