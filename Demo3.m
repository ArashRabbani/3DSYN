% Demo #3: Demonstrating specific surface area control in Berea sandstone reconstruction
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
params.TargetPoreSize = 0.5;

% Define specific surface area controller values (beta)
beta_values = [0.4, 0.5, 0.6];
reconstructions = cell(1, length(beta_values));

%% Generate reconstructions with different specific surface area controllers
fprintf('Generating reconstructions with different specific surface area controllers...\n');

for i = 1:length(beta_values)
    fprintf('\nGenerating reconstruction %d/%d (β = %.2f)\n', ...
        i, length(beta_values), beta_values(i));
    
    % Set specific surface area controller
    params.TargetSurfaceArea = beta_values(i);
    
    % Convert struct to name-value pairs
    paramCell = [fieldnames(params)'; struct2cell(params)'];
    paramCell = paramCell(:)';
    
    % Generate reconstruction
    reconstructions{i} = TDSYN('InputVolume', tissue_volume, paramCell{:});
    
    % Calculate specific surface area using perimeter
    perim = bwperim(reconstructions{i});
    surface_area = sum(perim(:)) / numel(reconstructions{i});
    fprintf('Surface area value: %.3f\n', surface_area);
end

%% Visualization
figure('Position', [100 100 1200 400]);
sgtitle('Berea Sandstone Reconstructions with Different Surface Area Controllers (β)', ...
    'FontSize', 14);

% Reconstructions
for i = 1:length(beta_values)
    subplot(1, 3, i);
    surfacex(reconstructions{i});
    perim = bwperim(reconstructions{i});
    surface_area = sum(perim(:)) / numel(reconstructions{i});
    title(sprintf('β = %.1f\nSurface Area = %.3f', beta_values(i), surface_area));
    xlabel('X'); ylabel('Y'); zlabel('Z');
end
print('img/Specific_surface.png','-dpng','-r200');
