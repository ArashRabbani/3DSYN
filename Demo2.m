% Demo #2: Demonstrating porosity feature control in Berea sandstone reconstruction
% Raw data from Figshare: Micro-CT image of Berea sandstone and extracted networks 
% (Dong and Blunt, 2009, DOI: http://dx.doi.org/10.1103/PhysRevE.80.036307)

clear; close all;
addpath('funcs');

% Data Loading and Parameter Setup
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
params.TargetSurfaceArea = 0.5;
params.TargetPoreSize = 0.5;

% Define porosity controller values (alpha)
alpha_values = [0.4, 0.5, 0.6];
reconstructions = cell(1, length(alpha_values));

% Generate reconstructions with different porosity controllers
fprintf('Generating reconstructions with different porosity controllers...\n');

for i = 1:length(alpha_values)
    fprintf('\nGenerating reconstruction %d/%d (α = %.2f)\n', ...
        i, length(alpha_values), alpha_values(i));
    
    % Set porosity controller
    params.TargetPorosity = alpha_values(i);
    
    % Convert struct to name-value pairs
    paramCell = [fieldnames(params)'; struct2cell(params)'];
    paramCell = paramCell(:)';
    
    % Generate reconstruction
    reconstructions{i} = TDSYN('InputVolume', tissue_volume, paramCell{:});
    
    % Calculate porosity value
    porosity = 1-sum(reconstructions{i}(:)) / numel(reconstructions{i});
    fprintf('Porosity value: %.3f\n', porosity);
end

% Visualization
figure('Position', [100 100 1200 400]);
sgtitle('Berea Sandstone Reconstructions with Different Porosity Controllers (α)', ...
    'FontSize', 14);

for i = 1:length(alpha_values)
    subplot(1, 3, i);
    surfacex(reconstructions{i});
    porosity = 1-sum(reconstructions{i}(:)) / numel(reconstructions{i});
    title(sprintf('α = %.1f\nPorosity = %.3f', alpha_values(i), porosity));
    xlabel('X'); ylabel('Y'); zlabel('Z');
end
print('img/Porosity.png','-dpng','-r200');
