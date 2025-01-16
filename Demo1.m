% Demo #1  Demonstrating differenr relaziations of the meniscal tissue
% reconstruction
clear; close all;
addpath('funcs');

%% Data Loading and Parameter Setup
% Load sample data
load('A.mat');
tissue_volume = logical(A > 0);

% Define parameters
params = struct();
params.TargetPorosity = 0.5;
params.TargetSurfaceArea = 0.5;
params.TargetPoreSize = 0.5;
params.SampleSize = [48, 48, 48];
params.VolumeSize = [115, 115, 200];
params.RotationAngle = 8;
params.MemoryAllowed = 0.5;
params.DownsampleFactor = 8;

% Create parameter sets with different seeds
params1 = params;
params1.RandomSeed = 42;
params2 = params;
params2.RandomSeed = 123;

% Convert structs to name-value pairs
paramCell1 = [fieldnames(params1)'; struct2cell(params1)'];
paramCell1 = paramCell1(:)';
paramCell2 = [fieldnames(params2)'; struct2cell(params2)'];
paramCell2 = paramCell2(:)';

%% Perform Calculations
% Calculate original porosity
original_porosity = sum(tissue_volume(:)) / numel(tissue_volume);

% Run first reconstruction
fprintf('Running first reconstruction...\n');
final_vol1 = TDSYN('InputVolume', tissue_volume, paramCell1{:});
synth1_porosity = sum(final_vol1(:)) / numel(final_vol1);
fprintf('First reconstruction complete.\n');

% Run second reconstruction
fprintf('Running second reconstruction...\n');
final_vol2 = TDSYN('InputVolume', tissue_volume, paramCell2{:});
synth2_porosity = sum(final_vol2(:)) / numel(final_vol2);
fprintf('Second reconstruction complete.\n');

%% Print Results
fprintf('\nPorosity Comparison:\n');
fprintf('Original: %.3f\n', original_porosity);
fprintf('Reconstruction 1: %.3f\n', synth1_porosity);
fprintf('Reconstruction 2: %.3f\n', synth2_porosity);

%% Visualization
% Create figure with subplots
figure('Position', [100 100 1200 400]);

% Plot original volume
subplot(1,3,1);
surfacex(tissue_volume);
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Original Volume\nPorosity: %.3f', original_porosity));

% Plot first reconstruction
subplot(1,3,2);
surfacex(final_vol1);
title(sprintf('Reconstruction 1 (Seed = 42)\nPorosity: %.3f', synth1_porosity));
xlabel('X'); ylabel('Y'); zlabel('Z');

% Plot second reconstruction
subplot(1,3,3);
surfacex(final_vol2);
title(sprintf('Reconstruction 2 (Seed = 123)\nPorosity: %.3f', synth2_porosity));
xlabel('X'); ylabel('Y'); zlabel('Z');

% Adjust figure properties
colormap('parula');
set(gcf, 'Color', 'white');
sgtitle('Comparison of Original and Reconstructed Volumes', 'FontSize', 14);

print('img/Meniscal_tissue_realiztions.png','-dpng','-r200');