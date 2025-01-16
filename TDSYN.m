function [finalVolume, reconstructedVolume]=TDSYN(varargin)
    % TDSYN - 3D Microstructure image synthesis using patch-based sampling and reconstruction
    %
    % Syntax:
    %   [finalVolume, reconstructedVolume] = TDSYN()
    %   [finalVolume, reconstructedVolume] = TDSYN('ParameterName', Value, ...)
    %
    % Input Parameters (Name-Value Pairs):
    %   'TargetPorosity'     - Target porosity value (default: 0.5)
    %   'TargetSurfaceArea'  - Target surface area value (default: 0.5)
    %   'TargetPoreSize'     - Target pore size value (default: 0.5)
    %   'SampleSize'         - Size of sampling window [x,y,z] (default: [40,40,40])
    %   'VolumeSize'         - Size of output volume [x,y,z] (default: [100,100,200])
    %   'RotationAngle'      - Angle for rotation transformations (default: 8)
    %   'MemoryAllowed'      - Maximum memory allowed in GB (default: 0.5)
    %   'InputVolume'        - Input 3D binary volume (required)
    %   'DownsampleFactor'   - Factor for downsampling (default: 10)
    %   'RandomSeed'         - Seed for random number generator (default: 5)
    
    % Parse input parameters
    p = inputParser;
    addParameter(p, 'TargetPorosity', 0.5, @isnumeric);
    addParameter(p, 'TargetSurfaceArea', 0.5, @isnumeric);
    addParameter(p, 'TargetPoreSize', 0.5, @isnumeric);
    addParameter(p, 'SampleSize', [40, 40, 40], @isnumeric);
    addParameter(p, 'VolumeSize', [100, 100, 200], @isnumeric);
    addParameter(p, 'RotationAngle', 8, @isnumeric);
    addParameter(p, 'MemoryAllowed', 0.5, @isnumeric);
    addParameter(p, 'InputVolume', [], @(x) isnumeric(x) || islogical(x));
    addParameter(p, 'DownsampleFactor', 10, @isnumeric);
    addParameter(p, 'RandomSeed', 5, @isnumeric);
    
    parse(p, varargin{:});
    
    % Extract parameters
    targetPorosity = p.Results.TargetPorosity;
    targetSurfaceArea = p.Results.TargetSurfaceArea;
    targetPoreSize = p.Results.TargetPoreSize;
    sampleSize = round(p.Results.SampleSize); % Ensure integer
    volumeSize = round(p.Results.VolumeSize); % Ensure integer
    rotationAngle = p.Results.RotationAngle;
    Memory_allowed = p.Results.MemoryAllowed;
    % Convert input volume to logical and then to single
    originalVolume = single(logical(p.Results.InputVolume));
    downsampleFactor = round(p.Results.DownsampleFactor); % Ensure integer
    randomSeed = p.Results.RandomSeed;
    
    % Check if input volume is provided
    if isempty(originalVolume)
        error('Input volume is required.');
    end
    
    % Ensure sample size is divisible by downsample factor
    sampleSize = ceil(sampleSize / downsampleFactor) * downsampleFactor;
    
    % Calculate number of samples based on memory constraint
    numSamples = floor((Memory_allowed*1024^3)/(4*7*prod(sampleSize))); % N for single precision (4 bytes)
    
    % Set random seed
    rng(randomSeed ,'threefry')
    fprintf('Creating the patch dataset...\n');
    % Generate distance maps for different morphological operations
    baseDistanceMap = single(normalize(bwdist(originalVolume))-normalize(bwdist(1-originalVolume))+1);
    
    openedVolume = morphologicalOperation(originalVolume, 0.75, 'open');
    openedDistanceMap = single(normalize(bwdist(openedVolume))-normalize(bwdist(1-openedVolume))+1);
    
    closedVolume = morphologicalOperation(originalVolume, 0.75, 'close');
    closedDistanceMap = single(normalize(bwdist(closedVolume))-normalize(bwdist(1-closedVolume))+1);
    
    erodedVolume = morphologicalOperation(originalVolume, 0.75, 'ero');
    erodedDistanceMap = single(normalize(bwdist(erodedVolume))-normalize(bwdist(1-erodedVolume))+1);
    
    dilatedVolume = morphologicalOperation(originalVolume, 0.75, 'dil');
    dilatedDistanceMap = single(normalize(bwdist(dilatedVolume))-normalize(bwdist(1-dilatedVolume))+1);


    % Generate samples with different transformations
    samples1 = sampler(baseDistanceMap, numSamples, sampleSize);
    
    rotatedMap1 = imrotate3(baseDistanceMap, rotationAngle, [0 0 1], 'nearest', 'crop');
    rotatedMap1 = trimEdges(rotatedMap1, round(size(rotatedMap1,1)*sind(10)));
    samples2 = sampler(rotatedMap1, numSamples, sampleSize);
    
    rotatedMap2 = imrotate3(baseDistanceMap, -rotationAngle, [0 0 1], 'nearest', 'crop');
    rotatedMap2 = trimEdges(rotatedMap2, round(size(rotatedMap2,1)*sind(10)));
    samples3 = sampler(rotatedMap2, numSamples, sampleSize);
    
    rotatedMap3 = imrotate3(openedDistanceMap, -rotationAngle, [0 1 0], 'nearest', 'crop');
    rotatedMap3 = trimEdges(rotatedMap3, round(size(rotatedMap3,1)*sind(10)));
    samples4 = sampler(rotatedMap3, numSamples, sampleSize);
    
    rotatedMap4 = imrotate3(closedDistanceMap, rotationAngle, [0 1 0], 'nearest', 'crop');
    rotatedMap4 = trimEdges(rotatedMap4, round(size(rotatedMap4,1)*sind(10)));
    samples5 = sampler(rotatedMap4, numSamples, sampleSize);
    
    rotatedMap5 = imrotate3(erodedDistanceMap, rotationAngle, [1 0 0], 'nearest', 'crop');
    rotatedMap5 = trimEdges(rotatedMap5, round(size(rotatedMap5,1)*sind(10)));
    samples6 = sampler(rotatedMap5, numSamples, sampleSize);
    
    rotatedMap6 = imrotate3(dilatedDistanceMap, -rotationAngle, [1 0 0], 'nearest', 'crop');
    rotatedMap6 = trimEdges(rotatedMap6, round(size(rotatedMap6,1)*sind(10)));
    samples7 = sampler(rotatedMap6, numSamples, sampleSize);
    
    % Combine all samples
    allSamples = cat(4, samples1, samples2, samples3, samples4, samples5, samples6, samples7);
    
    % Apply porosity corrections
    porosityScale = single(((rand(1,1,1,size(allSamples,4))*2)-1).*0.1);
    porosityScale = repmat(porosityScale, sampleSize(1), sampleSize(2), sampleSize(3), 1);
    porosityScale = porosityScale + 1;
    allSamples = allSamples .* porosityScale;
    
    randomNoise = single(((rand(size(allSamples))*2)-1).*0.1);
    allSamples = allSamples + randomNoise;
    
    % Clear intermediate variables
    clear samples1 samples2 samples3 samples4 samples5 samples6 samples7
    
    % Calculate features
    porosity = mean(mean(mean(allSamples>1,1),2),3);
    surfaceArea = single(abs(allSamples-1)<.15);
    specificSurface = mean(mean(mean(surfaceArea,1),2),3);
    maxPoreSize = max(max(max(allSamples,[],1),[],2),[],3);
    
    % Normalize features
    porosity = normalize(porosity);
    specificSurface = normalize(specificSurface);
    maxPoreSize = normalize(maxPoreSize);
    
    % Calculate distances to target properties
    targetProps = repmat([targetPorosity, targetSurfaceArea, targetPoreSize], numel(porosity), 1);
    actualProps = [squeeze(porosity) squeeze(specificSurface) squeeze(maxPoreSize)];
    distances = sqrt(sum((targetProps-actualProps).^2, 2));
    [~, sortedIndices] = sort(distances);
    
    % Select best matching samples
    selectedIndices = sortedIndices(1:round(numel(porosity)/8));
    selectedSamples = allSamples(:,:,:,selectedIndices);
    
    % Visualization of the feature space --------------------------
    % randomIndices = randi(numel(porosity), 1000, 1);
    % figure;
    % scatter3(actualProps(randomIndices,1), actualProps(randomIndices,2), actualProps(randomIndices,3), 20, actualProps(randomIndices,3));
    % xlabel('Porosity');
    % ylabel('Specific surface');
    % zlabel('Pore size');
    % axis square tight;
    % colormap winter;
    % hold on;
    % scatter3(actualProps(selectedIndices,1), actualProps(selectedIndices,2), actualProps(selectedIndices,3), 'red', 'filled');
    % legend({'All datapoints', 'Preferred datapoints'}, 'Location', 'north');
    % colorbar_handle = colorbar;
    % colorbar_handle.Label.String = 'Normalized pore size';
    
    % Rebuild final volume
    [reconstructedVolume, ~, errors] = Rebuild(selectedSamples, volumeSize, sampleSize, downsampleFactor);
    
    % Post-processing
    finalVolume = reconstructedVolume;
    finalVolume = imgaussfilt3(finalVolume, 0.75, 'Padding', 'circular');
    finalVolume = finalVolume < 1;
    

end
function MM=sampler(AB,N,s)
P=[randi(size(AB,1)-s(1),N,1) randi(size(AB,2)-s(2),N,1) randi(size(AB,3)-s(3),N,1)];
MM=single(zeros(s(1),s(2),s(3),N));
for I=1:N
    MM(:,:,:,I)=AB(P(I,1)+1:P(I,1)+s(1),P(I,2)+1:P(I,2)+s(2),P(I,3)+1:P(I,3)+s(3));  
end
end
function [reconstructedVolume, weightMap, errors] = Rebuild(samplePatches, targetSize, patchSize, downsampleRatio)
    % Initialize empty list for tracking used patches
    usedPatches = [];
    
    % Margin size for overlap between patches
    marginSize = 4;
    
    % Effective size after removing margins
    effectivePatchSize = patchSize - 2*marginSize;
    
    % Store original patches
    originalPatches = samplePatches;

    % Create mask for original resolution patches
    fullResMask = zeros(patchSize);
    fullResMask(1+marginSize:end-marginSize, 1+marginSize:end-marginSize, 1+marginSize:end-marginSize) = 1;
    fullResMask = 1 - fullResMask;
    fullResMask = fullResMask(:);
    
    % Initialize patch comparison matrix for full resolution
    fullResComparison = zeros(size(originalPatches, 4), sum(fullResMask==1));
    for patchIdx = 1:size(originalPatches, 4)
        currentPatch = originalPatches(:,:,:,patchIdx);
        currentPatch = currentPatch(:);
        currentPatch(fullResMask==0) = [];
        fullResComparison(patchIdx,:) = currentPatch;
    end

    % Create mask for downsampled patches
    downsampledMask = zeros(patchSize/downsampleRatio);
    downsampledMask(1+marginSize:end-marginSize, 1+marginSize:end-marginSize, 1+marginSize:end-marginSize) = 1;
    downsampledMask = 1 - downsampledMask;
    downsampledMask = downsampledMask(:);
    
    % Initialize downsampled patches and comparison matrix
    downsampledPatches = zeros(patchSize(1)/downsampleRatio, patchSize(2)/downsampleRatio, patchSize(3)/downsampleRatio, size(originalPatches, 4));
    downsampledComparison = zeros(size(originalPatches, 4), sum(downsampledMask==1));
    
    for patchIdx = 1:size(originalPatches, 4)
        currentPatch = originalPatches(:,:,:,patchIdx);
        currentPatch = imresize3(currentPatch, 1/downsampleRatio, 'nearest');
        downsampledPatches(:,:,:,patchIdx) = currentPatch;
        currentPatch = currentPatch(:);
        currentPatch(downsampledMask==0) = [];
        downsampledComparison(patchIdx,:) = currentPatch;
    end

    % Initialize output volume and weight map
    reconstructedVolume = zeros(targetSize + patchSize);
    weightMap = reconstructedVolume;
    
    % Calculate patch placement positions
    basePos1 = marginSize:effectivePatchSize(1):size(reconstructedVolume,1)-effectivePatchSize(1)-marginSize;
    basePos2 = marginSize:effectivePatchSize(2):size(reconstructedVolume,2)-effectivePatchSize(2)-marginSize;
    basePos3 = marginSize:effectivePatchSize(3):size(reconstructedVolume,3)-effectivePatchSize(3)-marginSize;
    
    % Generate all possible positions
    positionIdx = 1;
    for pos1 = 1:numel(basePos1)
        for pos2 = 1:numel(basePos2)
            for pos3 = 1:numel(basePos3)
                patchPositions(positionIdx,:) = [basePos1(pos1) basePos2(pos2) basePos3(pos3)];
                positionIdx = positionIdx + 1;
            end
        end
    end
    
    % Randomize placement order
    placementOrder = shuffleData([1:size(patchPositions,1)]');
    
    % Initialize progress tracking
    startTime = tic;
    barWidth = 50;  % Width of the progress bar
    lastPrintTime = 0;
    updateInterval = 1;  % Update every 0.1 seconds
    fprintf('Starting reconstruction process...\n');
    
    % Place patches
    for positionIdx = 1:size(patchPositions,1)
        if mod(positionIdx,1)==0
            usedPatches = [];
        end
        
        % Update progress display
        currentTime = toc(startTime);
        if currentTime - lastPrintTime > updateInterval || positionIdx == 1 || positionIdx == size(patchPositions,1)
            progress = positionIdx / size(patchPositions,1);
            filledWidth = round(barWidth * progress);
            bar = ['[' repmat('|', 1, filledWidth) repmat(' ', 1, barWidth - filledWidth) ']'];
            
            % Calculate estimated time remaining
            if positionIdx > 1
                timePerIteration = currentTime / positionIdx;
                remainingIterations = size(patchPositions,1) - positionIdx;
                estimatedTimeRemaining = timePerIteration * remainingIterations;
                timeStr = sprintf('ETA: %02d:%02d', floor(estimatedTimeRemaining/60), floor(mod(estimatedTimeRemaining,60)));
            else
                timeStr = 'ETA: calculating...';
            end
            
            % Print progress
            fprintf('\rProgress: %s %3.1f%% | %s', bar, progress*100, timeStr);
            lastPrintTime = currentTime;
        end
        
        currentPos = patchPositions(placementOrder(positionIdx),:);
        
        % Extract current region
        currentRegion = reconstructedVolume(currentPos(1)-marginSize+1:currentPos(1)+effectivePatchSize(1)+marginSize,...
                                         currentPos(2)-marginSize+1:currentPos(2)+effectivePatchSize(2)+marginSize,...
                                         currentPos(3)-marginSize+1:currentPos(3)+effectivePatchSize(3)+marginSize);
        occupiedMask = currentRegion > 0;

        % Find best matching patch
        [bestPatch, matchError, patchId, usedPatches] = FindSim(currentRegion, originalPatches, occupiedMask,...
                                                               downsampledPatches, downsampledMask, downsampledComparison,...
                                                               usedPatches, downsampleRatio);
        errors(positionIdx) = matchError;
        
        % Calculate blending weights
        blendWeights = bwdist(1-occupiedMask);
        blendWeights = blendWeights./max(blendWeights(:));
        blendWeights(isnan(blendWeights)) = 0;
        
        % Blend patches
        blendedRegion = (currentRegion).*blendWeights + bestPatch.*(1-blendWeights);
        
        % Update reconstructed volume
        reconstructedVolume(currentPos(1)-marginSize+1:currentPos(1)+effectivePatchSize(1)+marginSize,...
                          currentPos(2)-marginSize+1:currentPos(2)+effectivePatchSize(2)+marginSize,...
                          currentPos(3)-marginSize+1:currentPos(3)+effectivePatchSize(3)+marginSize) = blendedRegion;
        
        weightMap(currentPos(1)-marginSize+1:currentPos(1)+effectivePatchSize(1)+marginSize,...
                 currentPos(2)-marginSize+1:currentPos(2)+effectivePatchSize(2)+marginSize,...
                 currentPos(3)-marginSize+1:currentPos(3)+effectivePatchSize(3)+marginSize) = blendWeights;
        t=reconstructedVolume(1:targetSize(1), 1:targetSize(2), 1:targetSize(3)); 
        t(1,1,1)=.0001; 
        t(end,end,end)=2;  
    end
    
    % Clear the progress line and print final summary
    totalTime = toc(startTime);
    fprintf('\nReconstruction completed in %.1f seconds\n', totalTime);
    fprintf('Final average patch discrepancy error: %.4f\n', mean(errors));
    
    % Crop output to target size
    reconstructedVolume = reconstructedVolume(1:targetSize(1), 1:targetSize(2), 1:targetSize(3));
    weightMap = weightMap(1:targetSize(1), 1:targetSize(2), 1:targetSize(3));
end

function [bestMatchPatch, matchError, selectedPatchId, updatedUsedPatches] = FindSim(currentRegion, originalPatches, occupancyMask, downsampledPatches, downsampledMask, downsampledComparison, usedPatches, downsampleRatio)
    % If region is empty, return random patch
    if sum(occupancyMask(:)) == 0
        randomIndex = randi(size(originalPatches, 4), 1, 1);
        bestMatchPatch = originalPatches(:,:,:,randomIndex);
        matchError = 0;
        selectedPatchId = randomIndex;
        updatedUsedPatches = usedPatches;
        return;
    end

    % Downsample current region and mask
    downsampledRegion = imresize3(currentRegion, 1/downsampleRatio, 'nearest');
    downsampledOccupancy = imresize3(occupancyMask, 1/downsampleRatio, 'nearest');
    
    % Store original patches for final output
    fullResPatches = originalPatches;
    
    % Work with downsampled versions
    currentRegion = downsampledRegion;
    originalPatches = downsampledPatches;
    occupancyMask = downsampledOccupancy;
    
    % Prepare mask for comparison
    maskVector = occupancyMask(:);
    maskVector(downsampledMask == 0) = [];
    maskMatrix = repmat(maskVector, [1, size(originalPatches, 4)])';
    
    % Prepare region for comparison
    regionVector = currentRegion(:);
    regionVector(downsampledMask == 0) = [];
    regionMatrix = repmat(regionVector, [1, size(originalPatches, 4)])';
    
    % Calculate errors
    patchErrors = (abs(regionMatrix - downsampledComparison) .* maskMatrix);
    
    % Weight errors by mask
    maskWeights = occupancyMask(:);
    maskWeights(downsampledMask == 0) = [];
    maskWeightMatrix = repmat(maskWeights, [1, size(originalPatches, 4)])';
    
    % Calculate weighted average error
    averageErrors = sum(patchErrors .* maskWeightMatrix, 2) ./ sum(occupancyMask(:));
    
    % Sort patches by error
    [sortedErrors, sortedIndices] = sort(averageErrors);
    matchError = sortedErrors(1);
    
    % Remove already used patches
    usedMask = ismember(sortedIndices, usedPatches);
    sortedIndices(usedMask) = [];
    
    % Select best unused patch
    selectedPatchId = sortedIndices(1);
    updatedUsedPatches = [usedPatches selectedPatchId];
    
    % Handle case where no valid patches found
    if max(occupancyMask(:)) == 0
        selectedPatchId = randi(size(originalPatches, 4), 1, 1);
        matchError = 0;
    end
    
    % Return full resolution version of best patch
    bestMatchPatch = fullResPatches(:,:,:,selectedPatchId);
end
function shuffledArray = shuffleData(inputArray)
    % Randomly shuffles elements of an array while preserving its dimensions
    randomIndices = randperm(numel(inputArray));
    shuffledArray = reshape(inputArray(randomIndices), size(inputArray));
end

function [normalizedArray] = normalize(inputArray)
    % Normalizes array values to range [0,1]
    normalizedArray = double(inputArray);
    minVal = min(normalizedArray(:));
    maxVal = max(normalizedArray(:));
    
    % Return if array is constant
    if minVal == maxVal
        return;
    end
    
    normalizedArray = (normalizedArray - minVal) ./ (maxVal - minVal);
end

function [outputImage] = morphologicalOperation(inputImage, radius, operationType)
    % Applies morphological operations using distance transform
    % operationType can be: 'dilate', 'dil', 'erosion', 'ero', 'open', 'ope', 'close', 'clo'
    
    if strcmp(operationType, 'dil') || strcmp(operationType, 'dilate')
        outputImage = bwdist(inputImage);
        outputImage = outputImage <= (radius + 0.5);
    end
    
    if strcmp(operationType, 'ero') || strcmp(operationType, 'erosion')
        outputImage = bwdist(1-inputImage);
        outputImage = outputImage >= (radius + 0.5);
    end
    
    if strcmp(operationType, 'ope') || strcmp(operationType, 'open')
        outputImage = bwdist(1-inputImage);
        outputImage = outputImage >= (radius + 0.5);
        outputImage = bwdist(outputImage);
        outputImage = outputImage <= (radius + 0.5);
    end
    
    if strcmp(operationType, 'clo') || strcmp(operationType, 'close')
        outputImage = bwdist(inputImage);
        outputImage = outputImage <= (radius + 0.5);
        outputImage = bwdist(1-outputImage);
        outputImage = outputImage >= (radius + 0.5);
    end
end

function [trimmedArray] = trimEdges(inputArray, trimSize)
    % Trims specified number of elements from all edges of 2D or 3D array
    if ndims(inputArray) == 3
        trimmedArray = inputArray(trimSize+1:end-trimSize, ...
                                trimSize+1:end-trimSize, ...
                                trimSize+1:end-trimSize);
    end
    
    if ndims(inputArray) == 2
        trimmedArray = inputArray(trimSize+1:end-trimSize, ...
                                trimSize+1:end-trimSize);
    end
end

