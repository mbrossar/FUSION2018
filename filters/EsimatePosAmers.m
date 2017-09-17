function [points3d, errors] = EsimatePosAmers(pointTracks, ...
    camPoses, cameraParams)

numTracks = numel(pointTracks);
points3d = zeros(numTracks, 3);
numCameras = size(camPoses, 2);
cameraMatrices = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
for i = 1:numCameras
    id = camPoses(i).ViewId;
    R  = camPoses(i).Orientation;
    t  = camPoses(i).Location;
    size_t = size(t);
    if size_t(1) == 3
        t = t';
    end
    cameraMatrices(id) = cameraMatrix(cameraParams, R', -t*R');
end

for i = 1:numTracks
    track = pointTracks(i);
    points3d(i, :) = triangulateOnePoint(track, cameraMatrices);
end

if nargout > 1
    [~, errors] = reprojectionErrors(points3d, cameraMatrices, pointTracks);
end

%--------------------------------------------------------------------------
function point3d = triangulateOnePoint(track, cameraMatrices)

% do the triangulation
numViews = numel(track.ViewIds);
A = zeros(numViews * 2, 4);
for i = 1:numViews
    id = track.ViewIds(i);
    P = cameraMatrices(id)';
    A(2*i - 1, :) = track.Points(i, 1) * P(3,:) - P(1,:);
    A(2*i    , :) = track.Points(i, 2) * P(3,:) - P(2,:);
end

[~,~,V] = svd(A);
X = V(:, end);
X = X/X(end);
point3d = X(1:3)';

%--------------------------------------------------------------------------
function [errors, meanErrorsPerTrack] = reprojectionErrors(points3d, ...
    cameraMatrices, tracks)
numPoints = size(points3d, 1);
points3dh = [points3d, ones(numPoints, 1)];
meanErrorsPerTrack = zeros(numPoints, 1);
errors = [];
for i = 1:numPoints
    p3d = points3dh(i, :);
    reprojPoints2d = reprojectPoint(p3d, tracks(i).ViewIds, cameraMatrices);
    e = sqrt(sum((tracks(i).Points - reprojPoints2d).^2, 2));
    meanErrorsPerTrack(i) = mean(e);
    errors = [errors; e]; 
end

%--------------------------------------------------------------------------
function points2d = reprojectPoint(p3dh, viewIds, cameraMatrices)
numPoints = numel(viewIds);
points2d = zeros(numPoints, 2);
for i = 1:numPoints
    p2dh = p3dh * cameraMatrices(viewIds(i));
    points2d(i, :) = p2dh(1:2) ./ p2dh(3);
end

