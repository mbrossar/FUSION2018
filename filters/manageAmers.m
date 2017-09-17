function [S,PosAmers,ParamFilter,trackerBis,myTracks,PosAmersNew,...
    IdxAmersNew,trackCov,pointsMain,validityMain] = manageAmers(S,...
    PosAmers,ParamFilter,ParamGlobal,trackerBis,trajFilter,I,...
    pointsMain,validityMain,IdxImage,myTracks,pointsBis)

PosAmersNew = [];
IdxAmersNew = [];
trackCov = [];
MaxAmersNew = 10;

if sum(validityMain) < ParamFilter.NbAmersMin && I>120
    P = S'*S; %not computianny efficient
    NbAmersNew = min(sum(validityMain == 0),MaxAmersNew);
    [PosAmersNew,trackPoints,trackerBis,myTracks,trackCov] = ObserveAmersNew(ParamFilter,ParamGlobal,...
        trajFilter,I,IdxImage,trackerBis,myTracks,NbAmersNew,pointsMain,pointsBis,S);
    j = 1;
    IdxAmersNew = zeros(length(trackCov),1);
    IdxAmersOld = find(validityMain == 0);
    %if number of new landmarks is small
    IdxAmersOld = IdxAmersOld(1:length(trackCov));
    for ii = 1:length(IdxAmersOld)
        idxAmersOld = IdxAmersOld(ii);
        idxP = 15+(3*idxAmersOld-2:3*idxAmersOld);
        P(:,idxP) = 0;
        P(idxP,:) = 0;
        P(idxP,idxP) = trackCov{j};
        PosAmers(:,idxAmersOld) = PosAmersNew(j,:)';
        pointsMain(idxAmersOld,:) = trackPoints(:,j);
        validityMain(idxAmersOld) = 1;
        IdxAmersNew(j) = idxAmersOld;
        j = j+1;
    end
    if sum(pointsMain(:) < 0) > 0
        validityMain(pointsMain(:,1)<0) = 0;
        validityMain(pointsMain(:,2)<0) = 0;
        pointsMain(pointsMain(:)<0) = 1;
    end
    S = chol(P); %not computianny efficient
end
end

%--------------------------------------------------------------------------
function [posAmersNew,trackPoints,trackerBis,myTracks,trackCov] = ...
    ObserveAmersNew(ParamFilter,ParamGlobal,trajFilter,...
    I,IdxImage,trackerBis,myTracks,NbAmersNew,pointsMain,pointsBis,S)
% compute estimated locations for new landmarks
cameraParams = ParamFilter.cameraParams;
dirImage = ParamGlobal.dirImage;
fileImages = ParamGlobal.fileImages;
image = strcat(dirImage,int2str(fileImages(IdxImage)),'.png');
image = undistortImage(imread(image),cameraParams);

Newpoints = detectMinEigenFeatures(image);
Newpoints = selectUniform(Newpoints,NbAmersNew+150,size(image));

%number of views of candidate points
nbViews = zeros(1,length(myTracks));
for i = 1:length(myTracks)
    nbViews(i) = length(myTracks(i).ViewIds);
end

% tracking new points
posAmersNew = zeros(NbAmersNew,3);
trackPoints = ones(2,NbAmersNew);
trackCov = cell(NbAmersNew,1);
i = 1;
nbViewsMin = 7;
errorMax = 0.5;
PixelMin = 30;
while i <= NbAmersNew
    ok = 0;
    % find possible new point
    while ok == 0
        nbViews2 = find(nbViews>nbViewsMin);
        if isempty(nbViews2)
            if(nbViewsMin > 3)
                nbViewsMin = nbViewsMin - 1;
            end
            [trackerBis,myTracks] = razTrackerBis(ParamGlobal,ParamFilter,I,IdxImage,myTracks);
            for ii = 1:length(myTracks)
                nbViews(ii) = length(myTracks(ii).ViewIds);
            end
            errorMax = errorMax+1;
            nbViews2 = find(nbViews>nbViewsMin);
        end
        idx = randsample(nbViews2,1);
        ok = 1;
        pNew = myTracks(idx).Points(end,:);
        idxViewNew = myTracks(idx).ViewIds(end);
        for ii = 1:ParamFilter.NbAmers
            if norm(pNew-pointsMain(ii,:)) < PixelMin || idxViewNew < I
                ok = 0;
                break
            end
        end
        nbIdx = nbViews(idx);
        nbViews(idx) = 0;
    end
    
    % estimate location
    iReal = myTracks(idx).ViewIds(1);
    Rot = squeeze(trajFilter.Rot(:,:,iReal));
    x = trajFilter.x(:,iReal);
    camPoses = struct('ViewId',myTracks(idx).ViewIds(1),...
        'Orientation',Rot,'Location',x,'S',S);
    try
        for ii = 2:nbViewsMin% to be more time efficient (else use 2:nbIdx)
            iReal = myTracks(idx).ViewIds(round(ii/nbViewsMin*nbIdx));
            Rot = squeeze(trajFilter.Rot(:,:,iReal));
            x = trajFilter.x(:,iReal);
            camPoses(ii).ViewId = iReal;
            camPoses(ii).Orientation = Rot;
            camPoses(ii).Location = x;
            camPoses(ii).S = S;
        end
        myTracks(idx).ViewIds = myTracks(idx).ViewIds([1 round((2:nbViewsMin)*nbIdx/nbViewsMin)]);
        myTracks(idx).Points = myTracks(idx).Points([1 round((2:nbViewsMin)*nbIdx/nbViewsMin)],:);
        [xyzPoint,covariance] = myEsimatePosAmers(myTracks(idx),camPoses,cameraParams,ParamFilter);
        errors = sum(diag(covariance(end-2:end,end-2:end)));
    catch
        xyzPoint = ones(1,3);
        errorMax = errorMax*2;
        errors = errorMax+1;
        for iii = 1:length(myTracks)
            nbViews(iii) = length(myTracks(iii).ViewIds);
        end
        PixelMin = PixelMin/2;
    end
    
    % add is error is suficiently small
    if isempty(Newpoints)
        Newpoints = detectMinEigenFeatures(image);
        Newpoints = selectUniform(Newpoints,NbAmersNew+150,size(image));
    end
    idxNew = randi(length(Newpoints),1);
    if (errors < errorMax)
        posAmersNew(i,:) =  xyzPoint';
        trackPoints(:,i) = myTracks(idx).Points(end,:)';
        trackCov{i} = 3*10^-3*eye(3);%covariance(end-2:end,end-2:end);%2*10^-3*eye(3);
        myTracks(idx).ViewIds = I;
        myTracks(idx).Points = Newpoints(idxNew).Location;
        pointsBis(idx,:) = Newpoints(idxNew).Location;
        Newpoints(idxNew) = [];
        i = i+1;
    end
end
pointsBis(pointsBis(:)<=0) = 1;
trackerBis.setPoints(pointsBis);
end


