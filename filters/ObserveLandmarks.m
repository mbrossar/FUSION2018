function [y,yAmers,trackerMain,trackerBis,pointsMain,validityMain,myTracks,pointsBis] = ...
    ObserveLandmarks(trackerMain,trackerBis,dirImage,IdxImage,...
    fileImages,ParamFilter,Rot,x,PosAmers,I,S,myTracks)
%% trackerMain
cameraParams = ParamFilter.cameraParams;
chiC = ParamFilter.chiC;
RotC = chiC(1:3,1:3);
xC = chiC(1:3,4);
Pi = ParamFilter.Pi;

image = strcat(dirImage,int2str(fileImages(IdxImage)),'.png');
image = undistortImage(imread(image),cameraParams);

[pointsMain,validityMain] = trackerMain.step(image);
validityMain(pointsMain(:,1)<=0) = 0;
validityMain(pointsMain(:,2)<=0) = 0;
pointsMain(pointsMain(:)<=0) = 1;

y = zeros(2*sum(validityMain),1);
yAmers = zeros(sum(validityMain),1);
EcartPixelMax = ParamFilter.EcartPixelMax;
EcartPixelOut = EcartPixelMax;
j = 1;

for i = 1:length(validityMain)
    if validityMain(i) == 1
        PosAmers_i = PosAmers(:,i);
        pointsEst = Pi*( (Rot*RotC)' * (PosAmers_i-(x+Rot*RotC*xC)) );
        pixelEst = pointsEst(1:2)/pointsEst(3);
        if norm(pixelEst-pointsMain(i,:)') < EcartPixelMax %reject outlier
            y(2*j-1:2*j) = pointsMain(i,:)';
            yAmers(j) = i;
            j = j+1;
        elseif norm(pixelEst-pointsMain(i,:)') > EcartPixelOut
            validityMain(i) = 0;
        end
    end
end
y(2*j-1:end) = [];
yAmers(j:end) = [];

%for landmarks too close to each others
for ii = 1:2
    idx = randsample(find(validityMain == 1),1);
    pointIdx = pointsMain(idx,:);
    mMin = 10;
    mDis = (mMin+1)*ones(ParamFilter.NbAmers,1);
    for i = 1:length(validityMain)
        if validityMain(i) == 1 && i ~= idx
            pointI = pointsMain(i,:);
            mDis(i)= norm(pointIdx-pointI);
        end
    end
    [mDis,idx2] = min(mDis);
    if mDis < mMin % remove the most uncertain location
        idxS = 15+(3*idx-2:3*idx);
        idxS2 = 15+(3*idx2-2:3*idx2);
        S1 = sum(diag(S(idxS,idxS)));
        S2 = sum(diag(S(idxS2,idxS2)));
        if S2 > S1
            idx = idx2;
        end
        validityMain(idx) = 0;
    end
end

%% trackerBis
[pointsBis,validityBis] = trackerBis.step(image);
validityBis(pointsBis(:,1)<=0) = 0;
validityBis(pointsBis(:,2)<=0) = 0;
pointsBis(pointsBis(:)<=0) = 1;

mMin = 4;
for ii = 1:5
    idx = randsample(find(validityBis == 1),1);
    pointIdx = pointsBis(idx,:);
    mDis = (mMin+1)*ones(ParamFilter.NbAmers,1);
    for i = 1:length(validityBis)
        if validityBis(i) == 1 && i ~= idx
            pointI = pointsBis(i,:);
            mDis(i)= norm(pointIdx-pointI);
        end
    end
    [mDis,idx2] = min(mDis);
    if mDis < mMin
        nbView1 = length(myTracks(idx).ViewIds);
        nbView2 = length(myTracks(idx2).ViewIds);
        if nbView2 < nbView1
            idx = idx2;
        end
        validityBis(idx) = 0;
    end
end

Newpoints = detectMinEigenFeatures(image);
Newpoints = selectUniform(Newpoints,sum(validityBis == 0)+30,size(image));
for i = 1:length(validityBis)
    if validityBis(i) == 1
        myTracks(i).ViewIds = [myTracks(i).ViewIds I];
        myTracks(i).Points  =  [myTracks(i).Points;pointsBis(i,:)];
    elseif validityBis(i) == 0
        myTracks(i).ViewIds = I;
        if isempty(Newpoints)
            Newpoints = detectMinEigenFeatures(image);
            Newpoints = selectUniform(Newpoints,sum(validityBis == 0)+30,size(image));
        end
        j = randi(length(Newpoints),1);
        location = Newpoints(j).Location;
        Newpoints(j) = [];
        myTracks(i).Points = location;
        pointsBis(i,:) = location;
    end
end
trackerBis.setPoints(pointsBis);
end

