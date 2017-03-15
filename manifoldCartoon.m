%% makes a cartoon that shows a stimulation from intrinsic manifold

% figure 1: just decay back to threshold

% set up movie
shouldWrite = 1;
if shouldWrite
    saveLoc = 'D:\Presentations\2017 03 22\';
    movName = 'noDynamics';
    saveName = [saveLoc,movName, '.avi'];
    vidObj = VideoWriter(saveName);
    vidObj.FrameRate = 10;
    open(vidObj);
end
[x,y] = meshgrid(0:.05:2);
% make a random z surface that makes a pretty plot
z = (x.^2 - .4.*y.^3).^2;
% normalize to be between .5 and 2.5, don't want it to touch origin
z = z./max(z(:))*2+.5;

figH = figure('units', 'pix', 'outerposition', [0 50 1000 1000],...
    'Color', [1 1 1]);
surfH = surf(x,y,z);

direction = [0 0 1];
% rotate(surfH, direction, 30);

xlabel('FR 1', 'fontsize', 16, 'fontweight', 'bold');
ylabel('FR 2', 'fontsize', 16, 'fontweight', 'bold');
zlabel('FR 3', 'fontsize', 16, 'fontweight', 'bold');
axis([0 3 0 3 0 3]);

hold on;
scH = scatter3(x(1,1),y(1,1),z(1,1), 50,'filled');
scH.MarkerFaceColor = 'r';
if shouldWrite
    f = getframe(figH);
    writeVideo(vidObj,f);
end
N = size(x,1);
for i=1:50
    if i==1
        c1 = randi(N);
        c2 = randi(N);
    else
        % move by a few, can't be more than N or less than 0
        nToMove = 16;
        c1 = max([1 min([N c1+randi(nToMove)-nToMove/2])]);
        c2 = max([1 min([N c2+randi(nToMove)-nToMove/2])]);
        
        plotH = plot3([scH.XData, x(c1,c2)], [scH.YData, y(c1,c2)], ...
            [scH.ZData, z(c1,c2)], 'color', 'k', 'linewidth', 1.5);
        
    end
    scH.XData = x(c1,c2);
    scH.YData = y(c1,c2);
    scH.ZData = z(c1,c2);
    
    drawnow;
    pause(.1);
    if shouldWrite
        f = getframe(figH);
        writeVideo(vidObj,f);
    end
end
figH.Color = [0 0 1];
plotH = plot3([scH.XData, 3], [scH.YData, 3], ...
    [scH.ZData, 3], 'color', 'k', 'linewidth', 1.5);
scH.XData = 3;
scH.YData = 3;
scH.ZData = 3;

linearMat = [x(:) y(:) z(:)];
flat = cat(3,x,y,z);
euclid = @(X,v) (X(:,1)-v(1)).^2 + (X(:,2)-v(2)).^2 + (X(:,3)-v(3)).^2;

propToMove = .3;
distThresh = .01;
pause(.1);
if shouldWrite
    f = getframe(figH);
    writeVideo(vidObj,f);
end
figH.Color = [1 1 1];
for i=1:50
    if i==1
        curPoint = [3 3 3];
    end
    dists = euclid(linearMat, curPoint);
    [minDist,ind] = min(dists);
    vec = curPoint-linearMat(ind,:);
    if minDist>distThresh
        curPoint = curPoint-vec.*propToMove;
    else
        [c1,c2] = ind2sub(size(flat),ind);
        c1 = max([1 min([N c1+randi(nToMove)-nToMove/2])]);
        c2 = max([1 min([N c2+randi(nToMove)-nToMove/2])]);
        curPoint = linearMat(sub2ind(size(flat),c1,c2),:);
    end
    plotH = plot3([scH.XData, curPoint(1)], [scH.YData, curPoint(2)], ...
        [scH.ZData, curPoint(3)], 'color', 'k', 'linewidth', 1.5);
    scH.XData = curPoint(1);
    scH.YData = curPoint(2);
    scH.ZData = curPoint(3);
    drawnow;
    pause(.1);
    if shouldWrite
        f = getframe(figH);
        writeVideo(vidObj,f);
    end
end

if shouldWrite
    close(vidObj);
end

%% with some noise perturbations on its way back to manifold
shouldWrite = 1;
if shouldWrite
    movName = 'withDynamics';
    saveName = [saveLoc,movName, '.avi'];
    vidObj = VideoWriter(saveName);
    vidObj.FrameRate = 10;
    open(vidObj);
end
[x,y] = meshgrid(0:.05:2);
% make a random z surface that makes a pretty plot
z = (x.^2 - .4.*y.^3).^2;
% normalize to be between .5 and 2.5, don't want it to touch origin
z = z./max(z(:))*2+.5;

figH = figure('units', 'pix', 'outerposition', [0 50 1000 1000],...
    'Color', [1 1 1]);
surfH = surf(x,y,z);

direction = [0 0 1];
% rotate(surfH, direction, 30);

xlabel('FR 1', 'fontsize', 16, 'fontweight', 'bold');
ylabel('FR 2', 'fontsize', 16, 'fontweight', 'bold');
zlabel('FR 3', 'fontsize', 16, 'fontweight', 'bold');
axis([0 3 0 3 0 3]);

hold on;
scH = scatter3(x(1,1),y(1,1),z(1,1), 50,'filled');
scH.MarkerFaceColor = 'r';
if shouldWrite
    f = getframe(figH);
    writeVideo(vidObj,f);
end
N = size(x,1);
for i=1:50
    if i==1
        c1 = randi(N);
        c2 = randi(N);
    else
        % move by a few, can't be more than N or less than 0
        nToMove = 16;
        c1 = max([1 min([N c1+randi(nToMove)-nToMove/2])]);
        c2 = max([1 min([N c2+randi(nToMove)-nToMove/2])]);
        
        plotH = plot3([scH.XData, x(c1,c2)], [scH.YData, y(c1,c2)], ...
            [scH.ZData, z(c1,c2)], 'color', 'k', 'linewidth', 1.5);
        
    end
    scH.XData = x(c1,c2);
    scH.YData = y(c1,c2);
    scH.ZData = z(c1,c2);
    
    drawnow;
    pause(.1);
    if shouldWrite
        f = getframe(figH);
        writeVideo(vidObj,f);
    end
end
figH.Color = [0 0 1];
plotH = plot3([scH.XData, 3], [scH.YData, 3], ...
    [scH.ZData, 3], 'color', 'k', 'linewidth', 1.5);
scH.XData = 3;
scH.YData = 3;
scH.ZData = 3;

linearMat = [x(:) y(:) z(:)];
flat = cat(3,x,y,z);
euclid = @(X,v) (X(:,1)-v(1)).^2 + (X(:,2)-v(2)).^2 + (X(:,3)-v(3)).^2;

propToMove = .3;
distThresh = .01;
noiseMag = .2;
pause(.1);
if shouldWrite
    f = getframe(figH);
    writeVideo(vidObj,f);
end
figH.Color = [1 1 1];

for i=1:50
    if i==1
        curPoint = [3 3 3];
    end
    dists = euclid(linearMat, curPoint);
    [minDist,ind] = min(dists);
    vec = curPoint-linearMat(ind,:);
    if minDist>distThresh
        curPoint = curPoint-vec.*propToMove+noiseMag.*randn(1,3);
    else
        [c1,c2] = ind2sub(size(flat),ind);
        c1 = max([1 min([N c1+randi(nToMove)-nToMove/2])]);
        c2 = max([1 min([N c2+randi(nToMove)-nToMove/2])]);
        curPoint = linearMat(sub2ind(size(flat),c1,c2),:);
    end
    plotH = plot3([scH.XData, curPoint(1)], [scH.YData, curPoint(2)], ...
        [scH.ZData, curPoint(3)], 'color', 'k', 'linewidth', 1.5);
    scH.XData = curPoint(1);
    scH.YData = curPoint(2);
    scH.ZData = curPoint(3);
    drawnow;
    pause(.1);
    if shouldWrite
        f = getframe(figH);
        writeVideo(vidObj,f);
    end
end

if shouldWrite
    close(vidObj);
end
%% with some noise perturbations on its way back to manifold
shouldWrite = 1;
if shouldWrite
    movName = 'strongDynamics';
    saveName = [saveLoc,movName, '.avi'];
    vidObj = VideoWriter(saveName);
    vidObj.FrameRate = 10;
    open(vidObj);
end
[x,y] = meshgrid(0:.05:2);
% make a random z surface that makes a pretty plot
z = (x.^2 - .4.*y.^3).^2;
% normalize to be between .5 and 2.5, don't want it to touch origin
z = z./max(z(:))*2+.5;

figH = figure('units', 'pix', 'outerposition', [0 50 1000 1000],...
    'Color', [1 1 1]);
surfH = surf(x,y,z);

direction = [0 0 1];
% rotate(surfH, direction, 30);

xlabel('FR 1', 'fontsize', 16, 'fontweight', 'bold');
ylabel('FR 2', 'fontsize', 16, 'fontweight', 'bold');
zlabel('FR 3', 'fontsize', 16, 'fontweight', 'bold');
axis([0 3 0 3 0 3]);

hold on;
scH = scatter3(x(1,1),y(1,1),z(1,1), 50,'filled');
scH.MarkerFaceColor = 'r';
if shouldWrite
    f = getframe(figH);
    writeVideo(vidObj,f);
end
N = size(x,1);
for i=1:50
    if i==1
        c1 = randi(N);
        c2 = randi(N);
    else
        % move by a few, can't be more than N or less than 0
        nToMove = 16;
        c1 = max([1 min([N c1+randi(nToMove)-nToMove/2])]);
        c2 = max([1 min([N c2+randi(nToMove)-nToMove/2])]);
        
        plotH = plot3([scH.XData, x(c1,c2)], [scH.YData, y(c1,c2)], ...
            [scH.ZData, z(c1,c2)], 'color', 'k', 'linewidth', 1.5);
        
    end
    scH.XData = x(c1,c2);
    scH.YData = y(c1,c2);
    scH.ZData = z(c1,c2);
    
    drawnow;
    pause(.1);
    if shouldWrite
        f = getframe(figH);
        writeVideo(vidObj,f);
    end
end
figH.Color = [0 0 1];
plotH = plot3([scH.XData, 3], [scH.YData, 3], ...
    [scH.ZData, 3], 'color', 'k', 'linewidth', 1.5);
scH.XData = 3;
scH.YData = 3;
scH.ZData = 3;

linearMat = [x(:) y(:) z(:)];
flat = cat(3,x,y,z);
euclid = @(X,v) (X(:,1)-v(1)).^2 + (X(:,2)-v(2)).^2 + (X(:,3)-v(3)).^2;

propToMove = .3;
distThresh = .01;
noiseMag = .4;
pause(.1);
if shouldWrite
    f = getframe(figH);
    writeVideo(vidObj,f);
end
figH.Color = [1 1 1];

for i=1:50
    if i==1
        curPoint = [3 3 3];
    end
    dists = euclid(linearMat, curPoint);
    [minDist,ind] = min(dists);
    vec = curPoint-linearMat(ind,:);
    if minDist>distThresh
        curPoint = curPoint-vec.*propToMove+noiseMag.*randn(1,3);
    else
        [c1,c2] = ind2sub(size(flat),ind);
        c1 = max([1 min([N c1+randi(nToMove)-nToMove/2])]);
        c2 = max([1 min([N c2+randi(nToMove)-nToMove/2])]);
        curPoint = linearMat(sub2ind(size(flat),c1,c2),:);
    end
    plotH = plot3([scH.XData, curPoint(1)], [scH.YData, curPoint(2)], ...
        [scH.ZData, curPoint(3)], 'color', 'k', 'linewidth', 1.5);
    scH.XData = curPoint(1);
    scH.YData = curPoint(2);
    scH.ZData = curPoint(3);
    drawnow;
    pause(.1);
    if shouldWrite
        f = getframe(figH);
        writeVideo(vidObj,f);
    end
end

if shouldWrite
    close(vidObj);
end

%% what happens if perturbation silences everything
shouldWrite = 1;
if shouldWrite
    movName = 'silencing';
    saveName = [saveLoc,movName, '.avi'];
    vidObj = VideoWriter(saveName);
    vidObj.FrameRate = 10;
    open(vidObj);
end
[x,y] = meshgrid(0:.05:2);
% make a random z surface that makes a pretty plot
z = (x.^2 - .4.*y.^3).^2;
% normalize to be between .5 and 2.5, don't want it to touch origin
z = z./max(z(:))*2+.5;

figH = figure('units', 'pix', 'outerposition', [0 50 1000 1000],...
    'Color', [1 1 1]);
surfH = surf(x,y,z);

direction = [0 0 1];
% rotate(surfH, direction, 30);

xlabel('FR 1', 'fontsize', 16, 'fontweight', 'bold');
ylabel('FR 2', 'fontsize', 16, 'fontweight', 'bold');
zlabel('FR 3', 'fontsize', 16, 'fontweight', 'bold');
axis([0 3 0 3 0 3]);

hold on;
scH = scatter3(x(1,1),y(1,1),z(1,1), 50,'filled');
scH.MarkerFaceColor = 'r';
if shouldWrite
    f = getframe(figH);
    writeVideo(vidObj,f);
end
N = size(x,1);
for i=1:50
    if i==1
        c1 = randi(N);
        c2 = randi(N);
    else
        % move by a few, can't be more than N or less than 0
        nToMove = 16;
        c1 = max([1 min([N c1+randi(nToMove)-nToMove/2])]);
        c2 = max([1 min([N c2+randi(nToMove)-nToMove/2])]);
        
        plotH = plot3([scH.XData, x(c1,c2)], [scH.YData, y(c1,c2)], ...
            [scH.ZData, z(c1,c2)], 'color', 'k', 'linewidth', 1.5);
        
    end
    scH.XData = x(c1,c2);
    scH.YData = y(c1,c2);
    scH.ZData = z(c1,c2);
    
    drawnow;
    pause(.1);
    if shouldWrite
        f = getframe(figH);
        writeVideo(vidObj,f);
    end
end
figH.Color = [0 0 1];
plotH = plot3([scH.XData, 0], [scH.YData, 0], ...
    [scH.ZData, 0], 'color', 'k', 'linewidth', 1.5);
scH.XData = 0;
scH.YData = 0;
scH.ZData = 0;

linearMat = [x(:) y(:) z(:)];
flat = cat(3,x,y,z);

propToMove = .3;
distThresh = .01;
noiseMag = .0001;
pause(.1);
if shouldWrite
    f = getframe(figH);
    writeVideo(vidObj,f);
end
figH.Color = [1 1 1];

for i=1:50
    if i==1
        curPoint = [0 0 0];
    end
    dists = euclid(linearMat, curPoint);
    [minDist,ind] = min(dists);
    vec = curPoint-linearMat(ind,:);
    if minDist>distThresh
        curPoint = curPoint-vec.*propToMove+noiseMag.*randn(1,3);
    else
        [c1,c2] = ind2sub(size(flat),ind);
        c1 = max([1 min([N c1+randi(nToMove)-nToMove/2])]);
        c2 = max([1 min([N c2+randi(nToMove)-nToMove/2])]);
        curPoint = linearMat(sub2ind(size(flat),c1,c2),:);
    end
    plotH = plot3([scH.XData, curPoint(1)], [scH.YData, curPoint(2)], ...
        [scH.ZData, curPoint(3)], 'color', 'k', 'linewidth', 1.5);
    scH.XData = curPoint(1);
    scH.YData = curPoint(2);
    scH.ZData = curPoint(3);
    drawnow;
    pause(.1);
    if shouldWrite
        f = getframe(figH);
        writeVideo(vidObj,f);
    end
end

if shouldWrite
    close(vidObj);
end