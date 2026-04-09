function mutatedRnvec = aStarMutation(rnvec, model)
% aStarMutation  Mutates a 7×3 B-spline control‐point chromosome using A* on a reduced 2D grid.
%
%    mutatedRnvec = aStarMutation(rnvec, model)
%
% INPUTS:
%    rnvec : 7×3 matrix of control points [x, y, z].
%    model : struct with fields:
%            .X : 1×200 array of x‐coordinates for the grid
%            .Y : 200×1 array of y‐coordinates for the grid
%            .H : 200×200 matrix of terrain heights at each (Yidx, Xidx)
%
% OUTPUT:
%    mutatedRnvec : 7×3 matrix of control points (possibly mutated).

    mutatedRnvec = rnvec;  % Start with a copy

    % 1) Pick two random interior control points (indices 2..6)
    interiorIndices = 2:6;
    if numel(interiorIndices) < 2
        return;  % Nothing to do
    end
    pick = randperm(numel(interiorIndices), 2);
    i = interiorIndices(min(pick));
    j = interiorIndices(max(pick));
    if i == j
        return;  % No true sub‐path
    end

    startCP = rnvec(i, :);  % [x_i, y_i, z_i]
    goalCP  = rnvec(j, :);  % [x_j, y_j, z_j]

    % 2) Convert these (x, y) to global grid indices
    [~, startXIdx] = min(abs(model.X - startCP(1)));
    [~, startYIdx] = min(abs(model.Y - startCP(2)));
    [~, goalXIdx ] = min(abs(model.X - goalCP(1)));
    [~, goalYIdx ] = min(abs(model.Y - goalCP(2)));

    startCell = [startYIdx, startXIdx];  % format (row, col)
    goalCell  = [goalYIdx, goalXIdx];

    % 3) Define a small bounding box around start and goal
    R = 20;  % half‐width of the search region (you can tune this)
    nRows = numel(model.Y);
    nCols = numel(model.X);

    row_min = max(1, min(startYIdx, goalYIdx) - R);
    row_max = min(nRows, max(startYIdx, goalYIdx) + R);
    col_min = max(1, min(startXIdx, goalXIdx) - R);
    col_max = min(nCols, max(startXIdx, goalXIdx) + R);

    % 4) Call the new runAstar2D on the subgrid [row_min:row_max, col_min:col_max]
    pathCells = runAstar2D(model, startCell, goalCell, row_min, row_max, col_min, col_max);

    % If no path found, skip mutation
    if isempty(pathCells)
        return;
    end

    % 5) Convert the returned cell list (global indices) into continuous [x, y, z]
    %    Linear interpolation of z from startCP(3) to goalCP(3)
    nPoints = size(pathCells, 1);
    xyzPath = zeros(nPoints, 3);
    for k = 1:nPoints
        rGlobal = pathCells(k, 1);
        cGlobal = pathCells(k, 2);
        xVal = model.X(cGlobal);
        yVal = model.Y(rGlobal);
        t    = (k - 1) / (nPoints - 1);
        zVal = (1 - t)*startCP(3) + t*goalCP(3);
        xyzPath(k, :) = [xVal, yVal, zVal];
    end

    % 6) Overwrite interior control points between i and j
    nMid = j - i - 1;
    if nMid < 1 || nPoints <= 2
        return;  % Nothing to replace
    end

    % Sample exactly nMid waypoints from xyzPath(2:end‐1, :)
    sampleIdx = round(linspace(2, nPoints - 1, nMid));
    for m = 1:nMid
        mutatedRnvec(i + m, :) = xyzPath(sampleIdx(m), :);
    end
end