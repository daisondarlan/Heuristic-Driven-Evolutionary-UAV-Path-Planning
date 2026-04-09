function childRnvec = aStarCrossover(rnvecA, rnvecB, model)
% aStarCrossover_50  Performs crossover between two parent solutions (50x3) using A* as a "bridge."
%
%   childRnvec = aStarCrossover_50(rnvecA, rnvecB, model)
%
% INPUTS:
%   rnvecA, rnvecB : Each 50x3, the control points (parents).
%   model          : struct with fields:
%                    .X : 1x200 array of x-coordinates for the grid
%                    .Y : 200x1 array of y-coordinates for the grid
%                    .H : 200x200 matrix of terrain heights
%
% OUTPUT:
%   childRnvec : 50x3 matrix (the offspring). If no feasible crossover is found, returns [].

    % We'll return [] if no valid bridging path is found or splicing fails.
    childRnvec = [];

    % Indices for "interior" control points (excluding the first & last).
    interiorIndices = 2:6;

    % Randomly pick one interior index from Parent A and one from Parent B.
    i = interiorIndices(randi(numel(interiorIndices)));
    j = interiorIndices(randi(numel(interiorIndices)));

    % Extract the chosen waypoints from each parent.
    A_i = rnvecA(i,:);  % [xA_i, yA_i, zA_i]
    B_j = rnvecB(j,:);  % [xB_j, yB_j, zB_j]

    % 1) Attempt A* from A_i to B_j in the grid.
    pathCells = runAstar2D_bridge(model, A_i, B_j);

    if isempty(pathCells)
        % 2) If not feasible, try the reverse: B_j -> A_i
        pathCells = runAstar2D_bridge(model, B_j, A_i);

        if isempty(pathCells)
            % No feasible path in either direction => no crossover
            return;  % childRnvec stays []
        else
            % Found a path from B_j->A_i, so flip it to get A_i->B_j order
            newWaypoints = flipud( ...
                cellsToCoordinates_bridge(pathCells, model, B_j(3), A_i(3)) ...
            );
        end
    else
        % Found a path from A_i->B_j
        newWaypoints = cellsToCoordinates_bridge(pathCells, model, A_i(3), B_j(3));
    end

    % ---------------------------------------------------------------
    % Build the spliced path:
    %  Segment A1: from parent A(1) up to A(i) (including A_i)
    %  The new bridging sub-path (minus first & last duplicates)
    %  Segment B2: from parent B(j) to B(50) (including B_j)
    %
    % Then resample the combined path to 50 points.
    % ---------------------------------------------------------------

    segmentA1 = rnvecA(1:i, :);     % includes the point A_i
    segmentB2 = rnvecB(j:7, :);   % includes the point B_j

    % Exclude duplicates from newWaypoints
    % newWaypoints(1) = A_i, newWaypoints(end) = B_j
    if size(newWaypoints,1) > 2
        bridgeMid = newWaypoints(2:end-1, :);
    else
        bridgeMid = [];
    end

    combinedPath = [
        segmentA1;          % A(1..i) inclusive
        bridgeMid;          % bridging midpoints
        segmentB2(2:end,:)  % skip the first row of B2 (B_j) to avoid duplication
    ];

    % Now 'combinedPath' can have anywhere from 2 points up to many points.
    % We must downsample/upsample to exactly 50 points.
    nCombined = size(combinedPath,1);
    if nCombined < 2
        % Something went wrong; no valid path
        return;
    end

    % We'll uniformly sample 50 points in [1..nCombined].
    sampleIdx = linspace(1, nCombined, 7);
    childRnvec = zeros(7, 3);

    for k = 1:7
        idxFloat = sampleIdx(k);
        idxLow   = floor(idxFloat);
        idxHigh  = ceil(idxFloat);

        if idxLow == idxHigh
            % Exact integer index
            childRnvec(k,:) = combinedPath(idxLow,:);
        else
            % Linear interpolation between idxLow and idxHigh
            tFrac = idxFloat - idxLow;
            pLow  = combinedPath(idxLow, :);
            pHigh = combinedPath(idxHigh, :);
            childRnvec(k,:) = pLow + tFrac*(pHigh - pLow);
        end
    end

end % end of aStarCrossover_50


%% ------------------------------------------------------------------------
function pathCells = runAstar2D_bridge(model, ptA, ptB)
% runAstar2D_bridge   Wrapper to do A* from continuous [x,y,z] to [x,y,z].
%   pathCells = runAstar2D_bridge(model, ptA, ptB)
%
%   Converts (x,y) of ptA,ptB to nearest grid cells, calls runAstar2D,
%   returns discrete pathCells. If no path found, returns [].

    [~, startXIdx] = min(abs(model.X - ptA(1)));
    [~, startYIdx] = min(abs(model.Y - ptA(2)));
    [~, goalXIdx ] = min(abs(model.X - ptB(1)));
    [~, goalYIdx ] = min(abs(model.Y - ptB(2)));

    startCell = [startYIdx, startXIdx];  % (row, col)
    goalCell  = [goalYIdx, goalXIdx];

    % Use the discrete A* function
    pathCells = runAstar2D(model, startCell, goalCell);
end


%% ------------------------------------------------------------------------
function coords = cellsToCoordinates_bridge(pathCells, model, startZ, goalZ)
% cellsToCoordinates_bridge  Convert path of grid cells [row,col] to Nx3 [x,y,z].
%   We do linear interpolation of altitude from startZ to goalZ.
%
%   coords = cellsToCoordinates_bridge(pathCells, model, startZ, goalZ)

    coords = cellsToCoordinates(pathCells, model, startZ, goalZ);
end


%% ------------------------------------------------------------------------
function pathCells = runAstar2D(model, startCell, goalCell)
% runAstar2D  A simplified A* for a 2D grid using 4-direction moves.
%
%   pathCells = runAstar2D(model, startCell, goalCell)
%
%   model.X: 1x200
%   model.Y: 200x1
%   model.H: 200x200
%
%   startCell, goalCell: [row, col] in [1..200, 1..200]
%
%   Returns the path as Nx2 [row, col], or [] if none found.

    nRows = length(model.Y);
    nCols = length(model.X);

    if isequal(startCell, goalCell)
        pathCells = startCell;
        return;
    end

    % G, H, F cost arrays
    G = inf(nRows, nCols);
    H = inf(nRows, nCols);
    F = inf(nRows, nCols);

    startR = startCell(1); startC = startCell(2);
    goalR  = goalCell(1);  goalC  = goalCell(2);

    G(startR, startC) = 0;
    H(startR, startC) = euclidDist(startR, startC, goalR, goalC);
    F(startR, startC) = G(startR, startC) + H(startR, startC);

    % Parent for path reconstruction
    parent = zeros(nRows, nCols, 2, 'int32');

    % "open set" as a list of [row, col, F]
    openSet = [startR, startC, F(startR, startC)];
    CLOSED = false(nRows, nCols);

    % 4-direction moves
    moves = [ -1 0; 1 0; 0 1; 0 -1 ];

    pathFound = false;

    while ~isempty(openSet)
        % 1) Extract cell with lowest F
        [~, idxMin] = min(openSet(:,3));
        current = openSet(idxMin, 1:2);
        openSet(idxMin,:) = [];  % remove it
        cr = current(1); cc = current(2);

        CLOSED(cr, cc) = true;

        % Check if we reached goal
        if (cr == goalR && cc == goalC)
            pathFound = true;
            break;
        end

        % 2) Explore neighbors
        for m = 1:size(moves,1)
            nr = cr + moves(m,1);
            nc = cc + moves(m,2);

            % bounds check
            if nr<1 || nr>nRows || nc<1 || nc>nCols
                continue;
            end

            % if visited or not traversable, skip
            if CLOSED(nr, nc)
                continue;
            end

            if ~isAboveTerrain(nr, nc, model)
                continue;
            end

            % cost from current to neighbor
            stepCost = 1;  % in 4 directions
            newG = G(cr, cc) + stepCost;

            if newG < G(nr, nc)
                G(nr, nc) = newG;
                H(nr, nc) = euclidDist(nr, nc, goalR, goalC);
                F(nr, nc) = G(nr, nc) + H(nr, nc);

                parent(nr,nc,1) = cr;
                parent(nr,nc,2) = cc;

                openSet = [openSet; nr, nc, F(nr,nc)]; %#ok<AGROW>
            end
        end
    end

    if ~pathFound
        pathCells = [];
        return;
    else
        % Reconstruct path
        pathCells = [];
        r = goalR; c = goalC;

        while ~(r==0 && c==0)
            pathCells = [[r,c]; pathCells]; %#ok<AGROW>
            pr = parent(r,c,1);
            pc = parent(r,c,2);
            r = pr; c = pc;

            if (r == startR && c == startC)
                pathCells = [[r,c]; pathCells]; %#ok<AGROW>
                break;
            end
        end
    end
end


%% ------------------------------------------------------------------------
function yesno = isAboveTerrain(row, col, model)
% isAboveTerrain  Checks if UAV can occupy the cell [row,col].
% For demonstration, let's say it's traversable if H(row,col) < 999999
% or some large threshold. In a real scenario, you'd ensure altitude > terrain, etc.

    terrainZ = model.H(row, col);

    % Simple rule: pass if terrainZ < 1e5
    if terrainZ < 1e5
        yesno = true;
    else
        yesno = false;
    end
end


%% ------------------------------------------------------------------------
function d = euclidDist(r1, c1, r2, c2)
% euclidDist  Euclidean distance in grid coordinates
    d = sqrt((r1 - r2)^2 + (c1 - c2)^2);
end


%% ------------------------------------------------------------------------
function xyzPath = cellsToCoordinates(pathCells, model, startZ, goalZ)
% cellsToCoordinates  Convert Nx2 [row,col] to Nx3 [x,y,z],
%   linearly interpolating altitude from startZ to goalZ.
%
%   pathCells: Nx2 array of [row, col]
%   model.X   : 1x200
%   model.Y   : 200x1
%   model.H   : 200x200  (not strictly needed here unless you'd do another check)
%
%   We produce Nx3 [x, y, z].

    nPoints = size(pathCells,1);
    xyzPath = zeros(nPoints, 3);

    for k = 1:nPoints
        r = pathCells(k,1);
        c = pathCells(k,2);

        xVal = model.X(c);  % col -> X
        yVal = model.Y(r);  % row -> Y

        % Interpolate altitude from startZ to goalZ
        t = (k-1) / (nPoints-1);
        zVal = (1 - t)*startZ + t*goalZ;

        xyzPath(k,:) = [xVal, yVal, zVal];
    end
end
