%-----------------------------------------------------------------------------------------
function pathCells = runAstar2D(model, startCell, goalCell, row_min, row_max, col_min, col_max)
% runAstar2D  A* on a reduced 2D grid with a true min‐heap.
%
%   pathCells = runAstar2D(model, startCell, goalCell, row_min, row_max, col_min, col_max)
%
%   Searches only within rows [row_min:row_max] and cols [col_min:col_max].
%   Uses a binary min‐heap for the open set (O(log N) push/pop).
%
% INPUTS:
%   model     : struct with fields X (1×200), Y (200×1), H (200×200).
%   startCell : [startRow, startCol] in global coordinates.
%   goalCell  : [goalRow,  goalCol ] in global coordinates.
%   row_min, row_max, col_min, col_max : bounding indices for the subgrid.
%
% OUTPUT:
%   pathCells : K×2 list of [row, col] in global indices for the found path.
%               [] if no path is found.

    % Extract subgrid dimensions
    ROW0 = row_min;
    COL0 = col_min;
    nRowsRegion = row_max - row_min + 1;
    nColsRegion = col_max - col_min + 1;

    % Convert global start/goal to local indices
    sR = startCell(1) - ROW0 + 1;  sC = startCell(2) - COL0 + 1;
    gR = goalCell(1)  - ROW0 + 1;  gC = goalCell(2)  - COL0 + 1;

    % Bounds check: if start or goal outside the region, fail
    if sR < 1 || sR > nRowsRegion || sC < 1 || sC > nColsRegion || ...
       gR < 1 || gR > nRowsRegion || gC < 1 || gC > nColsRegion
        pathCells = [];
        return;
    end

    % Pre‐slice the traversability mask once
    subH = model.H(row_min:row_max, col_min:col_max);
    traversableMask = subH < 1e5;  % UAV can traverse if terrain < some threshold

    % Initialize cost arrays
    infVal = Inf;
    G = infVal * ones(nRowsRegion, nColsRegion);
    H = infVal * ones(nRowsRegion, nColsRegion);
    F = infVal * ones(nRowsRegion, nColsRegion);
    CLOSED = false(nRowsRegion, nColsRegion);
    parent = zeros(nRowsRegion, nColsRegion, 2, 'int32');

    % Heuristic: Euclidean distance in cell space
    for r = 1:nRowsRegion
        for c = 1:nColsRegion
            H(r, c) = euclidDist(r, c, gR, gC);
        end
    end

    % Initialize start node
    G(sR, sC) = 0;
    F(sR, sC) = H(sR, sC);

    % Prepare a min‐heap: fixed‐size arrays of length ≤ nRowsRegion*nColsRegion
    maxNodes = nRowsRegion * nColsRegion;
    heapSize = 0;
    heapF = zeros(maxNodes, 1);
    heapR = zeros(maxNodes, 1);
    heapC = zeros(maxNodes, 1);

    % Push function for the heap
    function push(rLoc, cLoc, fVal)
        heapSize = heapSize + 1;
        idx = heapSize;
        heapR(idx) = rLoc;
        heapC(idx) = cLoc;
        heapF(idx) = fVal;
        % Sift‐up
        while idx > 1
            parentIdx = floor(idx / 2);
            if heapF(parentIdx) <= heapF(idx)
                break;
            end
            % Swap nodes at idx and parentIdx
            [heapF(parentIdx), heapF(idx)] = deal(heapF(idx), heapF(parentIdx));
            [heapR(parentIdx), heapR(idx)] = deal(heapR(idx), heapR(parentIdx));
            [heapC(parentIdx), heapC(idx)] = deal(heapC(idx), heapC(parentIdx));
            idx = parentIdx;
        end
    end

    % Pop function for the heap
    function [rOut, cOut, fOut] = pop()
        if heapSize == 0
            rOut = []; cOut = []; fOut = [];
            return;
        end
        rOut = heapR(1);
        cOut = heapC(1);
        fOut = heapF(1);
        % Move last node to root
        heapF(1) = heapF(heapSize);
        heapR(1) = heapR(heapSize);
        heapC(1) = heapC(heapSize);
        heapSize = heapSize - 1;
        idx = 1;
        % Sift‐down
        while true
            left  = 2 * idx;
            right = left + 1;
            smallest = idx;
            if left <= heapSize && heapF(left) < heapF(smallest)
                smallest = left;
            end
            if right <= heapSize && heapF(right) < heapF(smallest)
                smallest = right;
            end
            if smallest == idx
                break;
            end
            [heapF(idx), heapF(smallest)] = deal(heapF(smallest), heapF(idx));
            [heapR(idx), heapR(smallest)] = deal(heapR(smallest), heapR(idx));
            [heapC(idx), heapC(smallest)] = deal(heapC(smallest), heapC(idx));
            idx = smallest;
        end
    end

    % Insert the start node
    push(sR, sC, F(sR, sC));

    % Movement offsets for 4‐connected grid (no diagonals)
    moves = [ -1, 0;  1, 0;  0, 1;  0, -1 ];

    pathFound = false;
    while heapSize > 0
        [rCurr, cCurr, fCurr] = pop();
        % If already closed (we saw a better path earlier), skip
        if CLOSED(rCurr, cCurr)
            continue;
        end
        % Mark as closed
        CLOSED(rCurr, cCurr) = true;
        % Check if reached goal
        if rCurr == gR && cCurr == gC
            pathFound = true;
            break;
        end
        % Explore neighbors
        for mi = 1:4
            nr = rCurr + moves(mi, 1);
            nc = cCurr + moves(mi, 2);
            % Bounds check
            if nr < 1 || nr > nRowsRegion || nc < 1 || nc > nColsRegion
                continue;
            end
            % Traversability check
            if ~traversableMask(nr, nc)
                continue;
            end
            if CLOSED(nr, nc)
                continue;
            end
            tentativeG = G(rCurr, cCurr) + 1;  % cost = 1 per move
            if tentativeG < G(nr, nc)
                G(nr, nc) = tentativeG;
                F(nr, nc) = tentativeG + H(nr, nc);
                parent(nr, nc, 1) = rCurr;
                parent(nr, nc, 2) = cCurr;
                push(nr, nc, F(nr, nc));
            end
        end
    end

    if ~pathFound
        pathCells = [];
        return;
    end

    % Reconstruct path in local indices
    localPath = [gR, gC];
    r = gR;  c = gC;
    while ~(r == sR && c == sC)
        pr = parent(r, c, 1);
        pc = parent(r, c, 2);
        r = pr;  c = pc;
        localPath = [[r, c]; localPath]; %#ok<AGROW>
    end

    % Convert local indices back to global (row, col)
    nPts = size(localPath, 1);
    pathCells = zeros(nPts, 2);
    for k = 1:nPts
        pathCells(k, 1) = localPath(k, 1) + ROW0 - 1;
        pathCells(k, 2) = localPath(k, 2) + COL0 - 1;
    end
end


%-----------------------------------------------------------------------------------------
function d = euclidDist(r1, c1, r2, c2)
    % Simple Euclidean distance in cell coordinates
    d = sqrt((r1 - r2)^2 + (c1 - c2)^2);
end