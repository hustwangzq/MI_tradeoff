function gammaGrid = SIM_build_boundary_aware_gamma_grid( ...
    gammaStart,gammaFormalMax,coarseStepDB,opts)
%SIM_BUILD_BOUNDARY_AWARE_GAMMA_GRID Dense grid near communication boundary.
%
% Far from the pure-communication ceiling, keep the original coarse dB step.
% Closer to the ceiling, progressively reduce the step so that a numerical
% failure does not skip the nearly-horizontal CRB tail.

if nargin < 4 || isempty(opts)
    opts = struct('nearBoundaryWidthDB',5);
end
if ~(isfinite(gammaStart) && gammaStart > 0 && ...
        isfinite(gammaFormalMax) && gammaFormalMax >= gammaStart)
    error('SIM_build_boundary_aware_gamma_grid:BadRange', ...
        'gammaStart/gammaFormalMax must define a finite positive interval.');
end

startDB = 10*log10(gammaStart);
endDB = 10*log10(gammaFormalMax);
if endDB <= startDB + 1e-12
    gammaGrid = gammaStart;
    return;
end

coarseStepDB = max(coarseStepDB,0.05);
nearWidth = get_option(opts,'nearBoundaryWidthDB',5.0);
midWidth = get_option(opts,'midBoundaryWidthDB',2.0);
veryNearWidth = get_option(opts,'veryNearBoundaryWidthDB',0.5);
nearStep = get_option(opts,'nearBoundaryStepDB',0.25);
midStep = get_option(opts,'midBoundaryStepDB',0.10);
veryNearStep = get_option(opts,'veryNearBoundaryStepDB',0.025);

pointsDB = startDB;
currentDB = startDB;
while currentDB < endDB-1e-12
    gapToEnd = endDB-currentDB;
    if gapToEnd <= veryNearWidth
        stepDB = min(coarseStepDB,veryNearStep);
    elseif gapToEnd <= midWidth
        stepDB = min(coarseStepDB,midStep);
    elseif gapToEnd <= nearWidth
        stepDB = min(coarseStepDB,nearStep);
    else
        stepDB = coarseStepDB;
    end
    nextDB = min(currentDB+stepDB,endDB);
    if nextDB <= currentDB + 1e-10
        break;
    end
    pointsDB(end+1,1) = nextDB; %#ok<AGROW>
    currentDB = nextDB;
end

if abs(pointsDB(end)-endDB) > 1e-9
    pointsDB(end+1,1) = endDB;
end
gammaGrid = 10.^(pointsDB/10);
gammaGrid = gammaGrid(:);
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end
