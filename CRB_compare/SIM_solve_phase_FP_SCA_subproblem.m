function result = SIM_solve_phase_FP_SCA_subproblem( ...
    gradF,cCurrent,gradC,rhoF,rhoC,trustRadius,constraintTol,alg)
%SIM_SOLVE_PHASE_FP_SCA_SUBPROBLEM Solve one small convex phase-step model.
%
% The real phase increment dTheta is obtained from
%
%   min  gradF(:)'*d + 0.5*rhoF*||d||_2^2
%
% subject to, for every user k,
%
%   c_k + gradC_k(:)'*d + 0.5*rhoC(k)*||d||_2^2 <= constraintTol,
%   |d_i| <= trustRadius.
%
% The current point d=0 is feasible whenever cCurrent<=constraintTol. The
% calling routine still evaluates the exact unregularized CRB and the exact
% SINR constraints before accepting the returned direction.

validateattributes(rhoF,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'rhoF');
validateattributes(trustRadius,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'trustRadius');
validateattributes(constraintTol,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'constraintTol');

phaseSize = size(gradF);
numVariables = numel(gradF);
numUsers = numel(cCurrent);

if ndims(gradC) ~= 3 || size(gradC,1) ~= phaseSize(1) || ...
        size(gradC,2) ~= phaseSize(2) || size(gradC,3) ~= numUsers
    error('SIM_solve_phase_FP_SCA_subproblem:GradientSizeMismatch', ...
        'gradC must have size [size(gradF), numel(cCurrent)].');
end

rhoC = rhoC(:);
if isscalar(rhoC)
    rhoC = repmat(rhoC,numUsers,1);
end
if numel(rhoC) ~= numUsers || any(~isfinite(rhoC)) || any(rhoC <= 0)
    error('SIM_solve_phase_FP_SCA_subproblem:InvalidConstraintCurvature', ...
        'rhoC must contain one positive finite value per user.');
end

objectiveGradient = real(gradF(:));
constraintGradient = zeros(numVariables,numUsers);
for k = 1:numUsers
    gradK = gradC(:,:,k);
    constraintGradient(:,k) = real(gradK(:));
end

result.success = false;
result.status = 'not-run';
result.cvxStatus = 'not-run';
result.cvxOptimalValue = Inf;
result.dTheta = zeros(phaseSize);
result.stepInfNorm = 0;
result.modelChange = 0;
result.predictedDecrease = 0;
result.maxModelConstraint = max(cCurrent(:));

try
    cvx_clear;
catch
end

try
    cvx_begin quiet
        variable d(numVariables)
        minimize( objectiveGradient.'*d + 0.5*rhoF*sum_square(d) )
        subject to
            d <= trustRadius;
            d >= -trustRadius;
            for k = 1:numUsers
                cCurrent(k) + constraintGradient(:,k).'*d + ...
                    0.5*rhoC(k)*sum_square(d) <= constraintTol;
            end
    cvx_end
catch ME
    result.status = ['CVX exception: ',ME.message];
    result.cvxStatus = 'exception';
    try
        cvx_clear;
    catch
    end
    return;
end

result.cvxStatus = char(string(cvx_status));
result.cvxOptimalValue = cvx_optval;
solved = contains(lower(result.cvxStatus),'solved');
if ~solved || any(~isfinite(d))
    result.status = ['phase subproblem failed: ',result.cvxStatus];
    return;
end

dTheta = reshape(real(d),phaseSize);
stepNormSquaredValue = sum(real(d).^2);
modelChange = objectiveGradient.'*real(d) + ...
    0.5*rhoF*stepNormSquaredValue;
modelConstraint = zeros(numUsers,1);
for k = 1:numUsers
    modelConstraint(k) = cCurrent(k) + ...
        constraintGradient(:,k).'*real(d) + ...
        0.5*rhoC(k)*stepNormSquaredValue;
end

result.success = true;
result.status = 'solved';
result.dTheta = dTheta;
result.stepInfNorm = norm(dTheta(:),Inf);
result.modelChange = modelChange;
result.predictedDecrease = max(0,-modelChange);
result.maxModelConstraint = max(modelConstraint);
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end
