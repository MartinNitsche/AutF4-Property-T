% This program solves the SDP created by MakeSDP in order to prove Property (T) for Aut(F_4).
% @author Martin Nitsche
% @version LapSym-MaxEq-01
% The mathematical details are explained in the accompanying article and the sources cited therein.
%
% The program also contains an ad-hoc routine for verifying our precomputed solution with integer
% arithmetic. The correctness of the computer proof relies only on the verification step and on
% the correctness of the problem data A, b, c.
% 
% To run this code you need to install MATLAB, and SeDuMi if you want to solve the SDP yourself.
% First, run MakeSDP (of the same version) in the same directory as this file, and make sure that
% it has written output into the 'data' subdirectory.
% Then open MATLAB, navigate to the directory that contains this script, and run it with 'run SolveSDP.m'.
% This program will need ca. 20 minutes to solve the SDP problem created by the unmodified MakeSDP.
% 
% Aside from basic manipulation of (sparse) matrices, this program uses only MATLAB's 'orth' and
% 'sqrtm' functions. All MATLAB documentation can be found at mathworks.com/help/matlab/
%
% If you want to experiment, you can modify the 'data/b.dat' file by hand before running this script,
% in order to prove Property (T) for the higher Aut(F_n) or to gain additional insight.
% You could also look at the dual solution y.


solverPrecision=1e-10;              % precision passed to SDP solver
target=int64(95);                   % to prove Property (T) the objective must be lower than this

disp('reading problem data A, b, c from disk...');
fid=fopen('data/A.dat','r');
sparseA=fscanf(fid,'%ld',[3,Inf])';
fclose(fid);
fid=fopen('data/b.dat','r');
sparseB=fscanf(fid,'%ld',[3,Inf])';
fclose(fid);
fid=fopen('data/c.dat','r');
sparseC=fscanf(fid,'%ld',[3,Inf])';
fclose(fid);

if input('use precomputed solution? y/n ','s')=='y'
  load('precomputedSolution');
  quickVerify(sparseA, sparseB, sparseC, precomputedSolution, target);
  if input('also verify with integer arithmetic? y/n ','s')=='y'
    integerVerify(sparseA, sparseB, sparseC, precomputedSolution, target);
  end
else
  if input('solving SDP with SeDuMi. continue? y/n ','s')=='y'
    solution=solveSDP(sparseA, sparseB, sparseC, solverPrecision);
    quickVerify(sparseA, sparseB, sparseC, solution, target);
    disp('note: verification with integer arithmetic is only available for the precomputed solution');
  else
    disp('abort. exiting.');
  end
end


% solve the SDP problem numerically with precision solverPrecision
% use the symmetry to reduce computation time
function solution = solveSDP(sparseA, sparseB, sparseC, solverPrecision)

  % The way we deal with the symmetry is a hardcoded special case of the calculation in
  % Marek Kaluba's Julia library (see the article referenced in the accompanying paper)
  % The full symmetry group G is an extension of the permutation group S_4 by the diagonal subgroup (Z_2)^4
  % We compute a minimal projection in R(G) for each central projection in R(G) as a product of projections
  % in the group algebras of the normal subgroup and the permutation group
  % To speed things up, the points of the semidefinite variable are sorted in blocks of G-invariant subsets.

  disp('creating orthonormal bases for minimal projections...');
  numProj = 20; % number of projections in R(G) to use
  load('data/permutationAction.dat'); perm=spconvert(permutationAction); % action of the permutation subgroup on the points
  load('data/diagonalAction.dat');    dgnl=spconvert(diagonalAction);    % action of the normal subgroup on the points
  load('data/blockSizes.dat');        sizes=full(spconvert(blockSizes)); % block sizes

  % coefficients for 16 times a rank-1 (central) projection of the normal subgroup
  % in order: +1/+1/+1/+1, -1/-1/-1/-1, -1/+1/+1/+1, +1/-1/-1/-1, -1/-1/+1/+1
  %        ++++  +++-  ++-+  ++--  +-++  +-+-  +--+  +---  -+++  -++-  -+-+  -+--  --++  --+-  ---+  ----
  cfDiag = [  1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1;
              1    -1    -1     1    -1     1     1    -1    -1     1     1    -1     1    -1    -1     1;
              1     1     1     1     1     1     1     1    -1    -1    -1    -1    -1    -1    -1    -1;
              1    -1    -1     1    -1     1     1    -1     1    -1    -1     1    -1     1     1    -1;
              1     1     1     1    -1    -1    -1    -1    -1    -1    -1    -1     1     1     1     1];

  % coefficients for 48 times a rank-1 projection of a subgroup of S_4, one for each irreducible representation
  % in order: whole group with representations 1/1s/2/3/3s, stabilizer of first element with representations 1/1s/2, and 4x <(12),(34)>
  %        1234  1243  1324  1342  1423  1432  2134  2143  2314  2341  2413  2431  3124  3142  3214  3241  3412  3421  4123  4132  4213  4231  4312  4321
  cfPerm = [  2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2     2    2;
              2    -2    -2     2     2    -2    -2     2     2    -2    -2     2     2    -2    -2     2     2    -2    -2     2     2    -2    -2    2;
              4     4    -2    -2    -2    -2     4     4    -2    -2    -2    -2    -2    -2    -2    -2     4     4    -2    -2    -2    -2     4    4;
              6    -6     3    -3    -3     3     6    -6     3    -3    -3     3     3    -3     3    -3     0     0    -3     3    -3     3     0    0;
              6     6    -3    -3    -3    -3    -6    -6     3     3     3     3     3     3    -3    -3     0     0     3     3    -3    -3     0    0;
              8     8     8     8     8     8     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    0;
              8    -8    -8     8     8    -8     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    0;
             16    16    -8    -8    -8    -8     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    0;
             12    12     0     0     0     0    12    12     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    0;
             12   -12     0     0     0     0    12   -12     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    0;
             12    12     0     0     0     0   -12   -12     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    0;
             12   -12     0     0     0     0   -12    12     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    0];
  whichDiag = [ 1 1 1 1 1 2 2 2 2 2 3 3 3 4 4 4 5  5  5  5];  % which diagonal-group projection to use for each minimal projection
  whichPerm = [ 1 2 3 4 5 1 2 3 4 5 6 7 8 6 7 8 9 10 11 12];  % which permutation-group projection to use for each minimal projection
  multiplicities = [1 1 2 3 3 1 1 2 3 3 4 4 8 4 4 8 6 6 6 6]; % dimensions of the corresponding irreducible representations
  numBlocks=size(sizes,1);                                    % number of blocks
  accSizes=sizes;                                             % accumulated block sizes, used later for indexing
  for i=2:size(accSizes,1)
    accSizes(i)=accSizes(i)+accSizes(i-1);
  end
  accSizes=accSizes-sizes;
  numPoints=sum(sizes);           % total number of points in the semidefinite variable
  projDiag=cfDiag*dgnl/16;        % rank-1 projections in group algebra of normal subgroup
  projPerm=cfPerm*perm/48;        % rank-1 projections in group algebra of permutation subgroup
  transf=cell(numProj,numBlocks); % orthonormal basis for image of i-th projection in j-th block
  dimens=cell(numProj,numBlocks); % dimension of image of i-th projection in j-th block

  % calculate transf for each minimal projection and each block
  for l=1:numProj
    i=whichDiag(l);
    j=whichPerm(l);
    fullDiag=reshape(projDiag(i,:),numPoints,numPoints);
    fullPerm=reshape(projPerm(j,:),numPoints,numPoints);
    for k=1:numBlocks
      partDiag=fullDiag(accSizes(k)+1:accSizes(k)+sizes(k),accSizes(k)+1:accSizes(k)+sizes(k));
      partPerm=fullPerm(accSizes(k)+1:accSizes(k)+sizes(k),accSizes(k)+1:accSizes(k)+sizes(k));
      proj=partDiag*partPerm; % obtain the final projections as products of projections for the diagonal and permutation subgroups
      transf{l,k}=orth(proj);
      dimens{l,k}=round(trace(proj));
      if size(transf{l,k},2) ~= dimens{l,k}
        disp('bad rank due to rounding error');
      end
    end
  end

  % for each projection piece together the orthonormal bases from different blocks
  for i=1:numProj
    accDimens{i,1}=0; % dimens accumulated over the blocks, used for indexing below
    for j=2:numBlocks
      accDimens{i,j}=accDimens{i,j-1}+dimens{i,j-1};
    end
    totalDims(i)=accDimens{i,numBlocks}+dimens{i,numBlocks}; % total dimension of i-th projection across all blocks
  end
  for i=1:numProj
    fullTransf{i}=sparse(numPoints,totalDims(i));
    for j=1:numBlocks
      fullTransf{i}(accSizes(j)+1:accSizes(j)+sizes(j),accDimens{i,j}+1:accDimens{i,j}+dimens{i,j})=sparse(transf{i,j});
    end
  end

  % convert objective and constraints
  % instead of a for-loop we use reshape and matrix multiplication to make use of MATLAB's strengths
  disp('converting objective and constraints to new basis...');
  AandC=[spconvert(double(sparseC))';spconvert(double(sparseA))];
  numRows=size(AandC,1);
  A0=reshape(AandC',numPoints,numRows*numPoints);
  for i=1:numProj
    A1=fullTransf{i}'*A0;
    A2=reshape(A1',numPoints,totalDims(i)*numRows)'*fullTransf{i};
    A3=reshape(A2',totalDims(i)*numRows,totalDims(i))';
    ATrans{i}=reshape(A3,totalDims(i)*totalDims(i),numRows)';
  end
  
  % bring data into SeDuMi compatible form and solve
  disp('solving with SeDuMi...');
  b=spconvert(double(sparseB));
  fullTrans=sparse(numRows,totalDims*totalDims');
  pos=1;
  for i=1:numProj
    fullTrans(:,pos:pos+totalDims(i)*totalDims(i)-1)=ATrans{i};
    pos = pos+totalDims(i)*totalDims(i);
  end
  cTrans=fullTrans(1,:)'; % the transformed objective c
  fullTrans(1,:)=[];      % the transformed constraint matrix A
  K.l=0;
  K.s=totalDims;
  pars.eps=solverPrecision;
  [x,y,info]=sedumi(fullTrans,-b,cTrans,K,pars); % the -b results from our and SeDuMi's sing conventions
  
  % convert back to solution of original problem by averaging over the group of symmetries
  disp('converting back to solution of original problem...');
  solPreAvg=zeros(numPoints);  % solution before averaging
  solAvgDgnl=zeros(numPoints); % solution averaged over normal subgroup
  solution=zeros(numPoints);   % solution to original problem
  % for each minimal projection transform back corresponding part of the solution
  pos=1;                       % keeping track of where the SDP blocks for the projections sit in the solution vector x
  for i=1:numProj
    partX=reshape(x(pos:pos+totalDims(i)*totalDims(i)-1),totalDims(i),totalDims(i));
    solPreAvg=solPreAvg+fullTransf{i}*partX*fullTransf{i}';
    pos=pos+totalDims(i)*totalDims(i);
  end
  % to make it invariant, first average it over normal subgroup...
  for j=1:size(dgnl,1)
    curDiag=reshape(dgnl(j,:),numPoints,numPoints);
    solAvgDgnl=solAvgDgnl+curDiag*solPreAvg*curDiag'/size(dgnl,1);
  end
  % and then average it over permutation subgroup
  for j=1:size(perm,1)
    curPerm=reshape(perm(j,:),numPoints,numPoints);
    solution=solution+curPerm*solAvgDgnl*curPerm'/size(perm,1);
  end
end


% Verify the solution, taking care of the fact that due to limited solver precision the
% solution is slightly infeasible. Rounding errors in floating-point arithmetic are estimated
% based on the IEEE 754 standard for the double-precision data type.
function quickVerify(sparseA, sparseB, sparseC, solution, target)
  disp(sprintf('Checking solution. To prove Property (T) we need a feasible point with objective <%d.', target));
  % make sure the primal solution lies in the cone of positive semidefinite matrices
  numPoints=size(solution,1);
  numConstraints=max(sparseA(:,1));
  solutionRoot = real(sqrtm(solution));
  solution = solutionRoot.'*solutionRoot;
  x= reshape(solution,numPoints*numPoints,1);
  % compute gap and constraint violations
  A=spconvert(double(sparseA));
  b=spconvert(double(sparseB));
  c=spconvert(double(sparseC));
  constraintViolation = A*x+b.';
  objective = c.'*x;
  gap = double(target)-objective;
  disp(sprintf('solution objective: %.4f => size of gap: %.4f', objective, gap));
  % account for non-zero constraint violations
  distanceErrorFactor = 36; % radius-squared of a ball around neutral element, containing all distances
  infeasCompensation = sum(max(constraintViolation,0))*distanceErrorFactor;
  disp(sprintf('To compensate for constraint violation, we have to reduce the gap by %e.', infeasCompensation));

  % account for rounding errors during the above computation
  accError=0;    % accumulated error: the sum of all error estimates propagated forward to the end result
  absRoot=abs(solutionRoot);
  errorXEntry=0; % error from matrix multiplication solutionRoot.'*solutionRoot
  errorXEntry=errorXEntry+rErr(max(max(absRoot))*max(max(absRoot)))*numPoints;   % scalar multiplication during matrix multiplication
  absSquare=absRoot.'*absRoot;
  errorXEntry=errorXEntry+rErr(max(max(absSquare)))*numPoints;                   % addition during matrix multiplication
  errorX=repmat(double(errorXEntry),numPoints*numPoints,1);
  accError=accError+abs(c.'*errorX)+distanceErrorFactor*sum(abs(A)*errorX);      % propagate error and add to accError
  if max(abs([sparseA(:,3);sparseC(:,3)]))>2 % if false, no rounding errors during scalar multication in A*x and c.'*x
    error('cannot handle constraint matrices with entries of absolute value >2');
  end
  errorViol=rErr(max(abs(A)*abs(x)+abs(b.')))*(size(sparseA,1)+size(sparseB,1)); % addition part of computation A*x, and summation
  accError=accError+distanceErrorFactor*errorViol;                               % propagate error and add to accError
  accError=accError+size(sparseC,1)*rErr(abs(c.')*abs(x));                       % add error from objective calculation
  accError=accError+rErr(double(target)-objective);                              % finally, error from gap and infeasCompensation computation
  accError=accError+distanceErrorFactor*numConstraints*rErr(sum(abs(constraintViolation)))+rErr(infeasCompensation);
  if accError>=2^52 % if false, it was justified to propagate rounding errors via floating point computation
    error('cannot handle errors of this magnitude');
  end
  approxError=double(accError)/(2^52);
  disp(sprintf('total rounding errors during verification: < %ld*2^{-52} (approx %.2e)', accError, approxError));
end

% give a bound for the rounding precision of the double x>0 in increments of 2^{-52}. used by quickVerify
function res = rErr(x)
  res = int64(ceil(2^floor(log2(abs(x)))));
end


% Multiply solution with a scaling factor, round and verify with integer arithmetic.
% This ad-hoc function is non-optimized and very sensitive to the specific parameters,
% and hence only intended for the precomputed solution.
function integerVerify(sparseA, sparseB, sparseC, solution, target)
  numPoints=size(solution,1);
  numConstraints=max(sparseA(:,1));
  verificationScale=int64(400000000);
  disp('computing matrix root, scaling, and rounding to integer...');
  intSqrt=int64(double(verificationScale)*real(sqrtm(solution)));
  disp('multiplying integer matrix with its adjoint...');
  intSqrtAbs=abs(intSqrt);
  intSol=zeros(numPoints,'int64');        % approximate solution, guaranteed semidefinite positive
  intBounds=zeros(numPoints,'int64');     % intSqrtAbs.'*intSqrtAbs, for checking against int64 overflows
  for i=1:numPoints                       % manual multiplication of dense integer matrices
    intSol(i,:)    = sum(bsxfun(@times, intSqrt(:,i),    intSqrt   ), 1);
    intBounds(i,:) = sum(bsxfun(@times, intSqrtAbs(:,i), intSqrtAbs), 1);
  end
  assert(max(max(intBounds))<intmax('int64')); % MATLAB saturates overflows
  x=reshape(intSol,numPoints*numPoints,1);
  disp('checking integer approximation to solution...');
  intObj=int64(0);                        % objective of the rounded solution
  objOverflow=0;                          % whenever intObj>verificationScale^2, we put the overflow here
  for i=1:size(sparseC,1)                 % compute objective
    intObj=intObj+(x(sparseC(i,1))*sparseC(i,3));
    assert(abs(intObj)<intmax('int64'));
    if intObj>verificationScale*verificationScale
      objOverflow=objOverflow+1;
      intObj=intObj-verificationScale*verificationScale;
    end
  end
  intViol=int64(zeros(numConstraints,1)); % extent of constraint violation
  assert(max(abs(x))*max(abs(sparseA(:,3)))<intmax('int64'));
  for i=1:size(sparseA,1)                 % compute A*x
    intViol(sparseA(i,1))=intViol(sparseA(i,1))+x(sparseA(i,2))*sparseA(i,3);
    assert(abs(intViol(sparseA(i,1)))<intmax('int64'));
  end
  for i=1:size(sparseB,1)                 % compute A*x+b
    intViol(sparseB(i,2))=intViol(sparseB(i,2))+sparseB(i,3)*verificationScale*verificationScale;
    assert(abs(intViol(sparseA(i,1)))<intmax('int64'));
  end
  compensation=int64(0);                  % penalty to objective to compensate for constraint violation
  distanceErrorFactor = int64(36);        % radius-squared of a ball around neutral element, containing all distances
  for i=1:size(intViol)
    compensation=compensation+max(0,intViol(i))*distanceErrorFactor;
  end
  if intObj+compensation<(target-objOverflow)*verificationScale*verificationScale
    disp('success! non-vanishing gap verified with integer arithmetic');
  else
    disp('the existence of the gap could not be verified');
  end
end
