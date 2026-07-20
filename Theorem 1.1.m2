-- The purpose of this code is to prove the border rank of the determinant tensor of order 4 is at least 12 which matches the current known upper bound.
-- For details, see the preprint 'The border rank of the $4 \times 4$ determinant tensor is twelve' by Jong In Han, Jeong-Hoon Ju, Yeongrak Kim.
-- This code also refers to the paper 'New lower bounds for matrix multiplication and $\det_3$' by A. Conner, A. Harper, J. M. Landsberg and their implementation.


kk = QQ;
r = 11; -- We want to show that no candidates pass the tests for brk=11. Note that this code is made only for the case r=11. A different value of r requires a minor modification.

listError = {};


-- In v_(i,j), 'i' represents the index of the factor
-- and 'j' represents the vector in the standard basis.
-- For example, v_(2,3) represents the element v_3 of the standard basis {v_1,v_2,v_3,v_4}
-- in the vector space V_2.
T = kk[v_(1, 1)..v_(4, 4)];


-- doPermute returns the result of the right action to a tensor by permu.
-- For example, {1,4,2,3} corresponds to the permutation
--  | 1 2 3 4 |
--  | 1 4 2 3 |
-- If we execute doPermute(v_(1,2)*v_(2,3)*v_(3,4)*v_(4,1),{4,3,1,2}),
-- then it returns v_(1,1)*v_(2,4)*v_(3,2)*v_(4,3).
needsPackage "Permutations";
rawPermute = permu -> (
    subList := {};
    inversePermu := toList inverse permutation permu;
    for i from 1 to #permu do (
        for j from 1 to 4 do (
            subList = append(subList, v_(i, j) => v_(inversePermu_(i-1), j));
        );
    );
    return subList;
);
doPermute = (x, permu) -> sub(x, rawPermute permu);


-- As $V_1\otimes V_2\otimes V_3
--    =S_{\lambda_1}V\oplus S_{\lambda_2}V\oplus S_{\lambda_3}V\oplus \Lambda^3 V$,
-- we first find all 7-dimensional B-fixed subspaces of
-- $S_{\lambda_1}V\oplus S_{\lambda_2}V\oplus S_{\lambda_3}V$.


-- Young symmetrizer for \lambda_i
symmetrizer = i -> (
    if i == 1 then (
        return j -> j + doPermute(j, {2, 1, 3}) + doPermute(j, {3, 2, 1}) + doPermute(j, {1, 3, 2}) + doPermute(j, {2, 3, 1}) + doPermute(j, {3, 1, 2});
    ) else if i == 2 then (
        return j -> j + doPermute(j, {2, 1, 3}) - doPermute(j, {3, 2, 1}) - doPermute(j, {3, 1, 2});
    ) else if i == 3 then (
        return j -> j + doPermute(j, {3, 2, 1}) - doPermute(j, {2, 1, 3}) - doPermute(j, {2, 3, 1});
    ) else if i == 4 then (
        return j -> j - doPermute(j, {2, 1, 3}) - doPermute(j, {3, 2, 1}) - doPermute(j, {1, 3, 2}) + doPermute(j, {2, 3, 1}) + doPermute(j, {3, 1, 2});
    ) else error "Error in symmetrizer.";
);


-- get u_i(a,b,c)
u = (i, a, b, c) -> (
    return (symmetrizer i)(v_(1, a)*v_(2, b)*v_(3, c));
);


-- (wedge^3)V
Slambda4 = {u(4, 1, 2, 3), u(4, 1, 2, 4), u(4, 1, 3, 4), u(4, 2, 3, 4)};


-- list of weights
weights = {(1, 1, 1), (1, 1, 2), (1, 1, 3), (1, 1, 4), (1, 2, 2), (1, 2, 3), (1, 2, 4), (1, 3, 3), (1, 3, 4), (1, 4, 4), (2, 2, 2), (2, 2, 3), (2, 2, 4), (2, 3, 3), (2, 3, 4), (2, 4, 4), (3, 3, 3), (3, 3, 4), (3, 4, 4), (4, 4, 4)};


-- weight spaces. Those are represented by their bases.
ws_(1, 1, 1) = {u(1, 1, 1, 1)}
ws_(1, 1, 2) = {u(1, 1, 1, 2), u(2, 1, 1, 2), u(3, 1, 1, 2)}
ws_(1, 1, 3) = {u(1, 1, 1, 3), u(2, 1, 1, 3), u(3, 1, 1, 3)}
ws_(1, 1, 4) = {u(1, 1, 1, 4), u(2, 1, 1, 4), u(3, 1, 1, 4)}
ws_(1, 2, 2) = {u(1, 1, 2, 2), u(2, 2, 2, 1), u(3, 2, 2, 1)}
ws_(1, 2, 3) = {u(1, 1, 2, 3), u(2, 1, 2, 3), u(2, 1, 2, 3)+u(2, 1, 3, 2), u(3, 1, 3, 2), u(3, 1, 2, 3)+u(3, 1, 3, 2)}
ws_(1, 2, 4) = {u(1, 1, 2, 4), u(2, 1, 2, 4), u(2, 1, 2, 4)+u(2, 1, 4, 2), u(3, 1, 4, 2), u(3, 1, 2, 4)+u(3, 1, 4, 2)}
ws_(1, 3, 3) = {u(1, 1, 3, 3), u(2, 3, 3, 1), u(3, 3, 3, 1)}
ws_(1, 3, 4) = {u(1, 1, 3, 4), u(2, 1, 3, 4), u(2, 1, 3, 4)+u(2, 1, 4, 3), u(3, 1, 4, 3), u(3, 1, 3, 4)+u(3, 1, 4, 3)}
ws_(1, 4, 4) = {u(1, 1, 4, 4), u(2, 4, 4, 1), u(3, 4, 4, 1)}
ws_(2, 2, 2) = {u(1, 2, 2, 2)}
ws_(2, 2, 3) = {u(1, 2, 2, 3), u(2, 2, 2, 3), u(3, 2, 2, 3)}
ws_(2, 2, 4) = {u(1, 2, 2, 4), u(2, 2, 2, 4), u(3, 2, 2, 4)}
ws_(2, 3, 3) = {u(1, 2, 3, 3), u(2, 3, 3, 2), u(3, 3, 3, 2)}
ws_(2, 3, 4) = {u(1, 2, 3, 4), u(2, 2, 3, 4), u(2, 2, 3, 4)+u(2, 2, 4, 3), u(3, 2, 4, 3), u(3, 2, 3, 4)+u(3, 2, 4, 3)}
ws_(2, 4, 4) = {u(1, 2, 4, 4), u(2, 4, 4, 2), u(3, 4, 4, 2)}
ws_(3, 3, 3) = {u(1, 3, 3, 3)}
ws_(3, 3, 4) = {u(1, 3, 3, 4), u(2, 3, 3, 4), u(3, 3, 3, 4)}
ws_(3, 4, 4) = {u(1, 3, 4, 4), u(2, 4, 4, 3), u(3, 4, 4, 3)}
ws_(4, 4, 4) = {u(1, 4, 4, 4)}

assert(#(flatten for i in weights list ws_i) == 60)


-- List of arrows
-- keys represent (source, target)
-- values represent the map
arrows = hashTable {{(1, 1, 2), (1, 1, 1)} => matrix{{1}}|matrix{{0}}|matrix{{0}},
    {(1, 1, 3), (1, 1, 2)} => matrix{{1}, {0}, {0}}|matrix{{0}, {1}, {0}}|matrix{{0}, {0}, {1}},
    {(1, 2, 2), (1, 1, 2)} => matrix{{2}, {0}, {0}}|matrix{{0}, {-1}, {0}}|matrix{{0}, {0}, {-1}},
    {(1, 1, 4), (1, 1, 3)} => matrix{{1}, {0}, {0}}|matrix{{0}, {1}, {0}}|matrix{{0}, {0}, {1}},
    {(1, 2, 3), (1, 1, 3)} => matrix{{1}, {0}, {0}}|matrix{{0}, {1}, {0}}|matrix{{0}, {1/2}, {0}}|matrix{{0}, {0}, {-2}}|matrix{{0}, {0}, {-1}},
    {(1, 2, 3), (1, 2, 2)} => matrix{{1}, {0}, {0}}|matrix{{0}, {-1/2}, {0}}|matrix{{0}, {-1}, {0}}|matrix{{0}, {0}, {1}}|matrix{{0}, {0}, {2}},
    {(2, 2, 2), (1, 2, 2)} => matrix{{3}, {0}, {0}},
    {(1, 2, 4), (1, 1, 4)} => matrix{{1}, {0}, {0}}|matrix{{0}, {1}, {0}}|matrix{{0}, {1/2}, {0}}|matrix{{0}, {0}, {-2}}|matrix{{0}, {0}, {-1}},
    {(1, 2, 4), (1, 2, 3)} => matrix{{1}, {0}, {0}, {0}, {0}}|matrix{{0}, {1}, {0}, {0}, {0}}|matrix{{0}, {0}, {1}, {0}, {0}}|matrix{{0}, {0}, {0}, {1}, {0}}|matrix{{0}, {0}, {0}, {0}, {1}},
    {(1, 3, 3), (1, 2, 3)} => matrix{{2}, {0}, {0}, {0}, {0}}|matrix{{0}, {0}, {-2}, {0}, {0}}|matrix{{0}, {0}, {0}, {0}, {1}},
    {(2, 2, 3), (1, 2, 3)} => matrix{{2}, {0}, {0}, {0}, {0}}|matrix{{0}, {2}, {0}, {0}, {0}}|matrix{{0}, {0}, {0}, {-1}, {0}},
    {(2, 2, 3), (2, 2, 2)} => matrix{{1}}|matrix{{0}}|matrix{{0}},
    {(1, 3, 4), (1, 2, 4)} => matrix{{1}, {0}, {0}, {0}, {0}}|matrix{{0}, {1}, {0}, {0}, {0}}|matrix{{0}, {0}, {1}, {0}, {0}}|matrix{{0}, {0}, {0}, {1}, {0}}|matrix{{0}, {0}, {0}, {0}, {1}},
    {(1, 3, 4), (1, 3, 3)} => matrix{{1}, {0}, {0}}|matrix{{0}, {-1/2}, {0}}|matrix{{0}, {-1}, {0}}|matrix{{0}, {0}, {1}}|matrix{{0}, {0}, {2}},
    {(2, 2, 4), (1, 2, 4)} => matrix{{2}, {0}, {0}, {0}, {0}}|matrix{{0}, {2}, {0}, {0}, {0}}|matrix{{0}, {0}, {0}, {-1}, {0}},
    {(2, 2, 4), (2, 2, 3)} => matrix{{1}, {0}, {0}}|matrix{{0}, {1}, {0}}|matrix{{0}, {0}, {1}},
    {(2, 3, 3), (1, 3, 3)} => matrix{{1}, {0}, {0}}|matrix{{0}, {1}, {0}}|matrix{{0}, {0}, {1}},
    {(2, 3, 3), (2, 2, 3)} => matrix{{2}, {0}, {0}}|matrix{{0}, {-1}, {0}}|matrix{{0}, {0}, {-1}},
    {(1, 4, 4), (1, 3, 4)} => matrix{{2}, {0}, {0}, {0}, {0}}|matrix{{0}, {0}, {-2}, {0}, {0}}|matrix{{0}, {0}, {0}, {0}, {1}},
    {(2, 3, 4), (1, 3, 4)} => matrix{{1}, {0}, {0}, {0}, {0}}|matrix{{0}, {1}, {0}, {0}, {0}}|matrix{{0}, {0}, {1}, {0}, {0}}|matrix{{0}, {0}, {0}, {1}, {0}}|matrix{{0}, {0}, {0}, {0}, {1}},
    {(2, 3, 4), (2, 2, 4)} => matrix{{1}, {0}, {0}}|matrix{{0}, {1}, {0}}|matrix{{0}, {1/2}, {0}}|matrix{{0}, {0}, {-2}}|matrix{{0}, {0}, {-1}},
    {(2, 3, 4), (2, 3, 3)} => matrix{{1}, {0}, {0}}|matrix{{0}, {-1/2}, {0}}|matrix{{0}, {-1}, {0}}|matrix{{0}, {0}, {1}}|matrix{{0}, {0}, {2}},
    {(3, 3, 3), (2, 3, 3)} => matrix{{3}, {0}, {0}},
    {(2, 4, 4), (1, 4, 4)} => matrix{{1}, {0}, {0}}|matrix{{0}, {1}, {0}}|matrix{{0}, {0}, {1}},
    {(2, 4, 4), (2, 3, 4)} => matrix{{2}, {0}, {0}, {0}, {0}}|matrix{{0}, {0}, {-2}, {0}, {0}}|matrix{{0}, {0}, {0}, {0}, {1}},
    {(3, 3, 4), (2, 3, 4)} => matrix{{2}, {0}, {0}, {0}, {0}}|matrix{{0}, {2}, {0}, {0}, {0}}|matrix{{0}, {0}, {0}, {-1}, {0}},
    {(3, 3, 4), (3, 3, 3)} => matrix{{1}}|matrix{{0}}|matrix{{0}},
    {(3, 4, 4), (2, 4, 4)} => matrix{{1}, {0}, {0}}|matrix{{0}, {1}, {0}}|matrix{{0}, {0}, {1}},
    {(3, 4, 4), (3, 3, 4)} => matrix{{2}, {0}, {0}}|matrix{{0}, {-1}, {0}}|matrix{{0}, {0}, {-1}},
    {(4, 4, 4), (3, 4, 4)} => matrix{{3}, {0}, {0}}}


-- Get the constraints that d_wt's must satisfy where d_wt is the dimension of the chosen subspace in the weight space of wt.
-- wt1posConstraints_wt1pos contains sequences (wt2posList, nullity).
-- wt1pos represents the index of wt1.
-- wt2posList contains a sublist of the indices of the weight spaces that wt1 maps into.
-- Nullity represents the nullity of the map wt1->(the direct sum of wt2's in the sublist)
-- These must satisfy: (sum of d_wt2 for wt2 in wt2posList) + nullity - d_wt1 >= 0.
for i from 0 to #weights-1 do (
    wt1pos = i;
    wt1posConstraints_wt1pos = {};
    listTarget = {};
    for j from 0 to #weights-1 do(
        wt2pos = j;
        wt1 = weights_i;
        wt2 = weights_j;
        if arrows#?{wt1, wt2} then (
            arrowMap = arrows#{wt1, wt2};
            listTarget = append(listTarget, (wt2, arrowMap, wt2pos));
        ) else continue;
    );

    listSubsets = subsets listTarget;
    for asubset in listSubsets do (
        if asubset == {} then continue;
        accMap = null;
        wt2posList = {};
        for pair in asubset do (
            wt2pos = pair_2;
            wt2posList = append(wt2posList, wt2pos);
            arrowMap = pair_1;
            if accMap === null then accMap = arrowMap else accMap = accMap || arrowMap;
        );
        nullity = numColumns accMap - rank accMap;
        wt1posConstraints_wt1pos = append(wt1posConstraints_wt1pos, (wt2posList, nullity));
    );
);


-- Using the constraints, compute the eligible combinations of d_wt using the recursive loop.
dimCombs = {};
recursiveLoop = (pos, remainingDim, L) -> (
    if pos == #weights then (
        if remainingDim == 0 then dimCombs = append(dimCombs, L);
        return;
    );

    wt := weights_pos;
    maxdim := min(remainingDim, #(ws_wt));
    isPassed := null;
    wt2posList := null;
    nullity := null;
    for i from 0 to maxdim do (
        isPassed = true;
        for constraint in wt1posConstraints_pos do (
            wt2posList = constraint_0;
            nullity = constraint_1;
            if (sum for j in wt2posList list L_j) + nullity - i < 0 then (
                isPassed = false;
                break;
            );
        );
        if isPassed then recursiveLoop(pos+1, remainingDim-i, append(L, i));
    );
);
recursiveLoop(0, r-4, {});


-- The maximum number of variables. We compute this to make the polynomial ring.
numvars = max for dimComb in dimCombs list (
    for i from 0 to #dimComb-1 do d_(weights_i) = dimComb_i;
    involvedWT = select(weights, wt -> d_wt != 0);
    sum for wt in involvedWT list d_wt * (#ws_wt-d_wt)
);


-- Make the polynomial ring
S = kk[x_0..x_(numvars-1)];
arrowsS = hashTable for i in keys arrows list i => sub(arrows#i, S);


-- Iterate for every possible combination of dimensions to get B-fixed subspaces
listBfixed = {};
count = 0;
for dimComb in dimCombs do (
    <<"Get B-fixed subspaces. In the loop: "<<count<<"/"<<#dimCombs<<endl;
    count = count+1;
    for i from 0 to #dimComb-1 do d_(weights_i) = dimComb_i;

    -- Now we have d_wt for every weight wt
    involvedWT = select(weights, wt -> d_wt != 0);

    for wt in involvedWT do (
        dimSubsp = d_wt;
        dimWS = #(ws_wt);
        pivot_wt = subsets(dimWS, dimSubsp);
        ind_wt = for i from 0 to #(pivot_wt)-1 list i;
    );

    listComb = for i in ind_(involvedWT_0) list {i};
    for i from 1 to #involvedWT-1 do(
        listComb = listComb ** ind_(involvedWT_i);
        listComb = for j in listComb list flatten toList j;
    );

    for comb in listComb do (
        getPivots = wt -> (
            posWT = position(involvedWT, i -> i == wt);
            posPivot = comb_posWT;
            return pivot_wt_posPivot;
        );

        L = {};

        numVar = 0;
        for wt in involvedWT do (
            dimSubsp = d_wt;
            dimWS = #(ws_wt);
            piv = getPivots wt;
            subsp_wt = mutableMatrix table(dimSubsp, dimWS, (i, j) -> 0_S);
            for i from 0 to #piv-1 do (
                pivRow = i;
                pivCol = piv_i;
                (subsp_wt)_(pivRow, pivCol) = 1_S;
                for j from pivCol+1 to dimWS-1 do(
                    if member(j, piv) then continue;
                    (subsp_wt)_(pivRow, j) = x_numVar;
                    numVar = numVar + 1;
                );
            );
        );

        for wt in involvedWT do(
            for wt2 in weights do (
                if d_wt2 != 0 then (
                    pivot2 = getPivots wt2;
                );
                if arrowsS#?{wt, wt2} then (
                    arrowMap = arrowsS#{wt, wt2};
                    img = matrix(subsp_wt) * transpose(arrowMap);
                    if d_wt2 != 0 then (
                        concatenated = mutableMatrix(img || matrix(subsp_wt2));
                        for j in pivot2 do(
                            for i from 0 to numRows img - 1 do (
                                if concatenated_(i, j) != 0 then rowAdd(concatenated, i, -concatenated_(i, j), numRows img + position(pivot2, jj -> jj == j));
                            );
                        );
                        for j from 0 to numColumns img - 1 do (
                            if member(j, pivot2) then continue else (
                                for i from 0 to numRows img - 1 do (
                                    if concatenated_(i, j) != 0 then L = append(L, concatenated_(i, j));
                                );
                            );
                        );
                    ) else (
                        L = join(L, select(flatten entries img, i -> i != 0));
                    );
                );
            );
        );

        if L == {} then I = ideal(0_S) else I = ideal L;
        if (1_S)%I == 0 then continue;
        listBfixed = append(listBfixed, (involvedWT, hashTable for wt in involvedWT list wt => matrix(subsp_wt), I, numVar));
    );
);


-- Get the inverse of x mod I
getInverseMod = (x, I) -> (
    if x == 0 then return (false, null);
    R := ring x;
    denom := (gens I) | (matrix{{x}});
    quot := (1_R) // denom;
    if denom * quot != matrix{{1_R}} then return (false, null) else return (true, quot_(-1, 0));
);


-- Reduce the matrix by doing row operations and column operations using units
reduction = (M, conditions, k) -> (
    nrow := numRows M;
    ncol := numColumns M;
    M = M % conditions;
    numPivot := 0;
    isUnitExist := false;
    inv := 0;
    while true do (
        if numPivot >= k then return (submatrix(M, toList(numPivot..nrow-1), toList(numPivot..ncol-1)), numPivot);
        isUnitExist = false;
        for i from numPivot to nrow-1 do(
            for j from numPivot to ncol-1 do(
                (isUnitExist, inv) = getInverseMod(M_(i, j), conditions);
                if isUnitExist then (
                    M = mutableMatrix M;
                    if numPivot != i then rowSwap(M, numPivot, i);
                    if numPivot != j then columnSwap(M, numPivot, j);
                    rowMult(M, numPivot, inv);
                    for ii from numPivot+1 to nrow-1 do rowAdd(M, ii, -M_(ii, numPivot), numPivot);
                    M = (matrix M) % conditions;
                    numPivot = numPivot + 1;
                    break;
                );
            );
            if isUnitExist then break;
        );
        -- If there is no unit
        if not isUnitExist then return (submatrix(M, toList(numPivot..nrow-1), toList(numPivot..ncol-1)), numPivot);
    );
);


-- Get a most common nonzero entry in the matrix
getMostCommonEntry = (M, I) -> (
    M = M % I;
    listNonzeroEntries := select(flatten entries M, i -> i != 0);
    if #listNonzeroEntries == 0 then return (false, null);
    return (true, first commonest listNonzeroEntries);
);


-- Get the locus of k x k minors of M under the given conditions
-- S represents the original ring without an additional variable
-- accUnit represents the multiplication of entries that are forced to be units by the additional variable
getConditions = (M, k, S, conditions, accUnit) -> (
    R := ring M;
    u := R_(-1);
    M = M % conditions;
    if k <= 0 then return ideal(1_S);
    if k == 1 then ( -- then every entry should be zero
        return sub(eliminate(u, conditions + ideal flatten entries M), S);
    );
    if M == 0 then return sub(eliminate(u, conditions), S); -- always satisfies. Return conditions
    numPivot := 0;
    (M, numPivot) = reduction(M, conditions, k);
    k = k - numPivot;
    if k <= 0 then return ideal(1_S);

    (isNonzeroExist, f) := getMostCommonEntry(M, conditions);
    if not isNonzeroExist then return sub(eliminate(u, conditions), S);

    d := degree(u, f);
    if d != 0 then (
        (C1, C2) := coefficients(f, Variables => u);
        L := {};
        for i from 0 to numColumns C1 - 1 do (
            L = append(L, C2_(i, 0) * accUnit^(d-degree(u, C1_(0, i))));
        );
        f = sum L;
    );

    -- In V(f)
    conditionsVf := trim(conditions + f);
    newM := 0;
    (newM, numPivot) = reduction(M, conditionsVf, k);
    kVf := k - numPivot;
    if kVf >= 1 then conditionsVf = getConditions(newM, kVf, S, conditionsVf, accUnit) else conditionsVf = ideal(1_S);

    -- In D(f)
    conditionsDf := trim(eliminate(u, conditions)+(u*f*accUnit-1));
    M = sub(M, u => u*f);
    (newM, numPivot) = reduction(M, conditionsDf, k);
    kDf := k - numPivot;
    if kDf >= 1 then conditionsDf = getConditions(newM, kDf, S, conditionsDf, f*accUnit) else conditionsDf = ideal(1_S);
    return intersect(conditionsVf, conditionsDf);
);


-- 2110 test via the map E1110 otimes V_1 to (wedge^2)V_1 otimes V_2 otimes V_3
survivor2110 = {};
count = 0;
for infoBfixed in listBfixed do (
    <<"2110 test. In the loop: "<< count << "/" << #listBfixed << endl;
    count = count+1;
    (involvedWT, hashSubsp, I, numVar) = infoBfixed;
    S2 = kk[x_0..x_(numVar-1)];
    S2u = kk[x_0..x_(numVar-1), uu];
    R2 = S2[v_(1, 1)..v_(4, 4), SkewCommutative => toList(v_(1, 1)..v_(4, 4))];

    hashSubspR2 = hashTable for i in keys hashSubsp list i => sub(hashSubsp#i, R2);

    E1110 = for i in Slambda4 list sub(i, R2); -- Every E1110 should contain (wedge^3)V
    for wt in involvedWT do (
        E1110 = join(E1110, flatten entries (hashSubspR2#wt * sub((matrix table(#(ws_wt), 1, (i, j) -> (ws_wt_i))), R2)));
    );

    basis2110 = flatten flatten for b1 in flatten entries sub(basis(2, kk[v_(1, 1)..v_(1, 4), SkewCommutative => toList(v_(1, 1)..v_(1, 4))]), R2) list for b2 in flatten entries sub(basis(1, kk[v_(2, 1)..v_(2, 4)]), R2) list for b3 in flatten entries sub(basis(1, kk[v_(3, 1)..v_(3, 4)]), R2) list sub(b1*b2*b3, R2);
    use R2;

    testMat = null;
    for i from 1 to 4 do (
        for j in E1110 do (
            matToAdd = matrix table(#basis2110, 1, (ii, jj) -> coefficient(basis2110_ii, j * v_(1, i)));
            if testMat === null then testMat = matToAdd else testMat = testMat | matToAdd;
        );
    );

    loc = getConditions(sub(testMat, S2u), numColumns testMat - r + 1, S2, sub(I, S2u), 1_S2u);
    isSurvived = ((1_S2)%loc != 0);
    if isSurvived then survivor2110 = append(survivor2110, (involvedWT, hashSubsp, loc, numVar));
);


-- 1210 test via the map E1110 otimes V_2 to V_1 otimes (wedge^2)V_2 otimes V_3
survivor1210 = {};
count = 0;
for infoBfixed in survivor2110 do (
    <<"1210 test. In the loop: "<< count << "/" << #survivor2110 << endl;
    count = count+1;
    (involvedWT, hashSubsp, I, numVar) = infoBfixed;
    S2 = kk[x_0..x_(numVar-1)];
    S2u = kk[x_0..x_(numVar-1), uu];
    R2 = S2[v_(1, 1)..v_(4, 4), SkewCommutative => toList(v_(1, 1)..v_(4, 4))];

    hashSubspR2 = hashTable for i in keys hashSubsp list i => sub(hashSubsp#i, R2);

    E1110 = for i in Slambda4 list sub(i, R2); -- Every E1110 should contain (wedge^3)V
    for wt in involvedWT do (
        E1110 = join(E1110, flatten entries (hashSubspR2#wt * sub((matrix table(#(ws_wt), 1, (i, j) -> (ws_wt_i))), R2)));
    );

    basis1210 = flatten flatten for b1 in flatten entries sub(basis(2, kk[v_(2, 1)..v_(2, 4), SkewCommutative => toList(v_(2, 1)..v_(2, 4))]), R2) list for b2 in flatten entries sub(basis(1, kk[v_(1, 1)..v_(1, 4)]), R2) list for b3 in flatten entries sub(basis(1, kk[v_(3, 1)..v_(3, 4)]), R2) list sub(b1*b2*b3, R2);
    use R2;

    testMat = null;
    for i from 1 to 4 do (
        for j in E1110 do (
            matToAdd = matrix table(#basis1210, 1, (ii, jj) -> coefficient(basis1210_ii, j * v_(2, i)));
            if testMat === null then testMat = matToAdd else testMat = testMat | matToAdd;
        );
    );

    loc = getConditions(sub(testMat, S2u), numColumns testMat - r + 1, S2, sub(I, S2u), 1_S2u);
    isSurvived = ((1_S2)%loc != 0);
    if isSurvived then survivor1210 = append(survivor1210, (involvedWT, hashSubsp, loc, numVar));
);


-- 1120 test via the map E1110 otimes V_3 to V_1 otimes V_2 otimes (wedge^2)V_3
survivor1120 = {};
count = 0;
for infoBfixed in survivor1210 do (
    <<"1120 test. In the loop: "<< count << "/" << #survivor1210 << endl;
    count = count+1;
    (involvedWT, hashSubsp, I, numVar) = infoBfixed;
    S2 = kk[x_0..x_(numVar-1)];
    S2u = kk[x_0..x_(numVar-1), uu];
    R2 = S2[v_(1, 1)..v_(4, 4), SkewCommutative => toList(v_(1, 1)..v_(4, 4))];

    hashSubspR2 = hashTable for i in keys hashSubsp list i => sub(hashSubsp#i, R2);

    E1110 = for i in Slambda4 list sub(i, R2); -- Every E1110 should contain (wedge^3)V
    for wt in involvedWT do (
        E1110 = join(E1110, flatten entries (hashSubspR2#wt * sub((matrix table(#(ws_wt), 1, (i, j) -> (ws_wt_i))), R2)));
    );

    basis1120 = flatten flatten for b1 in flatten entries sub(basis(2, kk[v_(3, 1)..v_(3, 4), SkewCommutative => toList(v_(3, 1)..v_(3, 4))]), R2) list for b2 in flatten entries sub(basis(1, kk[v_(1, 1)..v_(1, 4)]), R2) list for b3 in flatten entries sub(basis(1, kk[v_(2, 1)..v_(2, 4)]), R2) list sub(b1*b2*b3, R2);
    use R2;

    testMat = null;
    for i from 1 to 4 do (
        for j in E1110 do (
            matToAdd = matrix table(#basis1120, 1, (ii, jj) -> coefficient(basis1120_ii, j * v_(3, i)));
            if testMat === null then testMat = matToAdd else testMat = testMat | matToAdd;
        );
    );

    loc = getConditions(sub(testMat, S2u), numColumns testMat - r + 1, S2, sub(I, S2u), 1_S2u);
    isSurvived = ((1_S2)%loc != 0);
    if isSurvived then survivor1120 = append(survivor1120, (involvedWT, hashSubsp, loc, numVar));
);
-- There is a unique survivor.
if #survivor1120 != 1 then (
    listError = append(listError, "There is more than one survivor");
    error "There is more than one survivor";
);


-- 1111 test
<< "1111 test" << endl;
infoBfixed = first survivor1120; -- The unique survivor
(involvedWT, hashSubsp, I, numVar) = infoBfixed;

-- The survivor does not contain any parameter.
if numVar != 0 then (
    listError = append(listError, "The unique survivor contains parameters");
    error "The unique survivor contains parameters";
);

S2 = kk[x_0..x_(numVar-1)];
R2 = S2[v_(1, 1)..v_(4, 4)];

hashSubspR2 = hashTable for i in keys hashSubsp list i => sub(hashSubsp#i, R2);

E1110 = for i in Slambda4 list sub(i, R2);
for wt in involvedWT do (
    E1110 = join(E1110, flatten entries (hashSubspR2#wt * sub((matrix table(#(ws_wt), 1, (i, j) -> (ws_wt_i))), R2)));
);

-- The unique survivor is E1110 = <u(1,1,1,1),u(1,1,1,2),u(2,1,1,2),u(3,1,1,2),u(1,1,1,3),u(2,1,1,3),u(3,1,1,3)>+Slambda4 where Slambda4 = <u(4,1,2,3),u(4,1,2,4),u(4,1,3,4),u(4,2,3,4)>
use T;
assert(#E1110 == 11);
assert(image matrix{for i in E1110 list sub(i, T)} == image matrix {join({u(1, 1, 1, 1), u(1, 1, 1, 2), u(2, 1, 1, 2), u(3, 1, 1, 2), u(1, 1, 1, 3), u(2, 1, 1, 3), u(3, 1, 1, 3)}, Slambda4)});
use R2;

-- img1110 = E1110 otimes V_4
img1110 = flatten for i in E1110 list for j in toList(v_(4, 1)..v_(4, 4)) list i*j;

-- Get img0111, img1011, img1101 using the symmetry
img0111 = for i in img1110 list doPermute(i, {4, 3, 2, 1});
img1011 = for i in img1110 list doPermute(i, {3, 4, 1, 2});
img1101 = for i in img1110 list doPermute(i, {2, 1, 4, 3});

-- To get (E0111 otimes V_1) cap (E1011 otimes V_2) cap (E1101 otimes V_3) cap (E1110 otimes V_4)
J = intersect(ideal img1110, ideal img0111, ideal img1011, ideal img1101);
J = sub(J, T);

-- hilbertFunction(i,J) returns dim (T/J)_i
-- The following checks dim(J_1)=dim(J_2)=dim(J_3)=0.
for i from 1 to 3 do (
    if hilbertFunction(i, J) != hilbertFunction(i, T) then (
        listError = append(listError, "J has a degree <= 3 part");
        error "J has a degree <= 3 part";
    );
);

-- hilbertFunction(4,T) - hilbertFunction(4,J) gives dim(J_4)
isPassed1111 = (r <= hilbertFunction(4, T) - hilbertFunction(4, J));

<<"listError: " << listError << endl;  -- {}
<< "#dimCombs: " << #dimCombs << endl; -- 228
<< "#listBfixed: " << #listBfixed << endl;  -- 322
<< "#survivor2110: " << #survivor2110 << endl;  -- 128
<< "#survivor1210: " << #survivor1210 << endl;  -- 11
<< "#survivor1120: " << #survivor1120 << endl;  -- 1
<<"isPassed1111: " << isPassed1111 << endl;  -- false
if (listError == {}) and (#survivor1120==1) and (not isPassed1111) then << "The unique survivor does not pass the 1111 test." << endl;
