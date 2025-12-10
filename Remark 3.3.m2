-- This code is for Remark 3.3 in the paper 'The border rank of the $4 \times 4$ determinant tensor is twelve'
-- by Jong In Han, Jeong-Hoon Ju, Yeongrak Kim.

needsPackage "PieriMaps";

kk = QQ;

-- In v_(i,j), 'i' represents the order and 'j' represents the vector in the standard basis.
-- For example, v_(2,3) represents the element v_3 of the standard basis {v_1,v_2,v_3}
-- in the vector space V_2.
S=kk[v_(1,1)..v_(3,3)];
subS_1=kk[v_(1,1)..v_(1,3)];
subS_2=kk[v_(2,1)..v_(2,3)];
use S;
Sym2V3={v_(3,1)^2,v_(3,1)*v_(3,2),v_(3,1)*v_(3,3),v_(3,2)^2,v_(3,2)*v_(3,3),v_(3,3)^2};
--The order follows standardTableaux(3,{2})

shapeLambda = {     -- Corresponding Young tableaux
    null,           -- Indices of Young tableaux start from 1
    {6},            -- {1,2,3,4,5,6}
    {5,1},          -- {1,2,3,4,5},{6}
    {4,2},          -- {1,2,3,4},{5,6}
    {5,1},          -- {1,2,3,5,6},{4}
    {4,2},          -- {1,2,3,5},{4,6}
    {4,1,1},        -- {1,2,3,5},{4},{6}
    {3,3},          -- {1,2,3},{4,5,6}
    {3,2,1},        -- {1,2,3},{4,5},{6}
    {4,2},          -- {1,2,5,6},{3,4}
    {3,2,1},        -- {1,2,5},{3,4},{6}
    {2,2,2}         -- {1,2},{3,4},{5,6}
}

-- S_{{1,2,3,4,5,6}}V
pieriMap_(1,1)=pieri({6},{1,1},subS_1);
pieriMap_(1,2)=pieri({4},{1,1},subS_2);

-- S_{{1,2,3,4,5},{6}}V
pieriMap_(2,1)=pieri({5,1},{1,2},subS_1);
pieriMap_(2,2)=pieri({4},{1,1},subS_2);

-- S_{{1,2,3,4},{5,6}}V
pieriMap_(3,1)=pieri({4,2},{2,2},subS_1);
pieriMap_(3,2)=pieri({4},{1,1},subS_2);

-- S_{{1,2,3,5,6},{4}}V
pieriMap_(4,1)=pieri({5,1},{1,1},subS_1);
pieriMap_(4,2)=pieri({3,1},{1,2},subS_2);

--S_{{1,2,3,5},{4,6}}V
pieriMap_(5,1)=pieri({4,2},{1,2},subS_1);
pieriMap_(5,2)=pieri({3,1},{1,2},subS_2);

--S_{{1,2,3,5},{4},{6}}V
pieriMap_(6,1)=pieri({4,1,1},{1,3},subS_1);
pieriMap_(6,2)=pieri({3,1},{1,2},subS_2);

--S_{{1,2,3},{4,5,6}}V
pieriMap_(7,1)=pieri({3,3},{2,2},subS_1);
pieriMap_(7,2)=pieri({3,1},{1,2},subS_2);

--S_{{1,2,3},{4,5},{6}}V
pieriMap_(8,1)=pieri({3,2,1},{2,3},subS_1);
pieriMap_(8,2)=pieri({3,1},{1,2},subS_2);

--S_{{1,2,5,6},{3,4}}V
pieriMap_(9,1)=pieri({4,2},{1,1},subS_1);
pieriMap_(9,2)=pieri({2,2},{2,2},subS_2);

--S_{{1,2,5},{3,4},{6}}V
pieriMap_(10,1)=pieri({3,2,1},{1,3},subS_1);
pieriMap_(10,2)=pieri({2,2},{2,2},subS_2);

--S_{{1,2},{3,4},{5,6}}V
pieriMap_(11,1)=pieri({2,2,2},{3,3},subS_1);
pieriMap_(11,2)=pieri({2,2},{2,2},subS_2);

for l from 1 to 11 do (
    pieriMap_(l,1) = sub(pieriMap_(l,1),S); --  first Pieri map
    pieriMap_(l,2) = sub(pieriMap_(l,2),S); --  second Pieri map
    weightList = standardTableaux(3,shapeLambda_l);
    assert(#weightList == numColumns pieriMap_(l,1));
    dimSchur_l = #weightList;

    -- Get the element of S_(lambda_l)V embedded in Sym2V1 otimes Sym2V2 otimes Sym2V3
    -- by composing two Pieri maps.
    for i from 0 to dimSchur_l - 1 do (
        weight_(l,i) = sort(flatten(weightList_i));
        elt_(l,i) = 0;
        for j from 0 to numRows pieriMap_(l,1) - 1 do (
            for k from 0 to numRows pieriMap_(l,2) - 1 do (
                elt_(l,i) = elt_(l,i) + pieriMap_(l,1)_(j,i) * pieriMap_(l,2)_(k,j) * Sym2V3_k;
            );
        );
    );

    -- Get the weight diagram.
    -- The arrows of weight diagrams are saved in weightArrows.
    -- If the weight space has dimension greater than 1, then
    -- the indices with same weight is saved in indsWithSameWeight.
    -- In that case, the weight space is represented by the smallest
    -- corresponding index.
    indsToIgnore_l = {}; -- indices that does not represent the weight space
    for i from 0 to dimSchur_l - 1 do (
        weightArrows_(l,i) = {};
        indsWithSameWeight_(l,i) = {i};
        for j from 0 to dimSchur_l - 1 do (
            if i == j then continue;
            if sort(weight_(l,i) - weight_(l,j))=={0,0,0,0,0,1} then (
                weightArrows_(l,i) = append(weightArrows_(l,i),j);
            );
            if weight_(l,i) - weight_(l,j) == {0,0,0,0,0,0} then (
                if j > i then weightArrows_(l,i) = append(weightArrows_(l,i),j);
                indsWithSameWeight_(l,i) = append(indsWithSameWeight_(l,i),j);
            );
        );
        indsWithSameWeight_(l,i) = sort(indsWithSameWeight_(l,i));
        for j from 1 to #(indsWithSameWeight_(l,i))-1 do (
            indsToIgnore_l = append(indsToIgnore_l,indsWithSameWeight_(l,i)_j);
        );
    );
    indsToIgnore_l = sort toList set(indsToIgnore_l); -- remove redundancy

    validInds_l = for i from 0 to dimSchur_l-1 list if isMember(i,indsToIgnore_l)==false then i else continue;
);


-- pointingInds_(i,j) represents the smallest subspace of S_{\lambda_i}V containing
-- the j-th element of the basis of S_{\lambda_i}V, the elements pointed by the j-th element in the weight diagram,
-- and the elements pointed by those elements, and so on.
for l from 1 to 11 do (
    for i from 0 to dimSchur_l - 1 do (
        if isMember(i,indsToIgnore_l) then continue;
        pointingInds_(l,i) = set join({i},indsWithSameWeight_(l,i),weightArrows_(l,i));
        checkedInds = set join({i},indsWithSameWeight_(l,i));
        indsToGoUp = set(weightArrows_(l,i)) - checkedInds;
        while #indsToGoUp != 0 do (
            pickedInd = first elements indsToGoUp;
            indsToGoUp = indsToGoUp - set{pickedInd};
            checkedInds = checkedInds + set{pickedInd};
            pointingInds_(l,i) = pointingInds_(l,i) + set(weightArrows_(l,pickedInd));
            indsToGoUp = indsToGoUp + set(weightArrows_(l,pickedInd)) - checkedInds;
        );
    );
);


-- Ignore if pointingInds_(l,i) has more than 5 elts since it represents
-- the subspace of dimension > 5 in that case.
-- For each l, indsValidSubsp_l collects the indices i such that pointingInds_(l,i) has less or equal to 5 elts.
for l from 1 to 11 do (
    indsValidSubsp_l = {};
    for i in validInds_l do (
        if #(pointingInds_(l,i)) <= 5 then indsValidSubsp_l = append(indsValidSubsp_l, i);
    );
);


-- bFixedSubspOfDim_(i,dimSubsp) represents the list of B-fixed subspaces of dimension 'dimSubsp'.
for l from 1 to 11 do (
    indsComb_l = join(
        subsets(indsValidSubsp_l, 1),
        subsets(indsValidSubsp_l, 2),
        subsets(indsValidSubsp_l, 3),
        subsets(indsValidSubsp_l, 4),
        subsets(indsValidSubsp_l, 5)
    );
    bFixedSubspOfDim_(l,0)={set{}}; -- only the zero vector space when dimSubsp=0.
    for dimSubsp from 1 to 5 do (
        bFixedSubspOfDim_(l,dimSubsp) = {};
        for inds in indsComb_l do (
            if #inds == 0 then continue;
            sumset = sum for i in inds list pointingInds_(l,i);
            if #sumset == dimSubsp then (
                bFixedSubspOfDim_(l,dimSubsp) = append(bFixedSubspOfDim_(l,dimSubsp),sumset);

                -- remove redundant subspaces
                bFixedSubspOfDim_(l,dimSubsp) = toList set (bFixedSubspOfDim_(l,dimSubsp));
            );
        );
    );
);


-- Print the number of B-fixed subspace of dimension i in $S_{\lambda_l}V$ as (l,i,num of subspaces).
for l from 1 to 11 do (
    for i from 1 to 5 do print(l,i,#(bFixedSubspOfDim_(l,i)));
);


-- The set of whole entries in weight diagrams in the dual space.
-- There are 216 entries. 
-- In the code, we represent without the dual notation for the convenience.
wholeEntries = set flatten for l from 1 to 11 list for i from 0 to dimSchur_l-1 list elt_(l,i);

candidatesF222 = {};
indicesE222FromCandidateF222 = {};
for d1 from 0 to 5 do (
for d2 from 0 to 5-d1 do (
for d3 from 0 to 5-d1-d2 do (
for d4 from 0 to 5-d1-d2-d3 do (
for d5 from 0 to 5-d1-d2-d3-d4 do (
for d6 from 0 to 5-d1-d2-d3-d4-d5 do (
for d7 from 0 to 5-d1-d2-d3-d4-d5-d6 do (
for d8 from 0 to 5-d1-d2-d3-d4-d5-d6-d7 do (
for d9 from 0 to 5-d1-d2-d3-d4-d5-d6-d7-d8 do (
for d10 from 0 to 5-d1-d2-d3-d4-d5-d6-d7-d8-d9 do (
	d11 = 5-d1-d2-d3-d4-d5-d6-d7-d8-d9-d10;
    print(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11);
	for s1 in (bFixedSubspOfDim_(1,d1)) do(
	for s2 in (bFixedSubspOfDim_(2,d2)) do(
	for s3 in (bFixedSubspOfDim_(3,d3)) do(
    for s4 in (bFixedSubspOfDim_(4,d4)) do(
    for s5 in (bFixedSubspOfDim_(5,d5)) do(
    for s6 in (bFixedSubspOfDim_(6,d6)) do(
    for s7 in (bFixedSubspOfDim_(7,d7)) do(
    for s8 in (bFixedSubspOfDim_(8,d8)) do(
    for s9 in (bFixedSubspOfDim_(9,d9)) do(
    for s10 in (bFixedSubspOfDim_(10,d10)) do(
    for s11 in (bFixedSubspOfDim_(11,d11)) do(
        inds1= elements s1;
        inds2= elements s2;
        inds3= elements s3;
        inds4= elements s4;
        inds5= elements s5;
        inds6= elements s6;
        inds7= elements s7;
        inds8= elements s8;
        inds9= elements s9;
        inds10= elements s10;
        inds11= elements s11;
        entries1 = set for i in inds1 list elt_(1,i);
        entries2 = set for i in inds2 list elt_(2,i);
        entries3 = set for i in inds3 list elt_(3,i);
        entries4 = set for i in inds4 list elt_(4,i);
        entries5 = set for i in inds5 list elt_(5,i);
        entries6 = set for i in inds6 list elt_(6,i);
        entries7 = set for i in inds7 list elt_(7,i);
        entries8 = set for i in inds8 list elt_(8,i);
        entries9 = set for i in inds9 list elt_(9,i);
        entries10 = set for i in inds10 list elt_(10,i);
        entries11 = set for i in inds11 list elt_(11,i);
		candidatesF222 = append(candidatesF222, toList(wholeEntries - entries1 - entries2 - entries3 - 
                                                        entries4 - entries5 - entries6 - entries7 - 
                                                        entries8 - entries9 - entries10 - entries11));
	);
	);
	);
    );
    );
    );
	);
	);
    );
    );
    );
);
);
);
);
);
);
);
);
);
);


-- Do (322)-test to candidates of F222
survivor322test={};
count=0;
dimS7 = hilbertFunction(7, S); -- dimension of S_7
for candidateF222 in candidatesF222 do(
    print (count,#candidatesF222,#survivor322test);
    count=count+1;
    -- img represents the image of the multiplication map
	img = flatten for vec in candidateF222 list {vec * v_(1,1), vec * v_(1,2), vec * v_(1,3)};

    -- Get the dimension of img
    I=trim ideal img;
	dimImg = dimS7 - hilbertFunction(7, I);
    print dimImg;

	if dimImg <= 360-5 then (
		survivor322test = append(survivor322test, candidateF222);
	);
);

-- Do (232)-test to the survivors of the previous test
survivor232test={};
count=0;
for candidateF222 in survivor322test do (
    print (count,#survivor322test,#survivor232test);
    count=count+1;
    -- img represents the image of the multiplication map
	img = flatten for vec in candidateF222 list {vec * v_(2,1), vec * v_(2,2), vec * v_(2,3)};
	
	-- Get the dimension of img
    I=trim ideal img;
	dimImg = dimS7 - hilbertFunction(7, I);
    print dimImg;
	
	if dimImg <= 360-5 then (
		survivor232test = append(survivor232test, candidateF222);
	);
);

-- Do (223)-test to the survivors of the previous test
survivor223test={};
count=0;
for candidateF222 in survivor232test do (
	print (count,#survivor232test,#survivor223test);
    count=count+1;
    -- img represents the image of the multiplication map
	img = flatten for vec in candidateF222 list {vec * v_(3,1), vec * v_(3,2), vec * v_(3,3)};
	
	-- Get the dimension of img
    I=trim ideal img;
	dimImg = dimS7 - hilbertFunction(7, I);
    print dimImg;
	
	if dimImg <= 360-5 then (
		survivor223test = append(survivor223test, candidateF222);
	)
);
print(#survivor223test);

-- Get the vanishing locus of each F222.
needsPackage "MultiProjectiveVarieties";
L={}
for i in survivor223test do(
    I=ideal(i);
    X=PP_kk^{2,2,2};
    R=ring X;
    StoR = map(R,S,vars R);
    RI = StoR I;
    Z=projectiveVariety RI;
    L=append(L,decompose Z);
);
tally L -- result: every locus is a point.
