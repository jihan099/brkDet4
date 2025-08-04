-- The purpose of this code is to prove the border rank of the determinant tensor of order 4 is 12.
-- For details, see the preprint 'Border rank of the $4 \times 4$ determinant tensor is twelve'
-- by Jong In Han, Jeong-Hoon Ju, Yeongrak Kim.

kk = QQ;


-- In v_(i,j), 'i' represents the order and 'j' represents the vector in the standard basis.
-- For example, v_(2,3) represents the element v_3 of the standard basis {v_1,v_2,v_3,v_4}
-- in the vector space V_2.
S = kk[v_(1,1)..v_(4,4)];


-- doPermute returns the result of the right action to a tensor by permu.
-- Here, permu represents an element in the symmetric group of order 4.
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
        subList = append(subList, v_(i,j)=>v_(inversePermu_(i-1),j));
    );
    );
    return subList;
);
doPermute = (x, permu) -> sub(x, rawPermute permu);


-- As $V_1\otimes V_2\otimes V_3
--    =S_{\lambda_1}V\oplus S_{\lambda_2}V\oplus S_{\lambda_3}V\oplus \Lambda^3 V$,
-- we first find all 7-dimensional B-fixed subspaces of
-- $S_{\lambda_1}V\oplus S_{\lambda_2}V\oplus S_{\lambda_3}V$.


-- The following is the list of (a,b,c) such that u_1(a,b,c) appears in the weight
-- diagram of S_{\lambda_1}V.
indicesInWeightDiagramOfSV_1 = {
	(1,1,1),
	(1,1,2),
	(1,1,3),
	(1,2,2),
	(1,1,4),
	(1,2,3),
	(2,2,2),
	(1,2,4),
	(1,3,3),
	(2,2,3),
	(1,3,4),
	(2,2,4),
	(2,3,3),
	(1,4,4),
	(2,3,4),
	(3,3,3),
	(2,4,4),
	(3,3,4),
	(3,4,4),
	(4,4,4)
};


-- List of (a,b,c) such that u_2(a,b,c) appears in the weight diagram of S_{\lambda_2}V.
indicesInWeightDiagramOfSV_2 = {
	(1,1,2),
	(1,1,3),
	(2,2,1),
	(1,1,4),
	(1,2,3),
	(1,3,2),
	(1,2,4),
	(1,4,2),
	(3,3,1),
	(2,2,3),
	(1,3,4),
	(1,4,3),
	(2,2,4),
	(3,3,2),
	(4,4,1),
	(2,3,4),
	(2,4,3),
	(4,4,2),
	(3,3,4),
	(4,4,3)
};


-- List of (a,b,c) such that u_3(a,b,c) appears in the weight diagram of S_{\lambda_3}V.
indicesInWeightDiagramOfSV_3 = {
	(1,1,2),
	(1,1,3),
	(2,2,1),
	(1,1,4),
	(1,2,3),
	(1,3,2),
	(1,2,4),
	(1,4,2),
	(3,3,1),
	(2,2,3),
	(1,3,4),
	(1,4,3),
	(2,2,4),
	(3,3,2),
	(4,4,1),
	(2,3,4),
	(2,4,3),
	(4,4,2),
	(3,3,4),
	(4,4,3)
};


-- Construct the list consisting of $v_a\otimes v_b\otimes v_c$
-- for each (a,b,c) in indicesInWeightDiagramOfSV_i for each i=1,2,3.
for i from 1 to 3 do (
	eltsBeforeYoungSymmetrized_i = for j from 0 to length(indicesInWeightDiagramOfSV_i)-1
									list (v_(1,indicesInWeightDiagramOfSV_i_j_0)
										*v_(2,indicesInWeightDiagramOfSV_i_j_1)
										*v_(3,indicesInWeightDiagramOfSV_i_j_2));
);


-- Symmetrize the elements in eltsBeforeYoungSymmetrized_i by the Young symmetrizer $c_{\lambda_i}$.
-- This gives the list of entries in the weight diagram of S_{\lambda_i}V.
eltsSV_1 = for j in eltsBeforeYoungSymmetrized_1 list (j + doPermute(j,{2,1,3}) + doPermute(j,{3,2,1}) + doPermute(j,{1,3,2})
														+ doPermute(j,{2,3,1}) + doPermute(j,{3,1,2}));
eltsSV_2 = for j in eltsBeforeYoungSymmetrized_2 list (j + doPermute(j,{2,1,3}) - doPermute(j,{3,2,1}) - doPermute(j,{3,1,2}));
eltsSV_3 = for j in eltsBeforeYoungSymmetrized_3 list (j + doPermute(j,{3,2,1}) - doPermute(j,{2,1,3}) - doPermute(j,{2,3,1}));


-- Create a MutableHashTable whose key is an element in eltsSV_i and
-- the value is the position of that element in the list in eltsSV_i.
positionInEltsSV = new MutableHashTable;
for i from 1 to 3 do (
    for j from 1 to length(eltsSV_i) do (
        positionInEltsSV # (eltsSV_i_(j-1)) = j;
    );
);


-- The list weightArrows_i represents the weight diagram of S_{\lambda_i}V.
weightArrows_1 = {
	{},
	{1},
	{2},
	{2},
	{3},
	{3,4},
	{4},
	{5,6},
	{6},
	{6,7},
	{8,9},
	{8,10},
	{9,10},
	{11},
	{11,12,13},
	{13},
	{14,15},
	{15,16},
	{17,18},
	{19}
};

weightArrows_2 = {
	{},
	{1},
	{1},
	{2},
	{2,3,6},
	{},
	{4,5,8},
	{},
	{5},
	{5},
	{7,9,12},
	{},
	{7,10},
	{9,10},
	{11},
	{11,13,14,17},
	{},
	{15,16},
	{16},
	{18,19}
};

weightArrows_3 = weightArrows_2;


-- pointingInds_(i,j) represents the smallest subspace of S_{\lambda_i}V containing
-- the j-th element of eltsSV_i, the elements pointed by the j-th element in the weight diagram,
-- and the elements pointed by those elements, and so on.
for i from 1 to 3 do (
    for j from 1 to length(eltsSV_i) do (
        pointingInds_(i,j) = set{eltsSV_i_(j-1)};
    );
);
for i from 1 to 3 do(
    for j from 1 to length(weightArrows_i) do(
        for k in weightArrows_i_(j-1) do (
            pointingInds_(i,j) = pointingInds_(i,j) + pointingInds_(i,k);
        );
    );
);


-- Ignore if pointingInds_(i,j) has more than 7 elts since it represents
-- the subspace of dimension > 7 in that case.
-- For each i, indsValidSubsp_i collects the indices j such that pointingInds_(i,j) has less or equal to 7 elts.
indsList_1 = toList(1..20);
indsList_2 = {1,2,3,4,5,7,9,10,11,13,14,15,16,18,19,20};
indsList_3 = indsList_2;
for i from 1 to 3 do (
    indsValidSubsp_i = {};
    for j in indsList_i do (
        if #(pointingInds_(i,j)) <= 7 then indsValidSubsp_i = append(indsValidSubsp_i, j);
    );
);


-- bFixedSubspOfDim_(i,dimSubsp) represents the list of B-fixed subspaces of dimension 'dimSubsp'.
-- It is a list consists of the bases of such subspaces.
for i from 1 to 3 do (
    indsComb_i = join(
        subsets(indsValidSubsp_i, 1),
        subsets(indsValidSubsp_i, 2),
        subsets(indsValidSubsp_i, 3),
        subsets(indsValidSubsp_i, 4),
        subsets(indsValidSubsp_i, 5),
        subsets(indsValidSubsp_i, 6),
        subsets(indsValidSubsp_i, 7)
    );
    bFixedSubspOfDim_(i,0)={set{}}; -- set{} corresponds to the zero vector space
    for dimSubsp from 1 to 7 do (
        bFixedSubspOfDim_(i,dimSubsp) = {};
        for inds in indsComb_i do (
            if #inds == 0 then continue;
            sumset = sum for j in inds list pointingInds_(i,j);
            if #sumset == dimSubsp then (
                bFixedSubspOfDim_(i,dimSubsp) = append(bFixedSubspOfDim_(i,dimSubsp),sumset);

                -- remove redundant subspaces
                bFixedSubspOfDim_(i,dimSubsp) = toList set (bFixedSubspOfDim_(i,dimSubsp));
            );
        );
    );
);


-- Print the number of B-fixed subspace of dimension j in $S_{\lambda_i}V$ as (i,j,num of subspaces).
for i from 1 to 3 do (
    for j from 1 to 7 do (
		<< "Number of B-fixed subspace of dim " << j << " in S_{lambda_" << i 
		<< "}V: "  << #(bFixedSubspOfDim_(i,j)) << endl;
	);
);

-- The set of whole entries in weight diagrams in
-- $S_{\lambda_1}V^*\oplus S_{\lambda_2}V^*\oplus S_{\lambda_3}V^*$.
-- Hence there are 60 entries.
-- In the code, we borrow the list eltsSV_i to represent the entries.
-- Note that eltsSV_i originally represents elements in S_{\lambda_i}V.
wholeEntries = set join(eltsSV_1, eltsSV_2, eltsSV_3);


-- Eliminate 7 entries from 60 entries.
-- candidatesF1110 is the list of such 53 entries for each selection of 7 vectors obtained above.
-- indicesE1110 is a sequence of three list L1, L2, L3.
-- L1 corresponds to the B-fixed subspace of S_{\lambda_1}V contained in E1110,
-- and the same for L2 and L3.
candidatesF1110 = {};
indicesE1110 = {};
for dim1 from 0 to 7 do (
for dim2 from 0 to 7-dim1 do (
	dim3 = 7-dim1-dim2;

	-- entries1 corresponds to the basis of a B-fixed subspace of dimension dim1 in S_{\lambda_1}V
	-- and the same for entries2 and entries 3.
	for entries1 in (bFixedSubspOfDim_(1,dim1)) do(
	for entries2 in (bFixedSubspOfDim_(2,dim2)) do(
	for entries3 in (bFixedSubspOfDim_(3,dim3)) do(
		candidatesF1110 = append(candidatesF1110, toList(wholeEntries - entries1 - entries2 - entries3));
        
        indicesEntries_1 = for elt in toList(entries1) list (positionInEltsSV#elt);
        indicesEntries_2 = for elt in toList(entries2) list (positionInEltsSV#elt);
        indicesEntries_3 = for elt in toList(entries3) list (positionInEltsSV#elt);

        indicesE1110 = append(indicesE1110,toSequence for i from 1 to 3 list indicesEntries_i);
	);
	);
	);
);
);


-- A mutable hash table whose key is the subspace that is a candidate to F1110 and the value is
-- corresponding indices of E1110.
getIndicesE1110 = new MutableHashTable;
for k from 0 to length(candidatesF1110)-1 do (
    getIndicesE1110 # (candidatesF1110_k) = indicesE1110_k;
)


-- Do (2110)-test to candidates of F1110
survivor2110test={};
dimS4 = hilbertFunction(4, S); -- dimension of S_4
for candidateF1110 in candidatesF1110 do(
    -- img represents the image of the multiplication map
	img = flatten for vec in candidateF1110 list {vec * v_(1,1), vec * v_(1,2), vec * v_(1,3), vec * v_(1,4)};
	
    -- Get the dimension of img
	I=trim ideal img;
	dimImg = dimS4 - hilbertFunction(4, I);
	
	if dimImg <= 160-11 then (
		survivor2110test = append(survivor2110test, candidateF1110);
	);
);
<< "The number of candidates passed (2110)-test: " << #survivor2110test << endl;-- result : only one candidate passed

-- Do (1210)-test to the survivors of the previous test
survivor1210test={};
for candidateF1110 in survivor2110test do (
    -- img represents the image of the multiplication map
	img = flatten for vec in candidateF1110 list {vec * v_(2,1), vec * v_(2,2), vec * v_(2,3), vec * v_(2,4)};
	
	-- Get the dimension of img
	I=trim ideal img;
	dimImg = dimS4 - hilbertFunction(4, I);
	
	if dimImg <= 160-11 then (
		survivor1210test = append(survivor1210test, candidateF1110);
	);
);
<< "The number of candidates passed (1210)-test also: "<< #survivor1210test << endl;-- result : it passes again

-- Do (1120)-test to the survivors of the previous test
survivor1120test={};
for candidateF1110 in survivor1210test do (
	-- img represents the image of the multiplication map
	img = flatten for vec in candidateF1110 list {vec * v_(3,1), vec * v_(3,2), vec * v_(3,3), vec * v_(3,4)};
	
	-- Get the dimension of img
	I=trim ideal img;
	dimImg = dimS4 - hilbertFunction(4, I);
	
	if dimImg <= 160-11 then (
		survivor1120test = append(survivor1120test, candidateF1110);
	)
);
-- result : it passes again
<< "The number of candidates passed (1120)-test also: "<< #survivor1120test << endl; -- result: 1


-- Let candidateF1110 be the unique survivor.
candidateF1110 = first survivor1120test;
<< "Indices of subspaces forming E1110: " << (getIndicesE1110 # candidateF1110) << endl; --  result: ({1, 2, 3}, {1, 2}, {1, 2})

-- Hence the following assertion should be passed.
assert(#(sum{set candidateF1110, pointingInds_(1,3), pointingInds_(2,2), pointingInds_(3,2)})==60);
-- result: passed


-- Do (1111)-test

-- M1110 is the image of the multiplication map of F1110 and V_4^*
-- which is a subspace of $V_1^* \otimes V_2^* \otimes V_3^* \otimes V_4^*$.
M1110 = flatten for vec in candidateF1110 list {
        vec * v_(4,1), vec * v_(4,2), vec * v_(4,3), vec * v_(4,4)
};


-- Now we get another subspace M0111 that should be the image
-- of the multiplication map of F0111 and V_1^*.
-- To get M0111, we permutate the indices of elements in M1110.
M0111 = for entry in M1110 list doPermute(entry, {4,1,2,3});


-- In the same way, we get M1011.
M1011 = for entry in M1110 list doPermute(entry, {1,4,2,3});


-- In the same way, we get M1101.
M1101 = for entry in M1110 list doPermute(entry, {1,2,4,3});


-- The image of the multiplication map
img = join(M1110, M0111, M1011, M1101);
dimImg = hilbertFunction(4, image matrix{img});
<< "The result of (1111)-test: "<< dimImg << endl; -- result: 246 > 256-11. Hence it does not pass (1111)-test.
