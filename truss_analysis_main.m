clear;
clc;
close all;

%% 1. STRUCTURAL DEFINITIONS
fprintf('Defining 34-Joint, 75-Member stable lattice structure...\n');

Joints = [
    1,-24,-3; 2,-16,-2; 3,-8,-1; 4,0,0; 5,10,0; 6,20,0; 7,30,0; 8,40,0; 9,50,0;
    10,60,0; 11,70,0; 12,80,0; 13,90,0; 14,100,0; 15,108,-1; 16,116,-2; 17,124,-3;
    18,-16,4; 19,-8,6; 20,0,8; 21,10,9.44; 22,20,10.56; 23,30,11.36; 24,40,11.84;
    25,50,12; 26,60,11.84; 27,70,11.36; 28,80,10.56; 29,90,9.44; 30,100,8;
    31,108,6; 32,116,4
    33,0,-5;
    34,100,-5;
];

Members = [
    1,1,2; 2,2,3; 3,3,4; 4,4,5; 5,5,6; 6,6,7; 7,7,8; 8,8,9;
    9,9,10; 10,10,11; 11,11,12; 12,12,13; 13,13,14; 14,14,15;
    15,15,16; 16,16,17;
    17,18,19; 18,19,20; 19,20,21; 20,21,22; 21,22,23; 22,23,24;
    23,24,25; 24,25,26; 25,26,27; 26,27,28; 27,28,29; 28,29,30;
    29,30,31; 30,31,32;
    31,2,18; 32,3,19; 33,4,20; 34,5,21; 35,6,22; 36,7,23; 37,8,24;
    38,9,25; 39,10,26; 40,11,27; 41,12,28; 42,13,29; 43,14,30;
    44,15,31; 45,16,32;
    46,1,18; 47,2,19; 48,3,18; 49,3,20; 50,4,19;
    51,20,5; 52,21,6; 53,22,7;
    54,23,8; 55,24,7;
    56,24,9; 57,25,8;
    58,25,10; 59,26,9;
    60,26,11; 61,27,10;
    62,28,11; 63,29,12; 64,30,13;
    65,14,31; 66,15,30; 67,15,32; 68,16,31; 69,17,32
    70,3,33;
    71,4,33;
    72,5,33;
    73,13,34;
    74,14,34;
    75,15,34;
];

%% 2. ANALYSIS ASSUMPTIONS
E = 200e9;
A = 0.01;
numJoints = size(Joints, 1);
numMembers = size(Members, 1);
numDOFs = 2 * numJoints;

MemberLengths = zeros(numMembers, 1);
for i = 1:numMembers
    J1 = Members(i, 2); J2 = Members(i, 3);
    coord1 = Joints(J1, 2:3); coord2 = Joints(J2, 2:3);
    L = sqrt((coord2(1) - coord1(1))^2 + (coord2(2) - coord1(2))^2);
    MemberLengths(i) = L;
end

%% 3. DEFINE DEGREES OF FREEDOM (DOF)
DOF_Map = reshape(1:numDOFs, 2, numJoints)';

Supports = [
    1, 1, 1;
    33, 1, 1;
    34, 1, 1;
    17, 1, 1;
];

restrainedDOFs = [];
for i = 1:size(Supports, 1)
    jointID = Supports(i, 1);
    if Supports(i, 2) == 1, restrainedDOFs = [restrainedDOFs, DOF_Map(jointID, 1)]; end
    if Supports(i, 3) == 1, restrainedDOFs = [restrainedDOFs, DOF_Map(jointID, 2)]; end
end
allDOFs = 1:numDOFs;
freeDOFs = setdiff(allDOFs, restrainedDOFs);

%% 4. ASSEMBLE GLOBAL STIFFNESS MATRIX (K)
K_global = zeros(numDOFs, numDOFs);
for i = 1:numMembers
    J1 = Members(i, 2); J2 = Members(i, 3);
    L = MemberLengths(i);
    x1 = Joints(J1, 2); y1 = Joints(J1, 3);
    x2 = Joints(J2, 2); y2 = Joints(J2, 3);
    c = (x2-x1)/L; s = (y2-y1)/L;
    k_elem_global = (E*A/L) * [c^2, c*s, -c^2, -c*s;
                               c*s, s^2, -c*s, -s^2;
                              -c^2, -c*s, c^2, c*s;
                              -c*s, -s^2, c*s, s^2];
    dofs = [DOF_Map(J1,1), DOF_Map(J1,2), DOF_Map(J2,1), DOF_Map(J2,2)];
    K_global(dofs, dofs) = K_global(dofs, dofs) + k_elem_global;
end
%% 5. APPLY DEAD LOADS
fprintf('Applying dead load (self-weight)...\n');

rho = 7850;
g = 9.81;

MemberWeight = rho * A .* MemberLengths * g;

Loads = zeros(numDOFs, 1);

for i = 1:numMembers
    J1 = Members(i, 2); J2 = Members(i, 3);
    w_half = MemberWeight(i) / 2;
    Loads(DOF_Map(J1, 2)) = Loads(DOF_Map(J1, 2)) - w_half;
    Loads(DOF_Map(J2, 2)) = Loads(DOF_Map(J2, 2)) - w_half;
end

deckLoad = 10e3;
bottomChordJoints = 1:17;
for j = bottomChordJoints
    Loads(DOF_Map(j, 2)) = Loads(DOF_Map(j, 2)) - deckLoad;
end

fprintf('Dead load applied successfully.\n');
%% 6. APPLY LIVE LOADS (on left side of bridge)
fprintf('Applying live load on center-right of bridge...\n');

liveLoad = 50e3;
liveLoadJoints = [24, 25, 26, 27, 28];

for j = liveLoadJoints
    Loads(DOF_Map(j, 2)) = Loads(DOF_Map(j, 2)) - liveLoad;
end

fprintf('Live load applied successfully.\n');

%% 7. SOLVE FOR DISPLACEMENTS
fprintf('Solving for displacements...\n');

U = zeros(numDOFs, 1);
U(freeDOFs) = K_global(freeDOFs, freeDOFs) \ Loads(freeDOFs);

fprintf('Displacement solution complete.\n');

%% 8. POST-PROCESSING AND VISUALIZATION
fprintf('Plotting undeformed and deformed shape...\n');

scaleFactor = 0.05 * max(max(Joints(:, 2:3)) - min(Joints(:, 2:3))) / max(abs(U));

Deformed = Joints;
for j = 1:numJoints
    Deformed(j, 2) = Deformed(j, 2) + scaleFactor * U(DOF_Map(j, 1));
    Deformed(j, 3) = Deformed(j, 3) + scaleFactor * U(DOF_Map(j, 2));
end

figure('Name','Bridge Deformation','Color','w'); hold on; axis equal;
title('Bridge Truss: Undeformed vs Deformed (Scaled)');
xlabel('X-coordinate (m)'); ylabel('Y-coordinate (m)');
grid on;

for i = 1:numMembers
    J1 = Members(i, 2); J2 = Members(i, 3);
    plot([Joints(J1, 2), Joints(J2, 2)], [Joints(J1, 3), Joints(J2, 3)], 'k-', 'LineWidth', 1);
end

for i = 1:numMembers
    J1 = Members(i, 2); J2 = Members(i, 3);
    plot([Deformed(J1, 2), Deformed(J2, 2)], [Deformed(J1, 3), Deformed(J2, 3)], 'b--', 'LineWidth', 1.5);
end

plot(Joints(Supports(:,1),2), Joints(Supports(:,1),3), 'rs', 'MarkerFaceColor','r', 'DisplayName','Supports');
plot(Joints(liveLoadJoints,2), Joints(liveLoadJoints,3), 'gv', 'MarkerFaceColor','g', 'DisplayName','Live Load Joints');

legend('Undeformed','Deformed (scaled)','Supports','Live Load Joints','Location','best');
fprintf('Plot generated successfully.\n');

%% 9. OPTIONAL: Display maximum deflection
[maxDeflection, maxNode] = max(abs(U(2:2:end)));
fprintf('Maximum vertical deflection: %.4f m at Joint %d\n', U(2*maxNode), maxNode);