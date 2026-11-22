% BENG 110 HW7 - Homogeneous (3-marker) Cauchy strain analysis
% Local functions at the end of a script need MATLAB R2016b or later.
clear; clc;

A = [0, 0]; B = [20, 0]; C = [0, 20];
ref = [A; B; C];

fprintf('TEST 1: HW Prob 8 phantom\n');
analyze(ref, [0, 0; 21, 1; -1, 19], 0);

fprintf('\nTEST 2: Known strain exx=.002, eyy=-.001, exy=.0015, no rotation\n');
H = [0.002, 0.0015; 0.0015, -0.001];
analyze(ref, ref + ref*transpose(H), 0);

fprintf('\nTEST 3: same strain + 10 deg rigid rotation\n');
th = deg2rad(10);
R = [cos(th), -sin(th); sin(th), cos(th)];
analyze(ref, (ref + ref*transpose(H))*transpose(R), 0);

fprintf('\nTEST 4: Rotate axes 30 deg (hand calc: exx prime = 0.002549)\n');
analyze(ref, ref + ref*transpose(H), 30);

fprintf('\nTEST 5: 3D input, triangle tilted in space\n');
c = cos(0.6); s = sin(0.6);
Q = [1, 0, 0; 0, c, -s; 0, s, c];
ref3 = [0, 0, 0; 20, 0, 0; 0, 20, 0] * transpose(Q);
H3 = [0.002, 0.0015, 0; 0.0015, -0.001, 0; 0, 0, 0];
d3 = ref3 * Q;
d3 = (d3 + d3 * transpose(H3)) * transpose(Q);
analyze(ref3, d3, 0);

fprintf('\nTEST 6: collinear markers\n');
try
    analyze([0, 0; 1, 1; 2, 2], [0, 0; 1, 1; 2, 2], 0);
catch err
    fprintf('Error caught: %s\n', err.message);
end

%% ---------------- local functions ----------------
function analyze(ref, dfm, phi_deg)
Hd = dispGradient(ref, dfm);
omega = 0.5*(Hd(2,1) - Hd(1,2));
epsNaive = 0.5*(Hd + transpose(Hd));
[theta, Hc] = removeRotation(Hd);
epsD = 0.5*(Hc + transpose(Hc));
epsR = rotateAxes(epsD, phi_deg);
fprintf('Rotation (small-strain est.) = %8.4f deg\n', rad2deg(omega));
fprintf('Rotation (polar decomp.)     = %8.4f deg\n', rad2deg(theta));
fprintf('Naive Cauchy strain : exx=%9.6f eyy=%9.6f exy=%9.6f\n', epsNaive(1,1), epsNaive(2,2), epsNaive(1,2));
fprintf('De-rotated strain   : exx=%9.6f eyy=%9.6f exy=%9.6f\n', epsD(1,1), epsD(2,2), epsD(1,2));
if phi_deg ~= 0
    fprintf('Axes rotated %g deg : exx=%9.6f eyy=%9.6f exy=%9.6f\n', phi_deg, epsR(1,1), epsR(2,2), epsR(1,2));
end
end

function Hd = dispGradient(ref, dfm)
% 3D input: express both configs in the plane of the reference triangle
if size(ref, 2) == 3
    e1 = ref(2,:) - ref(1,:);
    e1 = e1 / norm(e1);
    n = cross(ref(2,:) - ref(1,:), ref(3,:) - ref(1,:));
    n = n / norm(n);
    e2 = cross(n, e1);
    P = [transpose(e1), transpose(e2)];
    ref2 = (ref - ref(1,:)) * P;
    dfm = (dfm - ref(1,:)) * P;
    ref = ref2;
end
u = dfm - ref;
M = [ones(3,1), ref(:,1), ref(:,2)];
if abs(det(M)) < 1e-12
    error('Markers are collinear: strain is undetermined.');
end
a = M \ u(:,1);
b = M \ u(:,2);
Hd = [a(2), a(3); b(2), b(3)];
end

function [theta, Hc] = removeRotation(Hd)
F = eye(2) + Hd;
theta = atan2(F(2,1) - F(1,2), F(1,1) + F(2,2));
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
Hc = transpose(R) * F - eye(2);
end

function e2 = rotateAxes(e, phi_deg)
p = deg2rad(phi_deg);
Q = [cos(p), sin(p); -sin(p), cos(p)];
e2 = Q * e * transpose(Q);
end

