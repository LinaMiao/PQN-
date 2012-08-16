function test_suite = test
% test for various Spot operators
%
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca
%
initTestSuite;

function testOpDFTR
% test opDFTR using internal test
% added June 14, 2011

n = 100;
A = opDFTR(n);
e = A.test;

assertTrue(e);

function testOpExtension
% test opExtension using internal test
% added June 14, 2011

n  = 100;
nb = 10;
flag = 1;

A = opExtension(n,nb,flag);
e = A.test;

assertTrue(e);

function testOpLInterp1D
% test opLInterp1D using internal test
% added June 14, 2011

x1 = [0:10:1000]';
x2 = [0:2:1000]';

A = opLInterp1D(x1,x2);
e = boolean(A.test);

assertTrue(e);

function testOpLInterp2D
% test opLInterp2D using internal test
% added June 14, 2011

x1 = [0:10:1000]';
y1 = [0:20:2000]';
X2 = [sort(200 + 500*rand(10,1)) sort(200 + 500*rand(10,1))];

A = opLInterp2D(x1,y1,X2);
e = boolean(A.test);

assertTrue(e);

function testOpSmooth
% test opSmooth using internal test
% added June 14, 2011

n = 100;
k = 10;

A = opSmooth(n,k);
e = A.test;

assertTrue(e);

function testOpSpline1D
% test opSpline1D using internal test
% added June 14, 2011

x1 = [0:10:1000]';
x2 = [0:2:1000]';
BC = [1 0];

A = opSpline1D(x1,x2,BC);
e = boolean(A.test);

assertTrue(e);

