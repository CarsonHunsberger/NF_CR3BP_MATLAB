%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example Script for CR3BP Normal Form Code
% 
% 
% Made by: Carson Hunsberger 03/11/2025
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear; clc; close all;

%% Load normal form data
%   Simply specify a value of the CR3BP mass parameter and normal form
%   truncation degree, and then call the loadNF function to load the data.
% 
%   If loadNF is not called prior to using the other NF functions, data
%   corresponding to the default values of mu and N (located within the
%   loadNF function) will be loaded in.
% 
%   The user is encouraged to edit these default values, thereby
%   eliminating the need to include the following three lines of code when
%   working with their preferred CR3BP system.
% 
%   Available values of mu:
%       0.012154 (Used in academic papers) (Default)
%       0.01215058560962404068939157752993196481839 (Earth-Moon)
%       3.0542e-6 (Sun-Earth)

mu = 0.012154;
N=11;

loadNF(mu,N);

%% Transform from action-angle state to CR3BP state
%   To transform from an action-angle state to the corresponding state in
%   the standard CR3BP frame, simply use the AAtoRTB function.
% 
%   The function takes three (optional) arguments in addition to the
%   action-angle state. These allow the user to specify the libration point
%   (Lpt: 1, 2, or 3), the normal form (nftype: 'Birkhoff' or 'Resonant'),
%   and the method of transformation (method: 'anl' or 'num').
% 
%   If one of these optional arguments is left out, AAtoRTB will use the
%   default values of Lpt=1, nftype='Birkhoff', and method='anl'.
% 
%   The following two function calls are thus equivalent to
%       disp(AAtoRTB(AA))
%       disp(AAtoRTB(AA,method='anl')),
%   however the unnecessary arguments have been left in for clarity.

AA = [0 0.3 0 0 0 0];

disp(AAtoRTB(AA,Lpt=1,nftype='Birkhoff',method='anl'))
disp(AAtoRTB(AA,Lpt=1,nftype='Birkhoff',method='num'))


%% Transform multiple states at once
%   If the user has a large number of action-angle states that they would
%   like to transform to the CR3BP frame, it would be tedious to have to
%   create a loop to transform each state individually.
% 
%   Luckily, AAtoRTB (and indeed all 5 of the other main transformations:
%   AAtoNF, NFtoAA, NFtoRTB, RTBtoNF, and RTBtoAA) are able to handle
%   multiple states at once
% 
%   Furthermore, the output will mimic the shape of the input. If the user
%   prefers to store states in a Nx6 matrix, the output will be Nx6. If the
%   user instead prefers a 6xN matrix, the output will be 6xN.

AA0 = [0 0 0.2 0 0 0];

phi_array = linspace(0,2*pi,50);

AA = zeros(50,6);
AA(:,2) = AA0(2);
AA(:,3) = AA0(3);
AA(:,6) = phi_array;

RTB = AAtoRTB(AA,Lpt=1,nftype='Birkhoff');

figure; plot3(RTB(:,1),RTB(:,2),RTB(:,3),LineWidth=2);
grid on; xlabel('x'); ylabel('y'); zlabel('z'); axis equal;


%% Action Angle Propagation (Birkhoff Lissajous Example)
%   To propagate an action-angle state within the action-angle space,
%   simply use the AAprop function.
% 
%   The required inputs are a vector of times at which the propagated state
%   should be evaluated, as well as the initial action-angle state.
%
%   The optional arguments are Lpt and nftype, which once again assume
%   their default values if left out.
% 
%   In this example, an initial state is defined on a Lissajous torus. The
%   action-angle state is then propagated for 120 time units (enough time
%   to fill out the torus somewhat nicely) before being transformed to the
%   restricted three-body (RTB) space with AAtoRTB.

AA0 = [0 0.1 0.007 0 0 0];

tspan = linspace(0,120,2000);

AA = AAprop(tspan,AA0,Lpt=2,nftype='Birkhoff');

RTB = AAtoRTB(AA,Lpt=2,nftype='Birkhoff');

figure; plot3(RTB(:,1),RTB(:,2),RTB(:,3));
axis equal; grid on; xlabel('x'); ylabel('y'); zlabel('z');

%% Action Angle Propagation (Resonant Lissajous Example)
%   All of the examples thus far have involved the Birkhoff normal form, so
%   let us consider an example that deals with the resonant normal form.
% 
%   Consider a scenario nearly identical to the previous example, with the
%   exception being that we would now like to obtain the resonant version
%   of the Lissajous torus in RTB space.
% 
%   The resonant actions can be obtained from the Birkhoff actions 
%   through the following transformation
%       I1hat = I1
%       I2hat = I2
%       I3hat = I2+I3,
%   Which yields a slightly different initial action-angle state.
% 
%   When working with the resonant normal form, the relation I2hat < I3hat
%   clearly must always hold. Violations of this condition are the most
%   common mistake, and will quickly result in errors being thrown.
% 
%   Additionally, make sure that the correct nftype is specified, since the
%   default of 'Birkhoff' will result in errors or nonsensical outputs.

AA0 = [0 0.1 0.107 0 0 0];

tspan = linspace(0,120,2000);

AA = AAprop(tspan,AA0,Lpt=2,nftype='Resonant');

RTB = AAtoRTB(AA,Lpt=2,nftype='Resonant',method='anl');

figure; plot3(RTB(:,1),RTB(:,2),RTB(:,3));
axis equal; grid on; xlabel('x [DU]'); ylabel('y [DU]'); zlabel('z [DU]');

%% Stable and Unstable Manifold Example
  % In order to highlight two of the "partial" transformations, AAtoNF and
  % NFtoRTB, let us consider a scenario in which a need for them naturally
  % arises: obtaining the stable and unstable manifolds of a Lyapunov orbit
  % 
  % To begin, an action-angle state on a Lyapunov orbit is defined and then
  % propagated for (roughly) one time period, sampling 50 points along the
  % orbit.
  % 
  % These 50 action-angle states are then transformed into the normal form
  % space, where small step-off distances are introduced in the pxtilde
  % direction for the stable manifold or in the xtilde direction for the
  % unstable manifold.
  % 
  % These normal form states are then transformed to their corresponding
  % restricted three-body states, where they are then propagated in the RTB
  % space using RTBprop to obtain the manifolds. Note that the stable
  % manifold states must be propagated backward in time, since they are
  % approaching the Lyapunov orbit.
  % 
  % Only one half of the manifold tubes are plotted in this example.

AA0 = [0 0.4 0 0 0 0];

tspan = linspace(0,3,50);

AA = AAprop(tspan,AA0,Lpt=1,nftype='Birkhoff');

NF0 = AAtoNF(AA,nftype='Birkhoff');
NF = NF0;
NF(:,4) = NF(:,4) + 0.001; 
%This "step-off" value determines how far along the manifold the state lies

RTB = NFtoRTB(NF,Lpt=1,nftype='Birkhoff',method='anl');


manifold_tspan = linspace(0,3.5,50);

% Propagate and plot stable manifold
figure; hold on;
for n=1:length(RTB)
    stable_manifold = RTBprop(-manifold_tspan,RTB(n,:));
    plot(stable_manifold(:,1),stable_manifold(:,2),Color=[0 1 0 0.5]);
end

NF = NF0;
NF(:,1) = NF(:,1) - 0.001;
% Sign of "step-off" value is selected so that we move away from the moon!

RTB = NFtoRTB(NF,Lpt=1,nftype='Birkhoff',method='anl');

% Propagate and plot unstable manifold
for n=1:length(RTB)
    unstable_manifold = RTBprop(manifold_tspan,RTB(n,:));
    plot(unstable_manifold(:,1),unstable_manifold(:,2),Color=[1 0 0 0.5]);
end

plot(RTB(:,1),RTB(:,2),Color=[0 0 1],LineWidth=2);
xlabel('x [DU]'); ylabel('y [DU]'); axis equal; grid on;

%% Inverse Transformation
  % So far, we have only examined the transformations that go in the
  % forward direction, from the action-angle/normal form spaces to
  % restricted three-body space. We will now cover the transformations that
  % go in the opposite direction.
  % 
  % Transforming a restricted three-body state to its corresponding
  % action-angle state can be done by calling the RTBtoAA function.
  % 
  % Just like with AAtoRTB, the libration point, normal form type, and
  % transformation type are all optional arguments that can be included.
  % 
  % When transforming from RTB to AA, the numerical transformation is much
  % more accurate, albeit quite a bit slower than the analytical
  % transformation.
  % 
  % If only a rough estimate of the action-angle state is needed, the
  % analytical transformation suffices. If sensitivity matrices or
  % Jacobians are needed, then the numerical transformation should be used.
  % 
  % This example shows that the composition of the numerical transformation
  % and its inverse is much closer to the identity transformation.

RTBanl = AAtoRTB([0 0.1 0.007 0 0 0],Lpt=1,nftype='Birkhoff',method='anl');
RTBnum = AAtoRTB([0 0.1 0.007 0 0 0],Lpt=1,nftype='Birkhoff',method='num');


format long;
disp(RTBtoAA(RTBanl,Lpt=1,nftype='Birkhoff',method='anl')');
disp(RTBtoAA(RTBnum,Lpt=1,nftype='Birkhoff',method='num')');
format short;


%% Inverse Transformation Example
  % Just like with the forward transformation, the inverse transformation
  % can be broken up into two stages, with functions RTBtoNF and NFtoAA.
  % 
  % Consider the following scenario. An initial action-angle state is
  % defined such that it lies on a somewhat large Lyapunov orbit, and we
  % would like to obtain its propagated state within the RTB frame. Rather
  % than propagating in the action-angle space with AAprop, we want to
  % transform the AA state into the RTB and then propagate it according to
  % the full CR3BP equations of motion.
  % 
  % Naturally, this trajectory eventually departs from the desired
  % Lyapunov, and we would like to quantify this departure by transforming
  % the propagated restricted three-body states back into the normal form
  % space to examine the unstable component, xtilde, and then finally back
  % to the action-angle space, so we can see if the actions truly remain
  % constant.
  % 
  % Run the example and notice how the magnitude of the unstable component
  % increases, while the actions stay fairly constant. Notice also that the
  % numerical transformation is much more accurate.
  % 
  % This example takes a few seconds to run due to the use of the slow but
  % accurate numerical transformation.


AA0 = [0 0.8 0 0 -pi/2 0];
tspan = linspace(0,3,50);

RTB = RTBprop(tspan,AAtoRTB(AA0,Lpt=1,nftype='Birkhoff',method='anl'));

figure; plot(RTB(:,1),RTB(:,2),LineWidth=2);
axis equal; grid on; xlabel('x [DU]'); ylabel('y [DU]');


NFanl = RTBtoNF(RTB,Lpt=1,nftype='Birkhoff',method='anl');
NFnum = RTBtoNF(RTB,Lpt=1,nftype='Birkhoff',method='num');

AAanl = NFtoAA(NFanl,nftype='Birkhoff');
AAnum = NFtoAA(NFnum,nftype='Birkhoff');


%The rest is just plotting stuff
fig = figure; set(fig, 'DefaultTextInterpreter', 'latex');
plot(tspan,NFanl(:,1),LineWidth=2); hold on;
plot(tspan,NFnum(:,1),LineWidth=2); grid on;
xlabel('t [TU]'); ylabel('$\tilde{x}$');
legend('Analytical','Numerical');
ax = gca; ax.FontSize = 14;

figure;
subplot(1,3,1);
plot(tspan,AAanl(:,1),LineWidth=2); hold on;
plot(tspan,AAnum(:,1),LineWidth=2);
grid on; xlabel('t [TU]'); ylabel('I_1');
legend('Analytical','Numerical');
subplot(1,3,2);
plot(tspan,AAanl(:,2),LineWidth=2); hold on;
plot(tspan,AAnum(:,2),LineWidth=2);
grid on; xlabel('t [TU]'); ylabel('I_2');
legend('Analytical','Numerical');
subplot(1,3,3);
plot(tspan,AAanl(:,3),LineWidth=2); hold on;
plot(tspan,AAnum(:,3),LineWidth=2);
grid on; xlabel('t [TU]'); ylabel('I_3');
legend('Analytical','Numerical');

%% Smaller Function Examples
%   We have arrived at the last two smaller (yet often useful) functions.
% 
%   The time derivate of the action-angle state can be found by calling the
%   AApartials function. The desired libration point and normal form type
%   can be passed as optional arguments.
% 
%   The value of the action-angle Hamiltonian can be found using the Heval
%   function. Once again, the desired libration point and normal form type
%   can be passed as optional arguments.

AABirkhoff = [0 0.1 0.007 0 0 0];
AAResonant = [0 0.1 0.107 0 0 0];

disp(AApartials(AABirkhoff,Lpt=2,nftype='Birkhoff'))
disp(AApartials(AAResonant,Lpt=2,nftype='Resonant'))

disp(Heval(AABirkhoff,Lpt=2,nftype='Birkhoff'))
disp(Heval(AAResonant,Lpt=2,nftype='Resonant'))






