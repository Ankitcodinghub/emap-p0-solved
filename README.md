# emap-p0-solved
**TO GET THIS SOLUTION VISIT:** [EMAP P0 Solved](https://www.ankitcodinghub.com/product/emap-prof-costas-constantinou-p0-solved/)


---

📩 **If you need this solution or have special requests:** **Email:** ankitcoding@gmail.com  
📱 **WhatsApp:** +1 419 877 7882  
📄 **Get a quote instantly using this form:** [Ask Homework Questions](https://www.ankitcodinghub.com/services/ask-homework-questions/)

*We deliver fast, professional, and affordable academic help.*

---

<h2>Description</h2>



<div class="kk-star-ratings kksr-auto kksr-align-center kksr-valign-top" data-payload="{&quot;align&quot;:&quot;center&quot;,&quot;id&quot;:&quot;124602&quot;,&quot;slug&quot;:&quot;default&quot;,&quot;valign&quot;:&quot;top&quot;,&quot;ignore&quot;:&quot;&quot;,&quot;reference&quot;:&quot;auto&quot;,&quot;class&quot;:&quot;&quot;,&quot;count&quot;:&quot;2&quot;,&quot;legendonly&quot;:&quot;&quot;,&quot;readonly&quot;:&quot;&quot;,&quot;score&quot;:&quot;5&quot;,&quot;starsonly&quot;:&quot;&quot;,&quot;best&quot;:&quot;5&quot;,&quot;gap&quot;:&quot;4&quot;,&quot;greet&quot;:&quot;Rate this product&quot;,&quot;legend&quot;:&quot;5\/5 - (2 votes)&quot;,&quot;size&quot;:&quot;24&quot;,&quot;title&quot;:&quot;EMAP P0 Solved&quot;,&quot;width&quot;:&quot;138&quot;,&quot;_legend&quot;:&quot;{score}\/{best} - ({count} {votes})&quot;,&quot;font_factor&quot;:&quot;1.25&quot;}">

<div class="kksr-stars">

<div class="kksr-stars-inactive">
            <div class="kksr-star" data-star="1" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" data-star="2" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" data-star="3" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" data-star="4" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" data-star="5" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
    </div>

<div class="kksr-stars-active" style="width: 138px;">
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
            <div class="kksr-star" style="padding-right: 4px">


<div class="kksr-icon" style="width: 24px; height: 24px;"></div>
        </div>
    </div>
</div>


<div class="kksr-legend" style="font-size: 19.2px;">
            5/5 - (2 votes)    </div>
    </div>
Contents

1 Introduction 2

2 The basics behind the FDTD algorithmAssignment Project Exam Help 2

2.1 Defining the Lattice . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 2

2.3 Spatial Step Size . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4

2.5 Absorbing Boundary Conditions . . . . . . . . . . . . . . . . . . . . . . . . . . . 5

3 Examples 5

3.1 example1 — Propagation in Free Space (tmax = 10 ns) . . . . . . . . . . . . . . . 5

3.2 example2 — Propagation in Free Space (tmax = 50 ns) . . . . . . . . . . . . . . . 5

3.3 example3 — Propagation in the Presence of a PEC Obstacle . . . . . . . . . . . 6

4 Assignment 6

5 Submission details 9

6 Statement of expectations 10

7 Code Listing — fdtd 1 10

1 Introduction

The Finite-Difference Time-Domain (FDTD) method is a computational electromagnetic technique for solving for the electric and magnetic fields in arbitrary spatial domains in the time domain. In contrast to techniques such as the Finite Element Method (FEM) and the Method of Moments (MoM), this technique is straightforward to understand and is simple to program. A rudimentary 2D TMz code is included in Section §7 and is used to illustrate the main features of the method.

The aim of this assignment is to use the provided FDTD code in a series of numerical investigations, and compare quantitatively its predictions against theory, which the student is expected to research independently after the completion of the taught part of the EMAP module. A formal report is not required, but your assignment report needs to answer all the assignment questions, in a self-contained manner.

2 The basics behind the FDTD algorithm

2.1 Defining the Lattice

Assignment Project Exam Help

The basic FDTD method (in Cartesian coordinates) makes use of a regular lattice of interleaved electric and magnetic field components as originally proposed by Yee [1]. In the case of a 2D TMz lattice , it is possible to derive the following from Maxwell’s equations:

.

Consider now the (i,j)th TMz lattice cell, as shown in Fig. 1. Using the above notation, it is possible to form the update equations [2, p73ff] for the various field components, given by

Ez|in−+00..55,j+0.5 = Ca|i−0.5,j+0.5Ez|in−−00..55,j+0.5 +

Cb|i−0.5,j+0.5hHy|ni,j+0.5 − Hy|ni−1,j+0.5+ (1)

source

Assignment Project Exam Help

(2)

and

. (3)

The current term Jsource can be used to excite the lattice. The coefficient matrices Ca(·), Cb(·), Da(·) and Db(·) are used to incorporate different materials within the lattice and are are given by

and .

It should be noted that the indices in the coefficient matrices correspond to the locations of the field components that are being updated. Although appearing cumbersome, these update equations can be programmed in a straightforward fashion. Initially, all field components are initialized to zero, and the field components updated in the order Ez (1) followed by Hx (2) and Hy (3). These calculations are then repeated in sequence until sufficient number of iterations have been performed .

2.2 Excitation

It is always necessary to excite the lattice in some fashion, and the specific method by which this is done is problem dependent. One way is to specify the term Jsource according to a predefined time sequence; alternatively a given field component can be directly specified in a similar fashion (in the TMz case it is usual to specify a single Ez component). As the simulation progresses, the field will be observed to propagate outwards from the source.

In many cases it is desirable to obtain the response of a system at a fixed frequency3. Although it would seem logical to use sinusoidal excitation to determine this, in practice is is usually better to estimate the impulse response of the system using a wideband pulse such as a Gaussian derivative pulse [3, p88] given by

(4) which provides a unit peak amplitude atAssignment Project Exam Helpt−m = ±s. Plotting the impulse response can provide

a very good qualitative understanding of how the propagating electromagnetic wave interacts with objects in the problem domain.

from the impulse response. For example, if the time-harmonic response at frequency f is required

for the Ez component, the time-harmonic electric field Ez(f) is given by [4, p169]

)] (5)

n=0

In practice, Ez(f) can be calculated by maintaining a separate complex field buffer which is incrementally determined by adding the new contribution at each time step.

2.3 Spatial Step Size

In the FDTD method it is necessary to select a spatial step size of approximately λ/20 in the most electromagnetically ‘dense’ material in the solution domain (i.e. in the region with the greatest value of ǫr).

2.4 Time Step Size

To ensure stability, it is necessary to select a time step size that is less than or equal to the Courant limit. In the case of a uniform mesh in 2D with cell size ∆, the Courant limit is given by [3, p70]

∆tCourant

where u is the speed of light in the most electromagnetically ‘dense’ material in the lattice. In practice a time step of ∆t = 0.95∆tCourant is used to ensure any finite precision rounding errors do not cause numerical instability [5, p31].

2.5 Absorbing Boundary Conditions

The FDTD lattice is, by default, terminated on its periphery by a perfect conductor which acts to reflect any outwardly propagating fields (can you figure out why?). However, this is not appropriate for problems which have open boundaries, in which any outwardly propagating fields should be absorbed. In these cases it is necessary to modify the material coefficient matrices and/or the update equations in the vicinity of the boundaries to minimize any reflections. The development of high performance absorbing boundaries has been an area of active research for some time, and boundaries such as the Uniaxial Perfectly Matched Layer (UPML) [4, p212] and Convolutional Perfectly Matched Layer (CPML) [4, p225] can achieve very high levels of performance. These boundaries can, however, be complex to implement.

A much simpler boundary condition is the Absorbing Boundary Condition (ABC) discussed in

[3, pp82-83]. In the case of the boundary at +x, the new value of a tangential field component is given by Assignment Project Exam Help

A rudimentary FDTD code (fdtd 1) has been written in MATLAB and is included in Section §7. Various examples using this code will be investigated in this section.

3.1 example1 — Propagation in Free Space (tmax = 10 ns)

This example is for the code included in Section §7. The source is located at (20,200), and the total simulation time is 10 ns. You must run the code as is and observe the excitation waveform, and the pulse response at 10 ns. Why do you think it is not meaningful to extract the timeharmonic response from this result? Can you observe any unwanted numerical reflections from the ABC?

3.2 example2 — Propagation in Free Space (tmax = 50 ns)

The magnitudes of the fields in Fig. 3 are noticeably smaller than those in the earlier example, as all propagating fields have encountered the ABC on the periphery of the computational domain at least once. However, the residual field is still of appreciable magnitude, and the only way to reduce these is to use a higher performance absorbing boundary such as the UPML or

CPML. The time-harmonic response in Fig. 4 shows a dominant cylindrically-propagating wave,

Assignment Project Exam Help t (ns)

Figure 2: example2 — Excitation waveform.

3.3 example3 — Propagation in the Presence of a PEC Obstacle

A PEC obstacle has been defined with vertices at (150,100) and (300,250). This is done by including the following definition for pecblocks:

pec_blocks = [150 100 300 250];

As in example2, tmax = 50 ns, and the pulse response at 50 ns is plotted in Fig. 5 and the timeharmonic response in Fig. 6 (the excitation waveform is the same as in example2). The pulse response in Fig. 5 is somewhat complex, as a result from waves reflecting from and diffracting around the PEC block. The effects of reflection can also be seen in Fig. 6 with the presence of a standing wave between the source and the box, and a significantly reduced field amplitude behind the box as a result of diffraction .

4 Assignment

Your report should be in four sections each providing an answer to the following questions. The assessment criteria for each section are, (a) a clear description of what you have done: 10%,

FDTD: example2: Pulse response x 10−4

−50 0 50 100 150 200 250 300 350 400 450

Assignment Project Exam Helpx sample index

z

−50 0 50 100 150 200 250 300 350 400 450 x sample index

Figure 4: example2 — Time-harmonic response.

FDTD: example3: Pulse response x 10−4

−50 0 50 100 150 200 250 300 350 400 450

Assignment Project Exam Helpx sample index

z

−50 0 50 100 150 200 250 300 350 400 450 x sample index

Figure 6: example3 — Time-harmonic response.

(b) presentation of your simulation results in an appropriate form for interpretation and discussion: 20%, brief summary of relevant theory researched (including citations of key references, but avoiding giving an unnecessary tutorial), implementation and calculation of corresponding theoretical predictions: 30%, (c) discussion of numerical and theoretical results: 30%, and (d) drawing conclusions: 10%.

3. Compare quantitatively the field distribution in the shadow region of two cascaded PECAssignment Project Exam Help

knife-edges furthest away from the source, to the corresponding shadow field distribution of a rectangular PEC block obtained by bridging the knife-edges. The knife-edges are

pec blocks = [150 1 151 150 249 1 250 150];

5 Submission details

You need to submit a short report, no more than 8 sides of A4 excluding figures, in 11 point Sans Serif font (e.g. Arial), single line spacing and 1.5 cm margins all round. The report should have a cover and feedback sheet which can be downloaded from the module’s Canvas page and completed with your student ID number clearly visible on all pages.

6 Statement of expectations

An excellent report will provide sufficient information to enable the assessor to reproduce its results independently, will have thoroughly researched the state-of-the-art literature for the appropriate theory to compare with numerical simulation results, will produce insightful discussions and will draw scientifically sound conclusions.

A failing report will generate numerical simulation results which cannot be verified independently for corectness, will simply cite relevant literature without justification of its appropriateness, and will limit its discussions to straight-forward observations.

7 Code Listing — fdtd1

% fdtd_1

%

% TMz FDTD Code for EE4D Assignment #2

% Code written by Dr. M. J. Neve

% clear all; % Clear all variables from memory

% Problem parametersAssignment Project Exam Help

id = ’example1’; % Experiment identification string

f = 1e9; % Frequency (Hz)

samples_per_wavelength = 20; % Samples per wavelength – 20 is good

want_abc = 1; % Should absorbing boundary condition be used? 1-yes, 0-no want_plot_excitation = 1; % Plot excitation waveform? want_plot_pulse = 1; % Plot pulse waveform at final iteration?

want_plot_geometry = 1; % Overplot geometry? want_plot_movie = 0; % Plot/create a movie/animation? ncontours = 100; % Number of contours

% Constants

c = 2.997925e8; % Speed of light in vacuum (m/s) eps_0 = 8.854e-12; % Permittivity of free space (F/m) mu_0 = 4*pi*1e-7; % Permeability of free space (H/m) eta_0 = sqrt(mu_0/eps_0); % Intrinsic impedance of free space (ohms)

% Excitation s = 3e-10; m = 4*s;

xs_idx = 20; % Location of source (Ez field indices) ys_idx = 200;

% Define PEC blocks

% This section can be used to define PEC blocks

% Each row defines a separate PEC block

% Format is [xmin_idx ymin_idx xmax_idx ymax_idx]

% Indices are specified relative to locations of the Ez component % if Ez included components on all sides of cell included

%

%pec_blocks = [200 1 201 200]; pec_blocks = [];

% Preliminary calculations lambda = c / f; % Wavelength in free space

lx = lattice_size_in_wavelengths * lambda; % x lattice size in m

ly = lx; % y lattice size in m

nx = samples_per_wavelength * lattice_size_in_wavelengths; % Number of samples in x direction ny = nx;

dx = lx / (nx – 1); % delta x dy = ly / (ny – 1); % delta y

dt = 0.95 * dx / (c * sqrt(2)); % Time step according to Courant limit nt = ceil(tmax / dt); % Number of time steps required

% Preliminary calculations are now complete

% Preallocate field storage

ez = zeros(nx – 1, ny – 1); % TMz field components

hx = zeros(nx – 1, ny); hy = zeros(nx, ny – 1);

ei_z = zeros(size(ez)); % Time harmonic buffer at frequency f

% Preallocation material constants to free space ca = ones(size(ez));

cb = ones(size(ez)) * dt / (eps_0 * dx);

dax = ones(size(hx)); day = ones(size(hy)); dbx = ones(size(hx)) * dt / (mu_0 * dx); dby = ones(size(hy)) * dt / (mu_0 * dx);

% Process any defined PEC blocks

[n_pec_blocks, temp] = size(pec_blocks); % n_pec_blocks is the number of PEC blocks for ii = 1:n_pec_blocks xmin_idx = pec_blocks(ii,1); ymin_idx = pec_blocks(ii,2); xmax_idx = pec_blocks(ii,3);

ca(xmin_idx:xmax_idx,ymin_idx:ymax_idx) = -1.0;Assignment Project Exam Help ymax_idx = pec_blocks(ii,4);

cb(xmin_idx:xmax_idx,ymin_idx:ymax_idx) = 0.0;

% dax,dbx,day,dby do not change because magnetic loss of PEC is still 0 end

t = [0:nt-1] * dt;

v = -exp(0.5)*(t – m) / s .* exp(-(t – m).^2/(2*s^2));

% Main time step loop

ez = ca .* ez + cb .* (hy(2:nx,:) – hy(1:(nx-1),:) + hx(:,2:ny) – hx(:,1:(ny-1))); ez(xs_idx,ys_idx) = v(ii);

if (want_abc)

hx(:,1) = hx(:,1) * (1 – c * dt / dx) + hx(:,2) * c * dt / dx; hx(:,ny) = hx(:,ny) * (1 – c * dt / dx) + hx(:,ny-1) * c * dt / dx; hy(1,:) = hy(1,:) * (1 – c * dt / dx) + hy(2,:) * c * dt / dx;

hy(nx,:) = hy(nx,:) * (1 – c * dt / dx) + hy(nx-1,:) * c * dt / dx;

end

hx(:,2:ny-1) = dax(:,2:ny-1) .* hx(:,2:ny-1) + dbx(:,2:ny-1) .* (ez(:,2:ny-1) – ez(:,1:(ny-2))); hy(2:nx-1,:) = day(2:nx-1,:) .* hy(2:nx-1,:) + dby(2:nx-1,:) .* (ez(2:nx-1,:) – ez(1:(nx-2),:));

% Update time harmonic buffer

ei_z = ei_z + ez * (cos(2*pi*f*(ii-1)*dt) – j * sin(2*pi*f*(ii-1)*dt)); end

% The remaining code is for plotting/visualisation purposes only

if (want_plot_excitation) f0 = figure; plot(t*1e9,v); xlabel(’t (ns)’); ylabel(’v (V)’);

title(sprintf(’FDTD: %s: Excitation waveform’, id)); print(f0, ’-depsc2’, sprintf(’%s_excitation.eps’, id))

end if (want_plot_pulse) % Transpose of data needed to get matrix in correct orientation f1 = figure;

[ch,ch]=contourf(ez’,ncontours); set(ch,’edgecolor’,’none’);

axis equal; %caxis([-0.2 0.2]); colorbar; xlabel(’x sample index’); ylabel(’y sample index’);

title(sprintf(’FDTD: %s: Pulse response’, id)); if (want_plot_geometry) hold on; for ii = 1:n_pec_blocks x1 = (pec_blocks(ii,1) – 1); y1 = (pec_blocks(ii,2) – 1); x2 = pec_blocks(ii,3); y2 = pec_blocks(ii,4);

plot([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],’w’); end

plot(xs_idx, ys_idx, ’wo’); hold off

end

print(f1, ’-depsc2’, sprintf(’%s_pulse.eps’, id))

end

if (want_plot_time_harmonic) f2 = figure;

[ch,ch]=contourf(real(ei_z)’,ncontours);

set(ch,’edgecolor’,’none’);

axis equal;Assignment Project Exam Help

%caxis([-2 2]); colorbar; xlabel(’x sample index’); ylabel(’y sample index’);

if (want_plot_geometry) hold on; for ii = 1:n_pec_blocks x1 = (pec_blocks(ii,1) – 1);

x2 = pec_blocks(ii,3); y2 = pec_blocks(ii,4);

plot([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],’w’); end

plot(xs_idx, ys_idx, ’wo’); hold off

end

print(f2, ’-depsc2’, sprintf(’%s_timeharmonic.eps’, id))

end

if (want_plot_lattice) f3 = figure; hold on; for ii = 1:nx-1 for jj = 1:ny-1 % hx x1 = (ii – 1) * dx; y1 = (jj – 1) * dx; x2 = ii * dx;

y2 = (jj – 1) * dx;

plot([x1 x2],[y1 y2],’b’);

end

end hold off

print(f3, ’-depsc2’, sprintf(’%s_lattice.eps’, id))

end

if (want_plot_movie) f4 = figure;

frame = 0;

mesh(real(ei_z));

set(gca,’nextplot’,’replacechildren’);

for ii = 0:20:340 frame = frame + 1; theta = ii * pi/180;

[ch,ch]=contourf(real(ei_z * exp(j * theta))’,ncontours); set(ch,’edgecolor’,’none’);

axis equal; caxis([-3 3]); colorbar;

xlabel(’x sample index’); ylabel(’y sample index’);

title(sprintf(’FDTD: %s: re(ei_z), f = %5.2f GHz’, id, f/1e9)); if (want_plot_geometry) hold on; for ii = 1:n_pec_blocks x1 = (pec_blocks(ii,1) – 1); y1 = (pec_blocks(ii,2) – 1); x2 = pec_blocks(ii,3); y2 = pec_blocks(ii,4);

plot([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],’w’); end plot(xs_idx, ys_idx, ’wo’); hold off

end drawnow; mov(frame) = getframe;

end

References

[1] K. Yee, “Numerical solution of initial boundary value problems involving Maxwell’s equations in isotropic media,” IEEE Trans. Antennas Propagat., vol. 14, no. 3, pp.302307, 1966.

[2] A. Taflove and S. C. Hagness, Computational electrodynamics: the finite-difference timedomain method. Boston: Artech House, 2005.

[3] D. B. Davidson, Computational Electromagnetics for RF and Microwave Engineering, 2nd ed. Cambridge, 2011.

[4] U. S. Inan and R. A. Marshall, Numerical Electromagnetics: The FDTD Method. Cambridge, 2011.
