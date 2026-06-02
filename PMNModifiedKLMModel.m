% Design and implement a modified KLM model to describe the behavior of
% PMN-38, an electrostrictive relaxor ceramic

clear
close all

%% Choose functions to run

% Plotting functions
plotCapData = true;     % Plot capacitance (C-V and C-f) data acquired from Keithley 4200-A (two plots)
plotMagPhaseZ = false;  % Plot impedance magnitude and phase data acquired from NanoVNA (two plots for positive bias voltage, two plots for negative)
plotReImZ = false;      % Plot real and imaginary impedance data acquired from NanoVNA (two plots for positive bias voltage, two plots for negative)
plotZpeaks = false;     % Overlay resonance peaks over real part of impedance (one plot per bias voltage)
plotZfit = true;        % Overlay fitted curves over real and imaginary impedance data (two plots per bias voltage)
plotRSS = false;        % Plot residual sum of squares for the fit at each bias voltage (one plot)
plotCoeffs = true;      % Overlay fitted curves over effective material parameter data (one plot, six subplots)
plotOffsets = false;    % Plot real and imaginary offsets from fitted curves (one plot, two subplots)
plotKLMParams = false;  % Overlay fitted curves over effective KLM parameters (one plot, five subplots)
plot1WayTime = true;    % Overlay simulated one-way temporal impulse response over hydrophone data (one plot per bias voltage)
plot1WayFreq = false;   % Overlay simulated one-way frequency impulse response over hydrophone data (one plot per bias voltage)
plot2WayTime = true;    % Overlay simulated two-way temporal impulse response over pulse-echo data (one plot per bias voltage)
plot2WayFreq = false;   % Overlay simulated two-way frequency impulse response over pulse-echo data (one plot per bias voltage)
plot1WayMax = true;     % Overlay maximum of simulated one-way temporal impulse response over maximum of hydrophone data (one plot)
plot2WayMax = true;     % Overlay maximum of simulated two-way temporal impulse response over maximum of pulse-echo data (one plot)
plotLast = true;        % For data that generates a new plot for each bias voltage, plot only the last bias voltage

% Saving functions
saveWorkspace = false;

%% Define simulation parameters, transducer properties, and filepaths to data

% Define constants and simulation parameters
eps0 = 8.854e-12;   % Permittivity of free space

% File paths and names
samplename = 'C2';
filepath = "Data\NanoVNA\";
filename = strcat(filepath,samplename,'_');

% Define transducer properties
vdc = 0:10:250;
thickness = 865e-6;
E = vdc / thickness;
rho = 7920;         % PMN-38 ceramic density
A = (10e-3)^2;      % Cross-sectional area of transducer

% Define fit types
evenFit = fittype(@(a,b,x) a + b*(x).^2,'independent','x','coefficients',{'a','b'});
oddFit = fittype(@(a,x) a*(x),'independent','x','coefficients',{'a'});
loglinFit = fittype(@(a,b,c,x) a*log10(x)+b*x+c,'independent','x','coefficients',{'a','b','c'});
scaleFit = fittype(@(a,x) a*x, 'independent','x','coefficients',{'a'});


%% Fit dielectric constant data

% Load capacitance data from Keithley
C2MHzfile = "Data\Keithley\C2_C-V_2MHz.CSV";
C0Vfile = "Data\Keithley\C2_C-f_0V.CSV";

% Convert capacitance data into relative permittivity
C2MHzdata = readmatrix(C2MHzfile);
C2MHzVDC = C2MHzdata(:,2);
C2MHz = C2MHzdata(:,1);
epsilon2MHz = C2MHz*thickness/ (A*eps0);
epsilonE = C2MHzVDC / thickness;

C0Vdata = readmatrix(C0Vfile);
C0Vf = C0Vdata(:,2);
C0V = C0Vdata(:,1);
[~,idxcutoff] = min(abs(C0Vf-6e6));
C0Vf = C0Vf(1:idxcutoff);
C0V = C0V(1:idxcutoff);
epsilon0V = C0V*thickness/ (A*eps0);

% Fit data
[epsilon2MHzFit,gof2MHz] = fit(epsilonE,epsilon2MHz,evenFit,'StartPoint',[12e3,0]);
[epsilon0VFit,gof0V] = fit(C0Vf,epsilon0V,loglinFit,'StartPoint',[-2800,0,3e4]);

%% Initialize variables and figures

% Material and fitting parameters

% Positive bias voltages
ktP = zeros(size(vdc));
cDP = zeros(size(vdc));
alphaP = zeros(size(vdc));
epsP = zeros(size(vdc));
RoffP = zeros(size(vdc));
XoffP = zeros(size(vdc));
RSSP = zeros(size(vdc));

% Negative bias voltages
ktN = zeros(size(vdc));
cDN = zeros(size(vdc));
alphaN = zeros(size(vdc));
epsN = zeros(size(vdc));
RoffN = zeros(size(vdc));
XoffN = zeros(size(vdc));
RSSN = zeros(size(vdc));

% NanoVNA Load
ZL = 50;

% Initialize figures
numfig = 0;

if plotCapData
    numfig = numfig + 1;
    fCap2MHz = figure(numfig);
    plt1 = plot(epsilon2MHzFit,'r',epsilonE,epsilon2MHz,'bo');
    set(plt1,'LineWidth',2)
    ylim([12e3 13e3])
    xlabel("Bias Field (V/m)")
    ylabel(strcat(char(949)," (",char(949),"_0)"))
    title("Relative Permittivity as a function of Bias Voltage")
    legend("Data","Fit")
    legend("boxoff")
    fontsize(18,"points")
    numfig = numfig + 1;
    fCap0V = figure(numfig);
    plt2 = plot(epsilon0VFit,'r',C0Vf,epsilon0V,'bo');
    set(plt2,'LineWidth',2)
    xlabel("Frequency (Hz)")
    ylabel(strcat(char(949)," (",char(949),"_0)"))
    title("Relative Permittivity as a function of Frequency")
    xscale log
    legend("Data","Fit")
    legend("boxoff")
    fontsize(18,"points")
end

if plotMagPhaseZ
    numfig = numfig + 1;
    fZmagP = figure(numfig);
    hold on
    numfig = numfig + 1;
    fZmagN = figure(numfig);
    hold on
    numfig = numfig + 1;
    fZphaseP = figure(numfig);
    hold on
    numfig = numfig + 1;
    fZphaseN = figure(numfig);
    hold on
end

if plotReImZ
    numfig = numfig + 1;
    fZrealP = figure(numfig);
    hold on
    numfig = numfig + 1;
    fZrealN = figure(numfig);
    hold on
    numfig = numfig + 1;
    fZimagP = figure(numfig);
    hold on
    numfig = numfig + 1;
    fZimagN = figure(numfig);
    hold on
end

%% Process and fit impedance curves

for i=1:numel(E)
    for j=1:2

        if j==1
            polarity = '+';
            file = strcat(filename,'P',num2str(vdc(i)),'VDC.S1P');
        elseif j==2
            polarity = '-';
            file = strcat(filename,'N',num2str(vdc(i)),'VDC.S1P');
        end

        % Read in NanoVNA data
        data = readtable(file,'FileType','text');
        f = data.Var1;
        S11R = data.Var2;
        S11X = data.Var3;
        Zin = ZL*(1+S11R+1i*S11X) ./ (1-S11R-1i*S11X);
    
        w = 2*pi*f;
    
        Zin = smoothdata(Zin,"gaussian",20);
        Zreal = real(Zin);
        Zimag = imag(Zin);
    
    
        % Identify peaks to find speed of sound - average fundamental and third
        % harmonic
    
        % Fundamental
        fmin = 1.5e6;
        fmax = 3.5e6;
        [~,idxmin] = min(abs(f-fmin));
        [~,idxmax] = min(abs(f-fmax));   
        [pks1,locs1,~,prom] = findpeaks(Zreal(idxmin:idxmax),f(idxmin:idxmax));
        [~,idxpk1] = max(prom);
        f01 = locs1(idxpk1);
    
        % Third Harmonic
        fmin = 7e6;
        fmax = 8.5e6;
        [~,idxmin] = min(abs(f-fmin));
        [~,idxmax] = min(abs(f-fmax));    
        [pks3,locs3,~,prom] = findpeaks(Zreal(idxmin:idxmax),f(idxmin:idxmax));
        [~,idxpk3] = max(prom);
        f03 = locs3(idxpk3)/3;
    
        % Average and calculate speed of sound
        f0 = (f01 + f03) / 2;
        c = 2*thickness*f0;
    
    
        % Fit KLM input impedance (resonance method)
    
        % Select frequency range right around fundamental resonance for fitting
        fmin = 2.05e6;
        fmax = 3.05e6;
        [~,idxmin] = min(abs(f-fmin));
        [~,idxmax] = min(abs(f-fmax));
    
        % Calculate epsilon from Keithley data
        epsilonoff = epsilon2MHzFit(E(i)) - epsilon0VFit(2e6);
    
        realFunc = @(w,x) (-thickness./(w.*(A*eps0*(epsilon0VFit(w/(2*pi))+epsilonoff)))).*(2*((x(1))^2)/thickness).*((2*x(2)*exp(x(2)*thickness)*sin((w/c)*thickness)-(w/c)*(exp(2*x(2)*thickness)-1))...
            ./ ((x(2)^2 + (w/c).^2).*(exp(2*x(2)*thickness)+2*exp(x(2)*thickness)*cos((w/c)*thickness)+1))) + x(3);
        imagFunc = @(w,x) (-thickness./(w.*(A*eps0*(epsilon0VFit(w/(2*pi))+epsilonoff)))).*(1 - (2*((x(1))^2)/thickness).*((2*(w/c).*exp(x(2)*thickness).*sin((w/c)*thickness)+x(2)*(exp(2*x(2)*thickness)-1))...
            ./ ((x(2)^2 + (w/c).^2).*(exp(2*x(2)*thickness)+2*exp(x(2)*thickness)*cos((w/c)*thickness)+1)))) + x(4);
        x0 = [0.2,80,0,0];
        Zfunc = @(x,w) [realFunc(w,x) imagFunc(w,x)];
        xdata = w(idxmin:idxmax);
        ydata = [Zreal(idxmin:idxmax) Zimag(idxmin:idxmax)];
        [x,~] = lsqcurvefit(Zfunc,x0,xdata,ydata,[0,0,-Inf,-Inf],[1,Inf,Inf,Inf]);

        ymin = min(ydata);
        ymax = max(ydata);
        ypred = Zfunc(x,xdata);

        % Residual sum of squares of normalized data
        RSS = sum(((ypred(:,1) - ydata(:,1)) / (ymax(1) - ymin(1))).^2 + ((ypred(:,2) - ydata(:,2)) / (ymax(2) - ymin(2))).^2);
    
        % Calcualte and store material parameters
        if j==1
            cDP(i) = rho*c^2;
            ktP(i) = abs(x(1));
            alphaP(i) = x(2);
            epsP(i) = epsilon0VFit(c/(2*thickness)) + epsilonoff;
            RoffP(i) = x(3);
            XoffP(i) = x(4);
            RSSP(i) = RSS;
        elseif j==2
            cDN(i) = rho*c^2;
            ktN(i) = -abs(x(1));
            alphaN(i) = x(2);
            epsN(i) = epsilon0VFit(c/(2*thickness)) + epsilonoff;
            RoffN(i) = x(3);
            XoffN(i) = x(4);
            RSSN(i) = RSS;
        end
    
        % Generate fitted curves
        Zfit = Zfunc(x,w);
        Zrealfit = Zfit(:,1);
        Zimagfit = Zfit(:,2);
        Zfit = Zrealfit + 1i*Zimagfit;
    
    
        % Plot figures
    
        if plotMagPhaseZ
            if j==1
                figure(fZmagP)
                hold on
                plot(f/1e6,abs(Zin))
                figure(fZphaseP)
                hold on
                plot(f/1e6,angle(Zin)*180/pi)
            elseif j==2
                figure(fZmagN)
                hold on
                plot(f/1e6,abs(Zin))
                figure(fZphaseN)
                hold on
                plot(f/1e6,angle(Zin)*180/pi)
            end
        end
    
        if plotReImZ
            if j==1
                figure(fZrealP)
                hold on
                plot(f/1e6,real(Zin))
                figure(fZimagP)
                hold on
                plot(f/1e6,imag(Zin))
            elseif j==2
                figure(fZrealN)
                hold on
                plot(f/1e6,real(Zin))
                figure(fZimagN)
                hold on
                plot(f/1e6,imag(Zin))
            end
        end
    
        if plotZpeaks
            if ~plotLast | (plotLast & i==numel(E))
                numfig = numfig + 1;
                figure(numfig)
                hold on
                plot(f/1e6,Zreal,'b')
                scatter(locs1(idxpk1)/1e6,pks1(idxpk1),'r','LineWidth',2)
                scatter(locs3(idxpk3)/1e6,pks3(idxpk3),'r','LineWidth',2)
                xlabel("Frequency (MHz)")
                ylabel("Re(Z_{in}) (\Omega)")
                title(strcat("Sample ",samplename," ",polarity,num2str(vdc(i))," V"))
                fontsize(18,'points')
            end
        end
    
        if plotZfit
            if ~plotLast | (plotLast & i==numel(E))
                numfig = numfig + 1;
                figure(numfig)
                hold on
                plot(f(idxmin:idxmax)/1e6,real(Zin(idxmin:idxmax)),'b','LineWidth',2)
                plot(f(idxmin:idxmax)/1e6,real(Zfit(idxmin:idxmax)),'r--','LineWidth',2)
                h1 = ylabel("Re(Z_{in}) (\Omega)");
                xlabel("Frequency (MHz)")
                title(strcat("Sample ",samplename," ",polarity,num2str(vdc(i)),"V Re(Z_{in})"))
                legend("Data","Fit")
                legend("boxoff")
                fontsize(18,'points')
                numfig = numfig + 1;
                figure(numfig)
                hold on
                plot(f(idxmin:idxmax)/1e6,imag(Zin(idxmin:idxmax)),'b','LineWidth',2)
                plot(f(idxmin:idxmax)/1e6,imag(Zfit(idxmin:idxmax)),'r--','LineWidth',2)
                h2 = ylabel("Im(Z_{in}) (\Omega)");
                xlabel("Frequency (MHz)")
                title(strcat("Sample ",samplename," ",polarity,num2str(vdc(i))," V Im(Z_{in})"))
                legend("Data","Fit")
                legend("boxoff")
                fontsize(18,'points')
            end
        end
    end
end
    
if plotMagPhaseZ
    figure(fZmagP)
    xlabel("Frequency (MHz)")
    ylabel("|Z_{in}| (\Omega)")
    title(strcat("Sample ",samplename," Positive Bias Magnitude"))
    figure(fZmagN)
    xlabel("Frequency (MHz)")
    ylabel("|Z_{in}| (\Omega)")
    title(strcat("Sample ",samplename," Negative Bias Magnitude"))
    figure(fZphaseP)
    xlabel("Frequency (MHz)")
    ylabel("<Z_{in} (deg.)")
    title(strcat("Sample ",samplename," Positive Bias Phase"))
    figure(fZphaseN)
    xlabel("Frequency (MHz)")
    ylabel("<Z_{in} (deg.)")
    title(strcat("Sample ",samplename," Negative Bias Phase"))
end

if plotReImZ
    figure(fZrealP)
    xlabel("Frequency (MHz)")
    ylabel("Re(Z_{in}) (\Omega)")
    title(strcat("Sample ",samplename," Positive Bias Resistance"))
    figure(fZrealN)
    xlabel("Frequency (MHz)")
    ylabel("Re(Z_{in}) (\Omega)")
    title(strcat("Sample ",samplename," Negative Bias Resistance"))
    figure(fZimagP)
    xlabel("Frequency (MHz)")
    ylabel("Im(Z_{in}) (\Omega)")
    title(strcat("Sample ",samplename," Positive Bias Reactance"))
    figure(fZimagN)
    xlabel("Frequency (MHz)")
    ylabel("Im(Z_{in}) (\Omega)")
    title(strcat("Sample ",samplename," Negative Bias Reactance"))
end

%% Calculate, fit, and plot parameters

% Combine data from positive and negative bias voltages into one vector
vdc = [-fliplr(vdc) vdc(2:end)];
E = vdc / thickness;
cD = [fliplr(cDN(2:end)) (cDN(1)+cDP(1))/2 cDP(2:end)];
kt = [fliplr(ktN(2:end)) (ktN(1)+ktP(1))/2 ktP(2:end)];
alpha = [fliplr(alphaN(2:end)) (alphaN(1)+alphaP(1))/2 alphaP(2:end)];
epsilon = [fliplr(epsN(2:end)) (epsN(1)+epsP(1))/2 epsP(2:end)];
Roff = [fliplr(RoffN(2:end)) (RoffN(1)+RoffP(1))/2 RoffP(2:end)];
Xoff = [fliplr(XoffN(2:end)) (XoffN(1)+XoffP(1))/2 XoffP(2:end)];
RSS = [fliplr(RSSN(2:end)) (RSSN(1)+RSSP(1))/2 RSSP(2:end)];

% Calculate other material parameters
cE = cD.*(1-kt.^2);
e = kt.*sqrt(cD.*epsilon*eps0);

c0 = sqrt(cD/rho);
Z0 = sqrt(cD*rho);

% Exclude data from fitting if the impedance fit was not good
selection = RSS < 1.5;

Esel = E(selection);
cE = cE(selection);
cD = cD(selection);
kt = kt(selection);
e = e(selection);
epsilon = epsilon(selection);
alpha = alpha(selection);
c0 = c0(selection);
Z0 = Z0(selection);
Roff = Roff(selection);
Xoff = Xoff(selection);

% Fit data to polynomial curves
[cEFit,cEgof] = fit(Esel',cE',evenFit,'StartPoint',[15.7e10,-0.01]);
[epsilonFit,epsgof] = fit(Esel',epsilon',evenFit,'StartPoint',[12e3,0]);
[eFit,egof] = fit(Esel',e',oddFit,'StartPoint',1.5e-4);
alphaFit = mean(alpha);

% Generate fitted curves for plotting
Edense = linspace(min(E),max(E),201);

cECurve = cEFit(Edense);
epsilonCurve = epsilonFit(Edense);
eCurve = eFit(Edense);
alphaCurve = ones(size(Edense))*alphaFit;

cDCurve = cECurve + (eCurve.^2)./(epsilonCurve*eps0);
ktCurve = (eCurve) ./ sqrt(cDCurve.*epsilonCurve*eps0);

c0Curve = sqrt(cDCurve/rho);
Z0Curve = sqrt(cDCurve*rho);

% Plot figures

if plotRSS
    numfig = numfig + 1;
    fRSS = figure(numfig);
    plot(E/1e5,RSS,'LineWidth',2)
    xlabel("Bias Field (kV/cm)")
    ylabel("Residual Sum of Squares after Normalization")
    fontsize(18,'points')
end

if plotCoeffs
    numfig = numfig + 1;

    fCoeffs = figure(numfig);

    subplot(3,2,1)
    hold on
    plot(Esel/1e5,cE/1e10,'bo','LineWidth',2)
    plot(Edense/1e5,cECurve./1e10,'r','LineWidth',2)
    ylabel("c^E (N/m^2 x10^{10})")
    
    subplot(3,2,2)
    hold on
    plot(Esel/1e5,cD/1e10,'bo','LineWidth',2)
    plot(Edense/1e5,cDCurve./1e10,'r','LineWidth',2)
    ylabel("c^D (N/m^2 x10^{10})")
    legend("Data","Fit","Location","north")
    legend('boxoff')
    
    subplot(3,2,3)
    hold on
    plot(Esel/1e5,e,'bo','LineWidth',2)
    plot(Edense/1e5,eCurve,'r','LineWidth',2)
    ylabel("e (C/m^2)")
    
    subplot(3,2,4)
    hold on
    plot(Esel/1e5,kt,'bo','LineWidth',2)
    plot(Edense/1e5,ktCurve,'r','LineWidth',2)
    ylabel("k^t")
    
    subplot(3,2,5)
    hold on
    plot(Esel/1e5,epsilon,'bo','LineWidth',2)
    plot(Edense/1e5,epsilonCurve,'r','LineWidth',2)
    xlabel("Bias Field (kV/cm)")
    ylabel(strcat(char(949),"(\omega=\omega_0) (",char(949),"_0)"))
    ylim([11.85e3 11.95e3])
    
    subplot(3,2,6)
    hold on
    plot(Esel/1e5,alpha,'bo','LineWidth',2)
    plot(Edense/1e5,alphaCurve,'r','LineWidth',2)
    xlabel("Bias Field (kV/cm)")
    ylabel("\alpha (Nepers/m)")
    ylim([40 120])

    fontsize(18,'points')
end

if plotOffsets
    numfig = numfig + 1;

    fOffset = figure(numfig);

    subplot(2,1,1)
    hold on
    plot(Esel/1e5,Roff,'bo','LineWidth',2)
    ylabel("R Offset (\Omega)")
    
    subplot(2,1,2)
    hold on
    plot(Esel/1e5,Xoff,'bo','LineWidth',2)
    xlabel("Bias Field (kV/cm)")
    ylabel("X Offset (\Omega)")

    fontsize(18,'points')
end

if plotKLMParams
    numfig = numfig + 1;

    fKLMParams = figure(numfig);

    subplot(3,2,1)
    hold on
    plot(Esel/1e5,c0,'bo','LineWidth',2)
    plot(Edense/1e5,c0Curve,'r','LineWidth',2)
    ylabel("c_0 (m/s)")
    legend("Data","Fit")
    legend('boxoff')
    
    subplot(3,2,2)
    hold on
    plot(Esel/1e5,Z0/1e6,'bo','LineWidth',2)
    plot(Edense/1e5,Z0Curve./1e6,'r','LineWidth',2)
    ylabel("Z_0 (MRayls/m^2)")
    
    subplot(3,2,3)
    hold on
    plot(Esel/1e5,kt,'bo','LineWidth',2)
    plot(Edense/1e5,ktCurve,'r','LineWidth',2)
    ylabel("k^t")
    
    subplot(3,2,4)
    hold on
    plot(Esel/1e5,epsilon,'bo','LineWidth',2)
    plot(Edense/1e5,epsilonCurve,'r','LineWidth',2)
    xlabel("Bias Field (kV/cm)")
    ylabel(strcat(char(949)," (",char(949),"_0)"))
    ylim([11.85e3 11.95e3])
    
    subplot(3,2,5)
    hold on
    plot(Esel/1e5,alpha,'bo','LineWidth',2)
    plot(Edense/1e5,alphaCurve,'r','LineWidth',2)
    xlabel("Bias Field (kV/cm)")
    ylabel("\alpha (Nepers/m)")
    ylim([40 120])

    fontsize(18,'points')
end

%% Set up parameters for modified KLM model

% Simulation parameters
fmin = 1e6;
fmax = 5e6;
npts = 501;
Finc = (fmax-fmin)/(npts-1);
f = linspace(fmin,fmax,npts);
w = 2*pi*f;

% Loaded transducer parameters
A = (10e-3)^2;              % Cross-sectional area [m^2]
lT = 770e-6;                % PMN38 thickness [m]
lG = 920e-9;                % Electrode thickness [m]
lB = 6e-3;                  % Backing layer thickness [m]
lM = 220e-6;                % Matching layer thickness [m]
lF = 0;                     % Front transmission medium thickness [m]

% Material parameters for passive media

% Gold (electrodes)
ZG = 62.6e6;                % Acoustic impedance [Rayls/m^2]
cG = 3240;                  % Speed of sound [m/s]
alphaG = 0;                 % Attenuation [Nepers/m]

% Epoxy 301 (Backing and matching layers)
ZB = 3.05e6;                % Acoustic impedance [Rayls/m^2]
cB = 2650;                  % Speed of sound [m/s]
alphaB = 500;               % Attenuation [Nepers/m]

% Water
ZW = 1.48e6;                % Acoustic impedance [Rayls/m^2]
cW = 1480;                  % Speed of sound [m/s]
alphaW = 0;                 % Attenuation [Nepers/m]

% Air
ZA = 415;                   % Acoustic impedance [Rayls/m^2]
cA = 343;                   % Speed of sound [m/s]
alphaA = 0;                 % Attenuation [Nepers/m]

% Biasing voltages and field strengths used for loaded transducer tests
vdc = -250:25:250;
E = vdc / lT;

% Interpolate material parameters for the specified field strengths
cE = cEFit(E);
e = eFit(E);
alphaT = ones(size(E))*alphaFit;

cD = cE + (e.^2)./(epsilonFit(E)*eps0);
kt = e ./ sqrt(cD.*epsilonFit(E)*eps0);
cT = sqrt(cD/rho);
ZT = sqrt(cD*rho);

% Initialize variables
hTxmax = zeros(size(E));
hTxRxmax = zeros(size(E));
HPmax = zeros(size(E));
PEmax = zeros(size(E));

%% Build modified KLM model and predict acoustic responses

for i=1:numel(E)

    % Estimate epsilon in the simulated frequency range
    epsilonoff = epsilon2MHzFit(E(i)) - epsilon0VFit(2e6);
    epsilon = epsilon0VFit(w/(2*pi)) + epsilonoff;

    % Propagation parameters
    w0 = (pi/lT)*cT(i);
    gT = alphaT(i) + 1i*(w / cT(i));
    gG = alphaG + 1i*(w / cG);
    gB = alphaB + 1i*(w / cB);
    gW = alphaW + 1i*(w / cW);
    gA = alphaA + 1i*(w / cA);
    
    % Electrical parameters
    C0 = A*epsilon*eps0 / lT;                                                   % Clamped capacitance [F]
    Cprime = (-C0 / kt(i)^2) .* (pi*w/w0 ./ (sin(pi*w/w0)));                    % Resonant capacitance [F]
    phi = kt(i)*sqrt(pi./(w0*C0*A*ZT(i))).*(sin(pi*w/(2*w0)) ./ (pi*w/(2*w0))); % Electromechanical turns ratio [N/V]
    Rs = 50;                                                                    % Tx source impedance [ohms]
    Rl = 1e6;                                                                   % Rx input impedance [ohms]
    
    % Initialize impulse response and input impedance
    HTx = zeros(size(w));
    HRx = zeros(size(w));
    HTxRx = zeros(size(w));
    Zin = zeros(size(w));
    
    for j=1:numel(w)
        % Electrical transfer matrices
        ARs = [1 Rs; 0 1];
        ARl = [1 0; 1/Rl 1];
        AC0 = [1 1/(1i*w(j)*C0(j)); 0 1];
        ACprime = [1 1/(1i*w(j)*Cprime(j)); 0 1];
    
        % Transduction transfer matrix
        Axf = [phi(j) 0; 0 1/phi(j)];
    
        % Acoustic transfer matrices
        AT = [cosh(gT(j)*lT/2) A*ZT(i)*sinh(gT(j)*lT/2); (1/(A*ZT(i)))*sinh(gT(j)*lT/2) cosh(gT(j)*lT/2)];
        AG = [cosh(gG(j)*lG) A*ZG*sinh(gG(j)*lG); (1/(A*ZG))*sinh(gG(j)*lG) cosh(gG(j)*lG)];
        AM = [cosh(gB(j)*lM) A*ZB*sinh(gB(j)*lM); (1/(A*ZB))*sinh(gB(j)*lM) cosh(gB(j)*lM)];
        AB = [cosh(gB(j)*lB) A*ZB*sinh(gB(j)*lB); (1/(A*ZB))*sinh(gB(j)*lB) cosh(gB(j)*lB)];
        AQ = AT*AG*AB*[1 0; 1/(A*ZA) 1];
        AP = [1 0; AQ(2,1)/AQ(1,1) 1];
        AF = [cosh(gW(j)*lF) A*ZW*sinh(gW(j)*lF); (1/(A*ZW))*sinh(gW(j)*lF) cosh(gW(j)*lF)];
    
        ATx = ARs*AC0*ACprime*Axf*AP*AT*AG*AM*AF*[1 0; 1/(A*ZW) 1];
        ARx = [1 A*ZW; 0 1]*AF*AM*AG*AT*AP*(Axf^-1)*ACprime*AC0*ARl;

        ATxRx = ATx*ARx;
    
        % Calculate and store impulse response and input impedance data
        HTx(j) = 1/ATx(1,1);
        HRx(j) = 1/ARx(1,1);
        HTxRx(j) = 1/ATxRx(1,1);
        Zin(j) = ATx(1,1)/ATx(2,1);
    end

    % Convert frequency domain predictions into time domain

    % Zero-pad HTx and HTxRx and generate two-sided FFTs
    Fs = 20e6;
    Ts = 1/(2*Fs);
    HTxPad = [zeros(1,fmin/Finc) HTx zeros(1,(Fs-fmax)/Finc)];
    HTxPad = [HTxPad(1) HTxPad(2:end)/2 fliplr(conj(HTxPad(2:end)))/2];
    HTxRxPad = [zeros(1,fmin/Finc) HTxRx zeros(1,(Fs-fmax)/Finc)];
    HTxRxPad = [HTxRxPad(1) HTxRxPad(2:end)/2 fliplr(conj(HTxRxPad(2:end)))/2];

    % Inverse FFTs
    hTx = ifft(HTxPad);
    hTx = ifftshift(hTx);
    hTxRx = ifft(HTxRxPad);
    hTxRx = ifftshift(hTxRx);

    tSim1 = (0:length(hTx)-1)*Ts;
    tSim2 = (0:length(hTxRx)-1)*Ts;


    % Load measured acoustic data for comparison

    % Load one-way impulse response (hydrophone) data
    if vdc(i) > 0
        hydrophoneData = readmatrix(strcat("Data\Hydrophone\PULSE_P",num2str(abs(vdc(i))),"VDC.CSV"));
    elseif vdc(i) < 0
        hydrophoneData = readmatrix(strcat("Data\Hydrophone\PULSE_N",num2str(abs(vdc(i))),"VDC.CSV"));
    elseif vdc(i) == 0
        hydrophoneData = readmatrix(strcat("Data\Hydrophone\PULSE_0VDC.CSV"));
    end
    tMeas1 = hydrophoneData(3:end,1);
    vMeas1 = hydrophoneData(3:end,2);

    % Remove any NaN values
    check1 = ~isnan(tMeas1);
    check2 = ~isnan(vMeas1);
    tMeas1 = tMeas1(check1 & check2);
    vMeas1 = vMeas1(check1 & check2);

    % Load two-way impulse response (pulse-echo) data
    if vdc(i) < 0
        pulseEchoData = readmatrix(strcat("Data\Pulse-Echo\P",num2str(abs(vdc(i))),"VDC.CSV"));
    elseif vdc(i) > 0
        pulseEchoData = readmatrix(strcat("Data\Pulse-Echo\N",num2str(abs(vdc(i))),"VDC.CSV"));
    elseif vdc(i) == 0
        pulseEchoData = readmatrix(strcat("Data\Pulse-Echo\0VDC.CSV"));
    end
    tMeas2 = pulseEchoData(3:end,1);
    vMeas2 = pulseEchoData(3:end,2);

    % Remove any NaN values
    check1 = ~isnan(tMeas2);
    check2 = ~isnan(vMeas2);
    tMeas2 = tMeas2(check1 & check2);
    vMeas2 = vMeas2(check1 & check2);

    % High-pass filter to remove low-frequency artifacts
    FsHP = 1/(tMeas1(2)-tMeas1(1));
    FsPE = 1/(tMeas2(2)-tMeas2(1));
    vMeas1 = highpass(vMeas1,1e6,FsHP);
    vMeas2 = highpass(vMeas2,1e6,FsPE);

    % Find maximum value of measured signal
    tMeasMin1 = 28.5e-6;
    tMeasMax1 = 33.5e-6;
    tMeasMin2 = 59e-6;
    tMeasMax2 = 64e-6;
    [~,idxmin1] = min(abs(tMeas1-tMeasMin1));
    [~,idxmax1] = min(abs(tMeas1-tMeasMax1));
    [~,idxmin2] = min(abs(tMeas2-tMeasMin2));
    [~,idxmax2] = min(abs(tMeas2-tMeasMax2));
    if E(i)<0
        HPmax(i) = -max(abs(vMeas1(idxmin1:idxmax1)));
    else
        HPmax(i) = max(abs(vMeas1(idxmin1:idxmax1)));
    end
    PEmax(i) = max(abs(vMeas2(idxmin2:idxmax2)));

    % Find maximum value of simulated signal
    tSimMin = 60e-6;
    tSimMax = 70e-6;
    [~,idxmin1] = min(abs(tSim1-tSimMin));
    [~,idxmax1] = min(abs(tSim1-tSimMax));
    [~,idxmin2] = min(abs(tSim2-tSimMin));
    [~,idxmax2] = min(abs(tSim2-tSimMax));
    if E(i)==0
        hTxmax(i) = 0;
    elseif E(i)<0
        hTxmax(i) = -max(abs(hTx(idxmin1:idxmax1)));
    elseif E(i)>0
        hTxmax(i) = max(abs(hTx(idxmin1:idxmax1)));
    end
    if E(i)==0
        hTxRxmax(i) = 0;
    else
        hTxRxmax(i) = max(abs(hTxRx(idxmin2:idxmax2)));
    end


    % Convert time domain measured signals into frequency domain

    % FFT of acoustic data
    HPfft = fft(vMeas1,numel(HTxPad));
    PEfft = fft(vMeas2,numel(HTxRxPad));

    % Convert frequency response to dB - normalize to the maximum value
    % within the relevant frequency range    
    L1 = length(HPfft);
    L2 = length(PEfft);
    fHP = (FsHP/L1)*(0:L1-1);
    fPE = (FsPE/L2)*(0:L2-1);

    [~,idxmin] = min(abs(fHP-fmin));
    [~,idxmax] = min(abs(fHP-fmax));
    HPdB = 20*log10(abs(HPfft)/max(abs(HPfft(idxmin:idxmax))));
    [~,idxmin] = min(abs(f-fmin));
    [~,idxmax] = min(abs(f-fmax));
    HTxdB = 20*log10(abs(HTx)/max(abs(HTx(idxmin:idxmax))));
    [~,idxmin] = min(abs(fPE-fmin));
    [~,idxmax] = min(abs(fPE-fmax));
    PEdB = 20*log10(abs(PEfft)/max(abs(PEfft(idxmin:idxmax))));
    [~,idxmin] = min(abs(f-fmin));
    [~,idxmax] = min(abs(f-fmax));
    HTxRxdB = 20*log10(abs(HTxRx)/max(abs(HTxRx(idxmin:idxmax))));


    % Plot figures

    % Shift and scale simulated data to align with measured data
    if E(i)~=0
        idx1 = find(abs(vMeas1)==abs(HPmax(i)),1);
        idx2 = find(abs(vMeas2)==abs(PEmax(i)),1);
        if vMeas1(idx1)>0
            [~,idx3] = max(hTx);
        elseif vMeas1(idx1)<0
            [~,idx3] = max(-hTx);
        end
        if vMeas2(idx2)>0
            [~,idx4] = max(hTxRx);
        elseif vMeas2(idx2)<0
            [~,idx4] = max(-hTxRx);
        end
        deltat1 = tSim1(idx3) - tMeas1(idx1);
        deltat2 = tSim2(idx4) - tMeas2(idx2);
        tSim1 = tSim1 - deltat1;
        tSim2 = tSim2 - deltat2;
        hTx = hTx * HPmax(i) / hTxmax(i);
        hTxRx = hTxRx * PEmax(i) / hTxRxmax(i);
    end

    % Plot one-way temporal impulse response
    if plot1WayTime
        if ~plotLast | (plotLast & i==numel(E))
            numfig = numfig+1;
            figure(numfig)
            hold on
            plot(tMeas1*1e6,vMeas1*1e3,'b','LineWidth',2)
            plot(tSim1*1e6,hTx*1e3,'r--','LineWidth',2)
            xlabel("Time (µs)")
            h3 = ylabel("Signal (mV)");
            xlim([tMeasMin1 tMeasMax1]*1e6)
            title(strcat(num2str(vdc(i))," VDC 1-Way Temporal Impulse Response"))
            legend("Measured","Simulated")
            legend("boxoff")
            fontsize(18,"points")
        end
    end

    % Plot one-way frequency impulse response
    if plot1WayFreq
        if ~plotLast | (plotLast & i==numel(E))
            numfig = numfig+1;
            figure(numfig)
            hold on
            plot(fHP/1e6,HPdB,'b','LineWidth',2)
            plot(f/1e6,HTxdB,'r--','LineWidth',2)
            xlabel("Frequency (MHz)")
            h4 = ylabel("Impulse Response (dB)");
            xlim([fmin/1e6 fmax/1e6])
            title(strcat(num2str(vdc(i))," VDC 1-Way Frequency Impulse Response"))
            legend("Measured","Simulated")
            legend("boxoff")
            fontsize(18,"points")
        end
    end

    % Plot two-way temporal impulse response
    if plot2WayTime
        if ~plotLast | (plotLast & i==numel(E))
            numfig = numfig+1;
            figure(numfig)
            hold on
            plot(tMeas2*1e6,vMeas2,'b','LineWidth',2)
            plot(tSim2*1e6,hTxRx,'r--','LineWidth',2)
            xlabel("Time (µs)")
            h4 = ylabel("Signal (V)");
            xlim([tMeasMin2 tMeasMax2]*1e6)
            title(strcat(num2str(vdc(i))," VDC 2-Way Temporal Impulse Response"))
            legend("Measured","Simulated")
            legend("boxoff")
            fontsize(18,"points")
        end
    end

    % Plot two-way frequency impulse response
    if plot2WayFreq
        if ~plotLast | (plotLast & i==numel(E))
            numfig = numfig+1;
            figure(numfig)
            hold on
            plot(fPE/1e6,PEdB,'b','LineWidth',2)
            plot(f/1e6,HTxRxdB,'r--','LineWidth',2)
            xlabel("Frequency (MHz)")
            ylabel("Impulse Response (dB)")
            xlim([fmin/1e6 fmax/1e6])
            title(strcat(num2str(vdc(i))," VDC 2-Way Frequency Impulse Response"))
            legend("Measured","Simulated")
            legend("boxoff")
            fontsize(18,"points")
        end
    end

end


% Scale data to match
HFit1 = fit(hTxmax',HPmax',scaleFit);
HFit2 = fit(hTxRxmax',PEmax',scaleFit);

% Calculate coefficients of determination (R-squared)
SStot1 = sum((HPmax-mean(HPmax)).^2);           % Total sum-of-squares
SSres1 = sum(((HFit1(hTxmax)' - HPmax)).^2);    % Residual sum-of-squares
Rsq1 = 1-SSres1/SStot1;                         % R-squared

SStot2 = sum((PEmax-mean(PEmax)).^2);           % Total sum-of-squares
SSres2 = sum(((HFit2(hTxRxmax)' - PEmax)).^2);  % Residual sum-of-squares
Rsq2 = 1-SSres2/SStot2;                         % R-squared

if plot1WayMax
    numfig = numfig + 1;
    figure(numfig)
    hold on
    plot(E/1e5,HPmax*1e3,'bo','LineWidth',2)
    plot(E/1e5,HFit1(hTxmax)*1e3,'r','LineWidth',2)
    text(1,-1,strcat("R^2=",num2str(round(Rsq1,4))),"HorizontalAlignment","center")
    xlabel("Bias Field (kV/cm)")
    h5 = ylabel("Maximum Signal (mV)");
    title("One-way Impulse Response")
    legend("Measured","Simulated","Location","northwest")
    legend("boxoff")
    fontsize(18,"points")
end

if plot2WayMax
    numfig = numfig + 1;
    figure(numfig)
    hold on
    plot(E/1e5,PEmax,'bo','LineWidth',2)
    plot(E/1e5,HFit2(hTxRxmax),'r','LineWidth',2)
    text(0,0.3,strcat("R^2=",num2str(round(Rsq2,4))),"HorizontalAlignment","center")
    xlabel("Bias Field (kV/cm)")
    h6 = ylabel("Maximum Signal (V)");
    title("Two-way Impulse Response")
    legend("Measured","Simulated","Location","north")
    legend("boxoff")
    fontsize(18,"points")
end

% Save workspace
if saveWorkspace
    save(strcat(samplename,"_DATA.MAT"))
end