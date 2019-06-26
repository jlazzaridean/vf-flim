%This function fits the data in decay to an unconstrained 2 exponential
%decay model described by:

% F(tau1, tau2, a1, a2, t) = a1*exp(-t/tau1) + a2*exp(-t/tau2)

%Code was originally adapted (although heavily modified) from open source code
%orginally documented in:
%   Enderlein, J; Erdmann, R, Optics Communications, 1997, 134, 371-378. 
%All edits were made by Julia Lazzari-Dean. 
%Last comments 6/25/2019; last edits approximately 6/10/2019.

%Input parameters:
% decay - time resolved fluorescence signal (one decay at a time), as an N x
%  1 array (N is the ADC resolution)
% irf - measured instrument response, also as an Nx1 array (I recommend
%  trimming only to time bins of interest, e.g. 26-36 out of 256 for the
%  current MaiTai setup). The fit guess is reconvolved with the irf 
%  after application of a color shift (below) before chiSq is evaluated. 
% cShift - starting value of the color shift parameter.
% shiftFixed - if 0, will optimize cShift. If 1, will fix the shift to
%  whatever was passed as cShift.

% Note: this code is hard-coded to set the offset to 0, although the original
% code could also support floating offset. In other words, it assumes
% there is no time-independent background in the FLIM decay. This is reasonable
% with a hybrid detector in a dark room. The way this
%  data was taken, that is effectively 0 (and the code doesn't currently
%  support any other values). Signal before the laser pulse will be
%  interpreted as arising from the previous laser cycle.

%The function returns the following:
% tm - amplitude weighted average decay time (tm) (a1*t1 + a2*t2, where a1 +
% a2 = 1
% aF - weights of each time constant, in the order of the time constants
% tF - fitted decay times, in ascending order.
% cF - color shift, either the fixed value or the fitted value
% offset - should be 0 in this implementation.
% chiSq - goodness of fit, see below.
% exitFlag - the stopping condition from the Matlab optimization routine.

% chiSq = sum((guess(st:fi)-decay(st:fi)).^2./abs(z(st:fi)))/(fi-st-nParam);  
%(Difference between the measured and the fit value for all time bins used,
% weighted by the number of photon counts in that time bin, as is typical
% for Poisson processes.

%Hard-coded default parameters - 
%Start point: all coefficients equal, decay constants: 1, 2.5 ns
%256 time bins of ADC resolution, fitting only to time bins 23-240, laser
%period 12.5 ns

%The function uses the fmincon built-in Matlab function for optimization
%and minimizes the reduced chi squared of the fit, specifically calculated
%as the following (where st and fi are the time bins of the start and
%finish of the fitting, respectively). Note: This is a LOCAL NOT A GLOBAL
%optimization, selected to decrease run time.

function [tm, aFs, tFs, cF, offset, chiSq, exitFlag] = floptimize2_2exp(decay,irf,cShift,shiftFixed)

%initialize parameters correctly for deckard usage
p = 12.5; %time range
n = length(irf);
t = (1:n)'; %time scale in time bin units
dt = p/n;
tp = (0:dt:(p-dt))'; %time scale 
st = 23;
fi = 240;
gStart = [0.5 0.5 1 2.5]; %starting point for a1, a2, t1, t2 (in that order)

%if cShift = -1, let shift float; if it is another value, fix it to that
if (shiftFixed == 0)
    %shift is a free parameter for the function now - start it at 0
    param = [gStart cShift]; 
    nParam = 5; %number of free params we are fitting to
elseif (shiftFixed == 1)
    c = cShift;
    param = gStart; %the only parameter is starting lifetime
    nParam = 4;
else
    error('shiftFixed boolean should be either 0 or 1');
end

%this will lead to an offset of 0 in practice
offFixed = 1;

%setup for fmincon based optimization
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(1,nParam);
ub = zeros(1,nParam);
if (shiftFixed)
    for i=1:nParam
        ub(1,i) = 10;
        lb(1,i) = 0;
    end
else
    for i=1:(nParam-1)
        ub(1,i) = 10;
        lb(1,i) = 0;
    end
    %bound shift between +/- 5 time bin units
    ub(1,nParam) = 5;
    lb(1,nParam) = -5;
end
nonlcon = [];
[pOpt,chiSq,exitFlag] = fmincon(@expScoring,param,A,b,Aeq,beq,lb,ub,nonlcon,options);

    function [chiSquared] = expScoring(param)
        %	This nested function calculates the scoring function       
        %this function has access to all of the values within floptimize.

        if (shiftFixed == 1)
            %evalulate the function again with the current guess
            aTotG = param(1) + param(2); %for re-normalizing coefficients
            a1G = param(1)/aTotG;
            a2G = param(2)/aTotG;
            t1G = param(3);
            t2G = param(4);
            x = a1G*exp(-tp/t1G) + a2G*exp(-tp/t2G);
            %shift the IRF by the designated amount and reconvolve with the
            %result
            irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
            z = convol(irs, x);
        elseif (shiftFixed == 0)
            %evalulate the function again with the current guess
            aTotG = param(1) + param(2); %for re-normalizing coefficients
            a1G = param(1)/aTotG;
            a2G = param(2)/aTotG;
            t1G = param(3);
            t2G = param(4);
            cG = param(5);
            x = a1G*exp(-tp/t1G) + a2G*exp(-tp/t2G);
            %shift the IRF by the current guess for the shift amount and 
            %reconvolve with the current exponential guess
            irs = (1-cG+floor(cG))*irf(rem(rem(t-floor(cG)-1, n)+n,n)+1) + (cG-floor(cG))*irf(rem(rem(t-ceil(cG)-1, n)+n,n)+1);
            z = convol(irs, x);            
        else
            error('Boolean for shift is not set correctly.');
        end
             
        %fix the offset coefficient or leave it floating
        if(offFixed)
            z = [zeros(size(z,1),1) z];
        else
            z = [ones(size(z,1),1) z];
        end
        
        %scale the hypothesized fit to the decay
        scale = lsqnonneg(z,decay);
        z = z*scale;
                
        %%recalculate the scoring function
        chiSquared = sum((z(st:fi)-decay(st:fi)).^2./abs(z(st:fi)))/(fi-st-nParam);
    end

%extract the final values from the output
if (shiftFixed == 1)
    %evalulate the function again with the current guess
    aTotF = pOpt(1) + pOpt(2); %for re-normalizing coefficients
    a1F = pOpt(1)/aTotF;
    a2F = pOpt(2)/aTotF;
    t1F = pOpt(3);
    t2F = pOpt(4);
    cF = c;
elseif (shiftFixed == 0)
    %evalulate the function again with the current guess
    aTotF = pOpt(1) + pOpt(2); %for re-normalizing coefficients
    a1F = pOpt(1)/aTotF;
    a2F = pOpt(2)/aTotF;
    t1F = pOpt(3);
    t2F = pOpt(4);
    cF = pOpt(5);
else
    error('Boolean for shift is not set correctly.');
end

%recalculate the final state of the function
x = a1F*exp(-tp/t1F) + a2F*exp(-tp/t2F);
%shift the IRF by the current guess for the shift amount and
%reconvolve with the current exponential guess
irs = (1-cF+floor(cF))*irf(rem(rem(t-floor(cF)-1, n)+n,n)+1) + (cF-floor(cF))*irf(rem(rem(t-ceil(cF)-1, n)+n,n)+1);
gC = convol(irs, x);
%fix the offset coefficient or leave it floating
if(offFixed)
    gC = [zeros(size(gC,1),1) gC];
else
    gC = [ones(size(gC,1),1) gC];
end
coeff = lsqnonneg(gC,decay);
gCW = gC*coeff;

%from lsqnonneg - coeff(1) is the offset and coeff(2)is the scale factor
%between the two - I THINK.
offset = coeff(1);

%calculate tm at  the final state
tm = a1F*t1F + a2F*t2F;

%put the taus and coefficients into arrays so they can be returned easily
aF = [a1F a2F];
tF = [t1F t2F];
%sort them such that the taus are in ascending order for ease of processing
%later.
[aFs,tFs] = sortATs(aF, tF);

%make a plot of the guesses versus the data
hold off
subplot('position',[0.1 0.4 0.8 0.5])
plot(t,log10(decay),t,log10(irs./max(irs).*max(decay)),t,log10(gCW));
v = axis;
v(1) = min(t);
v(2) = max(t);
axis(v);
xlabel('Time in ns');
ylabel('Log Count');
s = sprintf('Shift, Offset = %3.3f  %3.3f',cF,offset);
text(max(t)/4*3,v(4)-0.05*(v(4)-v(3)),s);
s = sprintf('a1 Weight = %3.3f  %3.3f',aFs);
text(max(t)/4*3,v(4)-0.12*(v(4)-v(3)),s);
s = sprintf('Lifetimes = %3.3f  %3.3f',tFs);
text(max(t)/4*3,v(4)-0.19*(v(4)-v(3)),s);
s = sprintf('Weighted Avg Tau = %3.3f',tm);
text(max(t)/4*3,v(4)-0.26*(v(4)-v(3)),s);
subplot('position',[0.1 0.1 0.8 0.2])
plot(t(st:fi),(decay(st:fi)-gCW(st:fi))./sqrt(abs(gCW(st:fi))));
v = axis;
v(1) = min(t);
v(2) = max(t);

axis(v);
xlabel('Time in ns');
ylabel('Residue');
s = sprintf('%3.3f', chiSq);
text(max(t)/2,v(4)-0.1*(v(4)-v(3)),['\chi^2 = ' s]);
set(gcf,'units','normalized','position',[0.01 0.05 0.98 0.83])


end

