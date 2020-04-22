function [species_name median_f eta beta H F min_veg_height]=wavespeeds_annotated_single_species
%Find the wave speed of a stage-structured species population model with an
% integral redistribution kernel using estimated MGF.

% Concised Codes as used for:
% Bullock et al. Journal of Ecology 2012, 100, 104–115
%
% Note that all functions are set within one file.
% median_f is the output parameter of interest corresponding to cstar in
% Equation 11 of Bullock et al. 2012 = Wavespeed in meters per year (calculated over 25 stochastic runs).

% Codes copyright of CEH:
% Danny Hooftman (dann1@ceh.ac.uk) and Steven White (smwhit@ceh.ac.uk)
% Centre for Ecology & Hydrology
% Benson Lane
% Crowmarsh Gifford
% Wallingford
% Oxon
% OX10 8BB
% UK
% Contact: +44 (0)1491 692699



%% 1. set all the species specific parameters
% H =  Wald parameter: seed release height in meters (needed in Eqns 4 & 5)
% F =  Wald parameter: terminal velocity (needed in Eqns 4 & 5)
% min_veg_height =  Wald parameter: minumum vegetation height (needed in Eqn 6)
% eta = Weibull scale parameter needed in Eq 1. The parameter that was
% changed to create Figure 4.
% beta = Weibull shape parameter needed in Eq 1.
% eta and beta are calcuated seperately from the windspeeds in Weibull.m
% and changed to different scenario's in the master loop
% A = Bo to feed into Equation 9

% set starting set of eta and beta of the Weibull function
% calculated via Weibull_calc.m, 
% Standard =
% Summer: 
%   eta: 4.36617107061204
%   beta: 1.89038522659865
% Winter:
%   eta: 5.04694306475024
%   beta: 1.75080418884777

eta = 5.04694306475024;
beta = 1.75080418884777;

% species is Cirsium
H = 0.1;  % seed release height in meters (Tremlova & Munzbergova 2007)
F=  0.39; %terminal velocity: m/s (Tremlova & Munzbergova 2007) ;
min_veg_height = 0.05; %5 centimeters (Tremlova & Munzbergova 2007) ;
species_name = 'Cirsium';

% 4 stage matrix for Cristine's Circium (based on Munzbergova 2005
A=[0.122 0 0 9.008;...
    0.12 0.338 0.027 1.228;...
    0 0.401 0.944 0.784;...
    0 0.017 0.035 0.169];

% check whether parameters are valid
if H <= 0|| F <= 0|| min_veg_height <= 0
    error('Invalid Wald parameters (<=0), program is stopped')
end

%% 2. set parameters to global = no transfer needed between functions
global H_global;
global F_global;
global veg_height_global;
global species_name_global
global rho_global

H_global = H;
F_global = F;
veg_height_global = min_veg_height;
species_name_global= species_name;


%% 3. Make the CDF
%number of draws from the CDF curve.
Nsamp=1e6; 

% Here CDF for the wald model x weibull model is calculated
% by Bullock et al. Eqns 7 and 8. Go to randgen below for actual
% calcutions

display (' making CDF and PDF')
display ('   ')
[CDF,j] = randgen(eta,beta); 


%% 4. Minimilasation function (Eqn 11 in Bullock et al. 2012)
display (' minimalisation')
display ('   ')

% Make 25 runs to deal with stochasticity (see page 109 of Bullock et al 2012).
for run = 1:1:1
    global x  %#ok<*TLEV>
    % Make the list of actual draws from the CDF of Weibull x Wald to be
    % used in Equation 10.
    y = interp1(CDF,x(j),rand(Nsamp,1));
    
    %set some opimistion tolerences for the mimimisation
    options=optimset('TolFun',1e-4,'TolX',1e-4);
    
    % this is the actual minimalisation function. Equation 11 in Bullock et al. 2012;
    % the function that is mimimised is "func" and is below.
    % 0.0001 = the starting value for s (s > 0); !! s is the parameter that
    % is changing within the minimisation function.
    % !! f = the actual wavespeed per year in meters !!
    
    [minval, f] = fminsearch(@func,0.01,options);
    % calculate the lambda's and sensitivities of the spread matrix
    % rho is the spread matrix x demographci matrix (see func below)
    % Equation 12 in Bullock et al. 2012
     [sens_itivy_full elas_ticity_full lam_bda stable_stage_distribution_full reproductive_vector_full] =sensetivty(rho_global);
    
    % collate the results per run
        if run == 1
           individual_results_min_val = minval;
            individual_results_f = f
            individual_results_sens_itivy  = sens_itivy_full;
            individual_results_elas_ticity = elas_ticity_full;
            individual_results_lam_bda = lam_bda;
              individual_results_stable_stage_distribution = stable_stage_distribution_full;
            individual_results_reproductive_vector = reproductive_vector_full;
        else
            individual_results_min_val = [individual_results_min_val;minval];  %#ok<*AGROW>
            individual_results_f = [individual_results_f;f]
            individual_results_sens_itivy  = [individual_results_sens_itivy ; sens_itivy_full];
            individual_results_elas_ticity = [individual_results_elas_ticity;elas_ticity_full];
             individual_results_lam_bda  = [individual_results_lam_bda; lam_bda];
             individual_results_stable_stage_distribution = [individual_results_stable_stage_distribution;stable_stage_distribution_full];
          individual_results_reproductive_vector = [individual_results_reproductive_vector;reproductive_vector_full];
            
        end
end

% Make a individual output file name
output_file = ['individual_results','_',num2str(eta),'_',num2str(beta),'_',species_name_global,'.mat'];

% Write the output away.
save(output_file,'individual_results_min_val','individual_results_f',...
    'individual_results_lam_bda','individual_results_sens_itivy',...
    'individual_results_elas_ticity', 'individual_results_stable_stage_distribution',...
    'individual_results_reproductive_vector', 'y')    
      
% set the median f (= cstar in Bullock et al. 2012) among runs.
median_f = median(individual_results_f)

%% 5. function func to be called for minimalisation.
    function f = func(s)
        %This is the actual function that is mimised; Bullock et al. Eqn 11,
        % with !!f = cstar = wavespeed per year!!; s = waveshape
        f = (1/s)*log(max(abs(eig(Hmat(s)))));
       
        % the function in here is Hmat calculating the A * M matrix (rho =
        % rho in Bullock et al. 2012 Eqn 11.
        
        % >>> calculate for s
        function rho=Hmat(s)
            %Calculate the H matrix (uses the estimated MGF matrix M).
            % Note that only the reproductive life stages of the matrix are
            % allowed to be spread.
            % Note: besseli is the Bessel function Io in
            % Bullock et al. 2012 eqn 10.
            % y is the whole range of 10 x e6 draws from the CDF (above).
            % (1/Nsamp)*sum(besseli(0,s.*y))= Eqn 10 in Bullock et al. 2012
            M=[1, 1, 1, (1/Nsamp)*sum(besseli(0,s.*y));...
                1,1,1,(1/Nsamp)*sum(besseli(0,s.*y));...
                1,1,1,1;...
                1,1,1,1];
            
            rho=A.*M; % = rho = A * M  in Eqn 11)
            rho_global = rho;
        end
    end
end
%This is end of program, below more functions to be called above
%%
%% 6. Randgen (CDF Walt x weibull)
function [CDF,j] = randgen(eta,beta)
% This function calculated the CDF for the wald model x weibull model as
% given by Bullock et al. Eqns 7 and 8.

% CDF = cumulative density function
% j, axis belonging to the CDF
% x is range to be intergrated over
% eps = a minimum number not being 0.

%set up the range of x values - note, this has to be both sufficiently large and
%fine.
global x
interval = 1;
x=eps:interval:10000;


%define the pdf - this is the intergral kernel for Bullock et al. Eqn 7 & 8.
global pdf
pdf=zeros(size(x));
tic 
% bullock et al.Eqn 8, integrating over the whole x-range
for i=1:length(x)
    pdf(i)=integral(x,i,eta,beta);   
end

%the output might need normalising, so do that anyway!
pdf=pdf/trapz(x,pdf);

%producing the pdf takes a long time, so save output to file
global species_name_global
global H_change_global
global F_change_global
global min_veg_height_change_global

% name an unique output file and save results of the PDF to hard drive
output_file = ['pdf','_',num2str(eta),'_',num2str(beta),'_',species_name_global,'_',num2str(H_change_global),...
    '_', num2str(F_change_global),'_', num2str(min_veg_height_change_global),'.mat'];
save(output_file,'x','pdf')

% Generate Cumulative DF from the pdf- need to make sure points are unique for inversion
[CDF,j] = unique(cumtrapz(x,pdf));

% name an unique output file and save results of the CDF to hard drive
output_file = ['cdf','_',num2str(eta),'_',num2str(beta),'_',species_name_global,'_',num2str(H_change_global),...
    '_', num2str(F_change_global),'_', num2str(min_veg_height_change_global),'.mat'];
save(output_file,'j','CDF')
end
%% 7. The Weibull x wald function pdf  
function I=integral(x,i,eta,beta)
%define the integral kernel
% min = 0 and max = 25 are the windspeeds range
% x is the whole range of distances (in meters)
% eta and beta are model parameters
% i is the sequence number within the x-range (so x(i) is a given distance)

I=quad(@integrand,0,25);
% This is Bullcok et al. Eqn 7, integrating per x(i) over the whole wind
% range
    function integ=integrand(windspeed)
        integ=weibull(windspeed).*wald(windspeed);
        function weifun=weibull(windspeed)
            % this is Eqn 1 in Bullock et al. 2012
            % windspeed = u in Bullock et al. 2012
            weifun=(beta/eta).*((windspeed./eta).^(beta-1)).*exp(-(windspeed./eta).^beta);
        end
        function walfun=wald(windspeed)
            %The PSF of the wald function follows 
            % This is Eqn 3 in Bullock et al. 2012 following Katul et al.
            % 2005; for Bullock et al. 2012: windspeed = u i; x(i) = distance
            walfun=((lambda(windspeed)./(2*pi*(x(i)^3))).^(1/2)).*exp(-(lambda(windspeed).*...
                (x(i)-mu(windspeed)).^2)./(2*(mu(windspeed).^2)*x(i)));
            % for lambda = Eqn 5 in Bullock et al. 2012, see below
            % for mu = eqn 4 in Bullock et al. 2012, see below
        end
    end
end
%% 8. lambda calculation feeding into the waldfunction calculation
% (Bullock et al. Eqns 4-6)
% u = windspeed
function lam=lambda(u)
% Sub Calculation of Wald parameters, Eqn 5 of Bullock et al. 2012, incl
% Eqn 6.

% set parameters
global H_global
global veg_height_global
Karman_constant = 0.4;
Kolmogorov_constant = 3.125;
zm = 10; % measurment height of wind.
Adoubleyou = 1.3;
d = 0.7* veg_height_global;
zo = 0.1* veg_height_global;

lam=zeros(size(u));
for j=1:length(u)
    % Bullock et al. Eqn 6. following Skarpaas adn Shea (2007)

    % Ustar according to Skarpaas & Shea 2007, A2 and A4 (equation not in Bullock et al. 2012)
    Ustar = Karman_constant*u(j)*((log((zm-d)/zo))^-1); % see Skarpaas & Shea 2007
    
    % this is the waldmodel calculations as intergral from the release
    % height to (d + zo) x vegetation height.
    windspeed_fun = @(z) (Ustar/Karman_constant)* (log((z-d)/zo));
    windspeed = ((1/H_global)*(quadl(windspeed_fun,(d+zo),H_global)));
    % Note here windpspeed is U(u) in Bullock et al. 2012
    
    % Eqn 5 from Bullock et al, with Sigma calculated according to Skarpaas
    % & Shea 2007, A2 and A4 (equation not in Bullock et al. 2012)
    sigma = (2*(Adoubleyou^2))* sqrt((Karman_constant*(H_global-d)*Ustar)/(Kolmogorov_constant*windspeed)); % see Skarpaas & Shea 2007
    lam(j) = (H_global/sigma)^2;
end
end
%% 9. Mu calculation feeding into the waldfunction calculation
% u = windspeed
function m=mu(u)
 % Sub Calculation of Wald parameters, Eqn 4 of Bullock et al. 2012, incl
% Eqn 6.

% set parameters
global H_global
global F_global
global veg_height_global
Karman_constant = 0.4;
zm = 10;
d = 0.7* veg_height_global;
zo = 0.1* veg_height_global;
m=zeros(size(u));
for j=1:length(u)
     % Bullock et al. Eqn 6. following Skarpaas adn Shea (2007)

    % Ustar according to Skarpaas & Shea 2007, A2 and A4 (equation not in Bullock et al. 2012)
    Ustar = Karman_constant*u(j)*((log((zm-d)/zo))^-1);% see Skarpaas & Shea 2007
    % this is the waldmodel calculations as intergral from the release
    % height to (d + zo) x vegetation height.
    windspeed_fun = @(z) (Ustar/Karman_constant)* (log((z-d)/zo));
    windspeed = ((1/H_global)*(quadl(windspeed_fun,(d+zo),H_global)));
    
    %Eqn 4 from Bullock et al,
    m(j) = (H_global*windspeed)/F_global;
end
end

%% 10. sensitivty and lambdas of the spread (A * M) matrix (Eqn 12)
function [sens_itivy_full elas_ticity_full lam_bda stable_stage_distribution_full reproductive_vector_full] =sensetivty(rho_global)
%find the left eigenvector
[matlvec,matval] = eig(rho_global');
[~, loc] = max(abs(diag(matval)));
left_vec=abs(matlvec(:,loc)); % = left eigenvector: reproductive vector (= avg number of offspring to expected from a plant in that class)

%find the right eigenvector
[matrvec,matval] = eig(rho_global);
[maxmatval, loc] = max(abs(diag(matval)));
right_vec=abs(matrvec(:,loc)); % = right eigenvector: stable stage distribution= percentage of plants in each class at equilibrium

%dom eigenvalue
lam_bda = maxmatval;

right_vec = right_vec./sum(right_vec); % rescales the rt eigenvector (uvec) to sum to 1
left_vec=left_vec./(left_vec'*right_vec); % rescales the lft eigenvector (vvec) to calc sens

sens_itivy = left_vec*right_vec'; % gets the sensitivity matrix
elas_ticity= sens_itivy.*(rho_global./lam_bda); %which scales each sij by each aij/lambda

stable_stage_distribution =right_vec; % already normalized to sum to 1
reproductive_vector=left_vec/left_vec(1);% normalize  so that newborns have rv of 1


sens_itivy_full = reshape(sens_itivy, 1, (length(sens_itivy))^2);
elas_ticity_full = reshape(elas_ticity, 1, (length(elas_ticity))^2);
stable_stage_distribution_full = reshape(stable_stage_distribution, 1, (length(stable_stage_distribution)));
reproductive_vector_full = reshape(reproductive_vector, 1, (length(reproductive_vector)));
end
