%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE SOME OF THIS AS ITS JOANAS AND WE WANT OURS  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load Data
clc, clear, close all

Mars = importdata('Mars.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Commonly useful things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Labels
gravPotLabel = "Gravity Potential [m^2/s^2]";
gravVectorLabel = "Gravity Vector, R-component [m/s^2]";
gravGradTensorLabel = "Gravity Gradient Tensor, Trr-component [1/s^2]";
longLabel = "Longitude [deg]";
latLabel = "Latitude [deg]";

% Paths
projectPath = pwd;
matlabPath = fullfile(projectPath, "PS2", "Ass2", "MATLAB_files");

% Booleans (often changing
titling = false;
saving = true;
showing = false;

% Simulation Parameters
lonLim = [-180 180 1]; % Plot all longs, with 1 degree precision
latLim = [-90 90 1]; % Plot all lats, with 1 degree precision
% SHbounds = [4 10]; % take out J0,1,2,3. Also truncate to 75
SHbounds = [4 75]; % take out J0,1,2,3. Also truncate to 75
Model = struct();
Model.Re = 3389500; % Mars radius, m
m = 6.4185e23; % Mars mass, kg
G = 6.67406e-11; % m3/kg/s^2
Model.GM = G*m ;%4.282837e13; %m^3/s^2

if showing == false
    set(gcf,'visible','off')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial set of things to calculate and plot, Q1a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% New parameters (ie height)
height = 200; % Height "ASL" in m

% Run simulation
data = model_SH_synthesis_PS(lonLim,latLim,height,SHbounds,Mars,Model);

%%%%%%% Do part 1a plot %%%%%%%%%%%%%%%

% Plot gravity potential
doPlot(data.pot, longLabel, latLabel, gravPotLabel, "Gravity Potential (200m)", titling, "gravPot_200m", matlabPath, saving, showing)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q1b plots - all things for different heights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% New parameters (ie height)
height = 0; % Height "ASL" in m

% Run simulation
data = model_SH_synthesis_PS(lonLim,latLim,height,SHbounds,Mars,Model);

%%%%%%% Do Q1b plots %%%%%%%%%%%%%%%

% Plot gravity potential
doPlot(data.pot, longLabel, latLabel, gravPotLabel, "Gravity Potential (ASL)", titling, "gravPot_ASL", matlabPath, saving, showing)

%Plot R-component
doPlot(data.vec.R, longLabel, latLabel, gravVectorLabel, "Gravity Vector, R-component (ASL)", titling, "gravVec_ASL", matlabPath, saving, showing)

%Plot Trr-component
doPlot(data.ten.Trr, longLabel, latLabel, gravGradTensorLabel, "Gravity Gradient Tensor, Trr-component (ASL)", titling, "gravTens_ASL", matlabPath, saving, showing)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q1c plots - all things for different heights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% New parameters (ie height)
heights = [0, 200000, 400000, 1000000]; % Height "ASL" in m

for height = heights
    % Run simulation
    data = model_SH_synthesis_PS(lonLim,latLim,height,SHbounds,Mars,Model);
    
    % Plot figure of R-vector component
    titleName = sprintf("Gravity Vector, R-component (%s km)", num2str(height/1000));
    saveName = sprintf("gravVec_%s", num2str(height/1000));
    doPlot(data.vec.R, longLabel, latLabel, gravVectorLabel, titleName, titling, saveName, matlabPath, saving, showing)
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doPlot(dataToPlot, xLabel, yLabel, colorBarLabel, Title, titling, saveName, savePath, saving, showing)
    if showing == false
        set(0,'DefaultFigureVisible','off')
    else
        set(0,'DefaultFigureVisible','on')
    end
    fig=figure();
    
    imagesc(dataToPlot);
    colormap;
    colorbar;
    
    if titling 
        title(Title); 
    end
    
    xlabel(xLabel);
    ylabel(yLabel);
    c = colorbar;
	c.Label.String = colorBarLabel;
    if saving
        fullSavePath = fullfile(savePath, saveName);
        saveas(fig, fullSavePath, "png");   
    end
end


function [data] = model_SH_synthesis_PS(lonLim,latLim,height,SHbounds,V,Model)
% 
% This function is responsible for the SH synthesis of a given model.
%
% input:
%           - lonLim: longitude limits [min-longitude max-longitude resolution] in degree
%           - latLim: latitude limits  [min-longitude max-longitude resolution] in degree
%           - height: height of computation surface in meters [scalar]
%           - SHbounds: domain of order and degree SH coeff. [nmin nmax]
%           - V: full set of SHcoefficients constructed by model_SH_analysis.m
%                  V =     [0 0 C00 S00]
%                          [1 0 C10 S10]
%                          [. . .   .  ]
%                          [. . .   .  ]
%                          [. . .   .  ]
%                          [1 1 C11 S11]
%                          [1 2 C12 S12]
%                          [1 3 C13 S13]
%                          [. . .   .  ]
%                          [. . .   .  ]
%                          [. . .   .  ]
%                          [2 1 C21 S21]
%                          [2 2 C22 S22]
%                          [. . .   .  ]
%                          [. . .   .  ]
%                          [nmax mmax Cnm Snm]
%
%           - Model: model structure constructed by inputModel.m contains Re and GM values
%                   Model = struct()
%                   Model.Re = radius of planet in meters
%                   Model.GM = gravitational constant of planet in m^3/s^2
%
% output:   - data structure with gravity potential [m^2/s^2], gravity vector [m/s^2 or 1e5*mGal] 
%             and gravity gradient tensor [1/s^2 or 1e9*Eotvos]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct the correct input variables

lon = lonLim(1):lonLim(3):lonLim(2);
lat = latLim(1):latLim(3):latLim(2);
    
Lon = repmat(lon,length(lat),1);
Lat = repmat(lat',1,length(lon));

r = Model.Re + height;

%% get potential field and others

if isscalar(r)
    [data] = gravityModule(Lat,Lon,r,SHbounds,V,Model.Re,Model.GM);
else
    [data] = gravityModule_full(Lat,Lon,r,SHbounds,V,Model.Re,Model.GM);
end

% saving the input values
data.latLim =    latLim; 
data.lonLim =    lonLim; 
data.height =    height;
data.SHbounds =  SHbounds;
data.Model = Model;

end

function [data] = gravityModule(Lat,Lon,r,SHbounds,V,Re,GM)
%
% Calculating the the gravity field for the EIGENGL04C model
%
% input: Lat: latitude in degree [matrix]
%        Lon: longitude in degree [marix]
%        r: radial distance of computation surface (r=Re+h) [scalar]
%        SHbounds: [minSH max SH] example [0 5] or [2 7]
%        V : SH coefficients [same size as setleg files]
%        V format:  degree; order; Cnm; Snm
%        Re: radius of the Earth [meters]
%        GM: gravitational parameter of Earth 
%
% output: data: data structure of the gravity field
%
% programs used by the routine: 
%       - bsxfun
%       - getLegendre
%           - visu2plm_ww
%       
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the legendre polynomials

[setLeg] = getLegendre(SHbounds,Lat);

% Create the correct domain in the coefficients
nmin = SHbounds(1);
nmax = SHbounds(2);
V(V(:,1)>nmax|V(:,1)<nmin,:) = [];

% check is setLeg and V have the same length

if (size(setLeg.PMatrix,1)~=size(V,1))
    error('Length of Legendre polynomials should have the same length as the SH-coefficients')
end

% set data matrices

Pot   = zeros(size(Lon));
VecR  = zeros(size(Lon));
VecT  = zeros(size(Lon));
VecL  = zeros(size(Lon));
TenRR = zeros(size(Lon));
TenTT = zeros(size(Lon));
TenLL = zeros(size(Lon));
TenRT = zeros(size(Lon));
TenRL = zeros(size(Lon));
TenTL = zeros(size(Lon));

Trans1 = zeros(size(Lon));

% Set some variables

fac  = GM/Re;
Sfac = GM/Re/Re;
Tfac = GM/Re/Re/Re;

lambda = deg2rad(Lon(1,:));
theta  = deg2rad(Lat(:,1))';
sint   = sin(pi/2-theta);
cost   = cos(pi/2-theta);

% start longitude for-loop, rangnge over all lon grid steps
lon = Lon(1,:);

for l = 1:length(lon)
    
    % Range related coefficients
    fac1 =  (Re./r).^(V(:,1)+1);
    Sfac1 = (Re./r).^(V(:,1)+2);
    Tfac1 = (Re./r).^(V(:,1)+3);

    %%%%%%%%%%%%%%%% Potential %%%%%%%%%%%%%%%%%%%%%%
    % Sumation of the Grace data potential

    pot = sum(bsxfun(@times,setLeg.PMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                           + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                           .*fac1)),1);    

    %%%%%%%%%%%%%% Gradient vector %%%%%%%%%%%%%%%%%%
    % For the GRACE data files

    vecR =  sum(bsxfun(@times,setLeg.PMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                             + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                          .*Sfac1.*(V(:,1)+1))),1);

    vecT =  sum(bsxfun(@times,setLeg.DMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                             + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                             .*Sfac1)),1);

    vecL =  sum(bsxfun(@times,setLeg.PMatrix,((-(V(:,2)).*sin((V(:,2))*lambda(1,l)).*V(:,3)...
                                               +(V(:,2)).*cos((V(:,2))*lambda(1,l)).*V(:,4))...
                                               .*Sfac1)),1);

    %%%%%%%%%%%%%% Gradient tensor %%%%%%%%%%%%%%%%%%
    % Summation of the individual gradient elements

    tenrr =   sum(bsxfun(@times,setLeg.PMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                               + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                            .*Tfac1.*(V(:,1)+1).*(V(:,1)+2))),1);

    tentt =   sum(bsxfun(@times,setLeg.SMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                               + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                               .*Tfac1)),1);                

    tenll =   sum(bsxfun(@times,setLeg.PMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                               + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                           .*Tfac1.*-(V(:,2)).^2)),1);

    tenrt =   sum(bsxfun(@times,setLeg.DMatrix,((cos((V(:,2))*lambda(1,l)).*V(:,3)...
                                               + sin((V(:,2))*lambda(1,l)).*V(:,4))...
                                               .*Tfac1.*(V(:,1)+1))),1);

    tenrl =   sum(bsxfun(@times,setLeg.PMatrix,((-(V(:,2)).*sin((V(:,2))*lambda(1,l)).*V(:,3)...
                                               +  (V(:,2)).*cos((V(:,2))*lambda(1,l)).*V(:,4))...
                                         .*Tfac1.*(V(:,1)+1))),1);

    tentl =   sum(bsxfun(@times,setLeg.DMatrix,((-(V(:,2)).*sin((V(:,2))*lambda(1,l)).*V(:,3)...
                                               +  (V(:,2)).*cos((V(:,2))*lambda(1,l)).*V(:,4))...
                                                 .*Tfac1)),1);   

    % Finalizing all the gravity data matrices

    Pot(:,l) = (fac*pot)';

    VecR(:,l) =  (Sfac*vecR)';
    VecT(:,l) = (-Sfac*vecT)';                    % added a minus sign
    VecL(:,l) =  (Sfac*vecL./sint)';

    TenRR(:,l) =  (Tfac*tenrr)';
    TenTT(:,l) =  (Tfac*tentt)';
    TenLL(:,l) =  (Tfac*tenll./sint./sint)';
    TenRT(:,l) =  (Tfac*tenrt)';                  % added a minus sign
    TenRL(:,l) = (-Tfac*tenrl./sint)';
    TenTL(:,l) = (-Tfac*tentl./sint)';            % added a minus sign
    
    % Transformation factors

    Trans1(:,l) = (cost./sint)';
    
end

% Construct data structure
data = struct();

% Gravity potential of topography
data.pot     = Pot;
data.grd.lon = Lon;
data.grd.lat = Lat;
data.grd.r   = r;

% Gravity vector of topography
data.vec.R = VecR;
data.vec.T = VecT;
data.vec.L = VecL;
data.vec.X = -VecT;
data.vec.Y = -VecL;
data.vec.Z = VecR;

% Gravity gradient tensor
data.ten.Trr = TenRR;
data.ten.Ttt = TenTT;
data.ten.Tll = TenLL;
data.ten.Trt = TenRT;
data.ten.Trl = TenRL;
data.ten.Ttl = TenTL;
end

function [setLeg] = getLegendre(SHbounds,Lat)
%
% This function computes the Legendre polynomial matrices used in the SH
% software package. Als the first and second derivative are computed for
% gradient computations
%
% input:  SHbounds = vector [nmin nmax], boundaries of the SH degree and order
%         Lat: latitude vector [degree] of the geogrid
%
% output: PMatrix: Legendre polynomials
%         DMatrix: first derivative
%         SMatrix: second derivative
%
% software routines used: visu2plm_ww
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmin = SHbounds(1);
nmax = SHbounds(2);

CoLat = (90 - Lat(:,1))';
axl = (nmax)*(nmax+1)/2+nmax+1 - ((nmin-1)*(nmin)/2+nmin);

PnmF = zeros(axl,length(CoLat));
DnmF = zeros(axl,length(CoLat));
SnmF = zeros(axl,length(CoLat));

% Generate the Legendre function matrix for a particular theta
e = 1;
k = 1;
for s = 0:nmax

    b = e;
    
    % using wouters software to calculate polynomials
    [P,DP,SP] = visu2plm_ww(nmin:nmax,s,CoLat);
    
    % Get the correct set (no zeros)
    newP = P';
    newDP = DP';
    newSP = SP';
    
    if s>nmin
        newP(1:k,:) = [];
        newDP(1:k,:) = [];
        newSP(1:k,:) = [];
        k = k + 1;
    end
    
    % file in temperal matrix
    PnmF(b:b + size(newP,1) - 1,:) =  newP;    
    DnmF(b:b + size(newDP,1) - 1,:) = newDP;
    SnmF(b:b + size(newSP,1) - 1,:) = newSP;
    
    e = b + size(newP,1);

end

% Making legendre setting
setLeg = struct();
setLeg.PMatrix = PnmF;
setLeg.DMatrix = DnmF;
setLeg.SMatrix = SnmF;
end

function [p, dp, ddp] = visu2plm_ww(l,m,th)

% PLM Fully normalized associated Legendre functions for a selected order M
%
% HOW p       = plm(l,th)			- assumes M=0
%     p       = plm(l,m,th)
%     [p,dp]  = plm(l,m,th)
%     [p,ddp] = plm(l,m,th)
%
%
% IN  l  - degree (vector). Integer, but not necessarily monotonic.
%          For l < m a vector of zeros will be returned.
%     m  - order (scalar). If absent, m=0 is assumed.
%     th - co-latitude [deg] (vector)
% OUT p  - Matrix with Legendre functions. The matrix has length(TH) rows
%          and length(L) columns, unless L or TH is scalar. Then the output
%          vector follows the shape of respectively L or TH. 
%    dp  - Matrix with first derivative of Legendre functions. The matrix 
%          has length(TH) rows and length(L) columns, unless L or TH is 
%          scalar. Then the output vector follows the shape of respectively 
%          L or TH. 
%    ddp - Matrix with second (colatitude) derivative of Legendre functions. The matrix 
%          has length(TH) rows and length(L) columns, unless L or TH is 
%          scalar. Then the output vector follows the shape of respectively 
%          L or TH. 
% 
% See also LEGPOL, YLM, IPLM

%-----------------------------------------------------------------------------
% Nico Sneeuw, IAPG, TU-Munich                                       08/08/94
%-----------------------------------------------------------------------------
% Uses none
%-----------------------------------------------------------------------------
% Revision history:
%  - NS09/06/97: help text brushed up
%  - NS13/07/98: Pmm non-recursive anymore 
%  - NS0299:     further help text brush-up
%  - MW13/08/04: extension for first derivative
%  - MW24/11/04: speed up calculation
%  - Wouter van der Wal, June 15, 2006: extension for second derivative
%    using recursion formulas of Novak and Grafarend (2006)
%    tested with analytical second derivatives (coefficient pairs 1,1; 2,1; 2,2; 3,1; 3,3; 4,1; 4,2; 4,4) 
%    see: 2ndDerivativeLegendre.pdf
%-----------------------------------------------------------------------------


% Some input checking.
if nargin == 2
   th = m;
   m  = 0;
end
if min(size(l)) ~= 1,  error('Degree l must be vector (or scalar)'), end
if any(rem(l,1) ~= 0), error('Vector l contains non-integers.'), end
if max(size(m)) ~= 1,  error('Order m must be scalar.'), end
if rem(m,1) ~=0,       error('Order m must be integer.'), end


% Preliminaries.
[lrow,lcol] = size(l);
[trow,tcol] = size(th);
lmax = max(l);
if lmax < m, error('Largest degree still smaller than order m.'), end
n    = length(th);				% number of latitudes
t    = th(:)*pi/180;
x    = cos(t);
y    = sin(t);
lvec = l(:)';					% l can be used now as running index.

% Recursive computation of the temporary matrix ptmp, containing the Legendre
% functions in its columns, with progressing degree l. The last column of
% ptmp will contain zeros, which is useful for assignments when l < m.
ptmp  = zeros(n,lmax-m+2);
if nargout >= 2, dptmp = zeros(n,lmax-m+2); end
if nargout == 3, ddptmp = zeros(n,lmax-m+2); end

%--------------------------------------------------------------------
% sectorial recursion: PM (non-recursive, though)
%--------------------------------------------------------------------
%WW: produces sqrt( (2n+1)/2n )
% Novak and Grafarend (2006), eq 64
if m == 0
   fac = 1;
else
   mm  = 2*(1:m);
   fac = sqrt(2*prod((mm+1)./mm));   % extra sqrt(2) because summation not over negative orders
end

ptmp(:,1) = fac*y.^m;                                      % The 1st column of ptmp.
% TEST: for m = 1, theta = 10, fac = sqrt(6/2)
% ptmp(1,1) = sqrt(6/2)*sind(10)*1
if nargout >= 2
    dptmp(:,1) = m*fac*(y.^(m-1).*x);
end     % The 1st column of dptmp.
% recursion is beta_n,n * beta_n-1,n-1 * ...
% sin(theta) * sin(theta) * ...
% * cos(theta) * P_0,0 (which is 1, dP_0,0 is 0)
% note that the term beta_n,n*cos(theta)*P_n-1,n-1 of Novak and
% Grafarend(2006) equation 72 is taken care of by the m.
if nargout == 3
    ddptmp(:,1) = -m*fac*(y.^m) + m*(m-1)*fac*(y.^(m-2).*x.^2);
end     % The 1st column of ddptmp.


%--------------------------------------------------------------------
% l-recursion: P
%--------------------------------------------------------------------
for l = m+1:lmax
   col   = l - m + 1;			% points to the next column of ptmp
   root1 = sqrt( (2*l+1)*(2*l-1)/((l-m)*(l+m)) ) ;                      % beta_n,m (65) 
   root2 = sqrt( (2*l+1)*(l+m-1)*(l-m-1) / ( (2*l-3)*(l-m)*(l+m) ) );   % beta_n,m (65) * gamma_n,m (66)

   % recursion
   if l == m+1
       ptmp(:,col) = root1 *x.*ptmp(:,col-1);
   else
       ptmp(:,col) = root1 *x.*ptmp(:,col-1) - root2 *ptmp(:,col-2);
   end
       
   if nargout >= 2, 
       if l == m+1
           dptmp(:,col) = root1 *(x.*dptmp(:,col-1)-y.*ptmp(:,col-1)); 
       else
           dptmp(:,col) = root1 *(x.*dptmp(:,col-1)-y.*ptmp(:,col-1)) - root2 *dptmp(:,col-2); 
       end
   end

   if nargout == 3,
       if l == m+1
           ddptmp(:,col) = root1 *(-x.*ptmp(:,col-1) -2*y.*dptmp(:,col-1) + x.*ddptmp(:,col-1) );
       else
           ddptmp(:,col) = root1 *(-x.*ptmp(:,col-1) -2*y.*dptmp(:,col-1) + x.*ddptmp(:,col-1) ) - root2*ddptmp(:,col-2); 
       end
   end

end


% The Legendre functions have been computed. What remains to be done, is to
% extract the proper columns from ptmp, corresponding to the vector lvec. 
% If l or theta is scalar the output matrix p reduces to a vector. It should
% have the shape of respectively theta or l in that case.

p          = zeros(n,length(lvec));         % size declaration.
lind       = find(lvec < m);			    % index into l < m
pcol       = lvec - m + 1;			        % index into columns of ptmp
pcol(lind) = (lmax-m+2)*ones(size(lind));	% Now l < m points to last col.
p          = ptmp(:,pcol);			        % proper column extraction 
if nargout >= 2, dp = dptmp(:,pcol); end    % proper column extraction
if nargout == 3, ddp = ddptmp(:,pcol); end  % proper column extraction

if max(size(lvec))==1  & min(size(th))==1 & (trow == 1), 
    p = p'; 
    if nargout >= 2, dp = dp'; end
    if nargout == 3, ddp = ddp'; end
end
if max(size(th))==1 & min(size(lvec))==1  & (lcol == 1), 
    p = p'; 
    if nargout >= 2, dp = dp'; end
    if nargout == 3, ddp = ddp'; end
end
end