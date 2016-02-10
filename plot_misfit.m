function plot_misfit(DIP_DIR,DIP,RAKE)

fid = fopen('test_homo.out');
MISFIT_TAB = textscan(fid,'%f %f %f %f %f %f %f\r\n');
MISFIT = MISFIT_TAB{4};

STRIKE = DIP_DIR - 90; STRIKE(STRIKE < 0) = STRIKE(STRIKE < 0) + 360;
% Get P and T axes
TTREND = zeros(length(STRIKE),1);  TPLUNGE = zeros(length(STRIKE),1);
PTREND = zeros(length(STRIKE),1);  PPLUNGE = zeros(length(STRIKE),1);
for j = 1:length(STRIKE)
    [TTREND(j),TPLUNGE(j),PTREND(j),PPLUNGE(j)] = dc2axes(STRIKE(j),DIP(j),RAKE(j));
end

P = [PTREND,PPLUNGE]; T = [TTREND,TPLUNGE];
drawstereonet('Projection','wullf');
hold on;
[X,Y] = drawstereonet(PTREND, 90 - PPLUNGE,'Projection','wullf');
Hp = scatter3(X,Y,-0.1.*ones(size(X)),25,MISFIT,'filled');
set(Hp,'MarKerEdgeColor','k');
c = colorbar; ylabel(c,'Misfit Angle');
%-------------------------------------------------------------------------
function [TS,TD,PS,PD] = dc2axes(S1,D1,R1)

[S2, D2, R2] = define_second_plane(S1,D1,R1);
pure_strike_slip = false;
if abs(sin(R1*pi/180)) > 0.0001
    IM = floor(R1./abs(R1));
elseif abs(sin(R2*pi/180)) > 0.0001
    IM = floor(R2./abs(R2));
else
    pure_strike_slip = true;
end

if pure_strike_slip
    if cos(R1*pi/180) < 0.0
        PS = S1 + 45; TS = S1 - 45;
    else
        PS = S1 - 45; TS = S1 + 45;
    end
    I = PS >= 360;
    PS(I) = PS(I) - 360;
    I = TS >= 360;
    TS(I) = TS(I) - 360;
    
    PD = 0; TD = 0;
else
    cd1 =  cos(D1 * pi/180) * sqrt(2);
    sd1 =  sin(D1 * pi/180) * sqrt(2);
    cd2 =  cos(D2 * pi/180) * sqrt(2);
    sd2 =  sin(D2 * pi/180) * sqrt(2);
    cp1 = -cos(S1 * pi/180) * sd1;
    sp1 =  sin(S1 * pi/180) * sd1;
    cp2 = -cos(S2 * pi/180) * sd2;
    sp2 =  sin(S2 * pi/180) * sd2;
    
    amz = - (cd1 + cd2);
    amx = - (sp1 + sp2);
    amy = cp1 + cp2;
    
    dx = atan2(sqrt(amx * amx + amy * amy), amz) - pi/2;
    px = atan2(amy, - amx);
    if px < 0.0
        px = px + 2*pi;
    end
    
    amz = cd1 - cd2;
    amx = sp1 - sp2;
    amy = - cp1 + cp2;
    dy = atan2(sqrt(amx * amx + amy * amy), - abs(amz)) - pi/2;
    py = atan2(amy, - amx);
    if amz > 0.0
        py = py - pi;
    end
    
    if py < 0.0
        py = py + 2*pi;
    end
    
    if IM == 1
        PD = dy;
        PS = py;
        TD = dx;
        TS = px;
    else
        PD = dx;
        PS = px;
        TD = dy;
        TS = py;
    end
    
    TS = TS * 180/pi;
    TD = TD * 180/pi;
    PS = PS * 180/pi;
    PD = PD * 180/pi;
end
% ======================