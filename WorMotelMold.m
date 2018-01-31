%
% 3D printed negative mold for WorMotel 
% Matt Churgin, Chris Fang-Yen et al.
% from Churgin et al. eLife 2017 paper
% clear all

% xwells = 24;
% ywells = 16;
docornercuts = 0;

xwells = 20;
ywells = 12;

% xwells = 8;
% ywells = 6;
 
% xwells = 4;
% ywells = 3;

% COMPACT BORDERS

textstring = '130305-240WELL-ARMS';  
% no spaces.  should be similar to filename of .m file

steps_per_mm = 10;

baseheight = 4*steps_per_mm;  % height of mold at top of wells
welldepth = 3*steps_per_mm;

moatwidth = .5*steps_per_mm;  % must be odd in # pixels
moatdepth = 2*steps_per_mm;

unitcell = 4.5*steps_per_mm;

% cornerdiagonal = 

numborders = 6;  % number of borders around wells from inside to out
borderx = zeros(numborders,1); % x width of borders
bordery = zeros(numborders,1); % y width of borders
borderheight= zeros(numborders,1); % border heights
outerwallheight = baseheight + 10 * steps_per_mm;

borderx(1) = .5*steps_per_mm;   % extra outer moat width
bordery(1) = .5*steps_per_mm;   % 
borderheight(1) =  baseheight + moatdepth;

borderx(2) = 2*steps_per_mm;   % baseheight border 
bordery(2) = 2*steps_per_mm;   %  
borderheight(2) =  baseheight;

borderx(3) = 2*steps_per_mm;   % hydration border 
bordery(3) = 2*steps_per_mm;   % 
borderheight(3) =  1*steps_per_mm;

borderx(4) = 0*steps_per_mm;   % notch for cutting off (NOT USED)
bordery(4) = 0*steps_per_mm;   % 
borderheight(4) =  baseheight  + 0;

borderx(5) = 8.8*steps_per_mm;   % arm / outer wall 
bordery(5) = 5.6*steps_per_mm;   % 
borderheight(5) =  outerwallheight-0.1;

borderx(6) = 2*steps_per_mm;   % outer wall
bordery(6) = 2*steps_per_mm;   % outer wall
borderheight(6) = outerwallheight;

% % approximate humidifying volume in mL
% v1 = 2 * 7 * 85 * 3 /1000
% 
% % approximate media volume
% v2= 384 * 20 / 1000
% 
% % approximate moat volume
% v3= 1.5 * 0.5 * (16*116 + 24 * 75) / 1000
% 
% % total liquid volume
% v1+v2+v3
% 
% borderx(4) = 2*steps_per_mm;   % humidification moat wall (also device outer wall)
% borderheight(4) =  baseheight;


alignmentmarkdepth = 1*steps_per_mm;

% bordery = borderx;  % symmetric

bordersumx = sum(borderx);  % total border width
bordersumy = sum(bordery);

R1 = 0.2*unitcell/2; % radius of inner flat portion
R2 = .8*unitcell/2; % radius of outer flat portion
x = [0:1:unitcell-1]; y =  [0:1:unitcell-1];
[X,Y] = meshgrid(x,y);
R = sqrt((X-.5*unitcell).^2 + (Y-.5*unitcell).^2);
Z = 1-1*(1-cos(pi*(R-R1)/(R2-R1)))/2;  % Function describing well shape
Z(R<R1) = 1;
Z(R>R2) = 0;

figure(1); 
surf(X,Y,Z);  % UNIT HEIGHT PROFILE
axis vis3d; 

% MAKE WELL ARRAY
Zarray1 = baseheight + welldepth * repmat(Z, ywells, xwells);
Zarray = baseheight*ones(2*bordersumy+ywells*unitcell, 2*bordersumx+xwells*unitcell);
Zarray((bordersumy+1):(bordersumy+ywells*unitcell),(bordersumx+1):(bordersumx+xwells*unitcell)) = Zarray1;

[Xarray,Yarray] = meshgrid([0:size(Zarray,2)-1], [0:size(Zarray,1)-1]);

% MAKE BORDERS
for j=1:numborders  % work from inside out
    Zarray(Yarray>(size(Zarray,1)-sum(bordery(j:end)))) = borderheight(j);
    Zarray(Yarray<sum(bordery(j:end))) = borderheight(j);
    Zarray(Xarray<sum(borderx(j:end))) = borderheight(j);
    Zarray(Xarray>(size(Zarray,2)-sum(borderx(j:end)))) = borderheight(j);
end

% Make moats and lines
shift1=1;
for x1 = 0:xwells
    Zarray(bordersumy:end-bordersumy+shift1,(bordersumx+x1*unitcell-(moatwidth-1)/2+shift1):(bordersumx+x1*unitcell+(moatwidth-1)/2)+shift1) = baseheight + moatdepth;
    
%     Zarray(sum(bordery(4:end)):sum(bordery(3:end)),(bordersumx+x1*unitcell-(moatwidth-1)/2+shift1):(bordersumx+x1*unitcell+(moatwidth-1)/2)+shift1) = borderheight(3) + alignmentmarkdepth;
%     Zarray(end-sum(bordery(3:end))+2:end-sum(bordery(4:end))+1,(bordersumx+x1*unitcell-(moatwidth-1)/2+shift1):(bordersumx+x1*unitcell+(moatwidth-1)/2)+shift1) = borderheight(3) + alignmentmarkdepth;

end
for y1 = 0:ywells
    Zarray((bordersumy+y1*unitcell-(moatwidth-1)/2+shift1):(bordersumy+y1*unitcell+(moatwidth-1)/2+shift1), bordersumx+shift1:end-bordersumx+shift1) = baseheight + moatdepth;
%     Zarray((bordersumy+y1*unitcell-(moatwidth-1)/2+shift1):(bordersumy+y1*unitcell+(moatwidth-1)/2+shift1), sum(borderx(4:end)):sum(borderx(3:end))) = borderheight(3) + alignmentmarkdepth;
%     Zarray((bordersumy+y1*unitcell-(moatwidth-1)/2+shift1):(bordersumy+y1*unitcell+(moatwidth-1)/2+shift1), end-sum(borderx(3:end))+2:end-sum(borderx(4:end))+1) = borderheight(3) + alignmentmarkdepth;
end



dotext =1;
donumbersletters = 1;
textbordernumber = 2;
textheight = .8 * steps_per_mm;

if dotext
    % add text at border
    % textstartx = size(Zarray,2)-sum(borderx((textbordernumber):end)) + 6; 
    textstartx = sum(borderx((textbordernumber+1):end)) + round(borderx(textbordernumber)/6); 
    %- round(.25 * borderx(textbordernumber)); ;
%     textstarty = size(Zarray,1)-sum(bordery((textbordernumber+1):end)) -3 ;
    textstarty = sum(bordery((textbordernumber+1):end)) +round(bordery(textbordernumber));
    % textsize = round(.5 * borderx(textbordernumber));
    textimg = imrotate(text2img(textstring, 16), -90);
%     textimg = textimg(:,end:-1:1);  
    Zarray((textstarty:textstarty+size(textimg,1)-1), (textstartx):(textstartx+size(textimg,2)-1)) = ...
       Zarray((textstarty:textstarty+size(textimg,1)-1), (textstartx):(textstartx+size(textimg,2)-1))+ textheight*textimg;
end

if donumbersletters
    % add numbers and lines along X
    for i=1:xwells
         numberimg = text2img(num2str(i),12); 
         
    %      numberimg =      numberimg(:,end:-1:1);
         xpos = round(sum(borderx)+ ((i-.5)*unitcell)-0.5*size(numberimg,2));
         ypos = sum(bordery((textbordernumber+1):end))+ round(bordery(textbordernumber)/4);
         Zarray(ypos:(ypos+size(numberimg,1)-1), xpos:(xpos+size(numberimg,2)-1)) =   Zarray(ypos:(ypos+size(numberimg,1)-1), xpos:(xpos+size(numberimg,2)-1))  + numberimg * textheight;

    %      Zarray(ypos, xpos) =   Zarray(ypos, xpos)+100;
    %      Zarray(

    end

    % add letters along Y
    for i=1:ywells
         numberimg = text2img(char(64+i),12); 
    %      numberimg =      numberimg(:,end:-1:1); % mirror image
         xpos = size(Zarray,2) - (sum(borderx((textbordernumber+1):end))+ round(borderx(textbordernumber)/1.7));
         ypos = round( sum(bordery)+ ((i-.5)*unitcell)-0.5*size(numberimg,1));
         Zarray(ypos:(ypos+size(numberimg,1)-1), xpos:(xpos+size(numberimg,2)-1)) =   Zarray(ypos:(ypos+size(numberimg,1)-1), xpos:(xpos+size(numberimg,2)-1))  + numberimg * textheight;
    end
end

doarms = 1;
armwidth = 5 * steps_per_mm;
armheight = baseheight;

if doarms
    armmask = zeros(size(Zarray));
    armmask(round(1/4*size(Zarray,1)-armwidth/2):round(1/4*size(Zarray,1)+armwidth/2), :)=1;
    armmask(round(3/4*size(Zarray,1)-armwidth/2):round(3/4*size(Zarray,1)+armwidth/2), :)=1;
    armmask(:,round(1/4*size(Zarray,2)-armwidth/2):round(1/4*size(Zarray,2)+armwidth/2))=1;
    armmask(:,round(3/4*size(Zarray,2)-armwidth/2):round(3/4*size(Zarray,2)+armwidth/2))=1;
    
    armregion = boolean((Zarray == outerwallheight-0.1).*armmask);
    Zarray(armregion) = armheight;
    
    
    mask1 = (Zarray < outerwallheight-0.1);
    se = strel('square',30);
    mask1d = imdilate(mask1,se);
    Zarray(mask1d==0)= 1*steps_per_mm;
end

if docornercuts      % MAKE CORNER CUTS
    for x=86:120
        y = 785-103-8+x;
        Zarray(y,x+3:x+10) = baseheight;
    end
    
    for x=43:125
        y = 680+x;
        Zarray(y,x+3:x+12) = borderheight(4);
    end
    
    % xi = [54.  122.   54.   54.]';
    % yi = [732.  802  802.  732.]';
    xi = [51 129 48 48.]';
    yi = [729 807 806 729]';
    xi = [50.1921   50.3825  129.5852   41.2437   42.0053   50.1921]';
    yi = [721.4640  727.9373  807.5208  807.5208  721.2736  721.4640]';
    BW=roipoly(Zarray,xi,yi);
    Zarray(BW) = borderheight(5);

%     Zarray(760:764,97:103) = baseheight;
%     Zarray(765:766,102:103) = baseheight;
%     Zarray(786:790,123:126) = baseheight;

    % Copy on RHS

    LHS = Zarray(711:810, 41:140);

    Zarray(711:711+size(LHS,1)-1, 1142:1142+size(LHS,2)-1) = LHS(:,end:-1:1);
end


% Border of zeros around entire object (eats into outermost border)
Zarray(1,:) = 0;
Zarray(end,:) = 0;
Zarray(:,1) = 0;
Zarray(:,end) = 0;

figure(3);clf;
subplot(211);
plot(Zarray(round(bordersumy+unitcell/2),:));
title('X section');
axis equal
subplot(212);
plot(Zarray(:,round(bordersumx+unitcell/2)));
title('Y section');
axis equal
% 
% figure(4); clf;
% % subplot(211);
% % surf(Xarray,Yarray,Zarray);
% surf(Xarray, Yarray, Zarray);
% axis vis3d; axis equal;
% shading interp;
% colormap gray;
% % imagesc(Zarray);
% axis image
% 
% figure(5); clf;
% % surf(Xarray,Yarray,Zarray);
% % subplot(212);
% surf(-Xarray, -Yarray, -Zarray);
% axis vis3d; axis equal;
% shading interp;
% colormap gray;
% title('Mold');

figure(6); clf;
imagesc(-Zarray);
axis image;
title('-Mold');
% colormap gray;

Zrange = max(max(Zarray)) - min(min(Zarray));
Xrange = max(max(Xarray)) - min(min(Xarray));
Yrange = max(max(Yarray)) - min(min(Yarray));

disp('Dimensions of final object (X, Y, Z) (mm) = ');
disp([Xrange Yrange Zrange]/steps_per_mm)

maxrange = max(max(Xrange,Yrange), Zrange)


% calculate size of molded object

Zmax = max(max(Zarray));
crosssectionx = Zarray(round(size(Zarray,1)/4), :);
crosssectiony = Zarray(:,round(size(Zarray,2)/4))';

idx1 = find(crosssectionx(1:100) == Zmax, 1, 'last');
idx2 = length(crosssectionx) - find(crosssectionx(end-1:-1:end-100) == Zmax, 1, 'last');
objectsize_x = idx2 - idx1;

idx1 = find(crosssectiony(1:100) == Zmax, 1, 'last');
idx2 = length(crosssectiony) - find(crosssectiony(end-1:-1:end-100) == Zmax, 1, 'last');
objectsize_y = idx2 - idx1;

disp('Dimensions of molded object (X, Y, Z) (mm) = ');
disp([objectsize_x objectsize_y 0]/steps_per_mm);

disp('Microplate interior dimensions  = 116.8 74.5')

disp('Shortfall = ');
disp([objectsize_x objectsize_y 0]/steps_per_mm -[116.8 74.5 0 ]);

%%
figure(7);
% convert to a solid
[X_out, Y_out, Z_out] = surf2solid(Xarray,Yarray,Zarray, maxrange/steps_per_mm, 1);
% print('-dpng', 'Plots/surf2solid_ex.png');

axis vis3d;
axis equal;

tic
% output a STL version of the solid
filename = ['\\psf\Home\Dropbox\PROJECTS\Microwell384\' textstring '.stl'];
mode = 'binary';
surf2stl(filename,X_out,Y_out,Z_out,mode);
toc
