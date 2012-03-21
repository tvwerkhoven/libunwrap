%
% Unwraps a phase using obsolete floodfilling unwrapping routine.
% OBSOLETE!
%
function unwrapped = unwrap_flood(phsm, pup, toplot)

borderlimitfrac=0.75;
pupdim = size(pup,1);
itcount = 1e3;
  
% Create a bit tighter pupil with the same support
v = (0.5:pupdim) - pupdim/2;
[m1 m2] = meshgrid(v);
rm = sqrt(m1.^2 + m2.^2.);
smpup = zeros(pupdim,pupdim);
smpup(rm < pupdim/2-2) = 1;

wrapped = phsm;
donemask = 0*pup;    
gloco=0;

po1=-1;
po2=-1;
mispo=0;
totborders = -1;
alltotborders = zeros(itcount, 1);

for i1=1:itcount
  prevtotborders = totborders;  
  borders = wrapped(1:127,1:127) - wrapped(2:128,2:128);
  borders(abs(borders) < borderlimitfrac*2*pi) = 0;
  borders = borders.*smpup(1:127,1:127);  
  totborders = sum(borders(:) ~= 0);
  alltotborders(i1) = totborders;
  
  %borders(donemask(1:127,1:127) == 1) = 0;
  if max(abs(borders(:))) == 0
    fprintf('no more borders\n');
    break;
  end
  
  minin = find(borders ~= 0, 1);
  po1last=po1; po2last=po2;
  [po1 po2] = ind2sub(size(borders), minin);
  donemask = 0*pup;
  
  if exist('toplot', 'var')
    fprintf('start: %d, %d (totborders %d)\n', po1, po2, totborders)
  end
  %unwrapflood(po1,po2, 0);
  
  % Start the code
  if mod(i1,10) == 0
    usedonemask = 0;
  else
    usedonemask = 1;
  end
  wrapped = unwrap_oct(wrapped, pup, po1-1, po2-1, usedonemask, borderlimitfrac);

  
  if po1last==po1 && po2last==po2
    mispo = mispo+1;
  else
    mispo = 0;
  end    
  if po1last==po1 && po2last==po2 && mispo > 20 && prevtotborders == totborders
    fprintf('Internal error: algo doesnt work!!\n');
    break
  end
  
  if exist('toplot', 'var')
    subplot(2,2,1)
    imagesc(phsm)
    subplot(2,2,2)
    imagesc(wrapped)
    subplot(2,2,4)
    imagesc(borders)
    drawnow      
    %pause
  end
end

wrapped(pup!=0) = wrapped(pup!=0) - mean(wrapped(pup!=0));
unwrapped = wrapped;

end
