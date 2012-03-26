
%
% This script can be used to test unwrapping routines.
%

% Put here your phase and quality
if 1==1
  ph  = phsms(:,:,7);
  qua = ampls(:,:,7) .*pup;
  phdim = size(ph,1);
  
  unwrapped2 = unwrap_qua(ph, qua);
  unwrapped2(pup~=0) = unwrapped2(pup~=0)-mean(unwrapped2(pup~=0));
end  


% Create easy and nice test phase
if 1==1
  ph = [0,   0.1,   0.2,        0.3; ...
	0.1, 2*pi,  2*pi+0.3,  0.4; ...
	0.3, 0,    0,          0.1];
  qua = 1+[0, 0, 0, 0; ...
	   0, 1, 0, 0; ...
	   2,.1, 0, 0];
  phdim1 = size(ph,1);
  phdim2 = size(ph,2);
end



% Create temp names
tmpbase  = 'tmp';
quaname  = [tmpbase '_quali.raw'];
phname   = [tmpbase '_towrap.raw'];
resuname = [tmpbase '_unwrapped.raw'];


% Write for C code
fid = fopen(quaname, 'wb');
fwrite(fid, qua, 'double');
fclose(fid);
fid = fopen(phname, 'wb');
fwrite(fid, ph, 'double');
fclose(fid);


cmd = ['../c/test/unwrap-test ' phname ' ' quaname ' ' resuname ' ' ...
       num2str(phdim1) ' ' num2str(phdim2)];
%cmd = ['../c/test/unwrap-test ' resuname];
system(cmd);

  
% Read the unwrapped phase
fid = fopen(resuname, 'rb');
wrapped = fread(fid, 'double');
fclose(fid);
wrapped = reshape(wrapped, size(ph));
if size(pup,1) == phdim1
  wrapped(pup!=0) = wrapped(pup!=0) - mean(wrapped(pup!=0));
end
unwrappedd = wrapped;

if 1==0
  system(['rm ' quaname]);
  system(['rm ' phname]);
  system(['rm ' resuname]);
end

if size(ph,1) > 64
  clf; hold on
  plot(unwrappedd(64,:))
  plot(unwrapped2(64,:), 'r')
  plot(ph(64,:), 'g')
end
