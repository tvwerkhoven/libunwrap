
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


cmd = ['../c/unwrap/test_unwrap ' phname ' ' quaname ' ' resuname];
system(cmd);

  
% Read the unwrapped phase
fid = fopen(resuname, 'rb');
wrapped = fread(fid, 'double');
fclose(fid);
wrapped = reshape(wrapped, size(ph));
wrapped(pup!=0) = wrapped(pup!=0) - mean(wrapped(pup!=0));
unwrapped = wrapped;

if 1==0
  system(['rm ' quaname]);
  system(['rm ' phname]);
  system(['rm ' resuname]);
end

clf; hold on
plot(unwrapped(64,:))
plot(unwrapped2(64,:), 'r')
plot(ph(64,:), 'g')

