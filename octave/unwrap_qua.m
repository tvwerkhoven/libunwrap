%
% Uses a disk interface of an unwrapping routine
%
function unwrapped = unwrap_qua(ph, qua)

% Create temp names
tmpbase  = 'tmp';
quaname  = [tmpbase '_quali.raw'];
phname   = [tmpbase '_towrap.raw'];
resuname = [tmpbase '_unwrapped.raw'];

phdim1 = size(ph, 1);
phdim2 = size(ph, 2);

% Write for C code
fid = fopen(quaname, 'wb');
fwrite(fid, qua, 'double');
fclose(fid);
fid = fopen(phname, 'wb');
fwrite(fid, ph, 'double');
fclose(fid);

cmd = ['unwrap-test ' phname ' ' quaname ' ' resuname ' ' ...
       num2str(phdim1) ' ' num2str(phdim2)];
system(cmd);

% Read the unwrapped phase
fid = fopen(resuname, 'rb');
unwrapped = fread(fid, [phdim1 phdim2], 'double');
fclose(fid);

% Remove tmp files
system(['rm ' quaname]);
system(['rm ' phname]);
system(['rm ' resuname]);

end
