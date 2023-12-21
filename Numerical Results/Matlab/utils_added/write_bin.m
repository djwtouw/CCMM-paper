% write_bin function
function write_bin(val, fname)
    fid = fopen(fname, 'wb');
    fwrite(fid, val, 'double');
    fclose(fid);
end
