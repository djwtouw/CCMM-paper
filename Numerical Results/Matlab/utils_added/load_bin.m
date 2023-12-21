% load_bin function
function val = load_bin(fname)
    fid = fopen(fname, 'rb');
    val = fread(fid, Inf, 'double');
    fclose(fid);
end
