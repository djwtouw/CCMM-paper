function [values, rownames, colnames] = load_uci_results(path)
    values = readtable(path,'TreatAsMissing', {'.','NA'});
    fid = fopen(path);
    colnames = string(strsplit(fgetl(fid), ','));
    fclose(fid);
    rownames = string(table2cell(values(:, 1)));
    values = table2array(values(:, 2:end));
end
