function save_uci_results(path, values, rownames, colnames)
    fileID = fopen(path,'w');

    fprintf(fileID, "%s,", colnames(1:end-1));
    fprintf(fileID, "%s\n", colnames(end));

    for i = 1:length(rownames)
       fprintf(fileID, "%s,", rownames(i));
       fprintf(fileID, "%d,", values(i, 1:end-1));
       fprintf(fileID, "%d\n", values(i, end));
    end

    fclose(fileID);
end