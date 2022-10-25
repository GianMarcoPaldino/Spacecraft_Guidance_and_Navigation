function [ satrec, longstr1, longstr2 ] = read_TLE( SAT_ID, whichconst )
% Reads TLE number TLENUM in file (or last TLE if TLEUM is larger than number of tle)
    
    if ~ischar(SAT_ID)
        SAT_ID = num2str(SAT_ID,'%05d');
    end

    finput = strcat('tle\',SAT_ID,'.tle');

    fid = fopen(finput);

    % Read first line
    longstr1 = fgetl(fid);
        
    % Read second line
    longstr2 = fgetl(fid);

    % Initialize
    typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
    opsmode    = 'a';  % afspc approach
    [satrec] = twoline2rv( longstr1, longstr2, typerun, 'e', opsmode, whichconst);

    longstr1 = fgetl(fid);

    fclose(fid);
    
end

