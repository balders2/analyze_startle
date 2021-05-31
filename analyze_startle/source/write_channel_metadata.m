function write_channel_metadata(fname)

if ~exist('fname')
fnames = dir('*.mat')
fname = fnames(1).name
end;


load(fname);
subject = fname(1:length(fname)-4);



for channelnum=1:length(channels)
    channelmeta{channelnum} = rmfield(channels{channelnum}, 'data');
    
    fields = fieldnames(channelmeta{channelnum})
    for f=1:length(fields)
        try 
            channelmeta{channelnum}.(fields{f}) = strrep(channelmeta{channelnum}.(fields{f}),',', '.');
        catch
        end;
    end;
        

        
    
    
    channeltbl = struct2table(channelmeta{channelnum});
    channelnumtbl = table(channelnum);
    tbl = [channelnumtbl, channeltbl];
    s(channelnum) = table2struct(tbl);
end;

channelmetadata = struct2table(s);
writetable(channelmetadata, strcat(subject, '.channel_metadata.csv'));
