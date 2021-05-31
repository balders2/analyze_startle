function analyze_startle(fname)



% version created 5/16/18


configstr = fileread('config.txt');
configcell = strsplit(configstr,';');
for i=1:length(configcell)
    tmpstr = configcell{i};
    eval(tmpstr);
end

cfg.datatype = fname(length(fname)-2:end);

if ~isfield(cfg, 'doblcentergraphraw')
    cfg.doblcenter = 'no';
end;

if ~isfield(cfg, 'downinterp')
    cfg.downinterp = 'no';
end;

if ~isfield(cfg, 'filtertype')
    cfg.filtertype = 'none';
end;

if ~isfield(cfg, 'dorectify')
    cfg.dorectify = 'yes';
end;

if ~isfield(cfg, 'dosmooth')
    cfg.dosmooth = 'yes';
end;


expdir = pwd;
resdir = strcat(expdir, '/results');
emgdir = strcat(expdir, '/startle');
pngdir = strcat(expdir, '/pngs');

mkdir(resdir);
mkdir(emgdir);
mkdir(pngdir);


%import edf
if strcmp(cfg.datatype, 'EDF')
    [hdr, timeseries] = edfread(fname);
    subject = fname(1:length(fname)-4);
    %set sampling rate
    samplingrate = hdr.samples(cfg.channel);
    channelname = num2str(cfg.channel);
    
end


if strcmp(cfg.datatype, 'mat')
    clear isi samples_per_second
    load(fname);
    subject = fname(1:length(fname)-4);
    samplingrate = samples_per_second;
    cfg.altttlchannels = [];

    for c=1:length(channels)
        chanlen(c) = length(channels{c}.data);
    end;
    
        maxchanlen = max(chanlen);

        clear gcf;
        close all;
        set(gcf, 'Visible', 'Off');
        hold on;
        
        for c=1:length(channels)
        
        if length(channels{c}.data) < maxchanlen
            long = linspace(1,maxchanlen,maxchanlen);
            interpindx = linspace(1,length(long),length(channels{c}.data));
            timeseries(:,c) = interp1(interpindx,channels{c}.data,long);
        else
            timeseries(:,c) = channels{c}.data;
        end;

        if strcmp(channels{c}.name, 'Digital input')
            cfg.altttlchannels = [cfg.altttlchannels c];
        end;
        
        if strcmp(channels{c}.name, 'EMG100C')
            cfg.altchannel= c;
        end;
        
        subplot(length(channels),1,c);
        plot(timeseries(:,c))
        
    end;
    
    saveas(gcf, strcat(pngdir, '/', 'all_channels.png'));
    
    if ~isfield(cfg, 'channel')
        cfg.channel = cfg.altchannel;
    end
    
    if ~isfield(cfg, 'ttlchannels')
        cfg.ttlchannels = cfg.altttlchannels;
    end
    
    channelname = channels{cfg.channel}.name;
    
    try
        channelname = strrep(channelname,',', '_');
        channelname = strrep(channelname,' ', '_');
    catch
    end;

    clear channels
    
    
    
    cfg.ttlchannel = size(timeseries,2) + 1;
    d = linspace(0,length(cfg.ttlchannels)-1,length(cfg.ttlchannels));
    bin = 2.^d;
    dig = timeseries(:,cfg.ttlchannels)/5;
    binmat = repmat(bin,length(dig),1);
    dbin = binmat .* dig;
    timeseries(:,cfg.ttlchannel) = sum(dbin,2);
    timeseries = timeseries';
end

if ~isfield(cfg, 'codes')
    un = unique(timeseries(cfg.ttlchannel,:));
    cfg.codes = un(un>0);
end;

indices = [];
for i=1:length(cfg.codes)
    
    deriv = [0, diff(timeseries(cfg.ttlchannel,:))];
    indices = [indices, find(deriv == cfg.codes(i))];
    
end;

trialonsets = sort(indices);
%setup sampling rate
nyquistrate = samplingrate/2;
baselinelength = cfg.baselinelengthms/1000*samplingrate;
trialstart = cfg.trialstartms/1000*samplingrate;
trialend = cfg.trialendms/1000*samplingrate;

if strcmp(cfg.downinterp,'yes')
    totalsamples = size(timeseries,2);
    wnsamples = cfg.wnlengthms/1000*samplingrate;
    onsetmarkers = zeros([totalsamples,1]);
    onsetmarkers(trialonsets) = 1;
    
    onsetvals = zeros([totalsamples,1]);
    onsetvals(trialonsets) = timeseries(cfg.channel,trialonsets-baselinelength);
    
    
    wnline = ones([wnsamples,1]);
    wnmarkers = conv(onsetmarkers,wnline);
    wnmarkers = abs(wnmarkers-1);
    wnmarkers = wnmarkers(1:totalsamples);
    
    
    wnvals = conv(onsetvals,wnline);
    wnvals = wnvals(1:totalsamples);
    
    wndat = timeseries(cfg.channel,:)';
    dinterp = (wndat .* wnmarkers) + wnvals;
else
    dinterp = timeseries(cfg.channel,:);
end;


%bandpass
if strcmp(cfg.filtertype,'lowpass')
    maxorder = min([max([round(cfg.upperfreq/2,0), 2]), 9]);
    [b,a] = butter(maxorder,cfg.upperfreq/nyquistrate');
    df = filter(b,a,dinterp);
elseif strcmp(cfg.filtertype,'highpass')
    maxorder = 9;
    [b,a] = butter(maxorder,cfg.lowerfreq/nyquistrate', 'high');
    df = filter(b,a,dinterp);
elseif strcmp(cfg.filtertype,'bandpass')
    maxorder = min([max([round(cfg.upperfreq/2,0), 2]), 9])
    [b,a] = butter(maxorder,[cfg.lowerfreq cfg.upperfreq]/nyquistrate);
    df = filter(b,a,dinterp);
else
    df = timeseries(cfg.channel,:);
end;

%rectify
if strcmp(cfg.dorectify,'yes')
    dabs = abs(df);
    dshilbert = abs(hilbert(df));
else
    dabs = df;
    dshilbert = hilbert(df);
end;

%smooth

if strcmp(cfg.dosmooth,'yes')
    windowsize = cfg.rectifyms/1000*samplingrate;
    b = (1/windowsize)*ones(1,windowsize);
    a = 1;
    dsmooth = filter(b,a,dabs);
    dhsmooth = filter(b,a,dshilbert);
else
    dsmooth = dabs;
    dhsmooth = dshilbert;
end;

%calculate diff
%calculate peak onset
%add pos neg to figname

if ~isempty(trialonsets)
    smoothedstd = std(dsmooth(trialonsets(1):trialonsets(length(trialonsets))));
    hsmoothedstd = std(dhsmooth(trialonsets(1):trialonsets(length(trialonsets))));
    
    for j=1:length(trialonsets)
        rawbaseline = dinterp(trialonsets(j) - baselinelength:trialonsets(j));
        smoothedbaseline = dsmooth(trialonsets(j) - baselinelength:trialonsets(j));
        hsmoothedbaseline = dhsmooth(trialonsets(j) - baselinelength:trialonsets(j));
        rawtrial = dinterp(trialonsets(j) + trialstart:trialonsets(j) + trialend);
        filteredtrial = df(trialonsets(j) + trialstart:trialonsets(j) + trialend);
        rectifiedtrial = dabs(trialonsets(j) + trialstart:trialonsets(j) + trialend);
        hilberttrial = dshilbert(trialonsets(j) + trialstart:trialonsets(j) + trialend);
        smoothedtrial = dsmooth(trialonsets(j) + trialstart:trialonsets(j) + trialend);
        hsmoothedtrial = dhsmooth(trialonsets(j) + trialstart:trialonsets(j) + trialend);
        
        if strcmp(cfg.doblcentergraphraw,'yes')
            raw = dinterp(trialonsets(j) - baselinelength:trialonsets(j) + trialend) - mean(rawbaseline);
        else
            raw = dinterp(trialonsets(j) - baselinelength:trialonsets(j) + trialend);
        end;
        
        filtered = df(trialonsets(j) - baselinelength:trialonsets(j) + trialend);
        rectified = dabs(trialonsets(j) - baselinelength:trialonsets(j) + trialend);
        hilbertts = dshilbert(trialonsets(j) - baselinelength:trialonsets(j) + trialend);
        smoothed = dsmooth(trialonsets(j) - baselinelength:trialonsets(j) + trialend);
        hsmoothed = dhsmooth(trialonsets(j) - baselinelength:trialonsets(j) + trialend);
        ttlsignal = timeseries(cfg.ttlchannel,trialonsets(j) - baselinelength:trialonsets(j) + trialend);
        
        smoothedstdrep(j) = smoothedstd;
        hsmoothedstdrep(j) = hsmoothedstd;
        trialonset(j) = trialonsets(j)/samplingrate;
        blmean(j) = mean(smoothedbaseline);
        blminmax(j) = max(smoothedbaseline) - min(smoothedbaseline);
        blstd(j) = std(smoothedbaseline);
        hblmean(j) = mean(hsmoothedbaseline);
        hblminmax(j) = max(hsmoothedbaseline) - min(hsmoothedbaseline);
        hblstd(j) = std(hsmoothedbaseline);
        rawtrialpeak(j) = max(rawtrial);
        filteredtrialpeak(j) = max(filteredtrial);
        rectifiedtrialpeak(j) = max(rectifiedtrial);
        hilberttrialpeak(j) = max(hilberttrial);
        smoothedtrialpeak(j) = max(smoothedtrial);
        hsmoothedtrialpeak(j) = max(hsmoothedtrial);
        ttlval(j) = timeseries(cfg.ttlchannel,trialonsets(j));
        channelnametrial{j} = channelname;
        
        trialcount(j) = j;
        trialsub{j} = subject;
        try
            peakonset(j) = (min(find(smoothedtrial == smoothedtrialpeak(j))) + trialstart)/samplingrate;
        catch
            peakonset(j) = NaN;
        end;
        
        smootheddiff(j) = smoothedtrialpeak(j) - blmean(j);
        hsmootheddiff(j) = hsmoothedtrialpeak(j) - hblmean(j);
        
        if smootheddiff(j) > blminmax(j)
            blinknoblink = 'blink';
            blinktrial(j) = 1;
        else
            blinknoblink = 'noblink';
            blinktrial(j) = 0;
        end;
        
        if smoothedstd*cfg.artifactthresh > blstd(j)
            goodbad = 'gd';
            keeptrial(j) = 1;
        else
            goodbad = 'bd';
            keeptrial(j) = 0;
        end;
        titlestr = strcat(subject,'.', channelnametrial{j}, '.', num2str(j), '.', num2str(ttlval(j)), '.', goodbad, '.', blinknoblink);
        figfname = strcat(titlestr, '.png');
        xunits = linspace(baselinelength*-1,trialend,baselinelength+trialend+1)/samplingrate*1000;
        
        disp('################################save figure')
        clear gcf;
        close all;
        set(gcf, 'Visible', 'Off');
        hold on;
        title(titlestr);
        a1 = plot(xunits,raw);
        a2 = plot(xunits,filtered);
        a3 = plot(xunits,rectified);
        a5 = plot(xunits,hilbertts);
        a4 = plot(xunits,smoothed);
        a4 = plot(xunits,hsmoothed);
        l1 = 'raw';
        l2 = 'filtered';
        l3 = 'rectified';
        l4 = 'hilbert';
        l5 = 'smoothed';
        l6 = 'hsmoothed';
        legend(l1, l2, l3, l4, l5, l6);
        legend('Location','northwest');
        
        saveas(gcf, strcat(pngdir, '/', figfname));
        
        
        
    end;
    
    java.lang.Runtime.getRuntime.gc;
    clear data
    trialsub = trialsub';
    channelnametrial = channelnametrial';
    trialcount = trialcount';
    trialonset = trialonset';
    peakonset = peakonset';
    ttlval = ttlval';
    keeptrial = keeptrial';
    blinktrial = blinktrial';
    smoothedstdrep = smoothedstdrep';
    hsmoothedstdrep = hsmoothedstdrep';
    blmean = blmean';
    hblmean = hblmean';
    blminmax = blminmax';
    hblminmax = hblminmax';
    blstd = blstd';
    hblstd = hblstd';
    rawtrialpeak = rawtrialpeak';
    filteredtrialpeak = filteredtrialpeak';
    rectifiedtrialpeak = rectifiedtrialpeak';
    hilberttrialpeak = hilberttrialpeak';
    smoothedtrialpeak = smoothedtrialpeak';
    hsmoothedtrialpeak = hsmoothedtrialpeak';
    smootheddiff = smootheddiff';
    hsmootheddiff = hsmootheddiff';
    
    cleanraw = smootheddiff;
    noblinkidx = find(blinktrial == 0);
    cleanraw(noblinkidx) = 0;
    badidx = find(keeptrial == 0);
    cleanraw(badidx) = nan;
    tscore = (cleanraw-nanmean(cleanraw))/nanstd(cleanraw)*10 + 50;
    
    data = table(trialsub, channelnametrial, trialcount, trialonset, peakonset, ttlval, keeptrial, blinktrial, smoothedstdrep, hsmoothedstdrep, blmean, hblmean, blminmax, hblminmax, blstd, hblstd, rawtrialpeak, filteredtrialpeak, rectifiedtrialpeak, hilberttrialpeak, smoothedtrialpeak, hsmoothedtrialpeak, smootheddiff, hsmootheddiff, cleanraw, tscore);
    writetable(data, strcat(resdir, '/', subject, '.', channelname, '.csv'));
else
    errfid = fopen(strcat(resdir, '/', subject, '.txt'),'wt');
    fprintf(errfid, 'no trials detected')
    fclose(errfid);
end;

movefile(fname, emgdir);

clear all
clear classes
clear java
close all
java.lang.Runtime.getRuntime.gc;
