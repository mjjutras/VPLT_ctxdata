datfil='TT091103.2'
setnum='SET051'
stmtyp='bmp'


if strmatch(datfil(find(double(datfil==46))+1:end),'nex') % for nex files
    
    ini=datfil(1:2);
    if strcmp(ini,'IW')==1 || strcmp(ini,'iw')==1
        datfil=['S:\NEX files\Irwin\' datfil];
    elseif strcmp(ini,'MP')==1 || strcmp(ini,'mp')==1
        datfil=['S:\NEX files\Peepers\' datfil];
    elseif strcmp(ini,'WR')==1 || strcmp(ini,'wr')==1
        datfil=['S:\NEX files\Wilbur\' datfil];
    elseif strcmp(ini,'TT')==1 || strcmp(ini,'tt')==1
        datfil=['S:\NEX files\Timmy\' datfil];
    elseif strcmp(ini,'JN')==1 || strcmp(ini,'jn')==1
        datfil=['S:\NEX files\Guiseppe\' datfil];
    elseif strcmp(ini,'TD')==1 || strcmp(ini,'td')==1
        datfil=['S:\NEX files\Theodore\' datfil];
    end
    
    if ~exist(datfil,'file')
        datfil=['S:\NEX files\SUA-VPC\aligned\' datfil(find(datfil=='\',1,'last')+1:find(datfil=='.',1,'last')-1) '-aligned.nex'];
    end

    header = getnexheader(datfil);
    smpfrq = header.varheader(end-2).wfrequency;
    mrk    = getnexmrk(header,smpfrq);
    perdef = 'VPCpermrk';
    [per]  = feval(perdef,mrk);

    numsmp=7000;

    dat=proeyedat(header,per, smpfrq, numsmp);

    [x, y] = VPCgetclrchg(datfil);
    meanxorigin = mean([x(6) x(2) x(1) x(5) x(9) ],2);
    xscale = mean([6/(x(8)-meanxorigin) 3/(x(4)-meanxorigin) 3/(abs(x(3)-meanxorigin)) 6/(abs(x(7)-meanxorigin))],2);
    meanyorigin = mean([y(7) y(3) y(1) y(4) y(8) ],2);
    yscale = mean([6/(y(6)-meanyorigin) 3/(y(2)-meanyorigin) 3/(abs(y(5)-meanyorigin)) 6/(abs(y(9)-meanyorigin))],2);
    % xorigin = x(1);
    % xscale = mean([6/(x(8)-meanxorigin) 6/(abs(x(7)-meanxorigin))],2);
    % yorigin = y(1);
    % yscale = mean([6/(y(6)-meanyorigin) 6/(abs(y(9)-meanyorigin))],2);


    eyedat = [];
    for trlop=1:size(dat.trial,2)
        if size(dat.trial{trlop},2)>=100
            eyedat{trlop}(1,:) = (dat.trial{trlop}(1,:)-mean(dat.trial{trlop}(1,1:100),2)) .* xscale;
            eyedat{trlop}(2,:) = (dat.trial{trlop}(2,:)-mean(dat.trial{trlop}(2,1:100),2)) .* yscale;
        else
            eyedat{trlop}(1,:) = (dat.trial{trlop}(1,:)-mean(dat.trial{trlop}(1,1:size(dat.trial{trlop},2)),2)) .* xscale;
            eyedat{trlop}(2,:) = (dat.trial{trlop}(2,:)-mean(dat.trial{trlop}(2,1:size(dat.trial{trlop},2)),2)) .* yscale;
        end
    end

    clear cnd lt
    numrpt = size(per,2);
    i=1;
    for rptlop = 1:numrpt
        lt(i)=per(rptlop).endsmpind-per(rptlop).begsmpind;
        cnd(i)=per(rptlop).cnd;
        i=i+1;
    end

    clear trind
    i=1;
    for cndnum=min(cnd):max(cnd)
        if length(find(cnd==cndnum))>=2
            trind(i,[1 2])=find(cnd==cndnum,2,'first');
            i=i+1;
        end
    end

    ltmat=lt(trind);
    
else % for ctx files
    
    ini=datfil(1:2);
    if strcmp(ini,'IW')==1 || strcmp(ini,'iw')==1
        datfil=['R:\Buffalo Lab\Cortex Data\Irwin\' datfil];
    elseif strcmp(ini,'MP')==1 || strcmp(ini,'mp')==1
        datfil=['R:\Buffalo Lab\Cortex Data\Peepers\' datfil];
    elseif strcmp(ini,'WR')==1 || strcmp(ini,'wr')==1
        datfil=['R:\Buffalo Lab\Cortex Data\Wilbur\' datfil];
    elseif strcmp(ini,'TT')==1 || strcmp(ini,'tt')==1
        datfil=['R:\Buffalo Lab\Cortex Data\Timmy\' datfil];
    elseif strcmp(ini,'JN')==1 || strcmp(ini,'jn')==1
        datfil=['R:\Buffalo Lab\Cortex Data\Guiseppe\' datfil];
    elseif strcmp(ini,'TD')==1 || strcmp(ini,'td')==1
        datfil=['R:\Buffalo Lab\Cortex Data\Theodore\' datfil];
    end

    [time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(datfil);
    
    numrpt = size(event_arr,2);
    valrptcnt = 0;
    clear per clrchgind
    for rptlop = 1:numrpt
        if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) < 1010)) ~=0
            if size(find(event_arr(:,rptlop) == 200)) ~=0
                perbegind = find(event_arr(:,rptlop) == 23);
                perendind = find(event_arr(:,rptlop) == 24);
                cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                begtimdum = time_arr(perbegind,rptlop);
                endtimdum = time_arr(perendind,rptlop);
                if endtimdum > begtimdum
                    valrptcnt = valrptcnt + 1;
                    clrchgind(valrptcnt)=rptlop;
                    per(valrptcnt).begsmpind = begtimdum;
                    per(valrptcnt).endsmpind = endtimdum;
                    per(valrptcnt).begpos = 1;
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                    per(valrptcnt).blk = event_arr(blknumind,rptlop);
                    per(valrptcnt).allval = event_arr(:,rptlop);
                    per(valrptcnt).alltim = time_arr(:,rptlop);
                end
            end
        end
    end

%     if strmatch('WR080422.2',datfil(find(double(datfil==92),1,'last')+1:end)) % Wilbur's first Expert
%         samprate=10;
%     else
%         samprate=5;
%     end

    % figure out sampling rate
    eyedur=[];
    for k=1:size(time_arr,2)
        eyedur=[eyedur; time_arr(event_arr(:,k)==101,k)-time_arr(event_arr(:,k)==100,k)];
    end
    eogdur10=[];
    for k=1:size(eog_arr,2)
        eogdur10=[eogdur10; length(find(~isnan(eog_arr(:,k))))/.2];
    end
    eogdur5=[];
    for k=1:size(eog_arr,2)
        eogdur5=[eogdur5; length(find(~isnan(eog_arr(:,k))))/.4];
    end
    if abs(mean(eogdur5-eyedur))>abs(mean(eogdur10-eyedur))
        samprate=10;
    else
        samprate=5;
    end


    clear cnd
    numrpt = size(per,2);
    for rptlop = 1:numrpt
        cnd(rptlop)=per(rptlop).cnd;
    end

    evnnmb=2:2:size(eog_arr,1);
    oddnmb=1:2:size(eog_arr,1);

    clear x y
    cndlst=unique(cnd);
    for k=1:length(cndlst)
        cndind=find(cnd==cndlst(k));
        allind=clrchgind(cndind);
        for l=1:length(allind)
            x{k}(l)=mean(eog_arr(intersect(floor(((per(cndind(l)).begsmpind-1000)/samprate)*2):(floor((per(cndind(l)).endsmpind-1000)/samprate))*2,oddnmb),allind(l)));
            y{k}(l)=mean(eog_arr(intersect(floor(((per(cndind(l)).begsmpind-1000)/samprate)*2):(floor((per(cndind(l)).endsmpind-1000)/samprate))*2,evnnmb),allind(l)));
        end
    end

    clear meanx meany
    for k=1:length(x)
        meanx(k)=mean(x{k});
    end
    for k=1:length(y)
        meany(k)=mean(y{k});
    end

    clear x y
    x=meanx; y=meany;
    
%     figure;scatter(x,y) % 9 calibration points, should look like a cross,

    meanxorigin = mean([x(6) x(2) x(1) x(5) x(9) ],2);
    xscale = mean([6/(x(8)-meanxorigin) 3/(x(4)-meanxorigin) 3/(abs(x(3)-meanxorigin)) 6/(abs(x(7)-meanxorigin))],2);
    meanyorigin = mean([y(7) y(3) y(1) y(4) y(8) ],2);
    yscale = mean([6/(y(6)-meanyorigin) 3/(y(2)-meanyorigin) 3/(abs(y(5)-meanyorigin)) 6/(abs(y(9)-meanyorigin))],2);

    numrpt = size(event_arr,2);
    valrptcnt = 0;
    clear per vpcind
    for rptlop = 1:numrpt
        if size(find(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop) >= 1010)) ~=0
            if size(find(event_arr(:,rptlop) == 200)) ~=0
                perbegind = find(event_arr(:,rptlop) == 23);
                perendind = find(event_arr(:,rptlop) == 24);
                cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                begtimdum = time_arr(perbegind,rptlop);
                endtimdum = time_arr(perendind,rptlop);
                if endtimdum > begtimdum
                    valrptcnt = valrptcnt + 1;
                    vpcind(valrptcnt)=rptlop;
                    per(valrptcnt).begsmpind = begtimdum - 500;
                    per(valrptcnt).endsmpind = endtimdum;
                    per(valrptcnt).begpos = 1;
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                    per(valrptcnt).blk = event_arr(blknumind,rptlop);
                    per(valrptcnt).allval = event_arr(:,rptlop);
                    per(valrptcnt).alltim = time_arr(:,rptlop);
                end
            end
        end
    end

    eyedat = [];
    for trlop=1:size(per,2)
        eyedat{trlop}(1,:) = (eog_arr(intersect(floor(((per(trlop).begsmpind+500-1000)/samprate)*2):(floor((per(trlop).endsmpind-1000)/samprate)*2),oddnmb),vpcind(trlop)) - mean(eog_arr(intersect(floor(((per(trlop).begsmpind-1000)/samprate)*2):(floor((per(trlop).begsmpind+500-1000)/samprate)*2),oddnmb),vpcind(trlop)))) .* xscale;
        eyedat{trlop}(2,:) = (eog_arr(intersect(floor(((per(trlop).begsmpind+500-1000)/samprate)*2):(floor((per(trlop).endsmpind-1000)/samprate)*2),oddnmb)+1,vpcind(trlop)) - mean(eog_arr(intersect(floor(((per(trlop).begsmpind-1000)/samprate)*2):(floor((per(trlop).begsmpind+500-1000)/samprate)*2),oddnmb)+1,vpcind(trlop)))) .* yscale;
    end

    clear lt cnd
    numrpt = size(per,2);
    for rptlop = 1:numrpt
        lt(rptlop)=per(rptlop).endsmpind-per(rptlop).begsmpind-500;
        cnd(rptlop)=per(rptlop).cnd;
    end

    clear trind
    i=1;
    for cndnum=min(cnd):max(cnd)
        if length(find(cnd==cndnum))>=2
            trind(i,[1 2])=find(cnd==cndnum,2,'first');
            i=i+1;
        end
    end
    ltmat=lt(trind);
    
end

%% plot image w/ eyescan

% best to change this line so you can view pictures in small groups
% (instead of cndlop=1:200)
for cndlop=1:size(trind,1)
    
    if ltmat(cndlop,1)>2000 %for limiting to pics with a minimum encoding LT
    
        if strmatch(stmtyp,'ctx') % ctx
            if ~isempty(strfind(setnum,'b')) || str2double(setnum(find(double(setnum==84))+1:end))>=118
                if ~strcmp(setnum,'SET112b')
                    [imgmtx, dmns, notes]=loadcx(strcat('R:\Buffalo Lab\eblab\VPC Stimuli\CTX\',setnum,'\',num2str(cndlop+99),'.ctx'));
                else
                    [imgmtx, dmns, notes]=loadcx(strcat('R:\Buffalo Lab\eblab\VPC Stimuli\CTX\',setnum,'\eee',num2str(cndlop+99),'.ctx'));
                end
            else
                [imgmtx, dmns, notes]=loadcx(strcat('R:\Buffalo Lab\eblab\VPC Stimuli\CTX\',setnum,'\eee',num2str(cndlop+99),'.ctx'));
            end
            newimg = imgmtx - 127;
            lut = loadlut('bg001.lut');
            dbllut = lut / 255;
            figure;
            imshow(newimg, dbllut);

        elseif strmatch(stmtyp,'bmp') % bmp
            [imgmtx,map] = imread(strcat('R:\Buffalo Lab\eblab\VPC Stimuli\BMP\',setnum,'\',num2str(cndlop),'.bmp'));
            figure;
            image(imgmtx);
        end

        h2=axes('Position',get(gca,'Position'),'Layer','top');
        hold on
        s1=scatter(eyedat{trind(cndlop,1)}(1,:),eyedat{trind(cndlop,1)}(2,:),3,'oy','MarkerFaceColor','y');
        s2=scatter(eyedat{trind(cndlop,2)}(1,:),eyedat{trind(cndlop,2)}(2,:),3,'oc','MarkerFaceColor','c');
        l1=plot(eyedat{trind(cndlop,1)}(1,:),eyedat{trind(cndlop,1)}(2,:),'y');
        l2=plot(eyedat{trind(cndlop,2)}(1,:),eyedat{trind(cndlop,2)}(2,:),'c');
        xlim([-5.5 5.5]);
        ylim([-5.5 5.5]);
        axis off;
        title(['Encoding: ' num2str(ltmat(cndlop,1)) ' msec; Recognition: ' num2str(ltmat(cndlop,2)) ' msec']);
        toptitle([datfil((find(double(datfil)==92,1,'last')+1):end) ', Cnd ' num2str(cndlop+9) ', Item ' num2str(cndlop+30)]);
    end
end


%% plot image w/ eyescan - circles proportional to fixation duration

% construct filter for eye velocity measure
fltord = 20;
lowpasfrq = 40;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);

% analyze eye movement to determine location and duration (start/end times)
% of fixation periods
for cndlop=1:size(trind,1)
    
    novtrl=trind(cndlop,1);
    reptrl=trind(cndlop,2);
    
    % calculate velocity
    x_vnov= diff(filtfilt(flt,1, eyedat{novtrl}(1,:))) .* 1000;
    y_vnov= diff(filtfilt(flt,1, eyedat{novtrl}(2,:))) .* 1000;
    x_vrep= diff(filtfilt(flt,1, eyedat{reptrl}(1,:))) .* 1000;
    y_vrep= diff(filtfilt(flt,1, eyedat{reptrl}(2,:))) .* 1000;
    
    % combine x- and y-velocity to get overall eye velocity
    velnov = abs(complex(x_vnov,y_vnov));
    velrep = abs(complex(x_vrep,y_vrep));
    % lim = threshold for detecting saccade
    lim = 25;
    sacbegnov = find(diff(velnov > lim) > 0);
    sacendnov = find(diff(velnov > lim) < 0);
    sacbegrep = find(diff(velrep > lim) > 0);
    sacendrep = find(diff(velrep > lim) < 0);
    
    if velnov(end)>=lim
        if velnov(1)<lim
            tempbegnov=[1 sacendnov];
            tempendnov=sacbegnov;
        else
            tempbegnov=sacendnov;
            tempendnov=sacbegnov;
        end
    else
        if velnov(1)<lim
            tempbegnov=[1 sacendnov];
            tempendnov=[sacbegnov length(velnov)];
        else
            tempbegnov=sacendnov;
            tempendnov=[sacbegnov length(velnov)];
        end
    end
    
    if velrep(end)>=lim
        if velrep(1)<lim
            tempbegrep=[1 sacendrep];
            tempendrep=sacbegrep;
        else
            tempbegrep=sacendrep;
            tempendrep=sacbegrep;
        end
    else
        if velrep(1)<lim
            tempbegrep=[1 sacendrep];
            tempendrep=[sacbegrep length(velrep)];
        else
            tempbegrep=sacendrep;
            tempendrep=[sacbegrep length(velrep)];
        end
    end
    
    tempnov=[tempbegnov; tempendnov];
    temprep=[tempbegrep; tempendrep];
    
    % limit to fixation periods at least 100 ms long
    fixtimnov=tempnov(:,diff(tempnov,1)>=100);
    fixtimrep=temprep(:,diff(temprep,1)>=100);
    
    fixhor_nov=nan(size(fixtimnov,2),1);
    fixvrt_nov=nan(size(fixtimnov,2),1);
    for fixlop=1:size(fixtimnov,2)
        fixhor_nov(fixlop)=mean(eyedat{novtrl}(1,fixtimnov(1,fixlop):fixtimnov(2,fixlop))); % fixation position, horizontal
        fixvrt_nov(fixlop)=mean(eyedat{novtrl}(2,fixtimnov(1,fixlop):fixtimnov(2,fixlop))); % fixation position, vertical
    end

    fixhor_rep=nan(size(fixtimrep,2),1);
    fixvrt_rep=nan(size(fixtimrep,2),1);
    for fixlop=1:size(fixtimrep,2)
        fixhor_rep(fixlop)=mean(eyedat{reptrl}(1,fixtimrep(1,fixlop):fixtimrep(2,fixlop))); % fixation position, horizontal
        fixvrt_rep(fixlop)=mean(eyedat{reptrl}(2,fixtimrep(1,fixlop):fixtimrep(2,fixlop))); % fixation position, vertical
    end

    fixdurnov=diff(fixtimnov,1);
    fixdurrep=diff(fixtimrep,1);
    
    % plot figure
    if strmatch(stmtyp,'ctx') % ctx
        if ~isempty(strfind(setnum,'b')) || str2double(setnum(find(double(setnum==84))+1:end))>=118
            if ~strcmp(setnum,'SET112b')
                [imgmtx, dmns, notes]=loadcx(strcat('S:\VPC Stimuli\CTX\',setnum,'\',num2str(cndlop+99),'.ctx'));
            else
                [imgmtx, dmns, notes]=loadcx(strcat('S:\VPC Stimuli\CTX\',setnum,'\eee',num2str(cndlop+99),'.ctx'));
            end
        else
            [imgmtx, dmns, notes]=loadcx(strcat('S:\VPC Stimuli\CTX\',setnum,'\eee',num2str(cndlop+99),'.ctx'));
        end
        newimg = imgmtx - 127;
        lut = loadlut('bg001.lut');
        dbllut = lut / 255;
        figure;
        imshow(newimg, dbllut);
        
    elseif strmatch(stmtyp,'bmp') % bmp
        [imgmtx,map] = imread(strcat('S:\VPC Stimuli\BMP\',setnum,'\',num2str(cndlop),'.bmp'));
        figure;
        image(imgmtx);
    end

    h2=axes('Position',get(gca,'Position'),'Layer','top');
    hold on
    h=line(eyedat{trind(cndlop,1)}(1,:),eyedat{trind(cndlop,1)}(2,:),'Color','b');
    hold on
    h=line(eyedat{trind(cndlop,2)}(1,:),eyedat{trind(cndlop,2)}(2,:),'Color','r');
    for fixlop=1:length(fixdurnov)
        hold on
        h=scatter(fixhor_nov(fixlop),fixvrt_nov(fixlop),fixdurnov(fixlop),'ob','MarkerFaceColor','b');
    end
    for fixlop=1:length(fixdurrep)
        hold on
        h=scatter(fixhor_rep(fixlop),fixvrt_rep(fixlop),fixdurrep(fixlop),'or','MarkerFaceColor','r');
    end
    xlim([-5.5 5.5]);
    ylim([-5.5 5.5]);
    axis off;
    title(['Encoding: ' num2str(ltmat(cndlop,1)) ' msec; Recognition: ' num2str(ltmat(cndlop,2)) ' msec']);
    toptitle([datfil((find(double(datfil)==92,1,'last')+1):end) ', Cnd ' num2str(cnd(novtrl)-1000) ', Item ' num2str(cnd(novtrl)-1000+21)]);
    
end

