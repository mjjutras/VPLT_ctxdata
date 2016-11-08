datfil='TT091102.2'
setnum='SET051'


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

% get calibration trials, set up per structure containing trial info
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

% figure out sampling rate based on trial durations and amount of eye data
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

%     figure;scatter(x,y) % 9 calibration points, should look like a cross

% simple eyedata calibration, determines offset and gain from eye position
% during 9-point calibration
meanxorigin = mean([x(6) x(2) x(1) x(5) x(9) ],2);
xscale = mean([6/(x(8)-meanxorigin) 3/(x(4)-meanxorigin) 3/(abs(x(3)-meanxorigin)) 6/(abs(x(7)-meanxorigin))],2);
meanyorigin = mean([y(7) y(3) y(1) y(4) y(8) ],2);
yscale = mean([6/(y(6)-meanyorigin) 3/(y(2)-meanyorigin) 3/(abs(y(5)-meanyorigin)) 6/(abs(y(9)-meanyorigin))],2);

% get VPLT trials, set up per structure containing trial info
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

% translate the eye data from voltage units to screen coordinates
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
    
%% plot image w/ eyescan

% best to change this line so you can view pictures in small groups
% (instead of cndlop=1:200)
for cndlop=1:size(trind,1)
    
    if ltmat(cndlop,1)>4000 %for limiting to pics with a minimum encoding LT
    
        [imgmtx,map] = imread(strcat('R:\Buffalo Lab\eblab\VPC Stimuli\BMP\',setnum,'\',num2str(cndlop),'.bmp'));
        figure;
        image(imgmtx);

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

