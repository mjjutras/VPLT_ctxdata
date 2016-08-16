% vpcanactx loads data for a single VPC session from a cortex file,and 
% plots looking times as well as absolute and percent change in looking
% time between presentations

% v3 070409 mjj added clear per line to prevent errors when running for
% multiple sessions
% 
% v2 070226 MJJ removed input for initials, now only need to type filename

filename=input('Type datafile name without quotes:  ', 's')
if strcmp(filename(1:2),'IW')==1 || strcmp(filename(1:2),'iw')==1
    monkey='Irwin';
end
if strcmp(filename(1:2),'MP')==1 || strcmp(filename(1:2),'mp')==1
    monkey='Peepers';
end
if strcmp(filename(1:2),'WR')==1 || strcmp(filename(1:2),'wr')==1
    monkey='Wilbur';
end
if strcmp(filename(1:2),'TD')==1 || strcmp(filename(1:2),'td')==1
    monkey='Theodore';
end
name=strcat('S:\Cortex Data\',monkey,'\',filename);
[time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(name);

numrpt = size(event_arr,2);
valrptcnt = 0;
clear per
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

if exist('per') == 0
    error('Not a VPC session');
end

%%Plot VPC Behavioral Data%%
clear ('cnd','lt','ltred','encind','recind','delind','delind1','ltenc',...
    'ltrec','ltenc1','ltrec1','ltred','ltred1','ltper','ltper1');

numrpt = size(per,2);
i=1;
 for rptlop = 1:numrpt
   lt(i)=per(rptlop).endsmpind-per(rptlop).begsmpind;
   cnd(i)=per(rptlop).cnd;
   %%blk(rptlop)=per(rptlop).blk;
   i=i+1;
 end

%%Sort conditions to pair first and second presentations of a stimulus
%%Create an index of first presentations - encind 
%%Create an index of second presentations - recind
%%Determine the delay between first and second presentations - delind
%%Calculate looking times for first and second presentations - ltenc and ltrec
[cndsrt,indx]=sort(cnd,2,'ascend');

%%Next two lines (4:405) for sessions 4/21/06-6/13/06
if str2num(name(end-7:end))>=604212 && str2num(name(end-7:end))<=606132
    cndsrt=cndsrt(4:405);
    indx=indx(4:405);
end

for cndlop=(cndsrt(1,1)):(cndsrt(1,end))
    if size(find(cndsrt==cndlop),2)~=2
        disp(strcat('Error in condition matrix:  ',num2str(cndlop)));
        if size(find(cndsrt==cndlop),2)==3
            h=find(cndsrt==cndlop);
            cndsrt=cndsrt([1:(h(1,end)-1) (h(1,end)+1):size(cndsrt,2)]);
            indx=indx([1:(h(1,end)-1) (h(1,end)+1):size(indx,2)]);
            disp('Removing third occurrence of condition')
        end
        if size(find(cndsrt==cndlop),2)==1
            h=find(cndsrt==cndlop);
            cndsrt=cndsrt([1:(h-1) (h+1):size(cndsrt,2)]);
            indx=indx([1:(h-1) (h+1):size(indx,2)]);
            disp('Removing first occurrence of condition')
        end
    end
end

cndmat=reshape(cndsrt,2,size(cndsrt,2)/2);
indmat=reshape(indx,2,size(indx,2)/2);
encind=indmat(1,:);
recind=indmat(2,:);
delind=recind-encind;
ltenc=lt(encind);
ltrec=lt(recind);
ltred=ltenc-ltrec;



%%Include only trials in which the encoding time was greater than 500ms
dum=find(ltenc>=500);
dum1=find(ltrec>=500);
ltenc1=ltenc(dum);
ltrec1=ltrec(dum);
ltred1=ltred(dum);
delind1=delind(dum);

%%Calculate the reduction in lt as a percentage of the encoding time
ltper=(ltenc-ltrec)./ltenc;
ltper1=ltper(dum);


%%Graph the results%%
figure;
subplot(2,2,1)
%%Scatterplot of Looking Times during Encoding vs Recognition 
scatter(ltenc1,ltrec1)
line(0:5000,0:5000)
xlabel('Encoding [ms]', 'FontSize',10,'FontName', 'Arial')
ylabel('Recognition [ms]', 'FontSize',10,'FontName', 'Arial')

subplot(2,2,2)
%%Bar Graph of average looking time
ltmat1=[ltenc1;ltrec1]';
bar(mean(ltmat1))
hold on
x=[1 2];
y=mean(ltmat1);
s=std(ltmat1)./sqrt(size(ltmat1,1));
errorbar(x,y,s,s,'.b','LineWidth',1)


subplot(2,2,3)
%%Distribution of Reduction in Looking Times 
%%(Only Trials with encoding>500ms)
hist(ltred1,30)
xlabel('Reduction in Looking Time [ms]', 'FontSize',10,'FontName', 'Arial')
ylabel('Count', 'FontSize',10,'FontName', 'Arial')
text(-3000,20,strcat('med=',num2str(median(ltred1)),'ms'))

subplot(2,2,4)
%%Distribution of Percentage Reduction in Looking Times
%%(Only Trials with encoding>500ms)
hist(ltper1,30)
xlabel('Reduction in Looking Time [%]', 'FontSize',10,'FontName', 'Arial')
ylabel('Count', 'FontSize',10,'FontName', 'Arial')
text(-5,30,strcat('med=',num2str(median(ltper1))))

toptitle(name)