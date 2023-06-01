%%
clc
clear
load('ArrayData.mat');
load('CleanTrials.mat');
ClrTr = Intersect_Clean_Trials;
clear Intersect_Clean_Trials;
% c = firpmord([.5 1.5],[0 1],[0.01 0.01],200,'cell');
% B = firpm(c{:});

for i = 1:48
    x=chan(i).lfp(:,ClrTr);
    [pxx,f]=pwelch(x,hamming(100),50,[],200);
    pxx = log(pxx);
    pxx2=mean(pxx')';
    p = polyfit(f,pxx2,1);
    pinknoise4 = p(1)*f + p(2);
    pxx2 = pxx2 - pinknoise4;
    pxx3 = exp(pxx2);
    pxx2 = 10*log10(pxx3);
    PSD1(i,:)=pxx2';
    [ma,I1]=max(pxx2);
    if f(I1)<10
        pxx2(I1)=0;
    elseif f(I1)>50
        pxx2(I1)=-10;
    elseif f(I1)==12.5
        if pxx2(I1+1)- pxx2(I1)<.05
            pxx2(I1)=-10;
        end
    end
    [ma,I1]=max(pxx2);
    if f(I1)<10
        pxx2(I1)=-10;
    end
    [ma,I1]=max(pxx2);
    if f(I1)<10
        pxx2(I1)=-10;
    end
    [ma,I1]=max(pxx2);
    frequency(i)=f(I1);
    
    
end



afreq = zeros(5,10);
for i=1:48
    [id,jd] = find(ChannelPosition==i) ;
    if frequency(i)>50
        frequency(i)=frequency(i+3);
    end
    c = firpmord([f(I1)-3 f(I1)-1.5 f(I1)+1.5 f(I1)+3],[0 1 0],[0.1 0.01 0.1],200,'cell');
    Bnotch(i,:) = firpm(c{:});
    afreq(id,jd) = frequency(i) ;
end


count1=1;
for i = 1:4
    figure('WindowState','maximized')
    for j=1:12
        subplot(4,3,j)
        plot(f,PSD1(count1,:))
        xlabel('Frequency(Hz)','Interpreter','latex')
        ylabel('PSD','Interpreter','latex')
        title([' Channel $ ' num2str(count1) ' $ Average PSD'],'Interpreter','latex')
        xline(frequency(count1),'--r',['Frequency = ' num2str(frequency(count1))],'Interpreter','latex')
        grid on
        grid minor
        count1=count1+1;
    end
    saveas(gcf,['Fig',num2str(i),'.png'])
end


figure('WindowState','maximized')
imagesc(afreq)
for i=1:9
    xline(i+.5)
end
for i=1:4
    yline(i+.5)
end
yticks(0:5)
xticks(1:10)
c = colorbar;
c.Label.String = 'Frequency(Hz)';
title('Clusters Of Neurons','Interpreter','latex')
saveas(gcf,'Fig5.png')

PSD2 = PSD1';

figure('WindowState','maximized')
imagesc(PSD2(129:-1:1,:))
c = colorbar;
c.Label.String = 'Magnitude (dB)';
yticks(1:4:129);
yticklabels(num2str(f(129:-4:1)))


xlabel('Channels','Interpreter','latex')
ylabel('Frequency (Hz)','Interpreter','latex')
saveas(gcf,'Fig6.png')
%%
clear pxx
clear pxx2
clear pinknoise4
clear f


%%

s1 = zeros(400,127);
for i = 1:48
    x=chan(i).lfp(:,ClrTr);
    for j=1:490
        x1 = x(:,j);
        x1 = highpass(x1,10,200);

        s = stft(x1,200,'Window',kaiser(10,5),'OverlapLength',5,'FFTLength',400);
        
        s1 = s1 + abs(s);
        
        
        
    end
    
end
%%
figure('WindowState','maximized')
s2 = s1(300:-1:201,:);
imagesc(s2/(490*48))
yticks(0:10:100)
f1 = 100:-1:0;

t2 = linspace(-1.2,2,127);
for i = 1:11
    f2{i} = (10-i+1)*10;
    f3{i} = num2str(t2(12*(i-1)+1));
end
yticklabels(f2)
xticks(1:12:127)
xticklabels(f3)
xlabel('Time(s)','Interpreter','latex')
ylabel('Frequency (Hz)','Interpreter','latex')
title('LFP Signal Power Spectrum over Time','Interpreter','latex')
c = colorbar;
c.Label.String = 'Magnitude (dB)';
saveas(gcf,'Fig7.png')


%%
phase1 = zeros(7,7,490,641);
ind1(1:6,1)=1:6;
ind1(7:13,1)=1:7;
ind1(14:20,1)=1:7;
ind1(21:27,1)=1:7;
ind1(28:34,1)=1:7;
ind1(35:41,1)=1:7;
ind1(42:48,1)=1:7;
ind1(1:6,2)=1;
ind1(7:13,2)=2;
ind1(14:20,2)=3;
ind1(21:27,2)=4;
ind1(28:34,2)=5;
ind1(35:41,2)=6;
ind1(42:48,2)=7;

ind2(1:3,1)=2:4;
ind2(4:8,1)=1:5;
ind2(9:13,1)=1:5;
ind2(14:18,1)=1:5;
ind2(19:23,1)=1:5;
ind2(24:28,1)=1:5;
ind2(29:33,1)=1:5;
ind2(34:38,1)=1:5;
ind2(39:43,1)=1:5;
ind2(44:48,1)=1:5;

ind2(1:3,2)=1;
ind2(4:8,2)=2;
ind2(9:13,2)=3;
ind2(14:18,2)=4;
ind2(19:23,2)=5;
ind2(24:28,2)=6;
ind2(29:33,2)=7;
ind2(34:38,2)=8;
ind2(39:43,2)=9;
ind2(44:48,2)=10;

L=641;
Fs=200;
f = Fs*(0:(L/2))/L;

for j = 1:48
    x=chan(j).lfp(:,ClrTr);
    [pxx,f]=pwelch(x,hamming(100),50,[],200);
    pxx = log(pxx);
    pxx2=mean(pxx')';   
    p = polyfit(f,pxx2,1);
    pinknoise4 = p(1)*f + p(2);
    pinknoise1 = exp(-pinknoise4);
    pinknoise2 = ifft(pinknoise1);
    for i=1:490
        P8=filtfilt(1,abs(pinknoise2),x(:,i));
        P9(:,i)=P8;
    end
    
    
    x2 = filtfilt(Bnotch(j,:),1,P9);
    x3 = hilbert(x2);
    xphase = unwrap(angle(x3));
%     xphase = rad2deg(xphase);
    q1 = ind1(j,1);
    r1 = ind1(j,2);
    
    q2 = ind2(j,1);
    r2 = ind2(j,2);
    
    phase1(q1,r1,:,:) = xphase';
    phase2(q1,r1,:,:) = mod(rad2deg(xphase),360)';
    
end
%%
clear xphase2
for i=1:490
    for j=1:641
        xphase2(:,:) = cos(phase1(:,:,i,j));
        imagesc(xphase2)
        caxis([-1 1])
        colorbar
        title(['The Time is ',num2str((j-1)/200 - 1.2),' s'],'Interpreter','latex')
        pause(.1);
        
    end
end


%%
for i=1:7
    xd(7*(i-1)+1:7*i)=i;
    yd(7*(i-1)+1:7*i)=1:7;
end
x11=ones(1,49);
x4=[x11;xd;yd];
b1=x4*x4';
b2=inv(b1);


c1 = 0;
for i=1:490
    for j=1:641
        ph(:,:) = phase2(:,:,i,j);
        ph(7,1)=(ph(6,1)+ph(7,2))/2;
        ph1 = [ph(1,:),ph(2,:),ph(3,:),ph(4,:),ph(5,:),ph(6,:),ph(7,:)];
        b3=x4*ph1';
        b=b2*b3;
        rs=0;
        for qr=1:5:360
            for zetha=0:.5:18
                zavieh=qr;
                b(2)=zetha*cosd(zavieh);
                b(3)=zetha*sind(zavieh);
                phh1 =mod(b(1) + b(2)*xd + b(3)*yd,360);
                n1 = length(ph1);
                sec1 = (1/n1)*(sum(cosd(ph1-phh1)));
                sec2 = (1/n1)*(sum(sind(ph1-phh1)));
                rs1 = sqrt(sec1^2 + sec2^2);
                if rs1>=rs
                    ZETHA122(i,j)=zetha;
                    rs=rs1;
                    PGD122(i,j)=rs;
                    alpha1(i,j) = atan2d(b(3),b(2));
                    spf(i,j)=zetha;
                    pcc1=sum(sind(ph1-mean(ph1)).*sin(phh1-mean(phh1)));
                    pcc2 = sqrt(sum((sind(ph1-mean(ph1))).^2) * sum((sin(phh1-mean(phh1))).^2));
                    pcc=pcc1/pcc2;
                    PGD(i,j) = 1 - ((1-pcc^2)*47)/44;
                end
            end
        end
            
        
    end
    
end


%%
for i=1:490
    for j=1:641
        xphase2(:,:) = cos(phase1(:,:,i,j));
        imagesc(xphase2)
        caxis([-1 1])
        colorbar
        hold on
        plot([4,cosd(90+alpha1(i,j))+4],[4,sind(90+alpha1(i,j))+4],'r->','LineWidth',2)
        title(['The Time is ',num2str((j-1)/200 - 1.2),' s'],'Interpreter','latex')
        
        pause(.01);
        
    end
end

%%
figure('WindowState','maximized')
t1 = -1.2:.005:2;

plot(t1,mean(PGD122))
xlim([-1.2 2])
grid on
grid minor

title('Goodness of Fit over Time','Interpreter','latex')
xlabel('Time(s)','Interpreter','latex')
saveas(gcf,'Fig8.png')

%%
c1 = 0;
f2 = diff(alpha1')';

for i = 1:490
    for j = 2:641
        c1 = c1 + 1;
        f = cos(phase1(:,:,i,j));
        [dfdx,dfdy] = gradient(f) ;
        M1 = norm(dfdx);
        M2 = norm(dfdy);
        AG = (M1 + M2)/2 ; 
        f1 =  cos(phase1(:,:,i,j-1));
        A2 = abs(f1-f);
        A3 = mean(A2,"all")*100;
        v3 = A3/AG;
        
        v1(i,j) = v3;

    end
end

v4 = v1(:,2:end);
%%
figure('WindowState','maximized')

plot(Time,mean(v1,1))
xlabel('Time(s)','Interpreter','latex')
ylabel('Speed(cm/s)','Interpreter','latex')
title('Speed Over Time','Interpreter','latex')
xlim([Time(1) Time(end)])
grid on
grid minor
saveas(gcf,'Fig9.png')


%%

alpha2=mod(alpha1+90,360);
for i=1:641
    xce=alpha2(:,i);
    polarhistogram(xce,36)
    title(['The Time is ',num2str((i-1)/200 - 1.2),' s'],'Interpreter','latex')
    pause(.1);
end
polarhistogram(alpha2,36)

%%
for i=1:490
    spf1((i-1)*641+1:i*641)=spf(i,:);
    PGD2((i-1)*641+1:i*641) = PGD(i,:);
end
figure
hist(spf1)
title('Spatial Frequency Histogram')
figure
hist(PGD2)
title("PGD Histogram")
