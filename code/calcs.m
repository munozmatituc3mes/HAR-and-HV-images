%%%% calcs %%%%


%%%%% plots images

C11=zeros(28,28);
C12=zeros(28,28);
C13=zeros(28,28);
C14=zeros(28,28);
C15=zeros(28,28);
for iter=1:length(labelsc1)
if (grp2idx(labelsc1(iter))==1)
C11=C11+C1(:,:,1,iter);
end
if (grp2idx(labelsc1(iter))==2)
C12=C12+C1(:,:,1,iter);
end
if (grp2idx(labelsc1(iter))==3)
C13=C13+C1(:,:,1,iter);
end
if (grp2idx(labelsc1(iter))==4)
C14=C14+C1(:,:,1,iter);
end
if (grp2idx(labelsc1(iter))==5)
C15=C15+C1(:,:,1,iter);
end
end

figure(1)
contour(C11')
title("Average contour plot for walking for participant 1")
xlabel("horizontal acceleration")
ylabel("vertical acceleration")
figure(2)
contour(C12')
title("Average contour plot for running for participant 1")
xlabel("horizontal acceleration")
ylabel("vertical acceleration")
figure(3)
contour(C13')
title("Average contour plot for jumping for participant 1")
xlabel("horizontal acceleration")
ylabel("vertical acceleration")
figure(4)
contour(C14')
title("Average contour plot for climbing up for participant 1")
xlabel("horizontal acceleration")
ylabel("vertical acceleration")
figure(5)
contour(C15')
title("Average contour plot for climbing down for participant 1")
xlabel("horizontal acceleration")
ylabel("vertical acceleration")



C21=zeros(28,28);
C22=zeros(28,28);
C23=zeros(28,28);
C24=zeros(28,28);
C25=zeros(28,28);
for iter=1:length(labelsc2)
if (grp2idx(labelsc2(iter))==1)
C21=C21+C2(:,:,1,iter);
end
if (grp2idx(labelsc2(iter))==2)
C22=C22+C2(:,:,1,iter);
end
if (grp2idx(labelsc2(iter))==3)
C23=C23+C2(:,:,1,iter);
end
if (grp2idx(labelsc2(iter))==4)
C24=C24+C2(:,:,1,iter);
end
if (grp2idx(labelsc2(iter))==5)
C25=C25+C2(:,:,1,iter);
end
end

figure(6)
contour(C21')
title("Average contour plot for walking for participant 2")
xlabel("horizontal acceleration")
ylabel("vertical acceleration")
figure(7)
contour(C22')
title("Average contour plot for running for participant 2")
xlabel("horizontal acceleration")
ylabel("vertical acceleration")
figure(8)
contour(C23')
title("Average contour plot for jumping for participant 2")
xlabel("horizontal acceleration")
ylabel("vertical acceleration")
figure(9)
contour(C24')
title("Average contour plot for climbing up for participant 2")
xlabel("horizontal acceleration")
ylabel("vertical acceleration")
figure(10)
contour(C25')
title("Average contour plot for climbing down for participant 2")
xlabel("horizontal acceleration")
ylabel("vertical acceleration")


figure(1)
mesh(C11')
figure(2)
mesh(C12')
figure(3)
mesh(C13')
figure(4)
mesh(C14')
figure(5)
mesh(C15')




%%%%%%%%%%%%%% init data

clear ga_x ga_y ga_z;
clear ca_x ca_y ca_z;
clear cga_x cga_y cga_z;
clear coga_x coga_y coga_z;
clear cga coga;

%%%%%% Compensating g %%%%%%

for iter=125:length(attr_x)-125
ga_x(iter-124) = mean(attr_x(iter-124:iter+124));
ga_y(iter-124) = mean(attr_y(iter-124:iter+124));
ga_z(iter-124) = mean(attr_z(iter-124:iter+124));
ca_x(iter-124) = attr_x(iter)-ga_x(iter-124);
ca_y(iter-124) = attr_y(iter)-ga_y(iter-124);
ca_z(iter-124) = attr_z(iter)-ga_z(iter-124);
gmod=sqrt(ga_x(iter-124)^2+ga_y(iter-124)^2+ga_z(iter-124)^2);
cga_x = ca_x(iter-124)*ga_x(iter-124)/gmod;
cga_y = ca_y(iter-124)*ga_y(iter-124)/gmod;
cga_z = ca_z(iter-124)*ga_z(iter-124)/gmod;
cga(iter-124)=cga_x+cga_y+cga_z;
coga(iter-124)=sqrt((ca_x(iter-124)^2+ca_y(iter-124)^2+ca_z(iter-124)^2)-cga(iter-124)^2);
end

ga = sqrt(ga_x.^2+ga_y.^2+ga_z.^2);
ga(isnan(ga)) = 0;
mga = mean(ga);

cga=cga/mga;
cga(isnan(cga))=0;
coga=coga/mga;
coga(isnan(coga))=0;
clear attr_x attr_y attr_z;

awv = cga; 
awh = coga; 

arv = cga; 
arh = coga; 

ajv = cga; 
ajh = coga; 

auv = cga; 
auh = coga; 

adv = cga; 
adh = coga; 


av = awv;
ah = awh;
lables = ones(1,length(awh));

av = [av arv];
ah = [ah arh];
lables = [lables 2*ones(1,length(arh))];

av = [av ajv];
ah = [ah ajh];
lables = [lables 3*ones(1,length(ajh))];

av = [av auv];
ah = [ah auh];
lables = [lables 4*ones(1,length(auh))];

av = [av adv];
ah = [ah adh];
lables = [lables 5*ones(1,length(adh))];


index=1;
for iter=26:length(av)-25
if (sum(ah(iter-25:iter+25))>10)
av(index)=av(iter);
ah(index)=ah(iter);
index=index+1;
end
end

clear C labels imgx imgy img;

index = 1;
for iter=1:125:length(lables)-250
imgx=ah(iter:iter+250); 
imgy=av(iter:iter+250);
imgx(imgx>2)=2;
imgy(imgy>2)=2;
imgy(imgy<-2)=-2;
imgx=floor(13.9*imgx)+1;
imgy=floor(6.9*imgy)+15;
img=zeros(28,28);
for it = 1:length(imgx)
img(imgx(it), imgy(it)) = img(imgx(it), imgy(it))+1;
end
C(:,:,1,index) =img;
labels(index)=lables(iter);
index=index+1;
end

labelsc = categorical(labels);

labelsc1 = labelsc;
C1 = C;


labelsc=cat(2,labelsc1,labelsc2);
C=cat(4,C1,C2);
labelsc=cat(2,labelsc,labelsc3);
C=cat(4,C,C3);
%%% labelsc=cat(2,labelsc,labelsc4);
%%% C=cat(4,C,C4);
labelsc=cat(2,labelsc,labelsc5);
C=cat(4,C,C5);
labelsc=cat(2,labelsc,labelsc6);
C=cat(4,C,C6);
labelsc=cat(2,labelsc,labelsc7);
C=cat(4,C,C7);
labelsc=cat(2,labelsc,labelsc8);
C=cat(4,C,C8);
labelsc=cat(2,labelsc,labelsc9);
C=cat(4,C,C9);
labelsc=cat(2,labelsc,labelsc10);
C=cat(4,C,C10);
labelsc=cat(2,labelsc,labelsc11);
C=cat(4,C,C11);
labelsc=cat(2,labelsc,labelsc12);
C=cat(4,C,C12);
labelsc=cat(2,labelsc,labelsc13);
C=cat(4,C,C13);
labelsc=cat(2,labelsc,labelsc14);
C=cat(4,C,C14);
%%% labelsc=cat(2,labelsc,labelsc15);
%%% C=cat(4,C,C15);


%%%%%%%%%%%%%% end init data

layers = [
    imageInputLayer([28 28 1])
    
    convolution2dLayer(3,8,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(5)
    softmaxLayer
    classificationLayer];

options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.2,...
    'LearnRateDropPeriod',25,...
    'MaxEpochs',60,...
    'MiniBatchSize',128,...
    'Plots','training-progress');

neth2t = trainNetwork(C, labelsc, layers, options);


[YPredH, scores] = classify(neth2t,C);
YValidation = labelsc;
confmat2t = confusionmat(YPredH, YValidation)
confmat2t = confmat2t';
for i=1:5
sens2t(i)=confmat2t(i,i)/sum(confmat2t(i,:)); 
espec2t(i)=confmat2t(i,i)/sum(confmat2t(:,i)); 
end;
sens2t
sens2s
sens2
espec2t
espec2s
espec2


clear Cfilt labelscfilt;
index=1;
for iter=2:length(labelsc)-1
if (grp2idx(labelsc(iter))==grp2idx(labelsc(iter+1))) & (grp2idx(labelsc(iter))==grp2idx(labelsc(iter-1)))
Cfilt(:,:,1,index)=C(:,:,1,iter);
labelscfilt(index)=labelsc(iter);
index=index+1;
end
end

%%%%%%% 10-fold cross validation

clear Cfilttest labelscfilttest Cfiltval labelscfiltval;
indextest=1;
indexval=1;
for iter=1:length(labelscfilt)
if (mod(iter,20) == 9 | mod(iter,20) == 0)
Cfiltval(:,:,1,indexval)=Cfilt(:,:,1,iter);
labelscfiltval(indexval)=labelscfilt(iter);
indexval=indexval+1;
else
Cfilttest(:,:,1,indextest)=Cfilt(:,:,1,iter);
labelscfilttest(indextest)=labelscfilt(iter);
indextest=indextest+1;
end
end

%%%%%% training 10-fold

nethfilt10fold5 = trainNetwork(Cfilttest, labelscfilttest, layers, options);


[YPredH, scores] = classify(nethfilt10fold5,Cfiltval);
YValidation = labelscfiltval;
confmat10fold5 = confusionmat(YPredH, YValidation);
confmat10fold5 = confmat10fold5';
for i=1:5
sens10fold5(i)=confmat10fold5(i,i)/sum(confmat10fold5(i,:)); 
espec10fold5(i)=confmat10fold5(i,i)/sum(confmat10fold5(:,i)); 
end;
confmat10fold5
sens10fold5
espec10fold5


confmat10fold=confmat10fold1+confmat10fold2+confmat10fold3+confmat10fold4+confmat10fold5;
for i=1:5
sens10fold(i)=confmat10fold(i,i)/sum(confmat10fold(i,:)); 
espec10fold(i)=confmat10fold(i,i)/sum(confmat10fold(:,i)); 
end;
confmat10fold
sens10fold
espec10fold


%%%%%% training leave one out 

nethfiltleaveoneout = trainNetwork(Cfilt, labelscfilt, layers, options);


[YPredH, scores] = classify(neth2t,Cfilt);
YValidation = labelscfilt;
confmat2tf = confusionmat(YPredH, YValidation)
confmat2tf = confmat2tf';
for i=1:5
sens2tf(i)=confmat2tf(i,i)/sum(confmat2tf(i,:)); 
espec2tf(i)=confmat2tf(i,i)/sum(confmat2tf(:,i)); 
end;
sens2tf
sens2t
espec2tf
espec2t


for iter=1:length(labelscfilt)
if (grp2idx(labelscfilt(iter))==1) 
walk(iter)=1;
else
walk(iter)=0;
end
end

figure(1)
subplot(2,1,1)       % add first plot in 2 x 1 grid
result=scores(:,1)'.*walk;
plot(result)
title('right')

subplot(2,1,2)       % add second plot in 2 x 1 grid
result=scores(:,1)'.*(1-walk);
plot(result)
title('wrong')

%%%% filtering
clear walk walking;
for iter=7:length(labelscfilt)-6
walking(iter)=sum(scores(iter-6:iter+6,1));
if (grp2idx(labelscfilt(iter))==1) 
walk(iter)=1;
else
walk(iter)=0;
end
end

figure(1)
subplot(2,1,1)       % add first plot in 2 x 1 grid
result=walking.*walk;
plot(result)
title('right-walking')
resultwalking=result;

subplot(2,1,2)       % add second plot in 2 x 1 grid
result=walking.*(1-walk);
plot(result)
title('wrong-walking')
resultbadwalking=result;

clear run running;
for iter=7:length(labelscfilt)-6
running(iter)=sum(scores(iter-6:iter+6,2));
if (grp2idx(labelscfilt(iter))==2) 
run(iter)=1;
else
run(iter)=0;
end
end

figure(2)
subplot(2,1,1)       % add first plot in 2 x 1 grid
result=running.*run;
plot(result)
title('right-running')
resultrunning=result;

subplot(2,1,2)       % add second plot in 2 x 1 grid
result=running.*(1-run);
plot(result)
title('wrong-running')
resultbadrunning=result;

clear jump jumping;
for iter=7:length(labelscfilt)-6
jumping(iter)=sum(scores(iter-6:iter+6,3));
if (grp2idx(labelscfilt(iter))==3) 
jump(iter)=1;
else
jump(iter)=0;
end
end

figure(3)
subplot(2,1,1)       % add first plot in 2 x 1 grid
result=jumping.*jump;
plot(result)
title('right-jumping')
resultjumping=result;

subplot(2,1,2)       % add second plot in 2 x 1 grid
result=jumping.*(1-jump);
plot(result)
title('wrong-jumping')
resultbadjumping=result;

clear up upping;
for iter=7:length(labelscfilt)-6
upping(iter)=sum(scores(iter-6:iter+6,4));
if (grp2idx(labelscfilt(iter))==4) 
up(iter)=1;
else
up(iter)=0;
end
end

figure(4)
subplot(2,1,1)       % add first plot in 2 x 1 grid
result=upping.*up;
plot(result)
title('right-upping')
resultupping=result;

subplot(2,1,2)       % add second plot in 2 x 1 grid
result=upping.*(1-up);
plot(result)
title('wrong-upping')
resultbadupping=result;

clear down downing;
for iter=7:length(labelscfilt)-6
downing(iter)=sum(scores(iter-6:iter+6,5));
if (grp2idx(labelscfilt(iter))==5) 
down(iter)=1;
else
down(iter)=0;
end
end

figure(5)
subplot(2,1,1)       % add first plot in 2 x 1 grid
result=downing.*down;
plot(result)
title('right-downing')
resultdowning=result;

subplot(2,1,2)       % add second plot in 2 x 1 grid
result=downing.*(1-down);
plot(result)
title('wrong-downing')
resultbaddowning=result;


%%%%%% one minute classification

clear cat;
[YPredH, scores] = classify(neth2t,Cfilt);
YValidation = labelscfilt;
for iter=7:length(labelscfilt)-6
cat(iter-6)=mode(labelscfilt(iter-6:iter+6));
catpred(iter-6)=mode(YPredH(iter-6:iter+6));
end

confmat3tf = confusionmat(catpred, cat);
confmat3tf = confmat3tf'
for i=1:5
sens3tf(i)=confmat3tf(i,i)/sum(confmat3tf(i,:)); 
espec3tf(i)=confmat3tf(i,i)/sum(confmat3tf(:,i)); 
end;

sens3tf
sens2tf
espec3tf
espec2tf

%%%% coulored acceleration images

%%%%%%%%%%%%%% init data

clear ga_x ga_y ga_z;
clear ca_x ca_y ca_z;
clear cga_x cga_y cga_z;
clear coga_x coga_y coga_z;
clear cga coga;

%%%%%% Compensating g %%%%%%

for iter=125:length(attr_x)-125
ga_x(iter-124) = mean(attr_x(iter-124:iter+124));
ga_y(iter-124) = mean(attr_y(iter-124:iter+124));
ga_z(iter-124) = mean(attr_z(iter-124:iter+124));
ca_x(iter-124) = attr_x(iter)-ga_x(iter-124);
ca_y(iter-124) = attr_y(iter)-ga_y(iter-124);
ca_z(iter-124) = attr_z(iter)-ga_z(iter-124);
gmod=sqrt(ga_x(iter-124)^2+ga_y(iter-124)^2+ga_z(iter-124)^2);
cga_x = ca_x(iter-124)*ga_x(iter-124)/gmod;
cga_y = ca_y(iter-124)*ga_y(iter-124)/gmod;
cga_z = ca_z(iter-124)*ga_z(iter-124)/gmod;
cga(iter-124)=cga_x+cga_y+cga_z;
coga(iter-124)=sqrt((ca_x(iter-124)^2+ca_y(iter-124)^2+ca_z(iter-124)^2)-cga(iter-124)^2);
end

ga = sqrt(ga_x.^2+ga_y.^2+ga_z.^2);
ga(isnan(ga)) = 0;
mga = mean(ga);

cga=cga/mga;
cga(isnan(cga))=0;
coga=coga/mga;
coga(isnan(coga))=0;
clear attr_x attr_y attr_z;

awv = cga; 
awh = coga; 

arv = cga; 
arh = coga; 

ajv = cga; 
ajh = coga; 

auv = cga; 
auh = coga; 

adv = cga; 
adh = coga; 


av = awv;
ah = awh;
lables = ones(1,length(awh));

av = [av arv];
ah = [ah arh];
lables = [lables 2*ones(1,length(arh))];

av = [av ajv];
ah = [ah ajh];
lables = [lables 3*ones(1,length(ajh))];

av = [av auv];
ah = [ah auh];
lables = [lables 4*ones(1,length(auh))];

av = [av adv];
ah = [ah adh];
lables = [lables 5*ones(1,length(adh))];


index=1;
for iter=26:length(av)-25
if (sum(ah(iter-25:iter+25))>7)
av(index)=av(iter);
ah(index)=ah(iter);
index=index+1;
end
end

clear Ccoulored labelscoulored imgx imgy img;

index = 1;
for iter=1:25:length(av)-250
imgx=ah(iter:iter+250); 
imgy=av(iter:iter+250);
imgx(imgx>2)=2;
imgy(imgy>2)=2;
imgy(imgy<-2)=-2;
imgx=floor(13.9*imgx)+1;
imgy=floor(6.9*imgy)+15;
img=zeros(28,28,2);
for it = 1:length(imgx)
img(imgx(it), imgy(it),1) = img(imgx(it), imgy(it),1)+1-((it-1)/length(imgx));
img(imgx(it), imgy(it),2) = img(imgx(it), imgy(it),2)+(it/length(imgx));
end
Ccoulored(:,:,:,index) =img;
labelscoulored(index)=lables(iter);
index=index+1;
end

labelsccoulored = categorical(labelscoulored);

labelsccoulored15 = labelsccoulored;
Ccoulored15 = Ccoulored;

labelsccoulored=cat(2,labelsccoulored1,labelsccoulored2);
Ccoulored=cat(4,Ccoulored1,Ccoulored2);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored3);
Ccoulored=cat(4,Ccoulored,Ccoulored3);
%%% labelsccoulored=cat(2,labelsccoulored,labelsccoulored4);
%%% Ccoulored=cat(4,Ccoulored,Ccoulored4);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored5);
Ccoulored=cat(4,Ccoulored,Ccoulored5);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored6);
Ccoulored=cat(4,Ccoulored,Ccoulored6);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored7);
Ccoulored=cat(4,Ccoulored,Ccoulored7);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored8);
Ccoulored=cat(4,Ccoulored,Ccoulored8);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored9);
Ccoulored=cat(4,Ccoulored,Ccoulored9);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored10);
Ccoulored=cat(4,Ccoulored,Ccoulored10);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored11);
Ccoulored=cat(4,Ccoulored,Ccoulored11);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored12);
Ccoulored=cat(4,Ccoulored,Ccoulored12);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored13);
Ccoulored=cat(4,Ccoulored,Ccoulored13);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored14);
Ccoulored=cat(4,Ccoulored,Ccoulored14);
%%% labelsccoulored=cat(2,labelsccoulored,labelsccoulored15);
%%% Ccoulored=cat(4,Ccoulored,Ccoulored15);


layers = [
    imageInputLayer([28 28 2])
    
    convolution2dLayer(3,8,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(5)
    softmaxLayer
    classificationLayer];

options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.2,...
    'LearnRateDropPeriod',25,...
    'MaxEpochs',70,...
    'MiniBatchSize',128,...
    'Plots','training-progress');

nethcoulored = trainNetwork(Ccoulored, labelsccoulored, layers, options);


[YPredHcoulored, scorescoulored] = classify(nethcoulored,Ccoulored);
YValidationcoulored = labelsccoulored;
confmatcoulored = confusionmat(YPredHcoulored, YValidationcoulored)
confmatcoulored = confmatcoulored';
for i=1:5
senscoloured(i)=confmatcoulored(i,i)/sum(confmatcoulored(i,:)); 
especcoloured(i)=confmatcoulored(i,i)/sum(confmatcoulored(:,i)); 
end;
senscoloured
especcoloured


%%%%%%% 10-fold cross validation

clear Ccouloredtest labelsccouloredtest Ccouloredval labelsccouloredval;
indextest=1;
indexval=1;
for iter=1:length(labelsccoulored)
if (mod(iter,20) == 9 | mod(iter,20) == 0)
%%%% if mod(iter,10) == 3
Ccouloredval(:,:,:,indexval)=Ccoulored(:,:,:,iter);
labelsccouloredval(indexval)=labelsccoulored(iter);
indexval=indexval+1;
else
Ccouloredtest(:,:,:,indextest)=Ccoulored(:,:,:,iter);
labelsccouloredtest(indextest)=labelsccoulored(iter);
indextest=indextest+1;
end
end

%%%%%% training 10-fold

nethcoulored10fold5 = trainNetwork(Ccouloredtest, labelsccouloredtest, layers, options);

[YPredH, scores] = classify(nethcoulored10fold5,Ccouloredval);
YValidation = labelsccouloredval;
confmatcoulored10fold5 = confusionmat(YPredH, YValidation);
confmatcoulored10fold5 = confmatcoulored10fold5';
for i=1:5
senscoulored10fold5(i)=confmatcoulored10fold5(i,i)/sum(confmatcoulored10fold5(i,:)); 
especcoulored10fold5(i)=confmatcoulored10fold5(i,i)/sum(confmatcoulored10fold5(:,i)); 
end;
confmatcoulored10fold5
senscoulored10fold5
especcoulored10fold5


confmatcoulored10fold=confmatcoulored10fold1+confmatcoulored10fold2+confmatcoulored10fold3+confmatcoulored10fold4+confmatcoulored10fold5;
for i=1:5
senscoulored10fold(i)=confmatcoulored10fold(i,i)/sum(confmatcoulored10fold(i,:)); 
especcoulored10fold(i)=confmatcoulored10fold(i,i)/sum(confmatcoulored10fold(:,i)); 
end;
confmatcoulored10fold
senscoulored10fold
especcoulored10fold



%%%%%%%%%%%%%%%%%%%%%%%%% leave one out

labelsc=cat(2,labelsc1,labelsc2);
C=cat(4,C1,C2);
labelsc=cat(2,labelsc,labelsc3);
C=cat(4,C,C3);
labelsc=cat(2,labelsc,labelsc4);
C=cat(4,C,C4);
labelsc=cat(2,labelsc,labelsc5);
C=cat(4,C,C5);
labelsc=cat(2,labelsc,labelsc6);
C=cat(4,C,C6);
labelsc=cat(2,labelsc,labelsc7);
C=cat(4,C,C7);
labelsc=cat(2,labelsc,labelsc8);
C=cat(4,C,C8);
labelsc=cat(2,labelsc,labelsc9);
C=cat(4,C,C9);
labelsc=cat(2,labelsc,labelsc10);
C=cat(4,C,C10);
labelsc=cat(2,labelsc,labelsc11);
C=cat(4,C,C11);
labelsc=cat(2,labelsc,labelsc12);
C=cat(4,C,C12);
labelsc=cat(2,labelsc,labelsc13);
C=cat(4,C,C13);
labelsc=cat(2,labelsc,labelsc14);
C=cat(4,C,C14);
%%% labelsc=cat(2,labelsc,labelsc15);
%%% C=cat(4,C,C15);


labelscval=labelsc15;
Cval=C15;


layers = [
    imageInputLayer([28 28 1])
    
    convolution2dLayer(3,8,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(5)
    softmaxLayer
    classificationLayer];

options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.2,...
    'LearnRateDropPeriod',15,...
    'MaxEpochs',50,...
    'MiniBatchSize',128,...
    'Plots','training-progress');

clear Cfilt labelscfilt;
index=1;
for iter=2:length(labelsc)-1
if (grp2idx(labelsc(iter))==grp2idx(labelsc(iter+1))) & (grp2idx(labelsc(iter))==grp2idx(labelsc(iter-1)))
Cfilt(:,:,1,index)=C(:,:,1,iter);
labelscfilt(index)=labelsc(iter);
index=index+1;
end
end

clear Cfiltval labelscfiltval;
index=1;
for iter=2:length(labelscval)-1
if (grp2idx(labelscval(iter))==grp2idx(labelscval(iter+1))) & (grp2idx(labelscval(iter))==grp2idx(labelscval(iter-1)))
Cfiltval(:,:,1,index)=Cval(:,:,1,iter);
labelscfiltval(index)=labelscval(iter);
index=index+1;
end
end

nethfiltleaveoneout = trainNetwork(Cfilt, labelscfilt, layers, options);


[YPredH, scores] = classify(nethfiltleaveoneout, Cfiltval);
YValidation = labelscfiltval;
confmatleaveoneout = confusionmat(YPredH, YValidation);
confmatleaveoneout = confmatleaveoneout';
for i=1:5
sensleaveoneout(i)=confmatleaveoneout(i,i)/sum(confmatleaveoneout(i,:)); 
especleaveoneout(i)=confmatleaveoneout(i,i)/sum(confmatleaveoneout(:,i)); 
end;

confmatleaveoneout15=confmatleaveoneout
sensleaveoneout15=sensleaveoneout
especleaveoneout15=especleaveoneout
nethfiltleaveoneout15=nethfiltleaveoneout;

%%%%% totals
confmatleaveoneout=confmatleaveoneout1+confmatleaveoneout2+confmatleaveoneout3+confmatleaveoneout4+confmatleaveoneout5+confmatleaveoneout6+confmatleaveoneout7+confmatleaveoneout8+confmatleaveoneout9+confmatleaveoneout10+confmatleaveoneout11+confmatleaveoneout12+confmatleaveoneout13+confmatleaveoneout14+confmatleaveoneout15
%%% confmatleaveoneout=confmatleaveoneout1+confmatleaveoneout2+confmatleaveoneout3+confmatleaveoneout5+confmatleaveoneout6+confmatleaveoneout7+confmatleaveoneout8+confmatleaveoneout9+confmatleaveoneout10+confmatleaveoneout11+confmatleaveoneout12+confmatleaveoneout13+confmatleaveoneout14

accuracy=sum(diag(confmatleaveoneout))/sum(sum((confmatleaveoneout)));
accuracy1=sum(diag(confmatleaveoneout1))/sum(sum((confmatleaveoneout1)));
accuracy2=sum(diag(confmatleaveoneout2))/sum(sum((confmatleaveoneout2)));
accuracy3=sum(diag(confmatleaveoneout3))/sum(sum((confmatleaveoneout3)));
accuracy4=sum(diag(confmatleaveoneout4))/sum(sum((confmatleaveoneout4)));
accuracy5=sum(diag(confmatleaveoneout5))/sum(sum((confmatleaveoneout5)));
accuracy6=sum(diag(confmatleaveoneout6))/sum(sum((confmatleaveoneout6)));
accuracy7=sum(diag(confmatleaveoneout7))/sum(sum((confmatleaveoneout7)));
accuracy8=sum(diag(confmatleaveoneout8))/sum(sum((confmatleaveoneout8)));
accuracy9=sum(diag(confmatleaveoneout9))/sum(sum((confmatleaveoneout9)));
accuracy10=sum(diag(confmatleaveoneout10))/sum(sum((confmatleaveoneout10)));
accuracy11=sum(diag(confmatleaveoneout11))/sum(sum((confmatleaveoneout11)));
accuracy12=sum(diag(confmatleaveoneout12))/sum(sum((confmatleaveoneout12)));
accuracy13=sum(diag(confmatleaveoneout13))/sum(sum((confmatleaveoneout13)));
accuracy14=sum(diag(confmatleaveoneout14))/sum(sum((confmatleaveoneout14)));
accuracy15=sum(diag(confmatleaveoneout15))/sum(sum((confmatleaveoneout15)));


accuracy
accuracy1
accuracy2
accuracy3
accuracy4
accuracy5
accuracy6
accuracy7
accuracy8
accuracy9
accuracy10
accuracy11
accuracy12
accuracy13
accuracy14
accuracy15

accuracyvector=[accuracy1
accuracy2
accuracy3
accuracy4
accuracy5
accuracy6
accuracy7
accuracy8
accuracy9
accuracy10
accuracy11
accuracy12
accuracy13
accuracy14
accuracy15]



%%%%%%%%%%%%%%%%%%%%%%%%% leave one out coulored -->

labelsccoulored=cat(2,labelsccoulored1,labelsccoulored2);
Ccoulored=cat(4,Ccoulored1,Ccoulored2);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored3);
Ccoulored=cat(4,Ccoulored,Ccoulored3);
%%% labelsccoulored=cat(2,labelsccoulored,labelsccoulored4);
%%% Ccoulored=cat(4,Ccoulored,Ccoulored4);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored5);
Ccoulored=cat(4,Ccoulored,Ccoulored5);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored6);
Ccoulored=cat(4,Ccoulored,Ccoulored6);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored7);
Ccoulored=cat(4,Ccoulored,Ccoulored7);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored8);
Ccoulored=cat(4,Ccoulored,Ccoulored8);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored9);
Ccoulored=cat(4,Ccoulored,Ccoulored9);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored10);
Ccoulored=cat(4,Ccoulored,Ccoulored10);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored11);
Ccoulored=cat(4,Ccoulored,Ccoulored11);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored12);
Ccoulored=cat(4,Ccoulored,Ccoulored12);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored13);
Ccoulored=cat(4,Ccoulored,Ccoulored13);
labelsccoulored=cat(2,labelsccoulored,labelsccoulored14);
Ccoulored=cat(4,Ccoulored,Ccoulored14);

%{
index=1;
for iter=1:length(labelsccoulored15)
if grp2idx(labelsccoulored15(iter))==1 | grp2idx(labelsccoulored15(iter))==2 | grp2idx(labelsccoulored15(iter))==3 %%% | grp2idx(labelsccoulored15(iter))==5
labelsccoulored15f(index)=labelsccoulored15(iter);
Ccoulored15f(:,:,:,index)=Ccoulored15(:,:,:,iter);
index=index+1;
end
end
%}

labelsccoulored=cat(2,labelsccoulored,labelsccoulored15f);
Ccoulored=cat(4,Ccoulored,Ccoulored15f);


labelsccouloredval=labelsccoulored4;
Ccouloredval=Ccoulored4;


layers = [
    imageInputLayer([28 28 2])
    
    convolution2dLayer(3,8,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(5)
    softmaxLayer
    classificationLayer];

options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.4,...
    'LearnRateDropPeriod',10,...
    'MaxEpochs',25,...
    'MiniBatchSize',128,...
    'Plots','training-progress');

clear Ccouloredfilt labelsccouloredfilt;
index=1;
for iter=2:length(labelsccoulored)-1
if (grp2idx(labelsccoulored(iter))==grp2idx(labelsccoulored(iter+1))) & (grp2idx(labelsccoulored(iter))==grp2idx(labelsccoulored(iter-1)))
Ccouloredfilt(:,:,:,index)=Ccoulored(:,:,:,iter);
labelsccouloredfilt(index)=labelsccoulored(iter);
index=index+1;
end
end

clear Ccouloredfiltval labelsccouloredfiltval;
index=1;
for iter=2:length(labelsccouloredval)-1
if (grp2idx(labelsccouloredval(iter))==grp2idx(labelsccouloredval(iter+1))) & (grp2idx(labelsccouloredval(iter))==grp2idx(labelsccouloredval(iter-1)))
Ccouloredfiltval(:,:,:,index)=Ccouloredval(:,:,:,iter);
labelsccouloredfiltval(index)=labelsccouloredval(iter);
index=index+1;
end
end

nethcouloredfiltleaveoneout = trainNetwork(Ccouloredfilt, labelsccouloredfilt, layers, options);


[YPredH, scores] = classify(nethcouloredfiltleaveoneout, Ccouloredfiltval);
YValidation = labelsccouloredfiltval;
confmatleaveoneout = confusionmat(YPredH, YValidation);
confmatleaveoneout = confmatleaveoneout';
for i=1:5
sensleaveoneout(i)=confmatleaveoneout(i,i)/sum(confmatleaveoneout(i,:)); 
especleaveoneout(i)=confmatleaveoneout(i,i)/sum(confmatleaveoneout(:,i)); 
end;

confmatcouloredleaveoneout4=confmatleaveoneout
senscouloredleaveoneout4=sensleaveoneout
especcouloredleaveoneout4=especleaveoneout
nethcouloredfiltleaveoneout4=nethcouloredfiltleaveoneout;
sum(diag(confmatleaveoneout))/sum(sum((confmatleaveoneout)))


%%%%% totals
confmatcouloredleaveoneout=confmatcouloredleaveoneout1+confmatcouloredleaveoneout2+confmatcouloredleaveoneout3+confmatcouloredleaveoneout4+confmatcouloredleaveoneout5+confmatcouloredleaveoneout6+confmatcouloredleaveoneout7+confmatcouloredleaveoneout8+confmatcouloredleaveoneout9+confmatcouloredleaveoneout10+confmatcouloredleaveoneout11+confmatcouloredleaveoneout12+confmatcouloredleaveoneout13+confmatcouloredleaveoneout14+confmatcouloredleaveoneout15
%%% confmatcouloredleaveoneout=confmatcouloredleaveoneout1+confmatcouloredleaveoneout2+confmatcouloredleaveoneout3+confmatcouloredleaveoneout5+confmatcouloredleaveoneout6+confmatcouloredleaveoneout7+confmatcouloredleaveoneout8+confmatcouloredleaveoneout9+confmatcouloredleaveoneout10+confmatcouloredleaveoneout11+confmatcouloredleaveoneout12+confmatcouloredleaveoneout13+confmatcouloredleaveoneout14

accuracycoulored=sum(diag(confmatcouloredleaveoneout))/sum(sum((confmatcouloredleaveoneout)));
accuracycoulored1=sum(diag(confmatcouloredleaveoneout1))/sum(sum((confmatcouloredleaveoneout1)));
accuracycoulored2=sum(diag(confmatcouloredleaveoneout2))/sum(sum((confmatcouloredleaveoneout2)));
accuracycoulored3=sum(diag(confmatcouloredleaveoneout3))/sum(sum((confmatcouloredleaveoneout3)));
accuracycoulored4=sum(diag(confmatcouloredleaveoneout4))/sum(sum((confmatcouloredleaveoneout4)));
accuracycoulored5=sum(diag(confmatcouloredleaveoneout5))/sum(sum((confmatcouloredleaveoneout5)));
accuracycoulored6=sum(diag(confmatcouloredleaveoneout6))/sum(sum((confmatcouloredleaveoneout6)));
accuracycoulored7=sum(diag(confmatcouloredleaveoneout7))/sum(sum((confmatcouloredleaveoneout7)));
accuracycoulored8=sum(diag(confmatcouloredleaveoneout8))/sum(sum((confmatcouloredleaveoneout8)));
accuracycoulored9=sum(diag(confmatcouloredleaveoneout9))/sum(sum((confmatcouloredleaveoneout9)));
accuracycoulored10=sum(diag(confmatcouloredleaveoneout10))/sum(sum((confmatcouloredleaveoneout10)));
accuracycoulored11=sum(diag(confmatcouloredleaveoneout11))/sum(sum((confmatcouloredleaveoneout11)));
accuracycoulored12=sum(diag(confmatcouloredleaveoneout12))/sum(sum((confmatcouloredleaveoneout12)));
accuracycoulored13=sum(diag(confmatcouloredleaveoneout13))/sum(sum((confmatcouloredleaveoneout13)));
accuracycoulored14=sum(diag(confmatcouloredleaveoneout14))/sum(sum((confmatcouloredleaveoneout14)));
accuracycoulored15=sum(diag(confmatcouloredleaveoneout15))/sum(sum((confmatcouloredleaveoneout15)));


accuracycoulored
accuracycoulored1
accuracycoulored2
accuracycoulored3
accuracycoulored4
accuracycoulored5
accuracycoulored6
accuracycoulored7
accuracycoulored8
accuracycoulored9
accuracycoulored10
accuracycoulored11
accuracycoulored12
accuracycoulored13
accuracycoulored14
accuracycoulored15

accuracycouloredvector=[accuracycoulored
accuracycoulored1
accuracycoulored2
accuracycoulored3
accuracycoulored4
accuracycoulored5
accuracycoulored6
accuracycoulored7
accuracycoulored8
accuracycoulored9
accuracycoulored10
accuracycoulored11
accuracycoulored12
accuracycoulored13
accuracycoulored14
accuracycoulored15]

figure;
plot(accuracycouloredvector);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistical preprocessing and outliers


clear mav mah vav vah;
for iter=26:length(av)-25
mav(iter-25)=mean(av(iter-25:iter+25));
mah(iter-25)=mean(ah(iter-25:iter+25));
vav(iter-25)=var(av(iter-25:iter+25));
vah(iter-25)=var(ah(iter-25:iter+25));
end

%% plot(mah)
%%figure
%%plot(mav)
%%figure
%%plot(vah)
%%figure
%%plot(vav)

i1=1;
i2=1;
i3=1;
i4=1;
i5=1;
clear mav1 mah1 vav1 vah1 mav2 mah2 vav2 vah2 mav3 mah3 vav3 vah3 mav4 mah4 vav4 vah4 mav5 mah5 vav5 vah5;
for iter=126:125:length(av)-125
if lables(iter)==1
mav1(i1)=mean(av(iter-125:iter+125));
mah1(i1)=mean(ah(iter-125:iter+125));
vav1(i1)=var(av(iter-125:iter+125));
vah1(i1)=var(ah(iter-125:iter+125));
i1=i1+1;
end
if lables(iter)==2
mav2(i2)=mean(av(iter-125:iter+125));
mah2(i2)=mean(ah(iter-125:iter+125));
vav2(i2)=var(av(iter-125:iter+125));
vah2(i2)=var(ah(iter-125:iter+125));
i2=i2+1;
end
if lables(iter)==3
mav3(i3)=mean(av(iter-125:iter+125));
mah3(i3)=mean(ah(iter-125:iter+125));
vav3(i3)=var(av(iter-125:iter+125));
vah3(i3)=var(ah(iter-125:iter+125));
i3=i3+1;
end
if lables(iter)==4
mav4(i4)=mean(av(iter-125:iter+125));
mah4(i4)=mean(ah(iter-125:iter+125));
vav4(i4)=var(av(iter-125:iter+125));
vah4(i4)=var(ah(iter-125:iter+125));
i4=i4+1;
end
if lables(iter)==5
mav5(i5)=mean(av(iter-125:iter+125));
mah5(i5)=mean(ah(iter-125:iter+125));
vav5(i5)=var(av(iter-125:iter+125));
vah5(i5)=var(ah(iter-125:iter+125));
i5=i5+1;
end
end

mvv1=mean(vav1);
mvv2=mean(vav2);
mvv3=mean(vav3);
mvv4=mean(vav4);
mvv5=mean(vav5);
mvh1=mean(vah1);
mvh2=mean(vah2);
mvh3=mean(vah3);
mvh4=mean(vah4);
mvh5=mean(vah5);
mmv1=mean(mav1);
mmv2=mean(mav2);
mmv3=mean(mav3);
mmv4=mean(mav4);
mmv5=mean(mav5);
mmh1=mean(mah1);
mmh2=mean(mah2);
mmh3=mean(mah3);
mmh4=mean(mah4);
mmh5=mean(mah5);
vvv1=sqrt(var(vav1));
vvv2=sqrt(var(vav2));
vvv3=sqrt(var(vav3));
vvv4=sqrt(var(vav4));
vvv5=sqrt(var(vav5));
vvh1=sqrt(var(vah1));
vvh2=sqrt(var(vah2));
vvh3=sqrt(var(vah3));
vvh4=sqrt(var(vah4));
vvh5=sqrt(var(vah5));
vmv1=sqrt(var(mav1));
vmv2=sqrt(var(mav2));
vmv3=sqrt(var(mav3));
vmv4=sqrt(var(mav4));
vmv5=sqrt(var(mav5));
vmh1=sqrt(var(mah1));
vmh2=sqrt(var(mah2));
vmh3=sqrt(var(mah3));
vmh4=sqrt(var(mah4));
vmh5=sqrt(var(mah5));

clear C labels imgx imgy img;
index=1;
for iter=126:125:length(av)-125
if lables(iter)==1
if abs(mean(av(iter-125:iter+125))-mmv1)<vmv1 & abs(mean(ah(iter-125:iter+125))-mmh1)<vmh1 & abs(var(av(iter-125:iter+125))-mvv1)<vvv1 & abs(var(ah(iter-125:iter+125))-mvh1)<vvh1

imgx=ah(iter-125:iter+125); 
imgy=av(iter-125:iter+125);
imgx(imgx>2)=2;
imgy(imgy>2)=2;
imgy(imgy<-2)=-2;
imgx=floor(13.9*imgx)+1;
imgy=floor(6.9*imgy)+15;
img=zeros(28,28);
for it = 1:length(imgx)
img(imgx(it), imgy(it)) = img(imgx(it), imgy(it))+1;
end
C(:,:,1,index) =img;
labels(index)=lables(iter);

index=index+1;
end
end

if lables(iter)==2
if abs(mean(av(iter-125:iter+125))-mmv2)<vmv2 & abs(mean(ah(iter-125:iter+125))-mmh2)<vmh2 & abs(var(av(iter-125:iter+125))-mvv2)<vvv2 & abs(var(ah(iter-125:iter+125))-mvh2)<vvh2

imgx=ah(iter-125:iter+125); 
imgy=av(iter-125:iter+125);
imgx(imgx>2)=2;
imgy(imgy>2)=2;
imgy(imgy<-2)=-2;
imgx=floor(13.9*imgx)+1;
imgy=floor(6.9*imgy)+15;
img=zeros(28,28);
for it = 1:length(imgx)
img(imgx(it), imgy(it)) = img(imgx(it), imgy(it))+1;
end
C(:,:,1,index) =img;
labels(index)=lables(iter);
index=index+1;
end
end

if lables(iter)==3
if abs(mean(av(iter-125:iter+125))-mmv3)<vmv3 & abs(mean(ah(iter-125:iter+125))-mmh3)<vmh3 & abs(var(av(iter-125:iter+125))-mvv3)<vvv3 & abs(var(ah(iter-125:iter+125))-mvh3)<vvh3

imgx=ah(iter-125:iter+125); 
imgy=av(iter-125:iter+125);
imgx(imgx>2)=2;
imgy(imgy>2)=2;
imgy(imgy<-2)=-2;
imgx=floor(13.9*imgx)+1;
imgy=floor(6.9*imgy)+15;
img=zeros(28,28);
for it = 1:length(imgx)
img(imgx(it), imgy(it)) = img(imgx(it), imgy(it))+1;
end
C(:,:,1,index) =img;
labels(index)=lables(iter);
index=index+1;
end
end

if lables(iter)==4
if abs(mean(av(iter-125:iter+125))-mmv4)<vmv4 & abs(mean(ah(iter-125:iter+125))-mmh4)<vmh4 & abs(var(av(iter-125:iter+125))-mvv4)<vvv4 & abs(var(ah(iter-125:iter+125))-mvh4)<vvh4
imgx=ah(iter-125:iter+125); 
imgy=av(iter-125:iter+125);
imgx(imgx>2)=2;
imgy(imgy>2)=2;
imgy(imgy<-2)=-2;
imgx=floor(13.9*imgx)+1;
imgy=floor(6.9*imgy)+15;
img=zeros(28,28);
for it = 1:length(imgx)
img(imgx(it), imgy(it)) = img(imgx(it), imgy(it))+1;
end
C(:,:,1,index) =img;
labels(index)=lables(iter);
index=index+1;
end
end

if lables(iter)==5
if abs(mean(av(iter-125:iter+125))-mmv5)<vmv5 & abs(mean(ah(iter-125:iter+125))-mmh5)<vmh5 & abs(var(av(iter-125:iter+125))-mvv5)<vvv5 & abs(var(ah(iter-125:iter+125))-mvh5)<vvh5
imgx=ah(iter-125:iter+125); 
imgy=av(iter-125:iter+125);
imgx(imgx>2)=2;
imgy(imgy>2)=2;
imgy(imgy<-2)=-2;
imgx=floor(13.9*imgx)+1;
imgy=floor(6.9*imgy)+15;
img=zeros(28,28);
for it = 1:length(imgx)
img(imgx(it), imgy(it)) = img(imgx(it), imgy(it))+1;
end
C(:,:,1,index) =img;
labels(index)=lables(iter);
index=index+1;
end
end
end

labelsc = categorical(labels);

labelsc1 = labelsc;
C1 = C;



%%%%%%%%%%%%%%%%% coulored filtered


clear mav mah vav vah;
for iter=26:length(av)-25
mav(iter-25)=mean(av(iter-25:iter+25));
mah(iter-25)=mean(ah(iter-25:iter+25));
vav(iter-25)=var(av(iter-25:iter+25));
vah(iter-25)=var(ah(iter-25:iter+25));
end

%% plot(mah)
%%figure
%%plot(mav)
%%figure
%%plot(vah)
%%figure
%%plot(vav)

i1=1;
i2=1;
i3=1;
i4=1;
i5=1;
clear mav1 mah1 vav1 vah1 mav2 mah2 vav2 vah2 mav3 mah3 vav3 vah3 mav4 mah4 vav4 vah4 mav5 mah5 vav5 vah5;
for iter=126:125:length(av)-125
if lables(iter)==1
mav1(i1)=mean(av(iter-125:iter+125));
mah1(i1)=mean(ah(iter-125:iter+125));
vav1(i1)=var(av(iter-125:iter+125));
vah1(i1)=var(ah(iter-125:iter+125));
i1=i1+1;
end
if lables(iter)==2
mav2(i2)=mean(av(iter-125:iter+125));
mah2(i2)=mean(ah(iter-125:iter+125));
vav2(i2)=var(av(iter-125:iter+125));
vah2(i2)=var(ah(iter-125:iter+125));
i2=i2+1;
end
if lables(iter)==3
mav3(i3)=mean(av(iter-125:iter+125));
mah3(i3)=mean(ah(iter-125:iter+125));
vav3(i3)=var(av(iter-125:iter+125));
vah3(i3)=var(ah(iter-125:iter+125));
i3=i3+1;
end
if lables(iter)==4
mav4(i4)=mean(av(iter-125:iter+125));
mah4(i4)=mean(ah(iter-125:iter+125));
vav4(i4)=var(av(iter-125:iter+125));
vah4(i4)=var(ah(iter-125:iter+125));
i4=i4+1;
end
if lables(iter)==5
mav5(i5)=mean(av(iter-125:iter+125));
mah5(i5)=mean(ah(iter-125:iter+125));
vav5(i5)=var(av(iter-125:iter+125));
vah5(i5)=var(ah(iter-125:iter+125));
i5=i5+1;
end
end

mvv1=mean(vav1);
mvv2=mean(vav2);
mvv3=mean(vav3);
mvv4=mean(vav4);
mvv5=mean(vav5);
mvh1=mean(vah1);
mvh2=mean(vah2);
mvh3=mean(vah3);
mvh4=mean(vah4);
mvh5=mean(vah5);
mmv1=mean(mav1);
mmv2=mean(mav2);
mmv3=mean(mav3);
mmv4=mean(mav4);
mmv5=mean(mav5);
mmh1=mean(mah1);
mmh2=mean(mah2);
mmh3=mean(mah3);
mmh4=mean(mah4);
mmh5=mean(mah5);
vvv1=sqrt(var(vav1));
vvv2=sqrt(var(vav2));
vvv3=sqrt(var(vav3));
vvv4=sqrt(var(vav4));
vvv5=sqrt(var(vav5));
vvh1=sqrt(var(vah1));
vvh2=sqrt(var(vah2));
vvh3=sqrt(var(vah3));
vvh4=sqrt(var(vah4));
vvh5=sqrt(var(vah5));
vmv1=sqrt(var(mav1));
vmv2=sqrt(var(mav2));
vmv3=sqrt(var(mav3));
vmv4=sqrt(var(mav4));
vmv5=sqrt(var(mav5));
vmh1=sqrt(var(mah1));
vmh2=sqrt(var(mah2));
vmh3=sqrt(var(mah3));
vmh4=sqrt(var(mah4));
vmh5=sqrt(var(mah5));


clear Ccoulored labelscoulored imgx imgy img;

index=1;
for iter=126:125:length(av)-125
if lables(iter)==1
if abs(mean(av(iter-125:iter+125))-mmv1)<vmv1 & abs(mean(ah(iter-125:iter+125))-mmh1)<vmh1 & abs(var(av(iter-125:iter+125))-mvv1)<vvv1 & abs(var(ah(iter-125:iter+125))-mvh1)<vvh1

imgx=ah(iter-125:iter+125); 
imgy=av(iter-125:iter+125);
imgx(imgx>2)=2;
imgy(imgy>2)=2;
imgy(imgy<-2)=-2;
imgx=floor(13.9*imgx)+1;
imgy=floor(6.9*imgy)+15;
img=zeros(28,28,2);
for it = 1:length(imgx)
img(imgx(it), imgy(it),1) = img(imgx(it), imgy(it),1)+1-((it-1)/length(imgx));
img(imgx(it), imgy(it),2) = img(imgx(it), imgy(it),2)+(it/length(imgx));
end
Ccoulored(:,:,:,index) =img;
labelscoulored(index)=lables(iter);
index=index+1;
end
end

if lables(iter)==2
if abs(mean(av(iter-125:iter+125))-mmv2)<vmv2 & abs(mean(ah(iter-125:iter+125))-mmh2)<vmh2 & abs(var(av(iter-125:iter+125))-mvv2)<vvv2 & abs(var(ah(iter-125:iter+125))-mvh2)<vvh2

imgx=ah(iter-125:iter+125); 
imgy=av(iter-125:iter+125);
imgx(imgx>2)=2;
imgy(imgy>2)=2;
imgy(imgy<-2)=-2;
imgx=floor(13.9*imgx)+1;
imgy=floor(6.9*imgy)+15;
img=zeros(28,28,2);
for it = 1:length(imgx)
img(imgx(it), imgy(it),1) = img(imgx(it), imgy(it),1)+1-((it-1)/length(imgx));
img(imgx(it), imgy(it),2) = img(imgx(it), imgy(it),2)+(it/length(imgx));
end
Ccoulored(:,:,:,index) =img;
labelscoulored(index)=lables(iter);
index=index+1;
end
end

if lables(iter)==3
if abs(mean(av(iter-125:iter+125))-mmv3)<vmv3 & abs(mean(ah(iter-125:iter+125))-mmh3)<vmh3 & abs(var(av(iter-125:iter+125))-mvv3)<vvv3 & abs(var(ah(iter-125:iter+125))-mvh3)<vvh3

imgx=ah(iter-125:iter+125); 
imgy=av(iter-125:iter+125);
imgx(imgx>2)=2;
imgy(imgy>2)=2;
imgy(imgy<-2)=-2;
imgx=floor(13.9*imgx)+1;
imgy=floor(6.9*imgy)+15;
img=zeros(28,28,2);
for it = 1:length(imgx)
img(imgx(it), imgy(it),1) = img(imgx(it), imgy(it),1)+1-((it-1)/length(imgx));
img(imgx(it), imgy(it),2) = img(imgx(it), imgy(it),2)+(it/length(imgx));
end
Ccoulored(:,:,:,index) =img;
labelscoulored(index)=lables(iter);
index=index+1;
end
end

if lables(iter)==4
if abs(mean(av(iter-125:iter+125))-mmv4)<vmv4 & abs(mean(ah(iter-125:iter+125))-mmh4)<vmh4 & abs(var(av(iter-125:iter+125))-mvv4)<vvv4 & abs(var(ah(iter-125:iter+125))-mvh4)<vvh4
imgx=ah(iter-125:iter+125); 
imgy=av(iter-125:iter+125);
imgx(imgx>2)=2;
imgy(imgy>2)=2;
imgy(imgy<-2)=-2;
imgx=floor(13.9*imgx)+1;
imgy=floor(6.9*imgy)+15;
img=zeros(28,28,2);
for it = 1:length(imgx)
img(imgx(it), imgy(it),1) = img(imgx(it), imgy(it),1)+1-((it-1)/length(imgx));
img(imgx(it), imgy(it),2) = img(imgx(it), imgy(it),2)+(it/length(imgx));
end
Ccoulored(:,:,:,index) =img;
labelscoulored(index)=lables(iter);
index=index+1;
end
end

if lables(iter)==5
if abs(mean(av(iter-125:iter+125))-mmv5)<vmv5 & abs(mean(ah(iter-125:iter+125))-mmh5)<vmh5 & abs(var(av(iter-125:iter+125))-mvv5)<vvv5 & abs(var(ah(iter-125:iter+125))-mvh5)<vvh5
imgx=ah(iter-125:iter+125); 
imgy=av(iter-125:iter+125);
imgx(imgx>2)=2;
imgy(imgy>2)=2;
imgy(imgy<-2)=-2;
imgx=floor(13.9*imgx)+1;
imgy=floor(6.9*imgy)+15;
img=zeros(28,28,2);
for it = 1:length(imgx)
img(imgx(it), imgy(it),1) = img(imgx(it), imgy(it),1)+1-((it-1)/length(imgx));
img(imgx(it), imgy(it),2) = img(imgx(it), imgy(it),2)+(it/length(imgx));
end
Ccoulored(:,:,:,index) =img;
labelscoulored(index)=lables(iter);
index=index+1;
end
end
end

labelsccoulored = categorical(labelscoulored);

labelsccoulored2 = labelsccoulored;
Ccoulored2 = Ccoulored;


labelsccoulored1bu = labelsccoulored1;
Ccoulored1bu = Ccoulored1;

labelsccoulored1=labelsccoulored1bu;
Ccoulored1=Ccoulored1bu;


%%%%%%%%%%%%%%%%%% Data quality

index=1;
for iter=1:length(lables)
if lables(iter)==1
ahplot(index)=ah(iter);
avplot(index)=av(iter);
index=index+1;
end
end

figure(1)
plot(ahplot)
figure(2)
plot(avplot)


index=1;
for iter=1:length(lables)
if lables(iter)==2
ahplot(index)=ah(iter);
avplot(index)=av(iter);
index=index+1;
end
end

figure(3)
plot(ahplot)
figure(4)
plot(avplot)


index=1;
for iter=1:length(lables)
if lables(iter)==3
ahplot(index)=ah(iter);
avplot(index)=av(iter);
index=index+1;
end
end

figure(5)
plot(ahplot)
figure(6)
plot(avplot)


index=1;
for iter=1:length(lables)
if lables(iter)==4
ahplot(index)=ah(iter);
avplot(index)=av(iter);
index=index+1;
end
end

figure(7)
plot(ahplot)
figure(8)
plot(avplot)


index=1;
for iter=1:length(lables)
if lables(iter)==5
ahplot(index)=ah(iter);
avplot(index)=av(iter);
index=index+1;
end
end

figure(9)
plot(ahplot)
figure(10)
plot(avplot)

