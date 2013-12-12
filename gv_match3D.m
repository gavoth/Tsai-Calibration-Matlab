

%read next frames

fid1 = fopen('C:\Documents and Settings\gvoth\My Documents\data\cornell\part026um\freq16Hz\ptv1_20kHz_26um_16Hz_0.cam0.dat', 'r');
fid2 = fopen('C:\Documents and Settings\gvoth\My Documents\data\cornell\part026um\freq16Hz\ptv1_20kHz_26um_16Hz_0.cam1.dat', 'r');
fid3 = fopen('C:\Documents and Settings\gvoth\My Documents\data\cornell\part026um\freq16Hz\ptv1_20kHz_26um_16Hz_0.cam2.dat', 'r');

     
c1 = fread(fid1, 7, 'int32')

c2 = fread(fid2, 7, 'int32')

fclose(fid);

2ddat=[2.1 4 0 0 0 0 1; 4. 6 0 0 0 0 2; 5. 7 0 0 0 0 2;4.0 3.2 0 0 0 0 3];

[uvalues uidx]=unique(2ddat);

tmax= max(uvalues); %include all three cameras here.

for t=1, tmax
    ncam1=uidx1
    

