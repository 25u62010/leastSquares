clear all;
p=2^10-1;
ms=idinput([p,1,5],'prbs',[0,1],[-1,1]);
mSequences(:,1)=0.002:0.002:2.046*5;
mSequences(:,2)=ms;
out=sim("simuModel",10);
simuOut=out.simuOut;

fid = fopen('dataGenerator.txt','w');  
[r,c]=size(simuOut);  
for i=1:r  
    for j=1:c  
        fprintf(fid,'%5f ',simuOut(i,j));  
    end  
    fprintf(fid,'\r\n');  
end  
fclose(fid); 