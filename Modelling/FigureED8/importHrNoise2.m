function [RAz,NameSt]               =...
    importHrNoise2(fileName,allName,Er,Nr) 

delimiterIn                         =   ' ';
headerlinesIn                       =   0;
A                                   =...
    importdata(fileName,delimiterIn,headerlinesIn);
RAz                                 =   A.data;
RAz(:,4:6)                          =   RAz(:,5:7);
            
NameSt                              =   A.textdata(1:end,1);

for i = 1:length(NameSt)
    for j = 1:length(allName)
        if isequal(allName(j),NameSt(i))
            RAz(i,7:8)              =   [Er(j),Nr(j)];
            break
        end
    end
end

