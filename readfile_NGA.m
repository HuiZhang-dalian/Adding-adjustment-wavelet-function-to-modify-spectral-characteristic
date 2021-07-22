function [info,Data] = readfile_NGA(FilePath)
fid = fopen(FilePath);
for i=1:4
    str0=fgetl(fid);
end
 str1=regexp(str0,'=','split');
 str2=cell2mat(str1(3));
 ns=strfind(str2,'S');
 dt=str2num(str2(1:ns-1));
 info.NPTS=str1;
 info.dt=dt;
 data0=textscan(fid,'%f','headerlines',0);
 Data=(cell2mat(data0));
fclose(fid);
 
 