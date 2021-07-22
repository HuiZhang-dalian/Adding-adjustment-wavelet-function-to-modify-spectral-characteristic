clc;clear;
fid='C:\Users\dell\Desktop\张辉\IDA\ATC63\近场考虑频谱\0.5g\'
fid1='C:\Users\dell\Desktop\张辉\IDA\ATC63\图片\'
filename=dir([fid,'*.txt'])

%% 设置参数
fid2='C:\Users\dell\Desktop\张辉\IDA\ATC63\近场\' %有数据的文件夹
filename1=dir([fid2,'*.AT2']);

for i=1:length(filename)
Data=importdata([fid,filename(i).name])
subplot(2,1,1)
plot(Data(9001:end))


[info,Data1] =readfile_NGA([fid2,filename1(i).name]);
subplot(2,1,2)
plot(Data1)
saveas(gcf,[fid1,filename(i).name,'.jpg'])
end