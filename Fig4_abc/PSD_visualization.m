clc;clear all;
A={};stim={};sham={};
load('C:\视频任务数据\论文数据整理\task2 data and code\Fig4_abc\Results\10hz_relative_psd_power.mat', 'psd_full')
A{2}=squeeze(mean(psd_full(2,2,:,:,:,:),4));
A{1}=squeeze(mean(psd_full(1,2,:,:,:,:),4));
a=squeeze(sum(A{2},2));b=squeeze(sum(A{1},2));
stim{1}=a;stim{2}=b;delta_stim=(a-b)./abs(b);

load('C:\视频任务数据\论文数据整理\task2 data and code\Fig4_abc\Results\sham_relative_psd_power.mat', 'psd_full')
A{2}=squeeze(mean(psd_full(2,2,:,:,:,:),4));
A{1}=squeeze(mean(psd_full(1,2,:,:,:,:),4));
a=squeeze(sum(A{2},2));b=squeeze(sum(A{1},2));
sham{1}=a;sham{2}=b;delta_sham=(a-b)./abs(b);

hz=0.5:0.5:90;
figure('Position',[10,10,800,400]),
shadedErrorBar(hz,delta_stim,{@mean,@(x) std(x)/sqrt(size(delta_stim,1))},'lineprops',{'r'});
hold on;
shadedErrorBar(hz,delta_sham,{@mean,@(x) std(x)/sqrt(size(delta_sham,1))},'lineprops',{'b'});