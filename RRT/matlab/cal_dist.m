function [dist] = cal_dist(p1,p2)
%CAL_DIST �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
dp = p1 - p2;
dist = sqrt(dp(1)^2+dp(2)^2);
end

