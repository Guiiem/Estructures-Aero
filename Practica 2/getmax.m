function [max_sig, min_sig] = getmax(sig_t)

max_sig = zeros(size(sig_t,3),1);
min_sig = zeros(size(sig_t,3),1);
for i = 1:size(sig_t,3)
    max_sig(i,1) = max(sig_t(:,1,i));
    min_sig(i,1) = min(sig_t(:,1,i));
end
