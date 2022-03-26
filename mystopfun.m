function stopnow = mystopfun(problem, x, info, last)
    stopnow = (last >= 10 && info(last-3).cost - info(last).cost < 1e-3) || (last>=10);
%     stopnow = (last>=10);
end
% function stopnow = mystopfun(problem, x, info, last)
%     stopnow = (last >= 100 && info(last-9).cost - info(last).cost < 1e-9) || (last>=100);
% %     stopnow = (last>=20);
% end