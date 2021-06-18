function event=rd_event_sacs(sac_path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read all sac files ended with .SAC from a file
%  event=rd_event_sacs(sac_path);
%  
%  caution: sac_path DOES NOT  end with /
%  it can be got from cd
%  
%  by Yin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SACfiles=dir([sac_path '/*.SAC'])
% NofSac=size(SACfiles,1);
% for i=1:NofSac
%     i
%     event(i)=readsac([sac_path '/' SACfiles(i).name]);
% end
event=readsac([sac_path '/*.SAC']);
%event=event';

end
