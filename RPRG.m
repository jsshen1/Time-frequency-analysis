function [findex,interset] = RPRG(findex,thrf)
% Ridge path regrouping (RPRG)

%%%%%%%%%%%%%%%%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% findex：detected frequency indexs of the ridge curves
% thrf: threshold for finding the intersection intervals,i.e., D_f in the paper

%%%%%%%%%%%%%%%%%%%%%%% output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% findex：the frequency indexs of the regrouped ridge curves
% interset：the detected intersection intervals for the ridge curves


%% find all the intersection intervals
[multi,N] = size(findex);%multin denotes the number of the components；N is the number of the samples

id = nchoosek(1:multi,2);%
[numcp,~] = size(id);

numiter = 1;
interset = {};
for i = 1:numcp 
   tempindex = find( abs(findex(id(i,1),:) - findex(id(i,2),:)) <= thrf) ;
   if (~isempty(tempindex)) 
      tempindex2 = find ((diff(tempindex)) >= N/70); % find if there exist multiple intersection intervals for two overlapped ridge curves
      if (~isempty(tempindex2)) %
       for j = 1:length(tempindex2)+1 %
           if j==1
               interset(numiter,1) = {[id(i,1),id(i,2)]};% ridge numbers of the two intersected ridge curves
               interset(numiter,2) = {[tempindex(1),tempindex(tempindex2(1))]}; % starting  and ending time instants for the intersection intervals
               numiter = numiter+1;%
           else if j==length(tempindex2)+1
               interset(numiter,1) = {[id(i,1),id(i,2)]};%
               interset(numiter,2) = {[tempindex(tempindex2(end)+1),tempindex(end)]}; %
               numiter = numiter+1;%
               else
               interset(numiter,1) = {[id(i,1),id(i,2)]};%
               interset(numiter,2) = {[tempindex(tempindex2(j-1)+1),tempindex(tempindex2(j))]}; %
               numiter = numiter+1;%
               end
           end
       end
      else 
          interset(numiter,1) = {[id(i,1),id(i,2)]};%
          interset(numiter,2) = {[tempindex(1),tempindex(end)]}; %
          numiter = numiter+1;%
      end
   end
end

%% merge intersection intervals which are close to each other
if (~isempty(interset)) %
    [multinter,~] = size(interset);%
    
    if multinter>2  % if two more intersection intervals, then merge close intervals
        
    id = nchoosek(1:multinter,2);%
    [numcp,~] = size(id);
    %%%%%%%%%%%%%%%%%%%%%
    for i = 1:numcp
        %交点1
        comp1 = interset{id(i,1),1};%
        interv1 = interset{id(i,1),2};%
        locafr1 = mean(findex(comp1,interv1(1)));%
        locati1 = interv1(1);%
        %交点2
        comp2 = interset{id(i,2),1};%
        interv2 = interset{id(i,2),2};%
        locafr2 = mean(findex(comp2,interv2(1)));%
        locati2 = interv2(1);%
        
        if ((locati1 - locati2)^2 + (locafr1 - locafr2)^2) <= 4*thrf^2 %compute the distance of two intersection intervals
            ucomp = union(comp1,comp2);% merge ridge numbers
            interset(id(i,1),1) = {ucomp};%
            uset = union(interv1,interv2);% merge time instants
            interset(id(i,1),2) = {[uset(1),uset(end)]};%
            interset(id(i,2),:) = interset(id(i,1),:);% assign the merging results to the two intersection intervals
        end
    end
    
    % delete the same interesection intervals in the set interset 
    douinter = zeros(multinter,100);%
    for i = 1:multinter
        ddd = cat(2,interset{i,:});  %
        douinter(i,1:length(ddd)) = ddd;
    end
    [~,m,~] = unique(douinter,'rows'); %merge the rows which have the same elements
   dd = setdiff(1:multinter,m);%
   interset(dd,:) = [];% delete
   %%%%% sort the intersection intervals according to the ascending order of the ending time of each interval
   [multinter,~] = size(interset);%
   rightend = zeros(1,multinter);%
   for i = 1:multinter
        rightend(i) = interset{i,2}(2);  %
   end
   [~,sorin] = sort(rightend);%
   for i = 1:multinter
       sortset(i,:) = interset(sorin(i),:);  %
   end
   interset = sortset;%
   
   
    end
   
%% regroup the ridge curves
   len = floor(N/70);%
   for i = 1:multinter
       comp = interset{i,1}; %
       numcomp = length(comp);%
       lefte = interset{i,2}(1); %
       righe = interset{i,2}(2); %
       
      if (lefte-len)<1 %
          slopr =  (findex(comp,righe + len) -  findex(comp,righe))/len;%
          for i = 1: numcomp
            findex(comp(i),1:righe) = round(findex(comp(i),righe) + slopr(i)*((1:righe)-righe));% linear prediction
          end
      end
      
      if (righe + len)>N %
          slopl =  (findex(comp,lefte) -  findex(comp,lefte-len))/len;%
          for i = 1: numcomp
            findex(comp(i),lefte:N) = round(findex(comp(i),lefte) + slopl(i)*((lefte:N)-lefte));% linear prediction
          end
      end
      
      if (lefte-len)>=1 && (righe + len)<=N
      slopl =  findex(comp,lefte) -  findex(comp,lefte-len);%slopes
      slopr =  findex(comp,righe + len) -  findex(comp,righe);%slopes
      slopl = slopl(:); slopr = slopr(:); %
      deltaslop = abs(bsxfun(@minus,slopr,slopl'));% connection matrix
      linearinterpol = zeros(numcomp,righe);%
      for i = 1:numcomp %
          [jj,kk] = find(deltaslop == min(deltaslop(:)));%
          jj = jj(1);
          kk = kk(1);%
          deltaslop(:,kk) = inf;%
          deltaslop(jj,:) = inf;%
          linearinterpol(jj,:) = [findex(comp(kk),1:(lefte-1)),round(linspace(findex(comp(kk),lefte),findex(comp(jj),righe),righe-lefte+1))];%linear interpolation
      end
      for jj = 1:numcomp 
          findex(comp(jj),1:righe) = linearinterpol(jj,:);%
      end
      end
      
   end         
end
   