function [Thr] = findoptRth(calcium, Opts)
%This code allows you to find the optimal threshold for a network that both
%optimized fit to a scalefree graph and constraints for the average number of
%connections per islet
%Input: calcium
%Opts: Method - either 'Scale-Free' or 'Degree'
%       if Scale-Free you also need Max and Min which correspond to max and
%       min average degree
%       if Degree you also need avDeg which the average degree of the islet


%Output: Threshold
% Jennifer Briggs, 2021

%check calcium direction
try
direction_hold = Opts.direction_hold;
if direction_hold ==0
  if size(calcium,1)<size(calcium,2)
     calcium = calcium';
  end
end
catch
  if size(calcium,1)<size(calcium,2)
     calcium = calcium';
  end
end

calcium = rescale(calcium, 0, 1);
[Rij, pval]=corr(calcium);
Rij = Rij-diag(diag(Rij)); %removes ones on diagonal
done = 0; %if the computer can't find a threshold with standard floating point, we increase the variable precision
increaseprecision = 0;

while done == 0
    if increaseprecision == 1
        Rij = vpa(Rij);
    end
        
switch Opts.Method  
    case 'Scale-Free'
        try 
            upper = Opts.Max;
        catch
            upper = 8 %set average number of connections for the top "upper" cells
        end


        fun = @(x)Lfunc(x,calcium,Rij); % Defines the cost function
        numcells = size(calcium,2); 
        Rijs = sort(Rij);
        av = mean(Rijs(end-upper-1:end,:)); %calculates the average correlation coeffiecient for top "upper" highest correlated cell pairs [9 x numcells]
        avs = sort(av); %Gives average correlation for top 9 highest correlated cellpairs
        ub = mean(avs(end-upper-1:end)); %Upper bound is set to be the average correlation for the top "upper"
        lb = mean(avs(1:2:end))
        x0 =  mean(avs(end-round(numcells/3):end)); %originial threshold


        %run optimization: fminsearchbnd can be downloaded here https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon
        try
        Thr = fminsearchbnd(fun,x0, lb, ub); %Find threshold that most fits power law while staying within upper and lower bounds
        catch
            disp('Threshold could not be found')
            Thr = mean([ub, x0]);
        end
    case 'Degree'
%     This function finds the threshold for a corrlation matrix (cor_mat) that gives the network a desired average degree k.
%     k = desired average degree
%     m = k*n/2 : where m is the total number of edges and n is the total number of nodes
%     p = 2m/(n(n-1)): the percent of edges we want compared to the number of possible edges
%     Once p is found, we sort all non-self correlations and find the threshold which returns p percent.
        k = Opts.avDeg;
        Rij = Rij - diag(diag(Rij)); %set diagonals to zeros
        all_cors = reshape(Rij, 1, []); % Make 1d list of all correlations
        all_cors_sort = sort(all_cors, 'descend');
        m = (k*length(Rij)); %desired number of edges
        Thr = vpa(all_cors_sort(ceil(m)))
end

if mean2(Rij>double(Thr)) > 0.75 %if most of the cells are not connected...
    increaseprecision = 1
    done = 0
else
    done = 1
end
end
end




function err = Lfunc(Threshold,calciumT,Rij)

numcells = size(calciumT,2);

Adj = Rij;
Adj(Adj >= Threshold) = 1;
Adj(Adj < Threshold) = 0;
% 
 Adj = Adj - diag(diag(Adj));        
 if mean(nonzeros(sum(Adj))) > 5 %here is where the lower bound is set!
%% 4. Determine number of "links" based on cov threshold
for i=1:numcells
    N (i,1) = nnz(Adj(:,i));  % N is matrix containing # of links for each cell (nnz returns number of nonzero elements)
end
% 


% %% 5. Creating a "probability of a cell to have k links" matrix

histArray=zeros(size(N))'; % prealocate
% a forloop to count how many times you encounter in a particular value:
for n=1:length(N)
    histArray(1,N(n)+1)=histArray(1,N(n)+1)+1; % every time you meet the particular value, you add 1 into to corresponding bin
end


histArrayPerc=histArray.*100/sum(histArray); % converting into % via dividing # of cell with specific # of links by total # of links 

m=find(histArray);    % returns all indexes of non-zero elements of the array
maxNonZeroLinks=max(m);   % number of links to consider when normalizing probabilty
k=1:1:maxNonZeroLinks;            % index of # of links (starting with 0 links)
kpercent=k.*100/(maxNonZeroLinks);  % convert # of links into % of limks
histArrayPercShort = histArrayPerc(1:maxNonZeroLinks);   % cropping the hisArray after the last non-0 value

histArrayShort = histArray(1:maxNonZeroLinks);
loghist = log(histArrayShort);
xlab = [1:length(histArrayPercShort)];
xlab(isinf(loghist))=[];
loghist(isinf(loghist))=[];


[s] = corrcoef(log(xlab),loghist)

err = abs(-1-s(2))
 else
     err = 2;
 end
end