% erfApprox_maxk_infN.m
%
% Uses infinite-N limits for distribution of top k of N choices from a
% Gaussian to assess whether top-k are above thresold 1 for on units to
% stay on while below threshold 2 for remaining off units to stay off. 
% If that is satisfied by an fraction f = k/N of on units then a stable
% state exists.
%
% Results are for a binary f(I) which is 0 for I<1 and 1 for I<=1.
%
% This code is used to produce Figure 9B in the paper 
% "Multistability in neural systems with random cross-connections"
%
% Paul Miller, June 5, 2023 pmiller@brandeis.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gvec = [0.01:0.01:5];       % grid of values of g
ds = 0.01;
svec = [ds/2:ds:1-ds/2];    % grid of values of s

Ng = length(gvec);
Ns = length(svec);

% Initialize arrays to store results
res_grid = zeros(Ns,Ng);
minf_grid = nan(Ns,Ng);
maxf_grid = nan(Ns,Ng);

% Grid for the fraction of active units to be tested for stability
df = 0.0001;
fvals = df/2:df:1-df/2;
off_OK_array = zeros(length(fvals),Ns,Ng);  % If inactive units stay inactive
on_OK_array = zeros(length(fvals),Ns,Ng);   % If active units stay active

for i_g = 1:Ng
    g = gvec(i_g);
    disp(g)
    
    % theta_2 is threshold of Gaussian selection from field, 
    % above which inactive units become active: the following lines
    % calculate then store the result of Eq. 12 in the paper
    theta_2 = (1/g)*sqrt(1./(2*fvals) );    
    off_OK = find(fvals > 0.5*erfc(theta_2) ); % indices where true
    off_val = zeros(size(fvals));
    off_val(off_OK) = 1;                        % set to 1 for just those indices
    off_OK_array(off_OK,:,i_g) = 1;             % store in array
    
    for i_s = 1:Ns;                             % Now loop through s
        s = svec(i_s);
        
        % theta_1 is threshold of Gaussian selection from field, 
        % above which active units remain active: the following lines
        % calculate then store the result of Eq. 12 in the paper 
        theta_1 = ((1-s)/g)*sqrt(1./(2*fvals) );        
        on_OK = find(fvals < 0.5*erfc(theta_1) ); % indices where true
        on_val = zeros(size(fvals));
        on_val(on_OK) = 1;                          % set to 1 for those indices
        on_OK_array(on_OK,i_s,i_g) = 1;             % store in array
        
        all_OK = find(off_val.*on_val);             % where both conditions are met
        range_on = length(all_OK)*df;               % range of f where stability is possible
        
        res_grid(i_s,i_g) = range_on;               % store range in array
        
        if ( range_on > 0 )                         % if there is a stable state
            minf_grid(i_s,i_g) = fvals(min(all_OK));        % minimum f for stability
            maxf_grid(i_s,i_g) = fvals(max(all_OK));        % maximum f for stability
        end
        
    end
    
end

% General results for the infinite-N binary f-I curve system
figure(1)
plot(gvec,res_grid)
figure(3)
plot(gvec,(0.5*(minf_grid+maxf_grid)))

figure(2)
imagesc(squeeze(on_OK_array(:,50,:)).* squeeze(off_OK_array(:,50,:)) )
set(gca,'YDir','normal')
set(gca,'YTickLabel',{0.2 0.4 0.6 0.8 1})
set(gca,'XTickLabel',{1 2 3 4 5})
xlabel('g')
ylabel('f')

% figure(4) is used for Figure 9B in the paper
figure(4)
imagesc(squeeze(on_OK_array(:,1,:))*2 + squeeze(off_OK_array(:,1,:)) )
set(gca,'YDir','normal')

save('erfApprox_infN_maxk_grid4_Jul2022.mat')


