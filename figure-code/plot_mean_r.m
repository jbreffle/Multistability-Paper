%plot_mean_r.m

if exist('target_s') && target_s==0.5
    load('inf_logistic_alt_alpha0.2_input0.mat')
    disp('Loading 1')
else
    disp('Loading 1')
	load('inf_logistic_alt_alpha0.2_s0.mat')
end

tol = 1e-6;

[amax bmax cmax dmax] = size(eta_std);

vals = zeros([amax,bmax,3]);


for a = 1:amax
    for b = 1:bmax

        thisval = [];
        for c = 1:cmax
            for d = 1:dmax
                eta = eta_std(a,b,c,d);
                if (eta > 0 )
                    new_eta = 1;
                    for i = 1:length(thisval)
                        if ( abs(eta-thisval(i)) < tol )
                            new_eta = 0;
                        end
                    end
                    if ( new_eta == 1 )
                        thisval(end+1) = eta;
                    end
                end
            end
        end
        
        if ( length(thisval) > 0 )
            newvals = sort(thisval);
            vals(a,b,1:length(newvals)) = thisval;
        end
        
    end
    
end

%{
figure; hold on
for i = 1:3
    figure; title(num2str(i))
    inds = find(vals(1,:,i) > 0 );
    
    if i ==1
        plot(gg(inds),vals(1,inds,i)./gg(inds),'k')
    else
        plot(gg(inds),vals(1,inds,i)./gg(inds),'k', 'HandleVisibility', 'off')
    end
    
end
% title('Inf-N, logistic, alpha=0.2, s=0.5'); assert(isequal(s, 0.5)); assert(isequal(alpha, 0.2));
xlabel('g')
ylabel('Mean rate')
keyboard
%}    


%% Switch to plotting as line

figure; hold on

tmp1 = vals(:); tmp2 = [gg(:), gg(:), gg(:)]; 
inds = vals>0.2 & vals>0.0;
[a,b] = sort(tmp1(inds));
y = tmp1(inds);
x = tmp2(inds);
plot(x(b), y(b)./x(b),'k')

tmp1 = vals(:); tmp2 = [gg(:), gg(:), gg(:)]; 
inds = vals<0.2 & vals>0.0;
[a,b] = sort(tmp1(inds));
y = tmp1(inds);
x = tmp2(inds);
plot(x(b), y(b)./x(b), 'k', 'HandleVisibility', 'off')

xlabel('g')
ylabel('Mean rate')



