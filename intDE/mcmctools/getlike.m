function getlike(varargin)

    plotdir = 'plots/';
    PriorsFileName = 'priors.txt';
    ploteachchain = 1;
    nbins = 40;
    
    RunDirID = char(varargin);
    RunDir = strcat('../chains/',RunDirID,'/');
    RunDirContents = dir(fullfile(RunDir,'*chain*'));
    dims = size(RunDirContents);
    nchains = dims(1);
    
    priorsinfo = fopen(strcat('../',PriorsFileName));
    C = textscan(priorsinfo, '%s %s %f %f %f');
    fclose(priorsinfo);
    paramsectn = C{:,1};
    paramnames = C{:,2};
    paramlower = C{:,3};
    paramupper = C{:,4};
    paramsigma = C{:,4};
    nparams = numel(paramnames);
    
    % setup spacing on the bins
    for p = 1:nparams
        min = paramlower(p);
        max = paramupper(p);
        dp = (max - min)/nbins;
        for b = 1:nbins
            val(p,b) = min+b*dp;
        end;
    end;
    
 
    
    for n = 1:nchains
        chaindata = load(strcat(RunDir,RunDirContents(n).name));
        chaindims = size(chaindata);
        nsamples = chaindims(1)-1;
        nparams = chaindims(2)-1;
        chainID = 1E4+n-1;
        
        % Get the likelihoods
        likedata = chaindata(:,nparams+1);
        
        % zero the bins
        for p=1:nparams
            for b=1:nbins
                bin(p,b)=0;
            end;
        end;
        
        % Get all parameters
        for p1=1:nparams
            
            % Get the samples data
            param1sample=chaindata(:,p1);
            
            for s=1:nsamples
                for b=1:nbins
                    if param1sample(s) >= val(p1,b) & param1sample(s) < val(p1,b+1)
                        bin(p1,b)=bin(p1,b)+1;
                    end;
                end;
            end;
            % If desired, plot the 1D likelihood from each chain
            if ploteachchain == 1
                plot(val(p1,:),bin(p1,:));
                xlim([paramlower(p1) paramupper(p1)]);
                xlabel(char(paramnames(p1)));
                outfile = strcat(plotdir,'chain_',num2str(chainID), '_',char(paramnames(p1)),'.eps');
                set(gcf, 'PaperUnits','inches');
                set(gcf, 'PaperPosition',[ 0 0 4 4]);
                print('-depsc2',outfile);
            end;
            
            for p2=p1+1:nparams
                 param2sample=chaindata(:,p2);
                 for b1=1:nbins
                     for b2=1:nbins
                         bin2d(b1,b2)=0;
                     end;
                 end;
                 for s=1:nsamples
                     for b1=1:nbins
                         for b2=1:nbins
                             if param1sample(s) >= val(p1,b1) & param1sample(s) < val(p1,b1+1) & param2sample(s) >= val(p2,b2) & param2sample(s) < val(p2,b2+1)
                                bin2d(b1,b2)=bin2d(b1,b2)+1;
                             end;
                         end;
                     end;
                 end;
                 if ploteachchain == 1
                     pcolor(val(p1,:),val(p2,:),bin2d);
                     xlim([paramlower(p1) paramupper(p1)]);
                     ylim([paramlower(p2) paramupper(p2)]);
                     xlabel(char(paramnames(p1)));
                     ylabel(char(paramnames(p2)));
                     outfile = strcat(plotdir,'surf_',num2str(chainID), '_',char(paramnames(p1)), '_',char(paramnames(p2)),'.eps');
                     set(gcf, 'PaperUnits','inches');
                     set(gcf, 'PaperPosition',[ 0 0 4 4]);
                     print('-depsc2',outfile);
                 end;
                 
            end;
            
        end; % finish looking at all parameters
        
      

        
    end; % finish reading in all the chains
    
    
    
    
    