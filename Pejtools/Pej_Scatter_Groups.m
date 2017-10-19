% This function scatters X and Y by different colors based on grouping G
% the important part is that it plots them in random order to avoid
% visualization artifacts

function Fig = Pej_Scatter_Groups(X,Y,G, Labels, Cbox, varargin)
% Fig = figure;
Fig =  gcf;
hold on

uG = sort(unique(G));
nGroups = length(uG);
RndIdx = randperm(length(G)); % break any orders in data
nBins = min(20, nGroups.^2);% plot the data in this many rounds
Bins = round(linspace(0, length(G), nBins+1));
if nargin<4 || isempty(Labels)
    if iscell(uG)
        Labels = uG;
    else
        for i = nGroups:-1:1
            Labels{i} = int2str(uG(i));
        end
    end
    
    if nargin<5
        if nGroups<=7
            Cbox = lines(nGroups);
        else
            Cbox = jet(nGroups);
        end
    end
    
    if length(varargin)==0
        varargin{1}='.';
    end
    for i = 1:nGroups
        plot(nan, nan, varargin{:}, 'DisplayName', Labels{i}, 'color', Cbox(i,:));
    end
    legend('show', 'location',  'westoutside');
    
    for R = 1:nBins
        Filt = RndIdx(Bins(R)+1:Bins(R+1)); % plot these guys
        tmpX = X(Filt);
        tmpY = Y(Filt);
        tmpG = G(Filt);
        
        if R==1
            Order = 1:nGroups;
        else
            Order = randperm(nGroups);
        end
        for i = Order
            Fmi = tmpG==uG(i);
            plot(tmpX(Fmi), tmpY(Fmi), varargin{:}, 'DisplayName', Labels{i}, 'color', Cbox(i,:));
        end
        %     if R ==1
        %         legend('show', 'location',  'westoutside');
        %     end
    end
    box on
end