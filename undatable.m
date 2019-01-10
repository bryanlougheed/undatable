function [udoutput, shadingmat] = undatable(inputfile,nsim,xfactor,bootpc,varargin)
% [udoutput, shadingmat] = undatable(inputfile,nsim,xfactor,bootpc)
%
% "Undatable" age-depth modelling software. 
% Version 1.0 (2018-12-27). For deailed description, see:
% Lougheed, B. C. and Obrochta, S. P. (2019), 
% "A rapid, deterministic age depth modelling routine for geological
% sequences with inherent depth uncertainty."
% Paleoceanography and Paleoclimatology. 
% Accepted Author Manuscript. doi:10.1029/2018PA003457
% https://doi.org/10.1029/2018PA003457
%
% ---- REQUIRED INPUT VARIABLES ----
%
% inputfile = string with location of Undatable input text
% file (e.g., 'inputfile.txt')
%
% nsim = number of iterations to run (e.g. 10^3, 10^4, 10^5)
%
% xfactor = Gaussian SAR uncertainty factor (e.g. 0.1, 0.2, etc)
%
% bootpc = Percent of age-depth constraints to bootstrap for 
% each age-depth model iteration (e.g. 10, 20, 40, etc.).
%
% ---- OUTPUT VARIABLES ----
%
% udoutput = n by 7 matrix containing age-depth model
% confidence intervals. Col 1 is depth value, Col 2 is median
% age, Col 3 is mean age, Col 4 is 2sigma lower age interval, 
% Col 5 1sigma lower age interval, Col 6 is 1sigma upper age
% interval, Col 7 is 2sigma upper age interval.
%
% shadingmat = n by 99 matrix containing the 1st to 99th
% percentiles of the age-depth model. Each row corresponds
% to the depth value of the corresponding row in the first
% column of udoutput.
%
% ---- OUTPUT TO HARD DRIVE ----
%
% The following will be saved to your working directory:
%
% An tabbed text file (yourinputname_admodel.txt), 
% containing the udoutput output variable.
%
% An Adobe PDF showing the age-depth plot will be saved,
% under the name yourinputname_admodel.pdf
%
% ---- OPTIONAL COMMANDS ----
%
% 'combine': Sum age PDFs with identical depth intervals
% interval. 1 = Yes, 0 = No. (e.g.: 'combine',1 ) Default = 1.
%
% 'savemat' = Save a .mat file (yourinputname_output.mat)
% containing all the variables produced. 1 = Yes, 0 = No.
% (e.g.: 'savemat',1 ) default = 0.
%
% 'writedir': Change the directory to which output files
% will by saved. (e.g.: 'writedir','myharddrive/somefolder/')
% Default is your working directory.
%
% 'debug' = Plot all age-depth model points used to make cloud
% 1 = Yes, 0 = No. (e.g.: 'debug',1') default = 0.
%
% 'plotme' = Enable/disable the plot window. 1 = Yes, 0 = No.
% (e.g.: 'plotme',0) default = 1.
%
% 'printme' = Enable/disable the saving of the Adobe PDF file.
% 1 = Yes, 0 = No. (e.g.: 'printme',0) default = 1. Will revert
% to 0 if the plot window is disabled.

%---INPUT PARSER
p=inputParser;
p.KeepUnmatched=true;
p.CaseSensitive=false;
defaultcombine = 1;
defaultplotme = 1;
defaultprintme = 1;
defaultsavemat = 0;
defaultdebug = 0;
defaultwritedir = '';
defaultguimode = 0;
defaultallowreversal = 0;
defaultrun1nsim = 2000;
if datenum(version('-date')) > datenum('May 19, 2013')
	addParameter(p,'combine',defaultcombine,@isnumeric);
	addParameter(p,'plotme',defaultplotme,@isnumeric);
	addParameter(p,'printme',defaultprintme,@isnumeric);
	addParameter(p,'savemat',defaultsavemat,@isnumeric);
	addParameter(p,'debug',defaultdebug,@isnumeric);
	addParameter(p,'writedir',defaultwritedir,@isstr);
	addParameter(p,'guimode',defaultguimode,@isnumeric);
	addParameter(p,'allowreversal',defaultallowreversal,@isnumeric);
	addParameter(p,'run1nsim',defaultrun1nsim,@isnumeric);
else
	addParamValue(p,'combine',defaultcombine,@isnumeric);
	addParamValue(p,'plotme',defaultplotme,@isnumeric);
	addParamValue(p,'printme',defaultprintme,@isnumeric);
	addParamValue(p,'savemat',defaultsavemat,@isnumeric);
	addParamValue(p,'debug',defaultdebug,@isnumeric);
	addParamValue(p,'writedir',defaultwritedir,@isstr);
	addParamValue(p,'guimode',defaultguimode,@isnumeric);
	addParamValue(p,'allowreversal',defaultallowreversal,@isnumeric);
	addParamValue(p,'run1nsim',defaultrun1nsim,@isnumeric);
end
parse(p,varargin{:});
depthcombine=p.Results.combine;
plotme=p.Results.plotme;
printme=p.Results.printme;
savemat=p.Results.savemat;
debugme=p.Results.debug;
writedir = p.Results.writedir;
guimode = p.Results.guimode;
allowreversal = p.Results.allowreversal;
run1nsim = p.Results.run1nsim;

% Check bootpc
if bootpc < 0
	bootpc = 0;
elseif bootpc >= 100
	error('It is not possible to bootstrap 100% or more, please select lower bootpc');
end

%---GET AND SORT INPUT DATA
[datelabel, depth1, depth2, depth, age, ageerr, datetype, calcurve, resage, reserr, dateboot] = udgetdata(inputfile);

%---MAKE AGE AND DEPTH PDFs
[medians, p68_2, p95_4, probtoplot, rundepth, rundepth1, rundepth2, rundepthpdf, runprob2sig, runboot, runncaldepth, udrunshuffle] = udmakepdfs(depth, depth1, depth2, age, ageerr, calcurve, resage, reserr, dateboot, depthcombine);

%---RUN THE AGE DEPTH LOOPS
if mean(depth2 - depth1) ~= 0
	% run1nsim = 2000; % default now set in input parser
	%message1 = ['Depth uncertainty detected: Running preliminary Monte Carlo age-depth loops'];
else
	run1nsim = nsim;
	%message1 = 'Running the Monte Carlo age-depth loops';
end

% run age depth loop for the first time (without anchors)
% disp(message1)
% disp(' ')
%save('testinput.mat','run1nsim', 'bootpc', 'xfactor', 'rundepth', 'rundepth1', 'rundepth2', 'rundepthpdf', 'runprob2sig', 'runboot', 'runncaldepth')
%error('stop here')
[agedepmat] = udrun(run1nsim, bootpc, xfactor, rundepth, rundepth1, rundepth2, rundepthpdf, runprob2sig, runboot, runncaldepth, udrunshuffle, allowreversal);

% summarise the data
interpinterval = 1;
depthstart = depth(1);
depthend = depth(end);
[summarymat, shadingmat, depthrange] = udsummary(depthstart, depthend, run1nsim, agedepmat, interpinterval, inputfile, writedir, bootpc, xfactor, depthcombine);

if mean(depth2 - depth1) ~= 0
	% make anchors based on first run
	% disp('Anchoring end points')
	% disp(' ')
	[rundepth, rundepth1, rundepth2, rundepthpdf, runprob2sig, runboot, runncaldepth] = udanchors(depthrange, depth, depth1, depth2, summarymat, rundepth, rundepth1, rundepth2, rundepthpdf, runprob2sig, runboot);
	% run age depth loop for the second time with anchors
	% disp('Running the final Monte Carlo age-depth loops with anchors')
	% disp(' ')
	
	[agedepmat] = udrun(nsim, bootpc, xfactor, rundepth, rundepth1, rundepth2, rundepthpdf, runprob2sig, runboot, runncaldepth, udrunshuffle, allowreversal);
	% summarise the data
	[summarymat, shadingmat, depthrange] = udsummary(depthstart, depthend, nsim, agedepmat, interpinterval, inputfile, writedir, bootpc, xfactor, depthcombine);
end

if sum(isnan(agedepmat(:,1,:))) == numel(agedepmat(:,1,:));
	warning('All age-depth runs failed to produce an age-depth relationship. Check run settings (bootpc too high/low?). Also, your data may not be suitable or contain typos.')
end

udoutput = [depthrange summarymat(:,1) summarymat(:,6) summarymat(:,2:5)]; % summarymat is: median, 2siglo, 1siglo, 1sighi, 2sighi, mean

% Save output (if savemat selected)
if savemat == 1 && guimode == 0
	[~,writename,~] = fileparts(inputfile);
	save([writedir,'/',writename,'.mat'])
elseif guimode == 1
	save([writedir,'/guitemp.mat'])
end


%---PLOT STUFF
if plotme == 1	
	figure(1)
	clf
	hold(gca,'on')
	
	% Plot density cloud
	for i = 1:49
		
		hi1sig=shadingmat(:,99-i);
		lo1sig=shadingmat(:,i);
		
		confx=[
			%right to left at top
			lo1sig(1); hi1sig(1);
			%down to bottom
			hi1sig(2:end);
			%left to right at bottom
			lo1sig(end);
			%up to top
			flipud(lo1sig(1:end-1))];
		
		confy=[
			%right to left at top
			depthrange(1,1); depthrange(1,1);
			%down to bottom
			depthrange(2:end,1);
			%left to right at bottom
			depthrange(end,1);
			%up to top
			flipud(depthrange(1:end-1,1))];
		
		patch(confx/1000,confy,[1-(i/49) 1-(i/49) 1-(i/49)],'edgecolor','none')
	end
	plot(summarymat(:,2)/1000,depthrange,'k--') % 95.4 range
	plot(summarymat(:,5)/1000,depthrange,'k--') % 95.4 range
	plot(summarymat(:,3)/1000,depthrange,'b--') % 68.2 range
	plot(summarymat(:,4)/1000,depthrange,'b--') % 68.2 range
	
	plot(summarymat(:,1)/1000,depthrange,'r')
	set(gca,'ydir','reverse')
	
	% PLOT PDFs
	% Colour schemes
	%                    1                   2                       3               4           5        6           7                8                9
	datetypes = {'14C marine fossil'; '14C terrestrial fossil'; '14C sediment'; 'Tephra'; 'Tie point'; 'Other'; 'Palaeomagnetism'; 'Paleomagnetism'; 'U/Th'};
	colours(:,:,1) = [41  128  185 ; 166 208 236]; % dark blue  ; light blue
	colours(:,:,2) = [34  153 85   ; 166 219 175]; % dark green ; light green
	colours(:,:,3) = [83 57  47    ; 191 156 145]; % dark brown ; light brown
	colours(:,:,4) = [192 57   43  ; 228 148 139]; % dark red   ; light red
	colours(:,:,5) = [254 194 0    ; 241 234 143]; % dark yellow; light yellow
	colours(:,:,6) = [64  64  64   ; 160 160 160]; % dark grey  ; light grey
	colours(:,:,7) = [201 106 18   ; 221 163 108]; % dark orange; light orange
	colours(:,:,8) = [201 106 18   ; 221 163 108]; % dark orange; light orange
	colours(:,:,9) = [131 39  147  ; 194 128 206]; % dark purple; light purple
	colours = colours/255; % RGB to Matlab
	
	d=ylim;
	probscale=.015;
	usedcolours = NaN(size(depth));
	
	% set up order in which to plot dates so that a clould of bulk dates with large uncertainty doesn't obsure tephras etc.
	plotorder = [];
	plotorder = [plotorder; find(strcmpi(datetype,'14C sediment'))];
	plotorder = [plotorder; find(strcmpi(datetype,'14C marine fossil'))];
	plotorder = [plotorder; find(strcmpi(datetype,'14C terrestrial fossil'))];
	plotorder = [plotorder; find(strcmpi(datetype,'tephra'))];
	plotorder = [plotorder; find(strcmpi(datetype,'tie point'))];
	plotorder = [plotorder; find(strcmpi(datetype,'other'))];
	plotorder = [plotorder; find(strcmpi(datetype,'palaeomagnetism'))];
	plotorder = [plotorder; find(strcmpi(datetype,'paleomagnetism'))];
	plotorder = [plotorder; find(strcmpi(datetype,'u/th'))];
	% in case user made a typo
	ndepth = 1:length(depth);
	typofields = ndepth(~ismember(ndepth,plotorder));
	plotorder = [plotorder; typofields'];
	
	for i = 1:length(depth)
		
		colourind = find(strcmpi(datetype{plotorder(i)},datetypes)==1);
		if isempty(colourind) == 1
			warning(['Date type "',datetype{plotorder(i)},'" unknown. Will use grey for plot colour.'])
			colourind = 6;
		end
		
		usedcolours(i) = colourind; % for legend later
		
		probnow = probtoplot{plotorder(i)};
		
		if ageerr(i) > 0
			
			% 2 sigma shading
			for j=1:size(p95_4{plotorder(i)},1)
				probindtemp=probnow(probnow(:,1)>p95_4{plotorder(i)}(j,2) & probnow(:,1)<p95_4{plotorder(i)}(j,1),:);
				probx=[probindtemp(:,1)/1000; flipud(probindtemp(:,1)/1000)];
				proby=[probindtemp(:,2)*((d(2)*probscale)/max(probnow(:,2)))+depth(plotorder(i),1); flipud(-1*probindtemp(:,2)*((d(2)*probscale)/max(probnow(:,2)))+depth(plotorder(i),1))];
				h1=patch(probx,proby,colours(2,:,colourind),'edgecolor','none');
			end
			
			% 1 sigma shading
			for j=1:size(p68_2{plotorder(i)},1)
				probindtemp=probnow(probnow(:,1)>p68_2{plotorder(i)}(j,2) & probnow(:,1)<p68_2{plotorder(i)}(j,1),:);
				probx=[probindtemp(:,1)/1000; flipud(probindtemp(:,1)/1000)];
				proby=[probindtemp(:,2)*((d(2)*probscale)/max(probnow(:,2)))+depth(plotorder(i),1); flipud(-1*probindtemp(:,2)*((d(2)*probscale)/max(probnow(:,2)))+depth(plotorder(i),1))];
				h1=patch(probx,proby,colours(1,:,colourind),'edgecolor','none');
			end
			
			% Outline
			probx=[probnow(:,1)/1000; flipud(probnow(:,1)/1000)];
			proby=[probnow(:,2)*((d(2)*probscale)/max(probnow(:,2)))+depth(plotorder(i),1); flipud(-1*probnow(:,2)*((d(2)*probscale)/max(probnow(:,2)))+depth(plotorder(i),1))];
			h1=patch(probx,proby,[1 1 1],'edgecolor','k','facecolor','none','linewidth',.2);
		
		elseif ageerr(i) == 0
			plot(age(i)/1000,depth(i),'kd','markeredgecolor',colours(1,:,colourind),'markerfacecolor',colours(1,:,colourind));
		end
		
	end
	plot(summarymat(:,1)/1000,depthrange,'r-') % median of age model runs (plot again on top of PDFs)
	set(gca,'ydir','reverse','tickdir','out','fontsize',12,'box','on')
	ylabel('Depth (cm)')
	xlabel('Calendar age (ka BP)')
	%title(strrep(inputfile,'.txt',''))
	grid on
	
	% plot the depth error bars
	for i = 1:length(depth)
		if depth1(i) <= depth2(i)
			plot( [medians(i)/1000 medians(i)/1000] , [depth1(i) depth2(i)], 'k-' )
		elseif depth1(i) > depth2(i)
			plot( [medians(i)/1000 medians(i)/1000] , [depth1(i)+depth2(i) depth1(i)-depth2(i)], 'k-' )
		end
	end
	
	set(gca, 'Layer', 'Top')
	
	% title
	if guimode == 0
		[~,NAME,~] = fileparts(inputfile); 
		NAME = strrep(NAME,'.txt','');
		NAME = strrep(NAME,'_udinput','');
		title(NAME);
	end
	

	% plot all the agedepth runs (debug mode)
	if debugme == 1
		for i = 1:size(agedepmat,3)
			plot(agedepmat(:,1,i)/1000,agedepmat(:,2,i),'r.','markersize',2)
			hold on
		end
		% disp(['Age reversals in agedepmat: ',num2str(sum(find(diff(agedepmat(:,1,:))) < 0))])
		% disp(['Depth reversals in agedepmat: ',num2str(sum(find(diff(agedepmat(:,2,:))) < 0))])
		% disp(['Median age reverals in summarymat: ',num2str(sum(find(diff(summarymat(:,1)) < 0)))])
	end
	
	udplotoptions % call udplotoptions.m to workspace
	
	% set paper size (cm)f
	set(gcf,'PaperUnits','centimeters')
	set(gcf, 'PaperSize',[plotwidth plotheight])
	% put figure in top left of paper
	set(gcf,'PaperPosition',[0 0 plotwidth plotheight])
	% make background white
	set(gcf,'InvertHardcopy','on');
	set(gcf,'color',[1 1 1]);
	
	% automatic legend (do last before printing so that it appears in the correct place)
	usedcolours = unique(usedcolours);
	txtxpos = min(xlim) + 0.97 * abs((max(xlim)-min(xlim)));
	txtypos = min(ylim) + 0.03 * abs((max(ylim)-min(ylim)));
	txtyinc = 0.04 * abs((max(ylim)-min(ylim)));
	if length(usedcolours) > 1
		for i = 1:length(usedcolours)
			if i > 1; txtypos = txtypos + txtyinc; end
			h = text(txtxpos, txtypos, strrep(datetypes{usedcolours(i)},'14C','^1^4C'),'horizontalalignment','right','color',colours(1,:,usedcolours(i)));
			set(h,'fontweight','bold')
		end
	end
	
	% set all fonts
	set(findall(gcf,'-property','FontSize'),'FontSize',textsize)
		
	% print the xfactor and bootpc to the bottom left corner
	txtxpos = min(xlim) + 0.03 * abs((max(xlim)-min(xlim)));
	txtypos = min(ylim) + 0.97 * abs((max(ylim)-min(ylim)));
	txtyinc = 0.04 * abs((max(ylim)-min(ylim)));
	text(txtxpos,txtypos-txtyinc,['xfactor = ',num2str(xfactor,'%.2g')],'FontSize',textsize)
	text(txtxpos,txtypos,['bootpc = ',num2str(bootpc,'%.2g')],'FontSize',textsize)
	
	
	
	% print
	if printme == 1		
		savename = strrep(inputfile,'.txt','_admodel.pdf');
		[~,NAME,EXT] = fileparts(savename);
		savename = [NAME,EXT];
		savename = [writedir,'/',savename];
		print(gcf, '-dpdf', '-painters', savename);
	end
end

end % end function

