function  [p95_4, p68_2, calprob, medage] = udmatcal(c14age, c14err, calcurve, yeartype, varargin)
% [p95_4, p68_2, calprob, medage] = matcal(c14age, c14err, calcurve, yeartype)
%
% Function for 14C age calibration using Bayesian higher posterior
% density analysis of a probability density function of calibrated age.
%
% Please see manuscript for more detailed information:
% Lougheed, B.C. & Obrochta, S.P. (2016). MatCal: Open Source Bayesian
% 14C Age Calibration in Matlab. Journal of Open Research Software. 4(1),
% p.e42. DOI: http://doi.org/10.5334/jors.130
%
% --- Required input parameters ---
%
% c14age:    Lab 14C determination in 14C yr BP.
%
% c14err:    Lab 14C determination uncertainty (1 sigma) in 14C yr.
%
% calcurve:  String specifying calibration curve to use, select from
%            the following (not case sensitive):
%            'IntCal13', 'Marine13', 'SHCal13, 'IntCal09', 'Marine09'
%            'IntCal04', 'Marine04', 'SHCal04, 'IntCal98', 'Marine98'
%
% yeartype:  String specifying how to report calibrated age.
%            Choices are 'CalBP' or 'BCE/CE'. (Not case sensitive)
%
% --- Optional input parameters ---
%
% resage:    Optional (parameter name and value). Specify reservoir
%            age in 14C yr. R(t) in the case of atmospheric calibration
%            curve, delta-R in the case of marine curve. (default = 0)
%            e.g. 'resage',320 for a reservoir age of 320
%
% reserr:    Optional (parameter name and value). Specify a 1 sigma
%            uncertainty for your chosen resage (default = 0)
%            e.g. 'reserr',50 for an uncertainty of 50
%
% plot:      Optional (parameter name and value). Return a calibration
%            plot of the to Figure 14. The plot displays the 1 and 2
%            sigma ranges of the calibration curve. The calibraiton
%            curve raw data is also shown if IntCal13 is selected.
%            Specify 1 to plot and 0 not to plot. (default = 1 for Matlab
%            users; default = 0 for Octave users)
%            e.g. 'plot',0 not to plot
%
% saveplot:  Optional (parameter name and value). Save an Adobe PDF of
%            the calibration plot to your working directory. Specify 1
%            to save and 0 not to save. (default = 0) Will be ignored
%            if plotting has been disabled.
%            e.g. 'saveplot',1 to save to your working directory.
%
% plotsize:  Optional (parameter name and value). Set the width and height
%            of the printed figure in cm. (default = 16).
%            e.g. 'plotsize',10 for 10 cm.
%
% fontsize:  Optional (parameter name and value). Set the value of the font
%            size in the output plot. (default = 8)
%            e.g. 'fontsize',12 for a font size of 12.
%
%
% --- Output data ---
%
% p95_4:     n by 3 matrix containing 95.45% calibrated age probability
%            range interval(s) calculated using highest posterior density.
%            Each row contains a probability range in Cols 1 and 2, and
%            the associated probability for that range in Col 3.
%            Probabilities are normalised to between zero and one.
%
% p68_2:     Same as p95_4, but for the 68.27% calibrated range.
%
% calprob:   n by 2 matrix containing an annualised calibrated age
%            probability density function for implementation in, e.g.,
%            age modelling. n is the annualised length of the chosen
%            calibration curve. Col 1 is a series of annual cal ages,
%            Col 2 contains their associated probability. All probabilities
%            are normalised such that they sum to 1.
%
% medage:    Median age calculated from calprob.
%
% --- Functional examples ---
%
% [p95_4 p68_2 prob] = matcal(1175, 30, 'IntCal13', 'BCE/CE')
% Calibrate a 14C age of 1175 ±30 14C yr BP using IntCal13 with output in BCE/CE.
%
%
% [p95_4 p68_2 prob] = matcal(23175, 60, 'Marine13', 'CalBP', 'resage', -50, 'reserr', 100, 'saveplot', 1)
% Calibrate a 14C age of 23175 ±60 14C yr BP using Marine13, with output in
% Cal BP, with delta-R of -50 ±100 14C yr and save a copy of the plot to
% your working directory as an Adobe PDF.
%
% [p95_4 p68_2 prob] = matcal(1175, 30, 'IntCal13', 'CalBP', 'plot', 0)
% Calibrate a 14C age of 1175±50 14C yr BP using IntCal13, with output in
% Cal BP and disable the plot window.
%
% ------------
%
% MatCal 2.4 (2018-06-25)
% Written using MatLab 2012a, compatible with 2017b.
% Please see manuscript for license information:
% http://doi.org/10.5334/jors.130

if nargin < 4
	error('Not enough input parameters (see help for instructions)')
end

matcalvers = 'MatCal 2.4 (Lougheed and Obrochta, 2016)';

udembedded = 1;

% Optional parameters input parser (parse varargin)

p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = false;
p.FunctionName='matcal';

defaultresage=0;
defaultreserr=0;
defaultplotme=1;
defaultprintme=0;
defaultplotsize=16;
defaultfontsize=8;

if exist('OCTAVE_VERSION', 'builtin') ~= 0
	defaultplotme = 0;
	addParamValue(p,'resage',defaultresage,@isnumeric);
	addParamValue(p,'reserr',defaultreserr,@isnumeric);
	addParamValue(p,'plot',defaultplotme,@isnumeric);
	addParamValue(p,'saveplot',defaultprintme,@isnumeric);
	addParamValue(p,'plotsize',defaultplotsize,@isnumeric);
	addParamValue(p,'fontsize',defaultfontsize,@isnumeric);
else
	if datenum(version('-date'))>datenum('May 19, 2013')
		addParameter(p,'resage',defaultresage,@isnumeric);
		addParameter(p,'reserr',defaultreserr,@isnumeric);
		addParameter(p,'plot',defaultplotme,@isnumeric);
		addParameter(p,'saveplot',defaultprintme,@isnumeric);
		addParameter(p,'plotsize',defaultplotsize,@isnumeric);
		addParameter(p,'fontsize',defaultfontsize,@isnumeric);
	else
		addParamValue(p,'resage',defaultresage,@isnumeric);
		addParamValue(p,'reserr',defaultreserr,@isnumeric);
		addParamValue(p,'plot',defaultplotme,@isnumeric);
		addParamValue(p,'saveplot',defaultprintme,@isnumeric);
		addParamValue(p,'plotsize',defaultplotsize,@isnumeric);
		addParamValue(p,'fontsize',defaultfontsize,@isnumeric);
	end
end

parse(p,varargin{:});
resage=p.Results.resage;
reserr=p.Results.reserr;
plotme=p.Results.plot;
printme=p.Results.saveplot;
plotsize=p.Results.plotsize;
fontsize=p.Results.fontsize;

% Cal curve case and symbols

headerlines = 11;

if strcmpi(calcurve, 'IntCal13') == 1
	calcurve = 'IntCal13';
	cite = '(Reimer et al., 2013)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine13') == 1
	calcurve = 'Marine13';
	cite = '(Reimer et al., 2013)';
	curvetype = 'mar';
elseif strcmpi(calcurve, 'SHCal13') == 1
	calcurve = 'SHCal13';
	cite = '(Hogg et al., 2013)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'IntCal09') == 1
	calcurve = 'IntCal09';
	cite = '(Reimer et al., 2009)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine09') == 1
	calcurve = 'Marine09';
	cite = '(Reimer et al., 2009)';
	curvetype = 'mar';
elseif strcmpi(calcurve, 'IntCal04') == 1
	calcurve = 'IntCal04';
	cite = '(Reimer et al., 2004)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine04') == 1
	calcurve = 'Marine04';
	cite = '(Hughen et al., 2004)';
	curvetype = 'mar';
elseif strcmpi(calcurve, 'SHCal04') == 1
	calcurve = 'SHCal04';
	cite = '(McCormac et al., 2004)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'IntCal98') == 1
	headerlines = 18;
	calcurve = 'IntCal98';
	cite = '(Stuiver et al., 1998)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'Marine98') == 1
	headerlines = 18;
	calcurve = 'Marine98';
	cite = '(Stuiver et al., 1998)';
	curvetype = 'mar';
elseif strcmpi(calcurve, 'IntCal13pCO2') == 1
	calcurve = 'IntCal13pCO2';
	cite = '(Reimer et al., 2013 + Galbraith et al. (2015) pCO2 effect)';
	curvetype = 'atm';
elseif strcmpi(calcurve, 'SHCal13pCO2') == 1
	calcurve = 'SHCal13pCO2';
	cite = '(Hogg et al., 2013 + Galbraith et al. (2015) pCO2 effect)';
	curvetype = 'atm';
else
	error(['Calibration curve "',calcurve,'" unknown. Please specify a valid calibration curve (see help for options)'])
end


if strcmp(curvetype, 'atm') == 1 && resage == 0
	extralabel = 0;
	if reserr > 0
		error('You have specified reserr without entering resage.')
	end
elseif strcmp(curvetype, 'atm') == 1 && resage > 0
	extralabel = 1;
	reslabel = 'R(t)';
	plot14Coriginal = 1;
elseif strcmp(curvetype, 'mar') == 1
	extralabel = 1;
	reslabel = '\DeltaR';
	plot14Coriginal = 1;
end


if strcmpi(yeartype, 'Cal BP') == 1 || strcmpi(yeartype, 'CalBP') == 1
	yearlabel = 'Cal yr BP';
elseif strcmpi(yeartype, 'BCE/CE') == 1
	yearlabel = 'BCE/CE';
else
	error(['Please specify a valid year type (see help for options)'])
end

% correct for reservoir age and store original ages in memory
if isnan(resage) == 1
	resage = 0;
end
if isnan(reserr) == 1;
	reserr = 0;
end

c14ageorig = c14age;
c14errorig = c14err;

c14age = c14age - resage;
c14err = sqrt(c14err^2 + reserr^2);

% 14C age in F14 space, see e.g. Bronk Ramsey (2008) doi:10.1111/j.1475-4754.2008.00394.x

f14age = exp(c14age/-8033);
f14err = f14age*c14err/8033;

% open cal curve data

File = fopen([calcurve,'.14c']);
Contents = textscan(File,'%f %f %f %f %f','headerlines',headerlines,'delimiter',',');
fclose(File);
curvecal = flipud(Contents{1});
curve14c = flipud(Contents{2});
curve14cerr = flipud(Contents{3});
curvef14 = exp(curve14c/-8033); % and also convert to F14 space (Libby half-life because IntCal in 14C years)
curvef14err = curvef14.*curve14cerr/8033; % and also convert to F14 space

% interpolate F14 cal curve to annual resolution

interpres = 1;
hicurvecal = curvecal(1):interpres:curvecal(end);
hicurvef14 = interp1(curvecal, curvef14, hicurvecal);
hicurvef14err = interp1(curvecal, curvef14err, hicurvecal);

% Calculate probability for every cal year in F14 space

calprob = NaN(length(hicurvecal),2);
calprob(:,1) = hicurvecal;

z = 0;
for i = 1:length(hicurvecal)
	
	z = z + 1;
	
	% equation from e.g. p.261 in Bronk Ramsey, 2008. doi:10.1111/j.1475-4754.2008.00394.x
	% split equation into parts for sanity's sake
	a = ( f14age - hicurvef14(i) )^2;
	b = 2 * (f14err^2 + hicurvef14err(i)^2);
	c = sqrt(f14err^2 + hicurvef14err(i)^2);
	calprob(z,2) = exp(-a/b) / c;
end

% normalise to 1
calprob(:,2) = calprob(:,2) / sum(calprob(:,2));

if udembedded ~= 1;

	% throw warning if 4sigma of 14C age exceeds 14C age limits in cal curve
	if (c14age + 4*c14err) > max(curve14c) || (c14age - 4*c14err) < min(curve14c)
		warning(['4sigma range of 14C age ',num2str(c14ageorig),char(177),num2str(c14errorig),' may exceed limits of calibration curve'])
	end
	
	% also throw warning if cal age PDF does not tail to zero at ends of cal curve (exceeds cal curve)
	if calprob(1,2) > 0.000001 || calprob(end,2) > 0.000001
		warning(['Calibrated age PDF for 14C age ',num2str(c14ageorig),char(177),num2str(c14errorig),' may exceed limits of calibration curve'])
	end

end

% find 1sig and 2sig intervals using highest posterior density (HPD)

hpd = calprob(:,1:2);

hpd = sortrows(hpd, 2);

hpd(:,3) = cumsum(hpd(:,2));

% 1 sig

hpd68_2 = hpd(hpd(:,3) >= 1-erf(1/sqrt(2)), :);

hpd68_2 = sortrows(hpd68_2,1);

ind1 = find(diff(hpd68_2(:,1)) > 1);

if isempty(ind1) == 1
	p68_2(1,1) = hpd68_2(end,1);
	p68_2(1,2) = hpd68_2(1,1);
	p68_2(1,3) = sum(hpd68_2(1:end,2));
else
	z = 0;
	for i = 1:length(ind1)
		z = z + 1;
		indy1(z,1) = ind1(i);
		z = z + 1;
		indy1(z,1) = ind1(i)+1;
	end
	indy1 = [ 1 ; indy1; length(hpd68_2(:,1)) ];
	z=0;
	for i = 2:2:length(indy1)
		z = z+1;
		p68_2(z,1) = hpd68_2(indy1(i),1);
		p68_2(z,2) = hpd68_2(indy1(i-1),1);
		p68_2(z,3) = sum(hpd68_2(indy1(i-1):indy1(i),2));
	end
	p68_2 = flipud(p68_2);
end


% 2 sig

hpd95_4 = hpd(hpd(:,3) >= 1-erf(2/sqrt(2)), :);

hpd95_4 = sortrows(hpd95_4,1);

ind2 = find(diff(hpd95_4(:,1)) > 1);

if isempty(ind2) == 1;
	p95_4(1,1) = hpd95_4(end,1);
	p95_4(1,2) = hpd95_4(1,1);
	p95_4(1,3) = sum(hpd95_4(1:end,2));
else
	z = 0;
	for i = 1:length(ind2)
		z = z + 1;
		indy2(z,1) = ind2(i);
		z = z + 1;
		indy2(z,1) = ind2(i)+1;
	end
	indy2 = [ 1 ; indy2; length(hpd95_4(:,1)) ];
	z=0;
	for i = 2:2:length(indy2)
		z = z+1;
		p95_4(z,1) = hpd95_4(indy2(i),1);
		p95_4(z,2) = hpd95_4(indy2(i-1),1);
		p95_4(z,3) = sum(hpd95_4(indy2(i-1):indy2(i),2));
	end
	p95_4 = flipud(p95_4);
end

% calculate median (can't use interp1 because of potential repeat values)
[~, median_ind] = min(abs(cumsum(calprob(:,2))-0.5));
medage = round(calprob(median_ind(1),1));

% Convert output to BC/BCE if necessary
if strcmpi(yeartype,'BCE/CE') == 1
	medage = (medage-1950) * -1;
	calprob(:,1) = (calprob(:,1)-1950) * -1;
	p95_4(:,1:2) = (p95_4(:,1:2)-1950) * -1;
	p68_2(:,1:2) = (p68_2(:,1:2)-1950) * -1;
end

%%%%% ----- Start of plotting module ----- %%%%%

if plotme==1
	
	if strcmpi(yeartype,'Cal BP') == 1 || strcmpi(yeartype,'CalBP') == 1
		calprob2 = calprob(cumsum(calprob(:,2)) > 0.001 & cumsum(calprob(:,2)) < 0.999);
		yrrng = (calprob2(end,1) - calprob2(1,1))/2;
		
		% round to nearest hundred for nice plot limits
		syr = (10^2) * round((calprob2(1,1)-yrrng) / (10^2));
		eyr = (10^2) * round((calprob2(end,1)+yrrng) / (10^2));
		
		ind = find(curvecal >= syr & curvecal <= eyr);
		curvecal = curvecal(ind);
		curve14c = curve14c(ind);
		curve14cerr = curve14cerr(ind);
		
	elseif strcmpi(yeartype,'BCE/CE') == 1
		
		calprob2 = calprob(cumsum(calprob(:,2)) > 0.001 & cumsum(calprob(:,2)) < 0.999);
		calprob2 = flipud(calprob2);
		yrrng = (calprob2(end,1) - calprob2(1,1))/2;
		
		% round to nearest hundred for nice plot limits
		syr = (10^2) * round((calprob2(1,1)-yrrng) / (10^2));
		eyr = (10^2) * round((calprob2(end,1)+yrrng) / (10^2));
		
		curvecal = (curvecal-1950) * -1;
		ind = find(curvecal >= syr & curvecal <= eyr);
		curvecal = curvecal(ind);
		curve14c = curve14c(ind);
		curve14cerr = curve14cerr(ind);
		
	end
	
	figure(14)
	clf
	
	%----- Plot ProbDistFunc
	axpdf = axes;
	axes(axpdf)
	area(calprob(:,1),calprob(:,2),'edgecolor','none')
	axpdfylims = ylim;
	axpdfxlims = xlim;
	area(calprob(:,1),calprob(:,2)*0.2,'edgecolor',[0 0 0],'facecolor',[0.9 0.9 0.9])
	hold on
	M = size(p95_4,1);
	for i = 1:M
		if strcmpi(yeartype,'Cal BP') == 1 || strcmpi(yeartype,'CalBP') == 1
			area( calprob(calprob(:,1) <= p95_4(i,1) & calprob(:,1) >= p95_4(i,2),1)  , calprob(calprob(:,1) <= p95_4(i,1) & calprob(:,1) >= p95_4(i,2),2)*0.2,'edgecolor','none','facecolor',[0.56 0.56 0.66])
		elseif strcmpi(yeartype,'BCE/CE') == 1
			area( calprob(calprob(:,1) >= p95_4(i,1) & calprob(:,1) <= p95_4(i,2),1)  , calprob(calprob(:,1) >= p95_4(i,1) & calprob(:,1) <= p95_4(i,2),2)*0.2,'edgecolor','none','facecolor',[0.56 0.56 0.66])
		end
	end
	M = size(p68_2,1);
	for i = 1:M
		if strcmpi(yeartype,'Cal BP') == 1 || strcmpi(yeartype,'CalBP') == 1
			area( calprob(calprob(:,1) <= p68_2(i,1) & calprob(:,1) >= p68_2(i,2),1)  , calprob(calprob(:,1) <= p68_2(i,1) & calprob(:,1) >= p68_2(i,2),2)*0.2,'edgecolor','none','facecolor',[0.5 0.5 0.6])
		elseif strcmpi(yeartype,'BCE/CE') == 1
			area( calprob(calprob(:,1) >= p68_2(i,1) & calprob(:,1) <= p68_2(i,2),1)  , calprob(calprob(:,1) >= p68_2(i,1) & calprob(:,1) <= p68_2(i,2),2)*0.2,'edgecolor','none','facecolor',[0.5 0.5 0.6])
		end
	end
	%----- Plot cal curve
	axcurve = axes;
	axes(axcurve)
	xdata = curvecal;
	ydata = curve14c;
	onesig = curve14cerr;
	fill([xdata' fliplr(xdata')],[ydata'+2*onesig' fliplr(ydata'-2*onesig')],[0.8 0.8 0.8],'edgecolor','none');
	hold on
	fill([xdata' fliplr(xdata')],[ydata'+onesig' fliplr(ydata'-onesig')],[0.6 0.6 0.6],'edgecolor','none');
	axcurveylims = ylim;
	axcurvexlims = [min(curvecal) max(curvecal)];
	
	%----- Plot raw data if intcal13 is selected
	if strcmpi('intcal13',calcurve) == 1
		
		axraw = axes;
		axes(axraw)
		
		rd = load('IntCal13 raw data.txt');
		
		rd_trees = rd(rd(:,1) >= 1 & rd(:,1) <= 8, :);
		rd_other = rd(rd(:,1) >= 9, :);
		
		rd_trees = rd_trees(rd_trees(:,3) <= 13900, :);
		rd_other = rd_other(rd_other(:,3) >= 13900, :);
		
		rd = [rd_trees; rd_other];
		
		if strcmpi(yeartype,'BCE/CE') == 1
			rd(:,3) = (rd(:,3)-1950) * -1;
			ind = find(rd(:,3) <= curvecal(1) & rd(:,3) >= curvecal(end));
		else
			ind = find(rd(:,3) >= curvecal(1) & rd(:,3) <= curvecal(end));
		end
		
		raw13_cal = rd(ind,3);
		raw13_calsigma = rd(ind,5);
		raw13_14c = rd(ind,6);
		raw13_14csigma = rd(ind,7);
		
		for i = 1:length(raw13_cal)
			
			% x error bars
			plot([raw13_cal(i)-raw13_calsigma(i) raw13_cal(i)+raw13_calsigma(i)],[raw13_14c(i) raw13_14c(i)],'-','color',[132/256 193/256 150/256])
			
			if i == 1
				hold on
			end
			
			% y error bars
			plot([raw13_cal(i) raw13_cal(i)],[raw13_14c(i)-raw13_14csigma(i) raw13_14c(i)+raw13_14csigma(i)],'-','color',[132/256 193/256 150/256])
			
			ymaxes(i) = raw13_14c(i)+raw13_14csigma(i);
			ymins(i) = raw13_14c(i)-raw13_14csigma(i);
			
		end
		
		axrawylims = [min(ymins) max(ymaxes)];
		axrawxlims = xlim;
		
	end
	
	%----- Plot 14C age normal distribution(s)
	axgauss = axes;
	axes(axgauss);
	
	gaussrange = [c14age-4*c14err:c14age+4*c14err];
	gauss = normpdf(gaussrange,c14age,c14err);
	
	gaussrangeorig = [c14ageorig-4*c14errorig:c14ageorig+4*c14errorig];
	gaussorig = normpdf(gaussrangeorig, c14ageorig, c14errorig);
	
	patch(gauss, gaussrange, 'blue');
	axgaussylims = ylim;
	axgaussxlims = xlim;
	axgaussxlims(2) = axgaussxlims(2)*5;
	if exist('plot14Coriginal','var') == 1
		a = patch(gaussorig, gaussrangeorig, 'blue');
		set(a,'edgecolor','none','facecolor',[0.8 0.5 0.5]);
		hold on
	end
	a = patch(gauss, gaussrange, 'blue');
	set(a,'edgecolor',[0 0 0],'facecolor',[0.7 0.4 0.4]);
	
	
	%----- set plot settings by axis, starting from back layer to front layer
	
	axes(axcurve)
	xlim(axcurvexlims)
	if strcmpi(yeartype, 'Cal BP') == 1 || strcmpi(yeartype, 'CalBP') == 1
		set(gca, 'XDir', 'reverse')
	end
	set(gca,'color','none')
	lab1 = ylabel('^1^4C yr BP');
	lab2 = xlabel(yearlabel);
	set( gca, 'TickDir', 'out' );
	if strcmpi('intcal13',calcurve) == 1
		ylim(axrawylims)
	else
		ylim(axcurveylims)
	end
	yt=get(gca,'ytick');
	ytl=textscan(sprintf('%1.0f \n',yt),'%s','delimiter','');
	set(gca,'yticklabel',ytl{1})
	xt=get(gca,'xtick');
	xtl=textscan(sprintf('%1.0f \n',xt),'%s','delimiter','');
	set(gca,'xticklabel',xtl{1})
	
	axes(axpdf)
	set(gca,'color','none')
	xlim(axcurvexlims)
	ylim(axpdfylims)
	set(gca,'xticklabel',[]);
	set(gca,'xtick',[]);
	set(gca,'yticklabel',[]);
	if strcmpi(yeartype, 'Cal BP') == 1 || strcmpi(yeartype, 'CalBP') == 1
		set(gca, 'XDir', 'reverse')
	end
	set(gca,'ytick',[]);
	
	axes(axgauss)
	set(gca,'color','none')
	xlim(axgaussxlims)
	set(gca,'xticklabel',[]);
	set(gca,'xtick',[]);
	set(gca,'yticklabel',[]);
	set(gca,'ytick',[]);
	if strcmpi('intcal13',calcurve) == 1
		ylim(axrawylims)
	else
		ylim(axcurveylims)
	end
	
	axes(axcurve) % bring curve forward
	
	if strcmpi('intcal13',calcurve) == 1
		axes(axraw)
		hold on
		if strcmpi(yeartype, 'Cal BP') == 1 || strcmpi(yeartype, 'CalBP') == 1
			set(gca, 'XDir', 'reverse')
		end
		set(gca,'color','none')
		xlim(axcurvexlims)
		ylim(axrawylims)
		set(gca,'xticklabel',[]);
		set(gca,'xtick',[]);
		set(gca,'yticklabel',[]);
		set(gca,'ytick',[]);
	else
		axes(axpdf)
	end
	
	%----- Plot some text on the final axis
	
	if strcmpi(yeartype, 'Cal BP') == 1 || strcmpi(yeartype, 'CalBP') == 1
		
		xlims = xlim;
		ylims = ylim;
		
		xaxrange = xlims(2) - xlims(1);
		yaxrange = ylims(2) - ylims(1);
		
		M = size(p95_4,1);
		
		text(xlims(2)-0.63*xaxrange, ylims(2)-0.03*yaxrange, ['^1^4C date: ',num2str(c14ageorig),' \pm ',num2str(c14errorig),' ^1^4C yr BP'])
		
		if extralabel == 1
			text(xlims(2)-0.63*xaxrange, ylims(2)-0.06*yaxrange, [reslabel,': ',num2str(resage),' \pm ',num2str(reserr), ' ^1^4C yr'])
			ypos = 0.11;
		else
			ypos = 0.08;
		end
		
		if M == 1
			text(xlims(2)-0.63*xaxrange, ylims(2)-ypos*yaxrange, ['Cal age 95.45% HPD interval:'])
		else
			text(xlims(2)-0.63*xaxrange, ylims(2)-ypos*yaxrange, ['Cal age 95.45% HPD intervals:'])
		end
		
		for i = 1:M
			ypos = ypos+0.03;
			text(xlims(2)-0.63*xaxrange, ylims(2)-ypos*yaxrange, [num2str(floor(p95_4(i,3)*1000)/10),'%: ',num2str(p95_4(i,1)),' to ',num2str(p95_4(i,2)),' cal yr BP'])
		end
		
		text(xlims(2)-0.02*xaxrange, ylims(2)-0.03*yaxrange, [matcalvers])
		text(xlims(2)-0.02*xaxrange, ylims(2)-0.06*yaxrange, [calcurve, ' ', cite])
		
	elseif strcmpi(yeartype,'BCE/CE')
		
		xlims = xlim;
		ylims = ylim;
		
		xaxrange = abs( xlims(1) - xlims(2) );
		yaxrange = ylims(2) - ylims(1);
		
		M = size(p95_4,1);
		
		text(xlims(1)+0.63*xaxrange, ylims(2)-0.03*yaxrange, ['^1^4C date: ',num2str(c14ageorig),' \pm ',num2str(c14errorig),' ^1^4C yr BP'])
		
		if extralabel == 1
			text(xlims(1)+0.63*xaxrange, ylims(2)-0.06*yaxrange, [reslabel,': ',num2str(resage),' \pm ',num2str(reserr), ' ^1^4C yr'])
			ypos = 0.11;
		else
			ypos = 0.08;
		end
		
		if M == 1
			text(xlims(1)+0.63*xaxrange, ylims(2)-ypos*yaxrange, ['Cal age 95.45% HPD interval:'])
		else
			text(xlims(1)+0.63*xaxrange, ylims(2)-ypos*yaxrange, ['Cal age 95.45% HPD intervals:'])
		end
		
		for i = 1:M
			ypos = ypos+0.03;
			text(xlims(1)+0.63*xaxrange, ylims(2)-ypos*yaxrange, [num2str(floor(p95_4(i,3)*1000)/10),'%: ',num2str(p95_4(i,1)),' to ',num2str(p95_4(i,2)),' ',yearlabel])
		end
		
		text(xlims(1)+0.02*xaxrange, ylims(2)-0.03*yaxrange, [matcalvers])
		text(xlims(1)+0.02*xaxrange, ylims(2)-0.06*yaxrange, [calcurve, ' ', cite])
		
	end
	
	if calprob(1,2) > 0.000001 || calprob(end,2) > 0.000001
		title('Warning! Age calibration may exceed limits of calibration curve.')
	elseif (c14age + 4*c14err) > max(curve14c) || (c14age - 4*c14err) < min(curve14c)
		title('Warning! Age calibration may exceed limits of calibration curve.')
	end
	
	%----- Fix fonts and appearance
	
	set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
	set(lab1,'FontWeight','bold')
	set(lab2,'FontWeight','bold')
	set(gcf,'color',[1 1 1]);
	
	%----- Prep plot for export
	if printme == 1
		
		% set figure size (cm)
		xSize = plotsize;
		ySize = plotsize;
		
		% set paper size (cm)f
		set(gcf,'PaperUnits','centimeters')
		Y = plotsize+2;
		X = plotsize+2;
		set(gcf, 'PaperSize',[X Y])
		
		% put figure in centre of paper
		xLeft = (X-xSize)/2;
		yBottom = (Y-ySize)/2;
		set(gcf,'PaperPosition',[xLeft yBottom xSize ySize])
		% make background white
		set(gcf,'InvertHardcopy','on');
		set(gcf,'color',[1 1 1]);
		
		print(figure(14), '-dpdf', '-painters', ['MatCal ',num2str(c14ageorig),char(177),num2str(c14errorig),'.pdf']);
	end
	
end


end % end function
