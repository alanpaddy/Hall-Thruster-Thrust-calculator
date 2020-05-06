%Matlab code to calculate the change in distance found using the laser sensor from
%   having the thruster turned on, then tuned of. 

% - Instuctions
%       - choose the path of the files you wish to analyse, including a backslash at the end (line 25)
%       - change the number of files you wish to analyse in line 31 (should be 3)
%       - run pragram (F5)
%       - choose 'change folder'
%       - a figure of the filtered raw data will pop up.
%       - choose the point at which the thuster was turned off by
%           right clicking on the chart
%       - two linear lines will be fitted to the data before and after the
%           point where you clicked.
%       - the y difference between the two linear lines will be calculated
%           at the x-point where the mouse was right clicked
%       - the distance calculated is plotted on a separate figure (you may
%           have to change the y-scale on figure to see it)
%       - the values of the difference calculated is stored in
%           "calculated_ydiff" file

clear all
close all
path='E:\Propulsion Lab\Data\CHT_Lab_7\';
cd(path)
%file to start at- may find low power data will not work if there is no step
firstFile=28;
%last file (this program will work out the change in distance for multiple
%files, so that you can get an average)
lastFile=30;
%laser sample rate
samplerate=576/30;%312.5;

for n=firstFile:1:lastFile;
    %load raw data--------------------------------------------
        number = num2str( n );

	file = strcat(path,'thrust_',number,'.csv'); %name of the file - name of the first file should be "thrust_1.csv"
	RawData = dlmread(file,',',3,0);
    file;
 
    %----------------------------------------------------

    %apply butterworth filter to current data set
    [b,a] = butter(1,0.1,'low'); % 
    signal_low = filtfilt(b, a, RawData(:,2));
    %plot(RawData(:,1),RawData(:,2))
    

    xTime = (linspace(0,length(signal_low)/samplerate,length(signal_low)))';  % creates x values from raw data

    %find minimum - i.e. where thruster has been turned off
    %[min_y min_x_s] = min(signal_low)

    % plot filtered data along with trendlines
    figure(n)
    plot(xTime,signal_low);
    hold on
    plot(xTime,RawData(:,2)) % plot raw data as well
    
    % maximizes figure
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1);

    %select manually points to take the thrust turned off as
    [min_x min_y] = getpts % right click to pick an x-value on each graph at which the thruster was turned off
    min_x_s = round(min_x*samplerate)
    plot(xTime,RawData(:,2), 'g')
    
    databefore = [xTime(1:min_x_s) signal_low(1:min_x_s)]; %+round(samplerate)
    dataafter  = [xTime(min_x_s:end) signal_low(min_x_s:end)]; %+round(samplerate)
    % plot(databefore(:,1), databefore(:,2))
    % plot(dataafter(:,1), dataafter(:,2))
    
    % fit first order polynomial to data before and after
    polyb4 = polyfit(databefore(:,1),databefore(:,2),1);
    polyafter = polyfit(dataafter(:,1), dataafter(:,2),1);

    % calculate value at x value corresponding to min y, for both trendline
    % before and after
    b4_valueatmin = (polyb4(:,1)*xTime(min_x_s)) + polyb4(:,2);
    after_valueatmin = (polyafter(:,1)*xTime(min_x_s)) + polyafter(:,2);
    
    % calculate difference between two y values from two trendline at minimum
    ydiff_trendlines(n) = after_valueatmin - b4_valueatmin;
    
    % Create two vectors of two trendlines
    b4_vector = [xTime (polyb4(:,1)*xTime) + polyb4(:,2)];
    after_vector = [xTime (polyafter(:,1)*xTime) + polyafter(:,2)];
    
    figure(n)
    plot(xTime,b4_vector(:,2), 'r')
    plot(xTime,after_vector(:,2), 'k')
    grid
    
end
figure(n+1)
plot(ydiff_trendlines,'o');
xlabel('data point no.');
ylabel('distance moved, mm');
%axis([0 7 0 6e-2]);

ydiff_trendlines = ydiff_trendlines';
csvwrite('calculated ydiff.csv',ydiff_trendlines)
av_ydifftrendlines = mean(ydiff_trendlines)
StDev_ydifftrendlines = std(ydiff_trendlines)
