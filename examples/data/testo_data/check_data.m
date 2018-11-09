%% Figure 1
figure(1)
image(imread('fig_27yo.png'))
hold on

% R Data
dat = load('Fig1_data_27yo_GnRH.csv');
x0 = 206; y0 = 666;
plot(x0, y1, 'bx')
plot(x0 + 1.818*dat(:,1), y0 - 395*dat(:,2), 'ro')


% L Data
dat = load('Fig1_data_27yo_LH.csv');
x0 = 206; y0 = 1495;
plot(x0, y1, 'bx')
plot(x0 + 1.818*dat(:,1), y0 - 100*(dat(:,2)-1), 'r')


%% Figure 2
figure(2)
image(imread('fig_40yo.png'))
hold on

% R Data
dat = load('Fig1_data_40yo_GnRH.csv');
x0 = 383; y0 = 909;
plot(x0, y0, 'bx')
plot(x0 + 2.282*dat(:,1), y0 - 1050*dat(:,2), 'ro')

% L Data
dat = load('Fig1_data_40yo_LH.csv');
x0 = 383; y0 = 1944;
plot(x0, y0, 'bx')
plot(x0 + 2.282*dat(:,1), y0 - 184*(dat(:,2)-1), 'r')