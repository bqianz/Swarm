%%
clear,close all;
load('cot_movie.mat');
load('sin_movie.mat');

%%
figure(1)
axis off;
title('cotangent potential');
movie(F,1)

%%
figure(2)
axis off;
title('sine potential');
movie(G,1)

